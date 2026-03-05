/*
*	Copyright (C) 2025 Gadi Lahav (gadi@rfwithcare.com)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Directional mode-matched waveguide absorber.
 *
 * Algorithm (runs co-located with Mur ABC at the port boundary):
 *   1. DoPostVoltageUpdates: compute a_E = overlap(E, mode_E)
 *   2. DoPostCurrentUpdates: compute a_H = overlap(H, mode_H)
 *   3. Apply2Current: a_bwd = (a_E - sign*Zw*a_H)/2, subtract from E and H
 *
 * The directional decomposition separates forward/backward waves using
 * the wave impedance relationship.  Only the backward (reflected) wave
 * is subtracted, so the forward (transmitted) wave passes through.
 *
 * Uses ENG_DISPATCH_ARGS macro for SSE/basic engine dispatch.
 */

#include "engine_ext_modeabsorb.h"
#include "operator_ext_modeabsorb.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"

Engine_Ext_ModeAbsorb::Engine_Ext_ModeAbsorb(Operator_Ext_ModeAbsorb* op_ext) :
	Engine_Extension(op_ext)
{
	m_Op_MA = op_ext;

	m_ny   = m_Op_MA->m_ny;
	m_nyP  = m_Op_MA->m_nyP;
	m_nyPP = m_Op_MA->m_nyPP;

	for (int i = 0; i < 3; ++i)
		m_posStart[i] = m_Op_MA->m_sheetX0[i];

	m_numLines_E[0] = m_Op_MA->m_numLines_E[0];
	m_numLines_E[1] = m_Op_MA->m_numLines_E[1];
	m_E_OverlapW[0]  = m_Op_MA->m_E_OverlapW[0];
	m_E_OverlapW[1]  = m_Op_MA->m_E_OverlapW[1];
	m_E_SubtractW[0] = m_Op_MA->m_E_SubtractW[0];
	m_E_SubtractW[1] = m_Op_MA->m_E_SubtractW[1];

	m_numLines_H[0] = m_Op_MA->m_numLines_H[0];
	m_numLines_H[1] = m_Op_MA->m_numLines_H[1];
	m_H_OverlapW[0]  = m_Op_MA->m_H_OverlapW[0];
	m_H_OverlapW[1]  = m_Op_MA->m_H_OverlapW[1];
	m_H_SubtractW[0] = m_Op_MA->m_H_SubtractW[0];
	m_H_SubtractW[1] = m_Op_MA->m_H_SubtractW[1];

	m_a_E = 0;
	m_a_H = 0;
	m_ZWave = m_Op_MA->m_ZWave;
	m_dirSign = m_Op_MA->m_dirSign;

	SetNumberOfThreads(1);
}

Engine_Ext_ModeAbsorb::~Engine_Ext_ModeAbsorb()
{
}

// ============================================================================
// DoPostVoltageUpdates: compute a_E from E after the FDTD voltage update.
// ============================================================================

template <typename EngType>
void Engine_Ext_ModeAbsorb::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	if (threadID != 0) return;

	unsigned int pos[3] = {0, 0, 0};
	pos[m_ny] = m_posStart[m_ny];

	m_a_E = 0;

	for (unsigned int posP = 0; posP < m_numLines_E[0]; ++posP)
	{
		pos[m_nyP] = m_posStart[m_nyP] + posP;
		for (unsigned int posPP = 0; posPP < m_numLines_E[1]; ++posPP)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + posPP;
			m_a_E += eng->EngType::GetVolt(m_nyP, pos)  * m_E_OverlapW[0][posP][posPP];
			m_a_E += eng->EngType::GetVolt(m_nyPP, pos) * m_E_OverlapW[1][posP][posPP];
		}
	}
}

void Engine_Ext_ModeAbsorb::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

// ============================================================================
// DoPostCurrentUpdates: compute a_H from H after the FDTD current update.
// ============================================================================

template <typename EngType>
void Engine_Ext_ModeAbsorb::DoPostCurrentUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	if (threadID != 0) return;

	unsigned int pos[3] = {0, 0, 0};
	pos[m_ny] = m_posStart[m_ny];

	m_a_H = 0;

	for (unsigned int posP = 0; posP < m_numLines_H[0]; ++posP)
	{
		pos[m_nyP] = m_posStart[m_nyP] + posP;
		for (unsigned int posPP = 0; posPP < m_numLines_H[1]; ++posPP)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + posPP;
			m_a_H += eng->EngType::GetCurr(m_nyP, pos)  * m_H_OverlapW[0][posP][posPP];
			m_a_H += eng->EngType::GetCurr(m_nyPP, pos) * m_H_OverlapW[1][posP][posPP];
		}
	}
}

void Engine_Ext_ModeAbsorb::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}

// ============================================================================
// Apply2Current: subtract backward wave from both E and H.
//
// Directional decomposition:
//   a_fwd = (a_E + sign*Zw*a_H) / 2   (forward wave)
//   a_bwd = (a_E - sign*Zw*a_H) / 2   (backward / reflected wave)
//
// Only a_bwd is subtracted, leaving the forward wave intact.
// ============================================================================

template <typename EngType>
void Engine_Ext_ModeAbsorb::Apply2CurrentImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	if (threadID != 0) return;

	double a_bwd = 0.5 * (m_a_E - m_dirSign * m_ZWave * m_a_H);

	unsigned int pos[3] = {0, 0, 0};
	pos[m_ny] = m_posStart[m_ny];

	// Subtract backward wave E-field component
	for (unsigned int posP = 0; posP < m_numLines_E[0]; ++posP)
	{
		pos[m_nyP] = m_posStart[m_nyP] + posP;
		for (unsigned int posPP = 0; posPP < m_numLines_E[1]; ++posPP)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + posPP;

			eng->EngType::SetVolt(m_nyP, pos,
				eng->EngType::GetVolt(m_nyP, pos) - a_bwd * m_E_SubtractW[0][posP][posPP]);
			eng->EngType::SetVolt(m_nyPP, pos,
				eng->EngType::GetVolt(m_nyPP, pos) - a_bwd * m_E_SubtractW[1][posP][posPP]);
		}
	}

	// Subtract backward wave H-field component
	// H_bwd = -dirSign * a_bwd / Zw * mode_H
	// Subtraction: H -= H_bwd = H + dirSign * a_bwd / Zw * mode_H
	double h_coeff = m_dirSign * a_bwd / m_ZWave;

	for (unsigned int posP = 0; posP < m_numLines_H[0]; ++posP)
	{
		pos[m_nyP] = m_posStart[m_nyP] + posP;
		for (unsigned int posPP = 0; posPP < m_numLines_H[1]; ++posPP)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + posPP;

			eng->EngType::SetCurr(m_nyP, pos,
				eng->EngType::GetCurr(m_nyP, pos) + h_coeff * m_H_SubtractW[0][posP][posPP]);
			eng->EngType::SetCurr(m_nyPP, pos,
				eng->EngType::GetCurr(m_nyPP, pos) + h_coeff * m_H_SubtractW[1][posP][posPP]);
		}
	}
}

void Engine_Ext_ModeAbsorb::Apply2Current(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2CurrentImpl, threadID);
}
