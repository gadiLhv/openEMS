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
 * Directional mode-matched waveguide absorber (CST-style port absorber).
 *
 * At each timestep:
 * 1. After E update: compute a_E = overlap integral of E-field with E-mode
 * 2. After H update: compute a_H = overlap integral of H-field with H-mode
 * 3. Decompose into forward/backward waves using wave impedance:
 *      a_bwd = (a_E - sign * Z_w * a_H) / 2
 * 4. Subtract ONLY the backward (reflected) wave from both E and H fields
 *
 * This gives a perfectly matched port absorber that works at all frequencies
 * including DC, without requiring Mur ABC or any external boundary condition.
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
// E-field overlap: compute a_E after voltage update
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

			// a_E += Volt * (mode_norm * area / edgeLen) = E * mode_norm * area
			m_a_E += eng->EngType::GetVolt(m_nyP, pos)  * m_E_OverlapW[0][posP][posPP];
			m_a_E += eng->EngType::GetVolt(m_nyPP, pos) * m_E_OverlapW[1][posP][posPP];
		}
	}
}

void Engine_Ext_ModeAbsorb::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

// Apply2Voltages: nothing to do here; E subtraction happens in Apply2Current
// after we have both a_E and a_H for the directional decomposition.
template <typename EngType>
void Engine_Ext_ModeAbsorb::Apply2VoltagesImpl(EngType* eng, int threadID)
{
	(void)eng;
	(void)threadID;
}

void Engine_Ext_ModeAbsorb::Apply2Voltages(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, threadID);
}

// ============================================================================
// H-field overlap: compute a_H after current update
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

			// a_H += Curr * (mode_norm * area / edgeLen) = H * mode_norm * area
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
// Apply2Current: subtract the backward (reflected) wave from BOTH H and E
// ============================================================================

template <typename EngType>
void Engine_Ext_ModeAbsorb::Apply2CurrentImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	if (threadID != 0) return;

	// Directional decomposition:
	//   a_fwd = (a_E + sign * Z_w * a_H) / 2   (forward/incident wave)
	//   a_bwd = (a_E - sign * Z_w * a_H) / 2   (backward/reflected wave)
	// We subtract only a_bwd from both fields.
	double a_bwd = 0.5 * (m_a_E - m_dirSign * m_ZWave * m_a_H);

	unsigned int pos[3] = {0, 0, 0};
	pos[m_ny] = m_posStart[m_ny];

	// Subtract backward wave E-field component (voltage correction)
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

	// Subtract backward wave H-field component (current correction)
	// For backward wave, H_bwd has opposite sign to forward wave's H.
	// So: Curr_new = Curr + sign * (a_bwd / Z_w) * H_SubtractW
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
