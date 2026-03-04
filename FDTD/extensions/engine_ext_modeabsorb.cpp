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
 * Mode-matched waveguide absorber engine extension.
 *
 * At each timestep:
 * 1. Compute the overlap integral of the field with the mode template
 * 2. Subtract mode_shape * overlap_coefficient from the field
 *
 * FDTD stores voltages (V = E * edge_length) and currents (I = H * edge_length),
 * not field values directly. The operator pre-computes two weight arrays:
 *   OverlapW = mode_norm * area / edgeLen  (converts Volt to field*mode*area)
 *   SubtractW = mode_norm * edgeLen        (converts field correction to Volt)
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
	m_dcBlockBeta = m_Op_MA->m_dcBlockBeta;
	m_a_H_dc = 0;

	SetNumberOfThreads(1);
}

Engine_Ext_ModeAbsorb::~Engine_Ext_ModeAbsorb()
{
}

// ============================================================================
// E-field absorption
// ============================================================================

// DoPostVoltageUpdates: E-field absorption disabled.
// The Mur ABC at the simulation boundaries already handles E-field absorption.
// Only H-field absorption is needed (the H-field doesn't decay with Mur ABC alone).
template <typename EngType>
void Engine_Ext_ModeAbsorb::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	(void)eng;
	(void)threadID;
}

void Engine_Ext_ModeAbsorb::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

// Apply2Voltages: E-field absorption disabled (handled by Mur ABC).
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
// H-field absorption
// ============================================================================

// DoPostCurrentUpdates: compute the H-field overlap integral with the mode
template <typename EngType>
void Engine_Ext_ModeAbsorb::DoPostCurrentUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	if (threadID != 0) return;

	unsigned int pos[3] = {0, 0, 0};
	pos[m_ny] = m_posStart[m_ny];

	double a_H_raw = 0;

	for (unsigned int posP = 0; posP < m_numLines_H[0]; ++posP)
	{
		pos[m_nyP] = m_posStart[m_nyP] + posP;
		for (unsigned int posPP = 0; posPP < m_numLines_H[1]; ++posPP)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + posPP;

			// a_H += Curr * (mode_norm * area / edgeLen) = H * mode_norm * area
			a_H_raw += eng->EngType::GetCurr(m_nyP, pos)  * m_H_OverlapW[0][posP][posPP];
			a_H_raw += eng->EngType::GetCurr(m_nyPP, pos) * m_H_OverlapW[1][posP][posPP];
		}
	}

	// DC-blocking filter: track the DC component and subtract it.
	// Only the AC (time-varying) part represents traveling waves to absorb.
	m_a_H_dc += m_dcBlockBeta * (a_H_raw - m_a_H_dc);
	m_a_H = a_H_raw - m_a_H_dc;
}

void Engine_Ext_ModeAbsorb::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}

// Apply2Current: subtract the mode-matched H-field component
template <typename EngType>
void Engine_Ext_ModeAbsorb::Apply2CurrentImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	if (threadID != 0) return;

	unsigned int pos[3] = {0, 0, 0};
	pos[m_ny] = m_posStart[m_ny];

	for (unsigned int posP = 0; posP < m_numLines_H[0]; ++posP)
	{
		pos[m_nyP] = m_posStart[m_nyP] + posP;
		for (unsigned int posPP = 0; posPP < m_numLines_H[1]; ++posPP)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + posPP;

			// Curr -= a_H * (mode_norm * edgeLen) = a_H * mode_norm * edgeLen
			// -> H_new = H - a_H * mode_norm (correct field subtraction)
			eng->EngType::SetCurr(m_nyP, pos,
				eng->EngType::GetCurr(m_nyP, pos) - m_a_H * m_H_SubtractW[0][posP][posPP]);
			eng->EngType::SetCurr(m_nyPP, pos,
				eng->EngType::GetCurr(m_nyPP, pos) - m_a_H * m_H_SubtractW[1][posP][posPP]);
		}
	}
}

void Engine_Ext_ModeAbsorb::Apply2Current(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2CurrentImpl, threadID);
}
