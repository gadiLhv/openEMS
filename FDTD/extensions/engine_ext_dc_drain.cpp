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

#include "engine_ext_dc_drain.h"
#include "operator_ext_dc_drain.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"

#include <cstdlib>
#include <iostream>

// Priority just below the absorbing_bc sheet (default 0) so that, in the
// voltage-apply phase (highest priority first), the damper runs AFTER the ABC
// has written its boundary voltages. Kept well above the excitation (-1000).
#define ENG_EXT_PRIO_DC_DRAIN (-100)

Engine_Ext_DC_Drain::Engine_Ext_DC_Drain(Operator_Ext_DC_Drain* op_ext) :
	Engine_Extension(op_ext)
{
	m_Op_drain = op_ext;
	m_Op       = op_ext->GetOperator();

	m_ny         = op_ext->m_ny;
	m_nyP        = op_ext->m_nyP;
	m_nyPP       = op_ext->m_nyPP;
	m_normalSign = op_ext->m_normalSign;
	m_alpha      = op_ext->m_alpha;
	m_cells      = op_ext->m_cells;

	for (int n = 0; n < 3; ++n)
	{
		m_posStart[n] = op_ext->m_posStart[n];
		m_posStop[n]  = op_ext->m_posStop[n];
	}
	m_numLines[0] = m_posStop[m_nyP]  - m_posStart[m_nyP]  + 1;
	m_numLines[1] = m_posStop[m_nyPP] - m_posStart[m_nyPP] + 1;

	SetPriority(ENG_EXT_PRIO_DC_DRAIN);

	// One boundary sheet -> single-threaded, matching the absorbing_bc sheet.
	SetNumberOfThreads(1);
}

Engine_Ext_DC_Drain::~Engine_Ext_DC_Drain()
{
}

template <typename EngType>
void Engine_Ext_DC_Drain::Apply2VoltagesImpl(EngType* eng, int threadID)
{
	if (m_Eng == NULL) return;
	// One boundary, one thread.
	if (threadID != 0) return;

	const int ny = m_ny, nyP = m_nyP, nyPP = m_nyPP;
	const int p0   = (int)m_posStart[ny];
	const int nLny = (int)m_Op->GetNumberOfLines(ny);
	const double keep = 1.0 - m_alpha;

	// Optional diagnostic: how much energy sits in the backing cells and how it
	// evolves. Enable with OPENEMS_DC_DRAIN_DEBUG.
	static const bool dbg = (getenv("OPENEMS_DC_DRAIN_DEBUG") != NULL);
	double dbg_E = 0.0;
	unsigned int dbg_n = 0;

	unsigned int pos[3] = {0,0,0};

	// Damp the field in the mesh cell(s) behind the sheet (opposite the normal),
	// i.e. at p0 - k*m_normalSign, over the sheet's transverse extent.
	for (int k = 1; k <= m_cells; ++k)
	{
		int q = p0 - k * m_normalSign;
		if (q < 0 || q >= nLny)
			continue;
		pos[ny] = (unsigned int)q;

		for (unsigned int i = 0; i < m_numLines[0]; ++i)
		{
			pos[nyP] = m_posStart[nyP] + i;
			for (unsigned int j = 0; j < m_numLines[1]; ++j)
			{
				pos[nyPP] = m_posStart[nyPP] + j;

				for (int d = 0; d < 3; ++d)
				{
					double v = eng->EngType::GetVolt(d, pos);
					if (dbg) { dbg_E += v * v; ++dbg_n; }
					eng->EngType::SetVolt(d, pos, v * keep);
				}
			}
		}
	}

	if (dbg)
	{
		// For comparison, sum |V|^2 in the mirror cell(s) in FRONT of the sheet
		// (interior / region-of-interest side, p0 + k*m_normalSign).
		double dbg_Efront = 0.0;
		for (int k = 1; k <= m_cells; ++k)
		{
			int qf = p0 + k * m_normalSign;
			if (qf < 0 || qf >= nLny) continue;
			pos[ny] = (unsigned int)qf;
			for (unsigned int i = 0; i < m_numLines[0]; ++i)
			{
				pos[nyP] = m_posStart[nyP] + i;
				for (unsigned int j = 0; j < m_numLines[1]; ++j)
				{
					pos[nyPP] = m_posStart[nyPP] + j;
					for (int d = 0; d < 3; ++d)
					{
						double v = eng->EngType::GetVolt(d, pos);
						dbg_Efront += v * v;
					}
				}
			}
		}
		unsigned int ts = m_Eng->GetNumberOfTimesteps();
		if (ts % 2000 == 0)
			std::cerr << "[DC_DRAIN] TS=" << ts
			          << " edges=" << dbg_n
			          << " behind|V|^2=" << dbg_E
			          << " front|V|^2=" << dbg_Efront << std::endl;
	}
}

void Engine_Ext_DC_Drain::Apply2Voltages(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, threadID);
}
