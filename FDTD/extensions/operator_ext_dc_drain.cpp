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

#include "operator_ext_dc_drain.h"
#include "operator_ext_absorbing_bc.h"
#include "engine_ext_dc_drain.h"

using std::endl;

Operator_Ext_DC_Drain::Operator_Ext_DC_Drain(Operator* op, Operator_Ext_Absorbing_BC* abc, double alpha, int cells)
	: Operator_Extension(op)
{
	m_alpha = alpha;
	m_cells = cells;

	m_ny   = abc->GetNy();
	// Tangential directions, same convention as the absorbing_bc sheet.
	m_nyP  = (m_ny + 1) % 3;
	m_nyPP = (m_ny + 2) % 3;

	// interior neighbour is at p0 + m_normalSign (matches the Mur stencil shift);
	// the backing cells are on the opposite side, p0 - k*m_normalSign.
	m_normalSign = abc->GetNormalSignPositive() ? +1 : -1;

	const unsigned int* x0 = abc->GetSheetX0();
	const unsigned int* x1 = abc->GetSheetX1();
	for (int n = 0; n < 3; ++n)
	{
		m_posStart[n] = x0[n];
		m_posStop[n]  = x1[n];
	}
}

Operator_Ext_DC_Drain::~Operator_Ext_DC_Drain()
{
}

bool Operator_Ext_DC_Drain::BuildExtension()
{
	// Pure geometric/dissipative extension; nothing to precompute.
	return true;
}

Engine_Extension* Operator_Ext_DC_Drain::CreateEngineExtention()
{
	m_Eng_Ext = new Engine_Ext_DC_Drain(this);
	return m_Eng_Ext;
}

void Operator_Ext_DC_Drain::ShowStat(std::ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
	ostr << " Backing damper alpha: " << m_alpha << " over " << m_cells << " cell(s)" << endl;
}
