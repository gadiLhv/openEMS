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

#ifndef ENGINE_EXT_DC_DRAIN_H
#define ENGINE_EXT_DC_DRAIN_H

#include "engine_extension.h"
#include "engine_extension_dispatcher.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

class Operator_Ext_DC_Drain;

//! Runs after the absorbing_bc sheet in the voltage-apply phase (lower
//! priority) and damps the field in the cell(s) behind the sheet, bleeding off
//! the trapped DC / near-DC energy that the ABC reflects.
class Engine_Ext_DC_Drain : public Engine_Extension
{
public:
	Engine_Ext_DC_Drain(Operator_Ext_DC_Drain* op_ext);
	virtual ~Engine_Ext_DC_Drain();

	// Only the voltage-apply phase is used; it is the only phase that may write
	// engine voltages, and by then the ABC has already written the boundary E.
	virtual void Apply2Voltages() {Engine_Ext_DC_Drain::Apply2Voltages(0);}
	virtual void Apply2Voltages(int threadID);

protected:
	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng, int threadID);

	Operator_Ext_DC_Drain*	m_Op_drain;
	Operator*				m_Op;			// cached operator for geometry / line counts

	int				m_ny, m_nyP, m_nyPP;
	int				m_normalSign;
	unsigned int	m_posStart[3];
	unsigned int	m_posStop[3];
	unsigned int	m_numLines[2];			// transverse cell counts of the sheet
	double			m_alpha;				// per-step voltage damping fraction
	int				m_cells;				// number of backing cells to damp
};

#endif // ENGINE_EXT_DC_DRAIN_H
