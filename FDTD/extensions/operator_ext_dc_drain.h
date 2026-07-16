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
 * Boundary backing damper ("DC removal") operator extension.
 *
 * Motivation
 * ----------
 * A wave/impedance absorbing boundary condition (Mur, SIBC) absorbs the
 * propagating band but reflects DC / near-DC content: a static field has zero
 * phase advance, so none of it passes and the low-frequency energy stalls the
 * run (a TEM line with DC in-band settles to ~-8 dB instead of ~-56 dB).
 *
 * A charge-projection ("divergence cleaning") drain was tried first and shown
 * to be ineffective here: the trapped DC on a TEM line is divergence-free at
 * the boundary (charge sits on the conductors, div D ~ 0 to machine precision),
 * so there is no spurious node charge to drain.
 *
 * Mechanism (this extension)
 * --------------------------
 * Instead of a charge drain, apply a plain dissipative damping to the field in
 * the mesh cell(s) directly BEHIND the absorbing sheet -- the side opposite the
 * sheet normal, i.e. away from the region of interest. Each timestep, after the
 * ABC has updated its boundary voltages, the voltages in those backing cells
 * are scaled by (1 - alpha). This bleeds off the static/near-DC energy that
 * piles up behind the absorber.
 *
 * Because the backing cells lie outside the region of interest (ports, probes
 * and the guide all sit on the interior side of the sheet), the damping does
 * not need to be frequency-selective: everything there can be dissipated
 * without corrupting the measured field. Any reflection the damped region
 * produces travels back into the absorber and is swallowed.
 *
 * Geometry: the sheet plane is p0 = m_posStart[m_ny]; the Mur stencil reads its
 * interior neighbour at p0 + m_normalSign (the region of interest), so the
 * backing cells are at p0 - k*m_normalSign, k = 1 .. m_cells.
 *
 * Wiring
 * ------
 * Created per Mur/SIBC sheet in openEMS::SetupAbsorbingSheets() when the
 * environment variable OPENEMS_DC_DRAIN is set. Tunables (all optional):
 *   OPENEMS_DC_DRAIN_ALPHA  per-step voltage damping fraction (default 0.05)
 *   OPENEMS_DC_DRAIN_CELLS  number of backing cells to damp   (default 1)
 *   OPENEMS_DC_DRAIN_DEBUG  print behind-cell energy diagnostics
 * No CSXCAD property and no Python setup are required; it reuses the geometry of
 * the absorbing_bc sheet it shadows.
 */

#ifndef OPERATOR_EXT_DC_DRAIN_H
#define OPERATOR_EXT_DC_DRAIN_H

#include "FDTD/operator.h"
#include "operator_extension.h"

class Operator_Ext_Absorbing_BC;

class Operator_Ext_DC_Drain : public Operator_Extension
{
	friend class Engine_Ext_DC_Drain;
public:
	//! Build a backing damper shadowing an already-initialized absorbing_bc sheet.
	Operator_Ext_DC_Drain(Operator* op, Operator_Ext_Absorbing_BC* abc, double alpha, int cells);
	~Operator_Ext_DC_Drain();

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual std::string GetExtensionName() const
	{
		return std::string("Boundary backing damper (DC removal)");
	}

	virtual void ShowStat(std::ostream &ostr) const;

	//! Accessor for the engine extension (m_Op is protected in the base).
	Operator* GetOperator() const { return m_Op; }

protected:
	// Sheet geometry copied from the shadowed absorbing_bc sheet.
	int				m_ny;			// normal direction of the sheet
	int				m_nyP;			// first  tangential direction
	int				m_nyPP;			// second tangential direction
	int				m_normalSign;	// +1 / -1 : interior neighbour is p0 + m_normalSign

	unsigned int	m_posStart[3];	// sheet start index (grid)
	unsigned int	m_posStop[3];	// sheet stop  index (grid)

	double			m_alpha;		// per-step voltage damping fraction (0..1]
	int				m_cells;		// number of backing cells to damp
};

#endif // OPERATOR_EXT_DC_DRAIN_H
