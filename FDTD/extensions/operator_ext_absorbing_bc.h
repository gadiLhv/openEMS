/*
*	Copyright (C) 2023-2025 Gadi Lahav (gadi@rfwithcare.com), Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
 * The localized boundary conditions currently implemented are denoted:
 * - First order Mur BC, single phased velocity
 * - First order Mur BC with super absorption, single phase velocity.
 * - Surface Impedance Absorbing Boundary Condition (Leontovich), SIBC.
 * The Mur variants were chosen after experimentation with the various boundary
 * conditions suggested in [1]. The amount of post processing and
 * absorption performance pointed towards the Single PV 1st order Mur.
 * Since this method also requires an extra mesh cell in the doublet grid,
 * the option to use the regular 1st order Mur was left available.
 * The SIBC follows the formulation in [2] -- only the tangential H on the
 * boundary plane is modified; tangential E on the boundary is assumed zero
 * (as it is for the outer face of a PEC-terminated grid).
 *
 * References:
 * [1] Betz, Vaughn Timothy, and R. Mittra. "Absorbing boundary conditions for the finite-difference time-domain analysis of guided-wave structures." Coordinated Science Laboratory Report no. UILU-ENG-93-2243 (1993).‏
 * [2] Y. Mao, A. Z. Elsherbeni, S. Li, T. Jiang, "Surface Impedance Absorbing Boundary for Terminating FDTD Simulations," ACES Journal, vol. 29, no. 12, pp. 1035-1046, 2014.
 */


#ifndef OPERATOR_EXT_ABSORBING_BC_H
#define OPERATOR_EXT_ABSORBING_BC_H

#include "FDTD/operator.h"
#include "operator_extension.h"
#include "tools/arraylib/array_ij.h"

#include "CSPropAbsorbingBC.h"

class Operator_Ext_Absorbing_BC : public Operator_Extension
{
	friend class Engine_Ext_Absorbing_BC;
public:

	// This should be a replica of the CSXCAD property, but can also be something
	// else, later
	enum ABCtype
	{
		UNDEFINED	= 0,
		MUR_1ST 	= 1,	// Mur's BC, 1st order
		MUR_1ST_SA 	= 2,	// Mur's BC, 1st order, with Super Absorption
		SIBC		= 3		// Surface Impedance Absorbing BC (Leontovich)
	};

	Operator_Ext_Absorbing_BC(Operator* op);
	~Operator_Ext_Absorbing_BC();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	//virtual bool IsMPISave() const {return true;}

	virtual std::string GetExtensionName() const
	{
		return std::string("Local absorbing boundary condition sheet");
	}

	virtual void ShowStat(std::ostream &ostr) const;

	//! Initialize all parameters, so the extension can be built later
	virtual bool SetInitParams(CSPrimitives* prim, CSPropAbsorbingBC* abc_prop);

protected:

	Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext);
	void Initialize();

	unsigned int			m_numCells;		// Number of cells in each primitive

	// Storage for the directions. ny is the normal direction to the sheet
	int				m_ny,
					m_nyP,
					m_nyPP;

	// Storage for sheet bounding box start and stop.
	unsigned int	m_sheetX0[3];
	unsigned int	m_sheetX1[3];

	unsigned int 	m_numLines[2];

	bool 			m_normalSignPositive;


	ABCtype			m_ABCtype;

	double			m_phaseVelocity;

	// SIBC: user-specified surface impedance (Ohm). 0 means "derive Z = sqrt(mu/eps)
	// per cell from the local material".
	double			m_surfaceImpedance;

	// Coefficients, to be initialized on-demand.
	// Mur/MurSA usage: K1 = (vp*dt - delta)/(vp*dt + delta), K2 = vp*dt/delta (SA only).
	// SIBC usage:      K1 = (1 - beta)/(1 + beta),
	//                  K2 = (2*dt/(mu*delta_ny)) / (1 + beta),
	//                  K3 = (dt/(mu*delta_tang)) / (1 + beta),  beta = dt*Z/(mu*delta_ny).
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K1_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K1_nyPP;
	ArrayLib::ArrayIJ<FDTD_FLOAT> 	m_K2_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K2_nyPP;
	ArrayLib::ArrayIJ<FDTD_FLOAT> 	m_K3_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K3_nyPP;



};

#endif // OPERATOR_EXT_ABSORBING_BC_H
