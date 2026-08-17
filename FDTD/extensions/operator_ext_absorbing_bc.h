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
 * These were chosen after experimentation with the various boundary
 * conditions suggested in [1]. The amount of post processing and
 * absorption performance pointed towards the Single PV 1st order Mur.
 * Since this method also requires an extra mesh cell in the doublet grid,
 * the option to use the regular 1st order Mur was left available.
 *
 * References:
 * [1] Betz, Vaughn Timothy, and R. Mittra. "Absorbing boundary conditions for the finite-difference time-domain analysis of guided-wave structures." Coordinated Science Laboratory Report no. UILU-ENG-93-2243 (1993).‏
 */


#ifndef OPERATOR_EXT_ABSORBING_BC_H
#define OPERATOR_EXT_ABSORBING_BC_H

#include <vector>

#include "FDTD/operator.h"
#include "operator_extension.h"
#include "tools/arraylib/array_ij.h"

#include "CSPropAbsorbingBC.h"

//! Starting FIR length of the Modal Mur per-cell delay, in taps.
/*!
  Deliberately a compile-time predef and NOT a user parameter. Its only
  effect is how faithfully the truncated filter reproduces exp(-j*beta*dz),
  which is an internal accuracy/cost trade with no physical meaning to a
  caller. Measured on the rectangular-waveguide testbed: 512 taps give
  0.49% band error and max|H| = 1.011, and S11 is essentially flat against
  tap count from 256 upward -- there is nothing here worth exposing.
  */
#define MODAL_MUR_NTAPS 512

//! Hard ceiling for the adaptive tap growth (see MODAL_MUR_MAX_ABS_H).
/*!
  512 taps is not universally enough. The kernel has to span the ONE-CELL
  TRANSIT TIME in timesteps, dz/(v*dt) = 1/nu, and nu is a property of the
  caller's mesh, not of the absorber: a line whose transverse mesh is much
  finer than its longitudinal one drives dt down while dz stays put. The
  PTFE coax test lands at nu = 0.0083, i.e. 120 timesteps per cell, and 512
  taps cannot represent that delay plus its dispersion tail -- measured
  max|H| = 1.066, which is NOT passive and will pump energy.
  */
#define MODAL_MUR_NTAPS_MAX 8192

//! Passivity target the adaptive tap growth aims for.
/*!
  |exp(-j*beta*dz)| <= 1 exactly, so any excess here is pure truncation
  error. Passivity is the one property the whole scheme rests on, so the
  filter is grown until the REALISED response respects it (or the ceiling is
  hit, which is then reported loudly).
  */
#define MODAL_MUR_MAX_ABS_H 1.02

class Operator_Ext_Absorbing_BC : public Operator_Extension
{
	friend class Engine_Ext_Absorbing_BC;
public:

	// This should be a replica of the CSXCAD property, but can also be something
	// else, later
	// Note: integer values must align with CSPropAbsorbingBC's absorbing type enum
	enum ABCtype
	{
		UNDEFINED	= 0,
		MUR_1ST 	= 1,	// Mur's BC, 1st order
		MUR_1ST_SA 	= 2,	// Mur's BC, 1st order, with Super Absorption
		MODAL		= 3,	// Modal absorption (Huygens injection form)
		MODAL_MUR	= 4		// Dispersive Modal Mur (one-way modal termination)
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

	ABCtype GetABCtype() const { return m_ABCtype; }
	int GetNy() const { return m_ny; }
	std::string GetEModeFileName() const { return m_EModeFileName; }
	std::string GetHModeFileName() const { return m_HModeFileName; }
	double GetZw() const { return m_Zw; }

	//! Copy the sheet bounding box (drawing units) into the caller's arrays.
	void GetSheetBoundingBox(double start[3], double stop[3], bool Efield = true) const;

	//! True once the Modal Mur taps and mode template were built successfully.
	bool IsModalMurReady() const { return m_MurReady; }

protected:

	Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext);
	void Initialize();

	//! Build everything the Modal Mur engine extension needs: the mode
	//! template at Yee edge positions, and the per-cell delay taps.
	bool BuildModalMur();

	unsigned int			m_numCells;		// Number of cells in each primitive

	// Storage for the directions. ny is the normal direction to the sheet
	int				m_ny,
					m_nyP,
					m_nyPP;

	// Storage for sheet bounding box start and stop (grid indices).
	unsigned int	m_sheetX0[3];
	unsigned int	m_sheetX1[3];

	unsigned int	m_sheetX0_h[3];
	unsigned int	m_sheetX1_h[3];

	// Physical coordinates of the sheet bounding box (drawing units), stored from SetInitParams.
	double			m_dSheetStart[3];
	double			m_dSheetStop[3];

	// H-Field PMM sheet coordinates
	double			m_hSheetStart[3];
	double			m_hSheetStop[3];

	// Mode file names for modal absorber (empty for non-modal types).
	std::string		m_EModeFileName;
	std::string		m_HModeFileName;

	// Wave impedance for modal absorber (0.0 if not set).
	double			m_Zw;

	unsigned int 	m_numLines[2];

	bool 			m_normalSignPositive;


	ABCtype			m_ABCtype;

	double			m_phaseVelocity;

	// ---- Dispersive Modal Mur (MODAL_MUR) --------------------------------
	// Modal cutoff. May be zero or negative; negative means kc^2 < 0, which is
	// how TEM / quasi-TEM lines are expressed. See CSPropAbsorbingBC.
	double			m_CutOffFrequency;
	bool			m_CutOffFrequencySet;

	// Index of the plane the modal amplitude is READ from: one cell into the
	// guide from the sheet. The sheet plane itself (m_sheetX0[m_ny]) is where
	// the delayed amplitude is written.
	unsigned int	m_MurReadPos;

	// Spacing between the read plane and the sheet plane [m]. This is the dz
	// the delay filter is designed for, so the mesh must be locally uniform
	// here for the lattice dispersion relation to hold.
	double			m_MurDz;

	// Per-cell delay FIR, newest sample first.
	std::vector<double>	m_MurTaps;

	// Mode template sampled at the Yee EDGE positions of each transverse
	// component, over the FULL sheet aperture, L2-normalised so that
	// sum (mP^2 + mPP^2) * dA = 1. That normalisation is what makes the read
	// integral and the deploy inverses of each other.
	//
	// Sampling at Yee edge positions (rather than at the node, as the mode
	// match does) matters: the two transverse components live half a cell
	// apart, and treating them as co-located is exactly the kind of half-cell
	// error that produced an edge-asymmetric absorber before.
	ArrayLib::ArrayIJ<double>	m_MurModeP;
	ArrayLib::ArrayIJ<double>	m_MurModePP;

	bool			m_MurReady;

	// Coefficients, to be initialized on-demand.
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K1_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K1_nyPP;
	ArrayLib::ArrayIJ<FDTD_FLOAT> 	m_K2_nyP;
	ArrayLib::ArrayIJ<FDTD_FLOAT>	m_K2_nyPP;



};

#endif // OPERATOR_EXT_ABSORBING_BC_H
