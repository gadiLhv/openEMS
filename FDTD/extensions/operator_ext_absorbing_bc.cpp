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

#include "operator_ext_absorbing_bc.h"
#include "engine_ext_absorbing_bc.h"
#include "modal_mur_taps.h"

#include "tools/array_ops.h"
#include "tools/global.h"

#include "CSPrimBox.h"
#include "CSModeFileParser.h"

#include <cmath>

using std::cerr;
using std::endl;

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op) : Operator_Extension(op)
{
	Initialize();
}

Operator_Ext_Absorbing_BC::~Operator_Ext_Absorbing_BC()
{
}

Operator_Ext_Absorbing_BC::Operator_Ext_Absorbing_BC(Operator* op, Operator_Ext_Absorbing_BC* op_ext) : Operator_Extension(op, op_ext)
{
	Initialize();
}

Operator_Extension* Operator_Ext_Absorbing_BC::Clone(Operator* op)
{
	return new Operator_Ext_Absorbing_BC(op, this);
}

void Operator_Ext_Absorbing_BC::Initialize()
{
	m_numCells = 0;

	m_ny = m_nyP = m_nyPP = -1;

	for (unsigned int dirIdx = 0 ; dirIdx < 3 ; m_sheetX0[dirIdx++] = 0);
	for (unsigned int dirIdx = 0 ; dirIdx < 3 ; m_sheetX1[dirIdx++] = 0);

	m_normalSignPositive = true;

	m_ABCtype = ABCtype::UNDEFINED;

	m_phaseVelocity = 0.0;
	m_Zw = 0.0;

	m_CutOffFrequency = 0.0;
	m_CutOffFrequencySet = false;
	m_MurReadPos = 0;
	m_MurDz = 0.0;
	m_MurTaps.clear();
	m_MurReady = false;
}

bool Operator_Ext_Absorbing_BC::SetInitParams(CSPrimitives* prim, CSPropAbsorbingBC* abc_prop)
{
	CSPrimBox* cSheet = dynamic_cast<CSPrimBox*>(prim);
	// If this just so happen to not be a sheet, ignore this
	if (!cSheet)
	{
		cerr << "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet validation failed, skipping. "
									<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;
		return false;
	}

	// Check that this is actually a sheet
	// Get box start and stop positions


	// Store physical bounding box coordinates before snapping, for use by SetupModalAbsorbProcessing.
	const double* dStart_phys = cSheet->GetStartCoord()->GetCoords(m_Op->m_MeshType);
	const double* dStop_phys  = cSheet->GetStopCoord()->GetCoords(m_Op->m_MeshType);
	for (int n = 0; n < 3; ++n)
	{
		m_dSheetStart[n] = dStart_phys[n];
		m_dSheetStop[n]  = dStop_phys[n];
	}

	// snap to the native coordinate system
	int Snap_Dimension =
		m_Op->SnapBox2Mesh(
				cSheet->GetStartCoord()->GetCoords(m_Op->m_MeshType),	// Start Coord
				cSheet->GetStopCoord()->GetCoords(m_Op->m_MeshType),	// Stop Coord
				m_sheetX0,	// Start Index
				m_sheetX1,	// Stop Index
				false,		// Dual H-mesh
				true);		// Full mesh?

	// Verify that snapped dimension is correct
	if (Snap_Dimension <= 0)
	{
		if (Snap_Dimension >= -1)
			cerr << "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet snapping failed! Dimension is: " << Snap_Dimension << " skipping. "
					<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;
		// Snap_Dimension == -2 means outside the simulation domain --> no special warning, but box probably marked as unused!
		return false;
	}

	// check that one of the dimensions is 1, and only one
	unsigned int sheetCheck = 0;
	m_ny = -1;
	for (int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
	{
		unsigned int Ncells = m_sheetX1[dimIdx] - m_sheetX0[dimIdx] + 1;

		// Count the number of dimensions that are single-cell length (start == stop)
		sheetCheck += Ncells == 1;

		// Update the normal direction
		m_ny = (Ncells == 1) ? dimIdx : m_ny;
	}

	if (sheetCheck != 1)
	{
		cerr 	<< "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet is not a sheet! Skipping. "
				<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;
		return false;
	}

	m_phaseVelocity = abc_prop->GetPhaseVelocity();
	// Currently not handling arbitrary phase velocity
	if (m_phaseVelocity == 0.0)
	{
		cerr 	<< "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Absorbing sheet Currently does not support per-material velocity. Setting to C0 "
				<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;

		m_phaseVelocity = __C0__;
	}

	// Copy all of the relevant data, so BuildExtension can work
	m_ABCtype = (ABCtype)(abc_prop->GetAbsorbingBoundaryType());
	m_normalSignPositive = abc_prop->GetNormalSignPositive();
	m_EModeFileName = abc_prop->GetEModeFileName();
	m_HModeFileName = abc_prop->GetHModeFileName();
	m_Zw = abc_prop->GetWaveImpedance();
	m_CutOffFrequency = abc_prop->GetCutOffFrequency();
	m_CutOffFrequencySet = abc_prop->IsCutOffFrequencySet();

	// Now the H-Field PMM coordinates
	// * Start by detecting the ePMM index

	for (unsigned int iy = 0 ; iy < 3 ; iy++)
	{
		m_sheetX0_h[iy] = m_sheetX0[iy];
		m_sheetX1_h[iy] = m_sheetX1[iy];
	}

	// Determine correct shift and notify user if it slides out of bounding box.
	// MODAL_MUR is exempt: it uses no H plane at all, and sitting ON the domain
	// edge is its RECOMMENDED placement (the sheet is meant to be the face of
	// the terminating PEC block), so the warning below would be backwards.
	unsigned int hShift = (unsigned int)m_normalSignPositive;
	// Shift one index back if this is catching negative direction propgating waves
	if ((m_ABCtype != ABCtype::MODAL_MUR)
	 && ((m_normalSignPositive && (m_sheetX0_h[m_ny] == 0)) || (!m_normalSignPositive && (m_sheetX0_h[m_ny] == (m_Op->GetNumberOfLines(m_ny) - 1)))))
	{
		cerr 	<< "Operator_Ext_Absorbing_BC::SetInitParams(): Warning: Trying to set local absorber on bonding box edge. Results will be erroneous"
				<< " ID: " << prim->GetID() << " @ Property: " << abc_prop->GetName() << endl;

		hShift = 0;
	}

	// Copy data
	for (unsigned int iy = 0 ; iy < 3 ; iy++)
	{
		m_sheetX0_h[iy] = m_sheetX0[iy];
		m_sheetX1_h[iy] = m_sheetX1[iy];

		m_hSheetStart[iy] = m_dSheetStart[iy];
		m_hSheetStop[iy]  = m_dSheetStop[iy];
	}
	// Update shifted coordinates for h-field
	m_sheetX0_h[m_ny] = m_sheetX1_h[m_ny] = m_sheetX0[m_ny] - hShift;
	// Update for the H-field sheet (PMM). The H mode-match runs on the DUAL mesh,
	// so hand it the dual-node coordinate of the intended plane: a primary-line
	// coordinate would be re-snapped onto the dual grid and can round to the
	// neighboring dual plane (observed: one cell off along the normal).
	m_hSheetStop[m_ny] = m_hSheetStart[m_ny] = m_Op->GetDiscLine(m_ny, m_sheetX1_h[m_ny], true);


	prim->SetPrimitiveUsed(true);

	return true;
}

bool Operator_Ext_Absorbing_BC::BuildExtension()
{

	double dT	= m_Op->GetTimestep();

	unsigned int	pos[] = {0,0,0};
	double			coord[] = {0.0,0.0,0.0};

	double 			delta;

	unsigned int	Ncells[] = {0,0,0},
					totCells = 1;

	// Count the number of cells in each direction
	for (unsigned int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
	{
		Ncells[dimIdx] = m_sheetX1[dimIdx] - m_sheetX0[dimIdx] + 1;
		totCells *= Ncells[dimIdx];
	}

	m_numCells = totCells;

	// Update remaining directions
	m_nyP 	= (m_ny + 1) % 3;
	m_nyPP 	= (m_ny + 2) % 3;

	// The position starts and stops in the same place.
	pos[m_ny] = m_sheetX0[m_ny];

	// The position is considered in the middle of the mesh cell.
	delta = fabs(m_Op->GetEdgeLength(m_ny,pos));

	// Initialize coefficients for this engine extension
	FDTD_FLOAT	vt_nyP  = m_phaseVelocity*dT,
				vt_nyPP = m_phaseVelocity*dT;

	// Initialize number of lines to be used in the arrayIJ
	m_numLines[0] = Ncells[m_nyP];
	m_numLines[1] = Ncells[m_nyPP];

	// The one-way form shares none of the Mur machinery below: no phase
	// velocity, no K1/K2 coefficients, no wave impedance. Build its own state
	// and return.
	if (m_ABCtype == ABCtype::MODAL_MUR)
		return BuildModalMur();

	// Initialize containers. If there are more BCs in the future, this needs to be updated with the respective conditions.
	m_K1_nyP.Init("K1_Coeff_nyP", m_numLines);
	m_K1_nyPP.Init("K1_Coeff_nyPP", m_numLines);

	if (m_ABCtype == ABCtype::MUR_1ST_SA)
	{
		m_K2_nyP.Init("K2_Coeff_nyP", m_numLines);
		m_K2_nyPP.Init("K2_Coeff_nyPP", m_numLines);
	}

	// Prepare containers for per-material assignment
	coord[m_ny] = m_Op->GetDiscLine(m_ny,pos[m_ny]);
	// Check if this is the first or last cell
	if (m_sheetX0[m_ny] == 0)
		coord[m_ny] +=  delta/2 / m_Op->GetGridDelta();
	if (m_sheetX0[m_ny] == (m_Op->GetNumberOfLines(m_ny,true)-1))
		coord[m_ny] += -delta/2 / m_Op->GetGridDelta();

	// Initialize array coefficients
	unsigned int	arrI = 0,
					arrJ;

	for (pos[m_nyP] = m_sheetX0[m_nyP] ; pos[m_nyP] <= m_sheetX1[m_nyP] ; ++pos[m_nyP])
	{
		coord[m_nyP] = m_Op->GetDiscLine(m_nyP,pos[m_nyP]);
		arrJ = 0;
		for (pos[m_nyPP] = m_sheetX0[m_nyPP] ; pos[m_nyPP] <= m_sheetX1[m_nyPP] ; ++pos[m_nyPP])
		{
			coord[m_nyPP] = m_Op->GetDiscLine(m_nyPP,pos[m_nyPP]);

			if (m_phaseVelocity == 0.0)
			{
				double eps,mue;
				CSProperties** prop = m_Op->GetGeometryCSX()->GetPropertiesByCoordsPriority(coord, CSProperties::MATERIAL, false);
				if(prop != NULL)
				{
					cerr << "Operator_Ext_Absorbing_BC::BuildExtension(): Warning: This shouldn't happen...";
					/*CSPropMaterial* mat = (*prop)->ToMaterial();

					eps = mat->GetEpsilonWeighted(int(m_dir[1]),coord);
					mue = mat->GetMueWeighted(int(m_dir[1]),coord);
					vt_nyP = __C0__ * dT / sqrt(eps*mue);

					eps = mat->GetEpsilonWeighted(int(m_dir[2]),coord);
					mue = mat->GetMueWeighted(int(m_dir[2]),coord);
					vt_nyPP = __C0__ * dT / sqrt(eps*mue);*/
				}
				else
				{
					eps = m_Op->GetBackgroundEpsR();
					mue = m_Op->GetBackgroundMueR();

					vt_nyP  = __C0__ * dT / sqrt(eps*mue);
					vt_nyPP = __C0__ * dT / sqrt(eps*mue);
				}

			}

			// If more boundary types of boundary conditions are added in the future,
			// this needs to be in a condition, as well.
			m_K1_nyP(arrI,arrJ)  = (vt_nyP  - delta) / (vt_nyP  + delta);
			m_K1_nyPP(arrI,arrJ) = (vt_nyPP - delta) / (vt_nyPP + delta);

			if (m_ABCtype == ABCtype::MUR_1ST_SA)
			{
				m_K2_nyP(arrI,arrJ)  = vt_nyP /delta;
				m_K2_nyPP(arrI,arrJ) = vt_nyPP/delta;
			}

			arrJ++;
		}
		arrI++;
	}

	return true;
}

bool Operator_Ext_Absorbing_BC::BuildModalMur()
{
	m_MurReady = false;

	if (m_EModeFileName.empty())
	{
		cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: MODAL_MUR needs an E mode file." << endl;
		return false;
	}
	if (!m_CutOffFrequencySet)
	{
		cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: MODAL_MUR needs a cutoff frequency"
		        " (fc, may be negative for TEM/QTEM). Set it via AddModalMurAbsorber(..., fc=<Hz>)." << endl;
		return false;
	}

	// ---- plane layout ------------------------------------------------------
	// The sheet plane is where the delayed amplitude is WRITTEN; in the
	// intended setup it is the face of a PEC block, so the E update re-zeroes
	// it every step and the write is a clean slate rather than an accumulation.
	// The READ plane is one cell into the guide.
	//
	// The shift follows the same convention the Mur path uses for its interior
	// neighbour (m_pos_ny0_shift_V): normalSignPositive means the interior lies
	// at HIGHER index. Direction is encoded entirely by which plane is read and
	// which is written -- there is no sign in the filter to get wrong.
	const unsigned int deployPos = m_sheetX0[m_ny];
	const int readShift = m_normalSignPositive ? +1 : -1;
	const long readPosL = (long)deployPos + readShift;

	if ((readPosL < 0) || (readPosL >= (long)m_Op->GetNumberOfLines(m_ny, true)))
	{
		cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: read plane falls outside the mesh."
		        " The absorber sheet needs at least one cell of guide on its inner side." << endl;
		return false;
	}
	m_MurReadPos = (unsigned int)readPosL;

	m_MurDz = fabs(m_Op->GetDiscLine(m_ny, m_MurReadPos, false)
	             - m_Op->GetDiscLine(m_ny, deployPos,    false)) * m_Op->GetGridDelta();
	if (m_MurDz <= 0.0)
	{
		cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: degenerate read/sheet plane spacing." << endl;
		return false;
	}

	// ---- mode template at the Yee EDGE positions ---------------------------
	CSModeFileParser modeFile;
	if (!modeFile.ParseFile(m_EModeFileName) || !modeFile.HasModeData())
	{
		cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: could not parse E mode file '"
		     << m_EModeFileName << "'." << endl;
		return false;
	}

	m_MurModeP.Init("modal_mur_mode_nyP", m_numLines);
	m_MurModePP.Init("modal_mur_mode_nyPP", m_numLines);

	// Mode-file coordinates are local to the sheet's physical start corner,
	// the same anchor ProcessModeMatch and the mode-file excitation both use.
	const double* origin = m_dSheetStart;

	unsigned int pos[3] = {0,0,0};
	pos[m_ny] = deployPos;

	double norm = 0.0;
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_sheetX0[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_sheetX0[m_nyPP] + j;

			// The two transverse E components do NOT live at the same point:
			// each sits at the centre of its own edge, half a cell apart. Sample
			// each where it actually lives -- the dual line along its own
			// direction, the primary line along the other.
			double vP, vPP, dummy;

			modeFile.LinInterp2(m_Op->GetDiscLine(m_nyP,  pos[m_nyP],  true)  - origin[m_nyP],
			                    m_Op->GetDiscLine(m_nyPP, pos[m_nyPP], false) - origin[m_nyPP],
			                    vP, dummy);

			modeFile.LinInterp2(m_Op->GetDiscLine(m_nyP,  pos[m_nyP],  false) - origin[m_nyP],
			                    m_Op->GetDiscLine(m_nyPP, pos[m_nyPP], true)  - origin[m_nyPP],
			                    dummy, vPP);

			if (std::isnan(vP)  || std::isinf(vP))  vP  = 0.0;
			if (std::isnan(vPP) || std::isinf(vPP)) vPP = 0.0;

			// Each component is half-shifted along ITS OWN direction, so it has
			// one FEWER valid Yee edge there: over N primary lines there are only
			// N-1 dual segments. Blanking that last line is not cosmetic -- the
			// mode file happily interpolates a non-zero value at the clamped
			// coordinate, which would both add a phantom row to the projection
			// and deposit charge on an edge lying on the aperture boundary. The
			// engine skips zero entries, so this removes them from both.
			if (i + 1 == m_numLines[0]) vP  = 0.0;
			if (j + 1 == m_numLines[1]) vPP = 0.0;

			m_MurModeP(i,j)  = vP;
			m_MurModePP(i,j) = vPP;

			double dA = m_Op->GetNodeArea(m_ny, pos, false);
			norm += (vP*vP + vPP*vPP) * dA;
		}
	}

	if (norm <= 0.0)
	{
		cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: mode template is identically zero over the sheet."
		        " Check the mode file's coordinate frame against the sheet's start corner." << endl;
		return false;
	}

	// L2-normalise so that sum (mP^2 + mPP^2) * dA = 1. This is what makes the
	// read integral a = int(E . m)dA and the deploy E = a*m exact inverses.
	norm = sqrt(norm);
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			m_MurModeP(i,j)  /= norm;
			m_MurModePP(i,j) /= norm;
		}

	// ---- per-cell delay taps -----------------------------------------------
	//  The kernel must span the one-cell transit time in timesteps,
	//  1/nu = dz/(v*dt), which is set by the CALLER's mesh: a line whose
	//  transverse mesh is much finer than its longitudinal one pushes dt down
	//  while dz stays put, and the delay stretches to hundreds of samples.
	//  So grow the filter until the REALISED response is passive, rather than
	//  trusting one fixed length. The length stays an internal matter -- it is
	//  never a user parameter -- and the growth is reported so a pathological
	//  mesh is visible rather than silent.
	//
	//  m_phaseVelocity is the medium's wave speed (C0 unless the caller set
	//  it). A dielectric-filled line MUST set it: c0/sqrt(eps_r) for a coax.
	double maxAbsH = 0.0, acausalFrac = 0.0;
	unsigned int nTaps = MODAL_MUR_NTAPS;
	const double nu = m_phaseVelocity * m_Op->GetTimestep() / m_MurDz;

	while (true)
	{
		if (!ModalMur::DesignDelayTaps(m_CutOffFrequency, m_phaseVelocity, m_MurDz,
		                               m_Op->GetTimestep(),
		                               nTaps, m_MurTaps, maxAbsH, acausalFrac))
		{
			cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Error: delay tap design failed." << endl;
			return false;
		}

		if (maxAbsH <= MODAL_MUR_MAX_ABS_H)
			break;

		if (nTaps >= MODAL_MUR_NTAPS_MAX)
		{
			// Out of room. Say exactly what is wrong and what to change: the
			// mesh is the lever here, not anything the absorber can fix.
			cerr << "Operator_Ext_Absorbing_BC::BuildModalMur(): Warning: realised max|H| = "
			     << maxAbsH << " > " << MODAL_MUR_MAX_ABS_H << " even at " << nTaps
			     << " taps, so the delay filter is NOT passive and this absorber may pump energy."
			     << " The mesh is the cause: nu = v*dt/dz = " << nu << " means "
			     << (1.0 / nu) << " timesteps of transit per cell. Coarsen the transverse"
			     << " mesh (which raises dt) or refine the mesh along the propagation"
			     << " direction near the sheet." << endl;
			break;
		}
		nTaps *= 2;
	}

	if (g_settings.GetVerboseLevel() > 0)
		cerr << "Operator_Ext_Absorbing_BC: Modal Mur, ny=" << m_ny
		     << " sheet=" << deployPos << " read=" << m_MurReadPos
		     << " dz=" << m_MurDz << " m, fc=" << m_CutOffFrequency
		     << " Hz, v=" << m_phaseVelocity << " m/s, nu=" << nu
		     << ", " << nTaps << " taps, max|H|=" << maxAbsH
		     << ", acausal=" << acausalFrac << endl;

	m_MurReady = true;
	return true;
}

Engine_Extension* Operator_Ext_Absorbing_BC::CreateEngineExtention()
{
	// Assigning to the base m_Eng_Ext is what makes GetEngineExtention() return
	// this instance later (matches the Operator_Ext_SteadyState pattern).
	// Without this, SetupModalAbsorbProcessing can't reach the engine extension.
	m_Eng_Ext = new Engine_Ext_Absorbing_BC(this);
	return m_Eng_Ext;
}

void Operator_Ext_Absorbing_BC::ShowStat(std::ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);

	ostr << " Total cells: " << m_numCells << endl;
}

void Operator_Ext_Absorbing_BC::GetSheetBoundingBox(double start[3], double stop[3], bool Efield) const
{
	const double* startArr = Efield ? m_dSheetStart : m_hSheetStart;
	const double* stopArr  = Efield ? m_dSheetStop  : m_hSheetStop ;

	for (int n = 0; n < 3; ++n)
	{
		start[n] = startArr[n];
		stop[n] = stopArr[n];
	}
}








