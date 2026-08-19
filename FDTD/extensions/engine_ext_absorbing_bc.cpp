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


/* This version of absorbing boundary conditions is based on this article:
 *
 * Betz, Vaughn Timothy, and R. Mittra. "Absorbing boundary conditions for the finite-difference time-domain analysis of guided-wave structures." Coordinated Science Laboratory Report no. UILU-ENG-93-2243 (1993).
 *
 * After some trial and error, it was discovered that the simplest and most efficient implementations
 * are:
 * 1. Mur first order boundary conditions
 * 2. First order Mur with "super-absorption".
 * Later I discovered that the latter is equivalent to the so-called "Surface impedance boundary
 * conditions" (SIBC).
 */


#include "engine_ext_absorbing_bc.h"
#include "operator_ext_absorbing_bc.h"
#include "FDTD/engine.h"
#include "FDTD/excitation.h"
#include "FDTD/extensions/operator_ext_excitation.h"
#include "FDTD/engine_sse.h"
#include "FDTD/engine_interface_fdtd.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
#include "tools/global.h"
#include "operator_ext_excitation.h"

Engine_Ext_Absorbing_BC::Engine_Ext_Absorbing_BC(Operator_Ext_Absorbing_BC* op_ext) :
	Engine_Extension(op_ext),
	m_K1_nyP(op_ext->m_K1_nyP),
	m_K1_nyPP(op_ext->m_K1_nyPP),
	m_K2_nyP(op_ext->m_K2_nyP),
	m_K2_nyPP(op_ext->m_K2_nyPP)
{

	m_Op_ABC = op_ext;
	m_ABCtype = int(m_Op_ABC->m_ABCtype);

	m_Eng_Interface = NULL;

	for (unsigned int dimIdx = 0 ; dimIdx < 3 ; dimIdx++)
	{
		m_posStart[dimIdx] = m_Op_ABC->m_sheetX0[dimIdx];
		m_posStop[dimIdx] = m_Op_ABC->m_sheetX1[dimIdx];
	}

	m_ny	= m_Op_ABC->m_ny;
	m_nyP 	= m_Op_ABC->m_nyP;
	m_nyPP 	= m_Op_ABC->m_nyPP;

	m_numLines[0] = m_Op_ABC->m_numLines[0];
	m_numLines[1] = m_Op_ABC->m_numLines[1];

	bool normalSignPositive = m_Op_ABC->m_normalSignPositive;

	m_Zw = m_Op_ABC->GetZw();
	m_normalSign = normalSignPositive ? +1 : -1;
	m_Hmm_prev = 0.0;
	m_corr_gain = -1.0;   // computed lazily on first use (needs the PMM grids)

	// One-way delay history. Starts at zero, which is the physically correct
	// initial condition: nothing has left through the sheet yet.
	m_MurHead = 0;
	if (m_ABCtype == int(Operator_Ext_Absorbing_BC::MODAL_MUR))
		m_MurHist.assign(m_Op_ABC->m_MurTaps.size(), 0.0);

	m_start_TS = 0;

	// On-the-fly modal correction: nothing learned, nothing seen, not frozen.
	m_otfcCount = 0;
	m_otfcSnsMax = 0.0;
	m_otfcFrozen = false;
	m_otfcHavePrev = false;

	// AMPLITUDE GATE REFERENCE -- the peak of the excitation waveform, exactly
	// as the prototype's maxEsrc. Do not learn until the sensed modal amplitude
	// is a set fraction of it.
	//
	// A running-maximum gate cannot do this job: on a rising edge the current
	// sample IS the running maximum, so the gate is trivially satisfied by the
	// numerical precursor that arrives long before the wave. Measured on the
	// coax with a running-max gate: sheet 1 learned at timestep 220 when the
	// wave needs 310 to reach it, sheet 2 at 2066 when it needs 3758, and the
	// termination went from -13.7 dB to -9.8 dB. The reference has to be the
	// EVENTUAL peak, which is what the excitation waveform supplies.
	m_otfcSrcPeak = 0.0;
	m_otfcSensePos = (unsigned int)-1;		// -1 = no usable sense plane, OTFC idle
	if (m_Op_ABC->m_Op != NULL)
	{
		Excitation* exc = m_Op_ABC->m_Op->GetExcitationSignal();
		if ((exc != NULL) && (exc->GetVoltageSignal() != NULL))
		{
			const FDTD_FLOAT* sig = exc->GetVoltageSignal();
			for (unsigned int n = 0; n < exc->GetLength(); ++n)
				if (fabs(sig[n]) > m_otfcSrcPeak)
					m_otfcSrcPeak = fabs(sig[n]);
		}

		// SENSE PLANE: ksrc + MODAL_OTFC_SENSE_OFFSET, shared by every sheet,
		// exactly as the prototype's single ksns serves both absorbers. Keyed
		// to the SOURCE and not to this sheet -- sampling next to the sheet
		// reads back the sheet's own write and iterates on it.
		Operator_Ext_Excitation* excExt = m_Op_ABC->m_Op->GetExcitationExtension();
		unsigned int srcPos = 0;
		if ((excExt != NULL) && excExt->GetExcitationPlane(m_ny, srcPos))
		{
			unsigned int cand = srcPos + MODAL_OTFC_SENSE_OFFSET;
			if (cand < m_Op_ABC->m_Op->GetNumberOfLines(m_ny))
				m_otfcSensePos = cand;
		}
		if ((m_otfcSensePos == (unsigned int)-1) && (g_settings.GetVerboseLevel() > 0))
			std::cerr << "Engine_Ext_Absorbing_BC: no single-plane excitation found along ny="
			          << m_ny << "; on-the-fly modal correction disabled for this sheet."
			          << std::endl;
	}

	// TIME GATE -- the guard the prototype gets from its look-ahead amplitude
	// threshold, and the one thing a running-maximum gate cannot provide. On a
	// rising edge the current sample IS the running maximum, so an
	// amplitude-relative gate is trivially satisfied by the numerical precursor
	// that arrives long before the wave does. Measured on the coax without this
	// gate: sheet 1 learned at timestep 220 when the wave needs 310 to reach it,
	// and sheet 2 at 2066 when it needs 3758. Both learned on dust, and the
	// termination went from -13.7 dB to -9.8 dB.
	//
	// The excitation's own peak is the natural reference: by then the source has
	// fully developed and, for any sane geometry, the wave is established at the
	// sheet. Same source of truth Engine_Ext_Mur_ABC uses for its start delay.
	m_otfcStartTS = 0;
	if (m_Op_ABC->m_Op != NULL)
	{
		Excitation* exc = m_Op_ABC->m_Op->GetExcitationSignal();
		if (exc != NULL)
			m_otfcStartTS = exc->GetMaxExcitationTimestep();
	}

	// Initialize shifted position for V
	m_pos_ny0_shift_V = m_posStart[m_ny] + (normalSignPositive  ? 1 : -1);

	// Initialize shifted position for I. Different for super-absorption
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MODAL)
	{
		// In this case, the H-field is absorbed
		m_pos_ny0_I = m_posStart[m_ny] + (normalSignPositive ? 0 : -1);
		m_pos_ny0_shift_I = m_posStart[m_ny] + (normalSignPositive ? 1 : -2);
	}
	else
	{
		// In this case, the H-field is absorbed
		m_pos_ny0_I = m_Op_ABC->m_sheetX0_h[m_ny];
		m_pos_ny0_shift_I = 0;
	}

	m_V_nyP.Init("volt_nyP",m_numLines);
	m_V_nyPP.Init("volt_nyPP",m_numLines);
	m_I_nyP.Init("curr_nyP",m_numLines);
	m_I_nyPP.Init("curr_nyPP",m_numLines);

	// One thread per boundary
	SetNumberOfThreads(1);
}

Engine_Ext_Absorbing_BC::~Engine_Ext_Absorbing_BC()
{
}

void Engine_Ext_Absorbing_BC::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	// This command assigns the number of jobs (primitives) handled by each thread
	m_linesPerThread = AssignJobs2Threads(m_numLines[0],m_NrThreads,false);

	// Basically cumulative sum. Starting point of each thread.
	m_threadStartLine.resize(m_NrThreads,0);
	m_threadStartLine.at(0) = 0;
	for (size_t n = 1; n < m_linesPerThread.size(); ++n)
		m_threadStartLine.at(n) = m_threadStartLine.at(n - 1) + m_linesPerThread.at(n - 1);
}

// The first order Mur boundary condition is based on the following step:
// The field at time step n, location i, is u
// E(i,n + 1) = E(i + s,n) + K1*[E(i,n + 1) - E(i,n)]
// where s is the +-1 shift, depending on the direction, and
// K1 = (vp*Dt - Dx)/(vp*Dt + Dx)

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPreVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	// Modal absorber: subtract the outgoing modal component from V and I at
	// the absorber plane. Runs on thread 0 only (no transverse parallelism for
	// the modal path). The Mur stencil below is skipped for MODAL.
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) == Operator_Ext_Absorbing_BC::MODAL)
	{
		if (threadID != 0) return;
		if (m_Eng_Interface == NULL) return;
		if (m_Zw <= 0.0) return;

		// Mode-match scalars, computed by the PMMs in PA at end of previous iter.
		double Emm      = m_Eng_Interface->GetLastModeMatchE();
		double Hmm_curr = m_Eng_Interface->GetLastModeMatchH();

		// Time-centre H against E (H sample lags by half a step).
		double Hmm  = 0.5 * (m_Hmm_prev + Hmm_curr);
		m_Hmm_prev  = Hmm_curr;

		// Modal amplitude of the wave leaving the domain through this absorber.
		// The absorbed propagation direction is p = -normalSign. Template sign
		// convention: a +ny-travelling wave measures Zw*Hmm = +Emm (c_h = +1);
		// verified against the working waveguide-port S-parameter math, which
		// reconstructs the incident wave as 0.5*(uf + Zref*if) from the same
		// mode-match probes. The outgoing-wave selector is therefore
		//   a = 0.5*(Emm + s*Zw*Hmm),  s = p*c_h = -normalSign.
		// NOTE: this sign and dH_factor below form a pair -- always flip together.
		double a = 0.5 * (Emm - m_normalSign * m_Zw * Hmm);

		// Diagnostic aid: measure-only mode. The PMMs keep recording (dumps stay
		// valid) but the absorber applies no correction, so the undisturbed field
		// can be compared against the mode templates. Enable by setting the
		// environment variable OPENEMS_ABC_MEASURE_ONLY.
		static const bool measure_only = (getenv("OPENEMS_ABC_MEASURE_ONLY") != NULL);
		if (measure_only)
			return;

		const Operator* op = m_Op_ABC->m_Op;

		// ------------------------------------------------------------------
		// DEPLOYMENT GRID. The corrections are painted on the operator's OWN
		// full-aperture sheet grid with templates sampled at the true Yee
		// component positions -- NOT on the mode-match (PMM) grids.
		//
		// The PMMs only MEASURE. Their grids deliberately drop the
		// domain-boundary lines, which is correct for an integral and fatal
		// for actuation: the node clipping is edge-asymmetric, so the
		// wall-adjacent E edges on one side were never painted while the
		// opposite side was. That asymmetry is a spurious source at the wall,
		// re-injected every timestep -- visible as hot spots hugging the
		// conductor at the sheet, and as a near-cutoff residue that stalls the
		// energy decay.
		//
		// Only the NORMAL PMM indices are still read below, and only to
		// calibrate the correction gain from the E/H plane separation.
		unsigned int startE[3], startH[3];
		for (int n = 0; n < 3; ++n)
		{
			startE[n] = m_Eng_Interface->GetModeMatchE_Start(n);
			startH[n] = m_Eng_Interface->GetModeMatchH_Start(n);
		}

		// Correction gain (discrete matching). Subtracting the FULL measured
		// amplitude every timestep over-corrects whenever the wave crosses less
		// than one cell per step: at nu = v*dt/dz ~ 0.1 the full-gain scheme
		// reflects ~80% of the field (measured, 1D and 3D). The matched per-step
		// gain for this staggered subtract-at-plane scheme is
		//   kappa = 2*nu/(1+nu),   nu = v*dt/dz_local
		// with dz_local = 2*|z_Eplane - z_Hplane| (the E/H planes are half a cell
		// apart) and v = Zw/mu0 (TEM in a non-magnetic medium). Verified in 1D:
		// R < 0.01 for nu in [0.1, 0.5], both absorber orientations.
		if (m_corr_gain < 0.0)
		{
			double zE = op->GetDiscLine(m_ny, startE[m_ny], false);
			double zH = op->GetDiscLine(m_ny, startH[m_ny], true);
			double dz = 2.0 * fabs(zE - zH) * op->GetGridDelta();
			double v  = m_Zw / (4e-7 * M_PI);
			double nu_loc = (dz > 0.0) ? (v * op->GetTimestep() / dz) : 1.0;
			m_corr_gain = 2.0 * nu_loc / (1.0 + nu_loc);
			if (g_settings.GetVerboseLevel() > 0)
				std::cerr << "Engine_Ext_Absorbing_BC: modal absorber local nu=" << nu_loc
				          << " -> correction gain " << m_corr_gain << std::endl;
		}
		a *= m_corr_gain;

		if (!m_Op_ABC->m_deployTemplatesValid)
			return;

		const unsigned int P  = m_numLines[0];
		const unsigned int PP = m_numLines[1];

		// V correction at the E plane. The templates are L2-normalized over the
		// sheet (sum m^2*dA = 1), so 'a' is the field-amplitude coefficient and
		// the per-edge voltage correction is  V -= a * m * EdgeLength (E = V/dl).
		unsigned int pos_v[] = {0,0,0};
		pos_v[m_ny] = m_posStart[m_ny];
		for (unsigned int i = 0; i < P; ++i)
		{
			pos_v[m_nyP] = m_posStart[m_nyP] + i;
			for (unsigned int j = 0; j < PP; ++j)
			{
				pos_v[m_nyPP] = m_posStart[m_nyPP] + j;
				// An E edge along nyP exists for i < P-1; along nyPP for j < PP-1.
				// Conductor cells and wall-normal components need no guard: true
				// -position sampling of the mode file gives them exactly 0.
				if (i < P - 1)
				{
					double mE = m_Op_ABC->m_dE_nyP(i, j);
					if (mE != 0.0)
						eng->EngType::SetVolt(m_nyP,  pos_v, eng->EngType::GetVolt(m_nyP,  pos_v) - a * mE * op->GetEdgeLength(m_nyP,  pos_v, false));
				}
				if (j < PP - 1)
				{
					double mE = m_Op_ABC->m_dE_nyPP(i, j);
					if (mE != 0.0)
						eng->EngType::SetVolt(m_nyPP, pos_v, eng->EngType::GetVolt(m_nyPP, pos_v) - a * mE * op->GetEdgeLength(m_nyPP, pos_v, false));
				}
			}
		}

		// I correction on the single dual plane this branch samples H on. The
		// selector above and this launcher sign are INDEPENDENT degrees of
		// freedom: the selector decides which wave is measured, this sign
		// decides which direction the canceling (E,H) pair radiates. Direction
		// verified by a one-shot kick experiment (arrival-time analysis at
		// neighboring probes): -normalSign launches OUTWARD through the
		// absorber for both orientations, as required.
		double dH_factor = -m_normalSign * a / m_Zw;
		unsigned int pos_i[] = {0,0,0};
		pos_i[m_ny] = m_Op_ABC->m_sheetX0_h[m_ny];	// operator's dual plane, not the PMM's
		for (unsigned int i = 0; i < P; ++i)
		{
			pos_i[m_nyP] = m_posStart[m_nyP] + i;
			for (unsigned int j = 0; j < PP; ++j)
			{
				pos_i[m_nyPP] = m_posStart[m_nyPP] + j;
				// I along nyP lives in the dual nyPP slot (j < PP-1); along
				// nyPP in the dual nyP slot (i < P-1). Mirror of the E loop.
				if (j < PP - 1)
				{
					double mH = m_Op_ABC->m_dH_nyP(i, j);
					if (mH != 0.0)
						eng->EngType::SetCurr(m_nyP,  pos_i, eng->EngType::GetCurr(m_nyP,  pos_i) - dH_factor * mH * op->GetEdgeLength(m_nyP,  pos_i, true));
				}
				if (i < P - 1)
				{
					double mH = m_Op_ABC->m_dH_nyPP(i, j);
					if (mH != 0.0)
						eng->EngType::SetCurr(m_nyPP, pos_i, eng->EngType::GetCurr(m_nyPP, pos_i) - dH_factor * mH * op->GetEdgeLength(m_nyPP, pos_i, true));
				}
			}
		}
		return;
	}

	// Bail out for any other non-Mur ABC type (defensive).
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST
	 && (Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA)
		return;

	unsigned int pos[] = {0,0,0};
	unsigned int pos_shift[] = {0,0,0};

	pos[m_ny] = m_posStart[m_ny];
	pos_shift[m_ny] = m_pos_ny0_shift_V;
	for (unsigned int i = m_threadStartLine.at(threadID) ; i < (m_threadStartLine.at(threadID) + m_linesPerThread.at(threadID)) ; i++)
	{
		// Store shifted location in this container
		pos_shift[m_nyP] = pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos_shift[m_nyPP] = pos[m_nyPP] = m_posStart[m_nyPP] + j;

			// E(i + s,n) - K1*E(i,n)
			m_V_nyP (i,j) = eng->EngType::GetVolt(m_nyP ,pos_shift) - m_Op_ABC->m_K1_nyP (i,j) * eng->EngType::GetVolt(m_nyP ,pos);
			m_V_nyPP(i,j) = eng->EngType::GetVolt(m_nyPP,pos_shift) - m_Op_ABC->m_K1_nyPP(i,j) * eng->EngType::GetVolt(m_nyPP,pos);
		}

	}

}

void Engine_Ext_Absorbing_BC::DoPreVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	// Mur-only step; modal absorber does not participate here.
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST
	 && (Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA)
		return;

	unsigned int pos_shift[] = {0,0,0};

	pos_shift[m_ny] = m_pos_ny0_shift_V;
	for (unsigned int i = m_threadStartLine.at(threadID) ; i < (m_threadStartLine.at(threadID) + m_linesPerThread.at(threadID)) ; i++)
	{
		// Store shifted location in this container
		pos_shift[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos_shift[m_nyPP] = m_posStart[m_nyPP] + j;

			// E(i + s,n) - K1*E(i,n) + K1*E(i + s,n) =
			// E(i + s,n) + [K1*E(i + s,n) - K1*E(i,n)]
			m_V_nyP (i,j) += m_K1_nyP (i,j) * eng->EngType::GetVolt(m_nyP ,pos_shift);
			m_V_nyPP(i,j) += m_K1_nyPP(i,j) * eng->EngType::GetVolt(m_nyPP,pos_shift);
		}

	}

}

void Engine_Ext_Absorbing_BC::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::Apply2VoltagesImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	// The one-way modal termination belongs exactly here: after the E update
	// and before the H update, so the value written on the sheet is what the
	// following H update sees -- the same ordering the Octave prototype uses.
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) == Operator_Ext_Absorbing_BC::MODAL_MUR)
	{
		if (threadID == 0)
			ApplyModalMur(eng);
		return;
	}

	// Mur-only step; modal absorber does not participate here.
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST
	 && (Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA)
		return;

	unsigned int pos[] = {0,0,0};

	pos[m_ny] = m_posStart[m_ny];
	for (unsigned int i = m_threadStartLine.at(threadID) ; i < (m_threadStartLine.at(threadID) + m_linesPerThread.at(threadID)) ; i++)
	{
		// Store shifted location in this container
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; j++)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;

			// E(i,n + 1) = E(i + s,n) + [K1*E(i + s,n) - K1*E(i,n)]
			eng->EngType::SetVolt(m_nyP ,pos, m_V_nyP (i,j));
			eng->EngType::SetVolt(m_nyPP,pos, m_V_nyPP(i,j));
		}

	}
}

void Engine_Ext_Absorbing_BC::Apply2Voltages(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2VoltagesImpl, threadID);
}

// =========================================================================
//  DISPERSIVE MODAL MUR  --  one-way modal termination
//
//  Impose, do not inject. The injection form has to MEASURE the outgoing
//  wave, which needs a direction test E = +-Zw*H, and Zw is singular at
//  cutoff and reactive below it -- so that test cannot work there, which is
//  what forced the caps, walls and guards the injection path carries.
//
//  A one-way condition needs no impedance at all. For a mode leaving through
//  the sheet plane, the modal amplitude there is just the DELAYED amplitude
//  one plane inside:
//
//      a_sheet(w) = a_inside(w) * exp(-j*beta(w)*dz)
//
//  so read the modal amplitude one cell in, run it through the delay filter,
//  and write the result onto the sheet. Three steps, no H field anywhere.
// =========================================================================
template <typename EngType, typename T>
bool Engine_Ext_Absorbing_BC::RelearnModalTemplate(EngType* eng,
                                                   ArrayLib::ArrayIJ<T>& modeP,
                                                   ArrayLib::ArrayIJ<T>& modePP,
                                                   unsigned int readPos)
{
	// Per-sheet opt-in, set by the caller through
	// CSPropAbsorbingBC::SetOTFC / AddModalAbsorber(otfc=True).
	if (!m_Op_ABC->m_OTFC)
		return false;
	if (m_otfcFrozen || !m_Op_ABC->m_otfcScratchValid)
		return false;
	if (m_Eng->GetNumberOfTimesteps() < m_otfcStartTS)
		return false;		// see m_otfcStartTS: do not learn on the precursor

	const Operator* op = m_Op_ABC->m_Op;
	unsigned int pos[3] = {0,0,0};
	pos[m_ny] = readPos;

	// LEARN ON THE EXISTING SUPPORT, not on the full aperture. Every cell where
	// the current template is exactly zero is zero for a reason the sensed
	// field cannot tell us: a conductor-owned edge (masked by GetVV == 0), or
	// the trailing Yee line a transverse component does not have along its own
	// direction. Re-deriving that support from a noisy field would resurrect
	// exactly the phantom edge row that the full-aperture fix removed. The
	// support is also demonstrably right where it matters -- on the coax the
	// template carries 0.00% of its energy anywhere the field is dead.
	//
	// Pass 1: the sensed modal amplitude, and the candidate's norm. The
	// candidate is (V/dl)/aSns, so aSns is the divisor the amplitude gate
	// below is protecting, and dividing by it is what keeps the template's
	// polarity locked to the field's instead of flipping with the pulse.
	double aSns = 0.0;
	double eNorm2 = 0.0;		// int E.E dA over the same support, for purity
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;
			double dA = op->GetNodeArea(m_ny, pos, false);
			if (modeP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyP, pos, false);
				if (dl != 0.0)
				{
					double e = eng->EngType::GetVolt(m_nyP, pos) / dl;
					aSns   += e * modeP(i,j) * dA;
					eNorm2 += e * e * dA;
				}
			}
			if (modePP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyPP, pos, false);
				if (dl != 0.0)
				{
					double e = eng->EngType::GetVolt(m_nyPP, pos) / dl;
					aSns   += e * modePP(i,j) * dA;
					eNorm2 += e * e * dA;
				}
			}
		}
	}
	if (aSns == 0.0)
		return false;
	// Purity of the field against the template currently in use. The template
	// is unit-norm, so this is (int E.m dA)^2 / (int E.E dA), the same quantity
	// ProcessModeMatch reports -- and the direct measure of whether a re-learn
	// is actually buying anything.
	const double purityBefore = (eNorm2 > 0.0) ? (aSns * aSns) / eNorm2 : 0.0;

	// AMPLITUDE GATE. The Octave prototype compares against a look-ahead peak
	// of the excitation waveform (0.005*maxEsrc); there is no unit relation in
	// openEMS between the excitation amplitude and a modal amplitude, so the
	// scale has to come from the signal itself. The running maximum is that
	// scale, and it makes the guarantee the prototype actually wanted: never
	// divide by a value far below the largest one seen.
	double absA = fabs(aSns);
	if (absA > m_otfcSnsMax)
		m_otfcSnsMax = absA;
	if ((m_otfcSrcPeak <= 0.0) || (absA < MODAL_OTFC_GATE_FRAC * m_otfcSrcPeak))
		return false;

	// Pass 2: build the normalised candidate, and measure how far it moved
	// from the previous one. Both are unit-norm, so the difference norm is a
	// direct relative shape change.
	double nrm2 = 0.0;
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;
			double dA = op->GetNodeArea(m_ny, pos, false);
			double cP = 0.0, cPP = 0.0;
			if (modeP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyP, pos, false);
				if (dl != 0.0) cP = (eng->EngType::GetVolt(m_nyP, pos) / dl) / aSns;
			}
			if (modePP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyPP, pos, false);
				if (dl != 0.0) cPP = (eng->EngType::GetVolt(m_nyPP, pos) / dl) / aSns;
			}
			nrm2 += (cP*cP + cPP*cPP) * dA;
		}
	}
	if (nrm2 <= 0.0)
		return false;
	const double invN = 1.0 / sqrt(nrm2);

	// Pass 3: normalise, and compare against the PREVIOUS candidate held in the
	// scratch -- not against the current template. Comparing against the
	// template would just re-measure the impurity we are here to remove (0.23
	// on the coax), which no sane tolerance would ever accept. What tells dust
	// from the mode is whether the candidate has stopped MOVING.
	double diff2 = 0.0;
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;
			double dA = op->GetNodeArea(m_ny, pos, false);
			double cP = 0.0, cPP = 0.0;
			if (modeP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyP, pos, false);
				if (dl != 0.0) cP = ((eng->EngType::GetVolt(m_nyP, pos) / dl) / aSns) * invN;
			}
			if (modePP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyPP, pos, false);
				if (dl != 0.0) cPP = ((eng->EngType::GetVolt(m_nyPP, pos) / dl) / aSns) * invN;
			}
			double dP  = cP  - m_Op_ABC->m_otfcCandP(i,j);
			double dPP = cPP - m_Op_ABC->m_otfcCandPP(i,j);
			diff2 += (dP*dP + dPP*dPP) * dA;
			m_Op_ABC->m_otfcCandP(i,j)  = cP;
			m_Op_ABC->m_otfcCandPP(i,j) = cPP;
		}
	}
	const double shapeChange = sqrt(diff2);

	// SHAPE GATE. Before the pulse body arrives, the sense plane holds only
	// numerical dust and its normalised shape jumps around; the prototype's
	// MSL variant guards this with a time gate and the comment that learning
	// on dust "detonates the run". Comparing successive candidates is the same
	// guard without needing to know when the pulse arrives, and it is backed by
	// measurement: the true transverse shape is constant to six decimals over
	// 48 cells, so a candidate that has stopped moving IS the mode.
	if (!m_otfcHavePrev)
	{
		m_otfcHavePrev = true;	// first candidate: nothing to compare against yet
		return false;
	}
	if (shapeChange >= MODAL_OTFC_SHAPE_TOL)
		return false;

	// Adopt: straight replacement, exactly as the prototype.
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			if (modeP(i,j)  != 0.0) modeP(i,j)  = (T)m_Op_ABC->m_otfcCandP(i,j);
			if (modePP(i,j) != 0.0) modePP(i,j) = (T)m_Op_ABC->m_otfcCandPP(i,j);
		}

	// Post-adopt audit: the new template must still be unit-norm, and the
	// field's projection onto it must not have collapsed -- that projection is
	// the whole signal the absorber acts on.
	double newNorm2 = 0.0, aNew = 0.0;
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;
			double dA = op->GetNodeArea(m_ny, pos, false);
			newNorm2 += ((double)modeP(i,j)*(double)modeP(i,j)
			           + (double)modePP(i,j)*(double)modePP(i,j)) * dA;
			if (modeP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyP, pos, false);
				if (dl != 0.0) aNew += (eng->EngType::GetVolt(m_nyP, pos)/dl) * modeP(i,j) * dA;
			}
			if (modePP(i,j) != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyPP, pos, false);
				if (dl != 0.0) aNew += (eng->EngType::GetVolt(m_nyPP, pos)/dl) * modePP(i,j) * dA;
			}
		}
	}

	++m_otfcCount;
	if (m_otfcCount >= MODAL_OTFC_MAX_UPDATES)
		m_otfcFrozen = true;

	if (g_settings.GetVerboseLevel() > 0)
		std::cerr << "Engine_Ext_Absorbing_BC: OTFC update " << m_otfcCount
		          << " on plane " << readPos << " @TS " << m_Eng->GetNumberOfTimesteps()
		          << ", |a|=" << absA << " (gate " << MODAL_OTFC_GATE_FRAC * m_otfcSrcPeak
		          << ", srcPeak " << m_otfcSrcPeak << ", |a|max so far " << m_otfcSnsMax << ")"
		          << ", purity before " << purityBefore
		          << ", newNorm " << sqrt(newNorm2) << ", a_old " << aSns << " -> a_new " << aNew
		          << ", shape change " << shapeChange
		          << (m_otfcFrozen ? " (frozen)" : "") << std::endl;
	return true;
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::ApplyModalMur(EngType* eng)
{
	if (!m_Op_ABC->IsModalMurReady())
		return;

	const Operator* op = m_Op_ABC->m_Op;
	const unsigned int nTaps = (unsigned int)m_Op_ABC->m_MurTaps.size();
	if (nTaps == 0)
		return;

	unsigned int pos[3] = {0,0,0};

	// ---- 0. on-the-fly modal correction -----------------------------------
	//  Re-learn the template BEFORE the read integral, so an update takes
	//  effect in the same timestep it was measured -- the prototype's ordering
	//  (its step 8 precedes the absorbers in step 9).
	//
	//  The sense plane is the read plane, optionally pushed further into the
	//  guide. m_normalSign gives the direction from the sheet into the guide,
	//  which is the same shift the read plane itself was built with.
	if (m_otfcSensePos != (unsigned int)-1)
		RelearnModalTemplate(eng, m_Op_ABC->m_MurModeP, m_Op_ABC->m_MurModePP,
		                     m_otfcSensePos);

	// ---- 1. read the modal amplitude one cell inside ----------------------
	//  a = integral( E . m ) dA, with E = V/dl on each edge. The template is
	//  L2-normalised over the aperture, so 'a' is the field amplitude of the
	//  mode and nothing else.
	pos[m_ny] = m_Op_ABC->m_MurReadPos;
	double aIn = 0.0;
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;

			double dA = op->GetNodeArea(m_ny, pos, false);

			// Skip cells the mode does not live in (PEC): they contribute
			// nothing to the projection and must not be written to later.
			double mP = m_Op_ABC->m_MurModeP(i,j);
			if (mP != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyP, pos, false);
				if (dl != 0.0)
					aIn += (eng->EngType::GetVolt(m_nyP, pos) / dl) * mP * dA;
			}

			double mPP = m_Op_ABC->m_MurModePP(i,j);
			if (mPP != 0.0)
			{
				double dl = op->GetEdgeLength(m_nyPP, pos, false);
				if (dl != 0.0)
					aIn += (eng->EngType::GetVolt(m_nyPP, pos) / dl) * mPP * dA;
			}
		}
	}

	// ---- 2. delay it by one cell ------------------------------------------
	//  Newest sample first: m_MurHead walks backwards so that history entry k
	//  is always the sample from k steps ago, without moving any data.
	m_MurHead = (m_MurHead == 0) ? (nTaps - 1) : (m_MurHead - 1);
	m_MurHist[m_MurHead] = aIn;

	double aOut = 0.0;
	unsigned int idx = m_MurHead;
	const std::vector<double>& h = m_Op_ABC->m_MurTaps;
	for (unsigned int k = 0; k < nTaps; ++k)
	{
		aOut += h[k] * m_MurHist[idx];
		idx = (idx + 1 == nTaps) ? 0 : (idx + 1);
	}

	// ---- 3. write it onto the sheet plane ---------------------------------
	//  A direct write, not a correction. The sheet is meant to sit on the face
	//  of a PEC block, whose operator coefficients re-zero these edges on every
	//  E update -- so the plane is a clean slate here and nothing accumulates.
	//  Non-modal content on the plane stays zero, i.e. it still sees PEC, which
	//  is exactly what an absorber for ONE mode should do with everything else.
	pos[m_ny] = m_posStart[m_ny];
	for (unsigned int i = 0; i < m_numLines[0]; ++i)
	{
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < m_numLines[1]; ++j)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;

			double mP = m_Op_ABC->m_MurModeP(i,j);
			if (mP != 0.0)
				eng->EngType::SetVolt(m_nyP, pos,
				                      aOut * mP * op->GetEdgeLength(m_nyP, pos, false));

			double mPP = m_Op_ABC->m_MurModePP(i,j);
			if (mPP != 0.0)
				eng->EngType::SetVolt(m_nyPP, pos,
				                      aOut * mPP * op->GetEdgeLength(m_nyPP, pos, false));
		}
	}
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPreCurrentUpdatesImpl(EngType* eng, int threadID)
{

	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	unsigned int 	pos[] = {0,0,0},
					pos_shift[] = {0,0,0};


	// If this isn't the appropriate boundary type, move on to the next primitive
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA)
		return;

	// For magnetic field, -1, due to dual grid
	unsigned int numLines_1 = std::min(
		m_threadStartLine.at(threadID) + m_linesPerThread.at(threadID),
		m_numLines[0] - 1
	);
	unsigned int numLine_0 = m_threadStartLine.at(threadID);

	pos[m_ny] = m_pos_ny0_I;
	pos_shift[m_ny] = m_pos_ny0_shift_I;
	for (unsigned int i = numLine_0 ; i < numLines_1 ; i++)
	{
		// Store shifted location in this container
		pos_shift[m_nyP] = pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < (m_numLines[1] - 1); j++)
		{
			pos_shift[m_nyPP] = pos[m_nyPP] = m_posStart[m_nyPP] + j;

			// H(i + s,n) - K1*H(i,n)
			m_I_nyP (i,j) = eng->EngType::GetCurr(m_nyP ,pos_shift) - m_K1_nyP (i,j)*eng->EngType::GetCurr(m_nyP ,pos);
			m_I_nyPP(i,j) = eng->EngType::GetCurr(m_nyPP,pos_shift) - m_K1_nyPP(i,j)*eng->EngType::GetCurr(m_nyPP,pos);
		}
	}

}

void Engine_Ext_Absorbing_BC::DoPreCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreCurrentUpdatesImpl, threadID);
}

// Super-absorption:
//
// 1. Re-iterate the Mur B.C.
// Hsa(i,n + 1) = H(i + s,n) + K1*[H(i,n + 1) - H(i,n)]
//
// 2. Update the H(i,n + 1) as such
// H(i,n + 1) = (K2*Hsa(i,n + 1) + Hc(i,n + 1))/(K2 + 1)
// Where Hsa(i,n + 1) is the field calculated by the boundary condition, and
// Hc(i,n + 1) is the field calculated by the FDTD step.
// and K2 = vp*Dt/Dx

template <typename EngType>
void Engine_Ext_Absorbing_BC::DoPostCurrentUpdatesImpl(EngType* eng, int threadID)
{

	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	unsigned int pos_shift[] = {0,0,0};

	// If this isn't the appropriate boundary type, move on to the next primitive
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA)
		return;

	// For magnetic field, -1, due to dual grid
	unsigned int numLine_1 = std::min<unsigned int>(
		m_threadStartLine.at(threadID) + m_linesPerThread.at(threadID),
		m_numLines[0] - 1
	);
	unsigned int numLine_0 = m_threadStartLine.at(threadID);

	pos_shift[m_ny] = m_pos_ny0_shift_I;
	for (unsigned int i = numLine_0 ; i < numLine_1 ; i++)
	{
		// Store shifted location in this container
		pos_shift[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < (m_numLines[1] - 1); j++)
		{
			pos_shift[m_nyPP] = m_posStart[m_nyPP] + j;

			// H(i + s,n) + K1*[H(i,n + 1) - H(i,n)]
			m_I_nyP (i,j) += m_K1_nyP (i,j)*eng->EngType::GetCurr(m_nyP ,pos_shift);
			m_I_nyPP(i,j) += m_K1_nyPP(i,j)*eng->EngType::GetCurr(m_nyPP,pos_shift);

		}
	}
}

void Engine_Ext_Absorbing_BC::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_Absorbing_BC::Apply2CurrentImpl(EngType* eng, int threadID)
{
	if (IsActive()==false) return;

	if (m_Eng==NULL) return;

	if (threadID >= m_NrThreads)
		return;

	unsigned int pos[] = {0,0,0};

	// Super-absorption is the only ABC type that touches currents here. Skip
	// for plain Mur and for the modal absorber (the latter does its work in
	// DoPreVoltageUpdatesImpl).
	if ((Operator_Ext_Absorbing_BC::ABCtype)(m_ABCtype) != Operator_Ext_Absorbing_BC::MUR_1ST_SA)
		return;

	// For magnetic field, -1, due to dual grid
	unsigned int numLine_1 = std::min<unsigned int>(
		m_threadStartLine.at(threadID) + m_linesPerThread.at(threadID),
		m_numLines[0] - 1
	);
	unsigned int numLine_0 = m_threadStartLine.at(threadID);

	pos[m_ny] = m_pos_ny0_I;
	for (unsigned int i = numLine_0 ; i < numLine_1 ; i++)
	{
		// Store shifted location in this container
		pos[m_nyP] = m_posStart[m_nyP] + i;
		for (unsigned int j = 0; j < (m_numLines[1] - 1); j++)
		{
			pos[m_nyPP] = m_posStart[m_nyPP] + j;

			// H(i + s,n) = (Hsa*K2 + Hc)/(1 + K2)
			eng->EngType::SetCurr(m_nyP ,pos, (m_I_nyP (i,j)*m_K2_nyP (i,j) + eng->EngType::GetCurr(m_nyP ,pos))/(m_K2_nyP (i,j) + 1.0));
			eng->EngType::SetCurr(m_nyPP,pos, (m_I_nyPP(i,j)*m_K2_nyPP(i,j) + eng->EngType::GetCurr(m_nyPP,pos))/(m_K2_nyPP(i,j) + 1.0));
		}
	}
}

void Engine_Ext_Absorbing_BC::Apply2Current(int threadID)
{
	ENG_DISPATCH_ARGS(Apply2CurrentImpl, threadID);
}
