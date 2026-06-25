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
#include "FDTD/engine_sse.h"
#include "FDTD/engine_interface_fdtd.h"
#include "tools/array_ops.h"
#include "tools/useful.h"
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

	m_start_TS = 0;

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

		// Modal amplitude of the wave leaving the domain through this absorber:
		//   normalSign = +1 (outward normal +ny): outgoing wave is a_pos = 0.5*(Emm + Zw*Hmm)
		//   normalSign = -1 (outward normal -ny): outgoing wave is a_neg = 0.5*(Emm - Zw*Hmm)
		// double a = 0.5 * (Emm + m_normalSign * m_Zw * Hmm);
		double a = 0.5 * (Emm - m_normalSign * m_Zw * Hmm);

		const Operator* op = m_Op_ABC->m_Op;

		unsigned int numLinesE_P  = m_Eng_Interface->GetModeMatchE_NumLines(0);
		unsigned int numLinesE_PP = m_Eng_Interface->GetModeMatchE_NumLines(1);
		unsigned int numLinesH_P  = m_Eng_Interface->GetModeMatchH_NumLines(0);
		unsigned int numLinesH_PP = m_Eng_Interface->GetModeMatchH_NumLines(1);

		// V correction at the E plane: V_comp -= a * modeE[comp][i][j] * EdgeLength_comp
		unsigned int pos_v[] = {0,0,0};
		pos_v[m_ny] = m_posStart[m_ny];
		for (unsigned int i = 0; i < numLinesE_P; ++i)
		{
			pos_v[m_nyP] = m_posStart[m_nyP] + i;
			for (unsigned int j = 0; j < numLinesE_PP; ++j)
			{
				pos_v[m_nyPP] = m_posStart[m_nyPP] + j;
				double el_nyP  = op->GetEdgeLength(m_nyP,  pos_v, false);
				double el_nyPP = op->GetEdgeLength(m_nyPP, pos_v, false);
				eng->EngType::SetVolt(m_nyP,  pos_v, eng->EngType::GetVolt(m_nyP,  pos_v) - a * m_Eng_Interface->GetModeDistE(0, i, j) * el_nyP);
				eng->EngType::SetVolt(m_nyPP, pos_v, eng->EngType::GetVolt(m_nyPP, pos_v) - a * m_Eng_Interface->GetModeDistE(1, i, j) * el_nyPP);
//				eng->EngType::SetVolt(m_nyP,  pos_v, eng->EngType::GetVolt(m_nyP,  pos_v) - a * m_Eng_Interface->GetModeDistE(0, i, j));
//				eng->EngType::SetVolt(m_nyPP, pos_v, eng->EngType::GetVolt(m_nyPP, pos_v) - a * m_Eng_Interface->GetModeDistE(1, i, j));

			}
		}

		// I correction at the H plane: I_comp -= normalSign * (a/Zw) * modeH[comp][i][j] * EdgeLength_comp_dual
		double dH_factor = -m_normalSign * a / m_Zw;
		unsigned int pos_i[] = {0,0,0};
		pos_i[m_ny] = m_pos_ny0_I;
		for (unsigned int i = 0; i < numLinesH_P; ++i)
		{
			pos_i[m_nyP] = m_posStart[m_nyP] + i;
			for (unsigned int j = 0; j < numLinesH_PP; ++j)
			{
				pos_i[m_nyPP] = m_posStart[m_nyPP] + j;
				double el_nyP_d  = op->GetEdgeLength(m_nyP,  pos_i, true);
				double el_nyPP_d = op->GetEdgeLength(m_nyPP, pos_i, true);
				eng->EngType::SetCurr(m_nyP,  pos_i, eng->EngType::GetCurr(m_nyP,  pos_i) - dH_factor * m_Eng_Interface->GetModeDistH(0, i, j) * el_nyP_d);
				eng->EngType::SetCurr(m_nyPP, pos_i, eng->EngType::GetCurr(m_nyPP, pos_i) - dH_factor * m_Eng_Interface->GetModeDistH(1, i, j) * el_nyPP_d);
//				eng->EngType::SetCurr(m_nyP,  pos_i, eng->EngType::GetCurr(m_nyP,  pos_i) - dH_factor * m_Eng_Interface->GetModeDistH(0, i, j));
//				eng->EngType::SetCurr(m_nyPP, pos_i, eng->EngType::GetCurr(m_nyPP, pos_i) - dH_factor * m_Eng_Interface->GetModeDistH(1, i, j));

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
