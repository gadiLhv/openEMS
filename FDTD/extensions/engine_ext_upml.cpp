/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_ext_upml.h"
#include "operator_ext_upml.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/useful.h"

Engine_Ext_UPML::Engine_Ext_UPML(Operator_Ext_UPML* op_ext) : Engine_Extension(op_ext)
{
	m_Op_UPML = op_ext;

	//this ABC extension should be executed first!
	m_Priority = ENG_EXT_PRIO_UPML;

	volt_flux.Init("volt_flux", m_Op_UPML->m_numLines);
	curr_flux.Init("curr_flux", m_Op_UPML->m_numLines);
	if (m_Op_UPML->m_CFS)
	{
		const char* nm[2][4] = {{"cfs_vG","cfs_vD","cfs_vF","cfs_vE"},{"cfs_iG","cfs_iD","cfs_iF","cfs_iE"}};
		for (int f=0; f<2; ++f)
		{
			cfs_G[f].Init(nm[f][0], m_Op_UPML->m_numLines);
			cfs_D[f].Init(nm[f][1], m_Op_UPML->m_numLines);
			cfs_F[f].Init(nm[f][2], m_Op_UPML->m_numLines);
			cfs_E[f].Init(nm[f][3], m_Op_UPML->m_numLines);
		}
	}

	SetNumberOfThreads(1);
}

Engine_Ext_UPML::~Engine_Ext_UPML()
{
}

void Engine_Ext_UPML::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_Op_UPML->m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}

template <typename EngType>
void Engine_Ext_UPML::DoPreVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng==NULL)
		return;

	if (threadID>=m_NrThreads)
		return;

	if (m_Op_UPML->m_CFS)
	{
		CFSPreImpl(eng, threadID, 0);
		return;
	}

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = m_Op_UPML->vv(0, loc_pos[0], loc_pos[1], loc_pos[2])   * eng->EngType::GetVolt(0,pos)
						 - m_Op_UPML->vvfo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetVolt(0,pos, volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
				volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

				f_help = m_Op_UPML->vv(1, loc_pos[0], loc_pos[1], loc_pos[2])   * eng->EngType::GetVolt(1,pos)
						 - m_Op_UPML->vvfo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetVolt(1,pos, volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
				volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

				f_help = m_Op_UPML->vv(2, loc_pos[0], loc_pos[1], loc_pos[2])   * eng->EngType::GetVolt(2,pos)
						 - m_Op_UPML->vvfo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetVolt(2,pos, volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
				volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;
			}
		}
	}
}

void Engine_Ext_UPML::DoPreVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_UPML::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	if (m_Op_UPML->m_CFS)
	{
		CFSPostImpl(eng, threadID, 0);
		return;
	}

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = eng->EngType::GetVolt(0,pos);
				eng->EngType::SetVolt(0,pos, f_help + m_Op_UPML->vvfn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

				f_help = volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = eng->EngType::GetVolt(1,pos);
				eng->EngType::SetVolt(1,pos, f_help + m_Op_UPML->vvfn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

				f_help = volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
				volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = eng->EngType::GetVolt(2,pos);
				eng->EngType::SetVolt(2,pos, f_help + m_Op_UPML->vvfn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * volt_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
			}
		}
	}
}

void Engine_Ext_UPML::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_UPML::DoPreCurrentUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	if (m_Op_UPML->m_CFS)
	{
		CFSPreImpl(eng, threadID, 1);
		return;
	}

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = m_Op_UPML->ii(0, loc_pos[0], loc_pos[1], loc_pos[2])   * eng->EngType::GetCurr(0,pos)
						 - m_Op_UPML->iifo(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetCurr(0,pos, curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));
				curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

				f_help = m_Op_UPML->ii(1, loc_pos[0], loc_pos[1], loc_pos[2])   * eng->EngType::GetCurr(1,pos)
						 - m_Op_UPML->iifo(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetCurr(1,pos, curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));
				curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

				f_help = m_Op_UPML->ii(2, loc_pos[0], loc_pos[1], loc_pos[2])   * eng->EngType::GetCurr(2,pos)
						 - m_Op_UPML->iifo(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetCurr(2,pos, curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
				curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = f_help;

			}
		}
	}
}

void Engine_Ext_UPML::DoPreCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPreCurrentUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_UPML::DoPostCurrentUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng==NULL)
		return;
	if (threadID>=m_NrThreads)
		return;

	if (m_Op_UPML->m_CFS)
	{
		CFSPostImpl(eng, threadID, 1);
		return;
	}

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT f_help;

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_UPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_Op_UPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_UPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_Op_UPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_UPML->m_StartPos[2];

				f_help = curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]) = eng->EngType::GetCurr(0,pos);
				eng->EngType::SetCurr(0,pos, f_help + m_Op_UPML->iifn(0, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(0, loc_pos[0], loc_pos[1], loc_pos[2]));

				f_help = curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]) = eng->EngType::GetCurr(1,pos);
				eng->EngType::SetCurr(1,pos, f_help + m_Op_UPML->iifn(1, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(1, loc_pos[0], loc_pos[1], loc_pos[2]));

				f_help = curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]);
				curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]) = eng->EngType::GetCurr(2,pos);
				eng->EngType::SetCurr(2,pos, f_help + m_Op_UPML->iifn(2, loc_pos[0], loc_pos[1], loc_pos[2]) * curr_flux(2, loc_pos[0], loc_pos[1], loc_pos[2]));
			}
		}
	}
}

void Engine_Ext_UPML::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_UPML::CFSPreImpl(EngType* eng, int threadID, int fld)
{
	// hand the main engine the integrated curl G; keep the old field
	const Operator_Ext_UPML* op = m_Op_UPML;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& on = fld ? op->i_on : op->v_on;
	unsigned int pos[3], l[3];
	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		l[0]=lineX+m_start.at(threadID);
		pos[0] = l[0] + op->m_StartPos[0];
		for (l[1]=0; l[1]<op->m_numLines[1]; ++l[1])
		{
			pos[1] = l[1] + op->m_StartPos[1];
			for (l[2]=0; l[2]<op->m_numLines[2]; ++l[2])
			{
				pos[2] = l[2] + op->m_StartPos[2];
				for (int n=0; n<3; ++n)
				{
					if (on(n,l[0],l[1],l[2])==0)
						continue;
					if (fld==0)
					{
						cfs_E[0](n,l[0],l[1],l[2]) = eng->EngType::GetVolt(n,pos);
						eng->EngType::SetVolt(n,pos, cfs_G[0](n,l[0],l[1],l[2]));
					}
					else
					{
						cfs_E[1](n,l[0],l[1],l[2]) = eng->EngType::GetCurr(n,pos);
						eng->EngType::SetCurr(n,pos, cfs_G[1](n,l[0],l[1],l[2]));
					}
				}
			}
		}
	}
}

template <typename EngType>
void Engine_Ext_UPML::CFSPostImpl(EngType* eng, int threadID, int fld)
{
	// G -> D (s_nP) -> F (s_n) -> E (s_nPP, eps/mu): three bilinear first-order filters
	const Operator_Ext_UPML* op = m_Op_UPML;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& on  = fld ? op->i_on  : op->v_on;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& dd  = fld ? op->i_dd  : op->v_dd;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& dgn = fld ? op->i_dgn : op->v_dgn;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& dgo = fld ? op->i_dgo : op->v_dgo;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& ff  = fld ? op->i_ff  : op->v_ff;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& fdn = fld ? op->i_fdn : op->v_fdn;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& fdo = fld ? op->i_fdo : op->v_fdo;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& ee  = fld ? op->i_ee  : op->v_ee;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& efn = fld ? op->i_efn : op->v_efn;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& efo = fld ? op->i_efo : op->v_efo;
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& G = cfs_G[fld];
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& D = cfs_D[fld];
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& F = cfs_F[fld];
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& Eo = cfs_E[fld];
	unsigned int pos[3], l[3];
	FDTD_FLOAT Gn, Dn, Fn;
	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		l[0]=lineX+m_start.at(threadID);
		pos[0] = l[0] + op->m_StartPos[0];
		for (l[1]=0; l[1]<op->m_numLines[1]; ++l[1])
		{
			pos[1] = l[1] + op->m_StartPos[1];
			for (l[2]=0; l[2]<op->m_numLines[2]; ++l[2])
			{
				pos[2] = l[2] + op->m_StartPos[2];
				for (int n=0; n<3; ++n)
				{
					const unsigned int i=l[0], j=l[1], k=l[2];
					if (on(n,i,j,k)==0)
						continue;
					Gn = fld ? eng->EngType::GetCurr(n,pos) : eng->EngType::GetVolt(n,pos);
					Dn = dd(n,i,j,k)*D(n,i,j,k) + dgn(n,i,j,k)*Gn - dgo(n,i,j,k)*G(n,i,j,k);
					Fn = ff(n,i,j,k)*F(n,i,j,k) + fdn(n,i,j,k)*Dn - fdo(n,i,j,k)*D(n,i,j,k);
					FDTD_FLOAT En = ee(n,i,j,k)*Eo(n,i,j,k) + efn(n,i,j,k)*Fn - efo(n,i,j,k)*F(n,i,j,k);
					G(n,i,j,k) = Gn;  D(n,i,j,k) = Dn;  F(n,i,j,k) = Fn;
					if (fld)
						eng->EngType::SetCurr(n,pos, En);
					else
						eng->EngType::SetVolt(n,pos, En);
				}
			}
		}
	}
}
