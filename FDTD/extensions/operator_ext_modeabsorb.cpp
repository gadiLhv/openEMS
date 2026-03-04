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

#include "operator_ext_modeabsorb.h"
#include "engine_ext_modeabsorb.h"

#include "tools/array_ops.h"

#include "CSPrimBox.h"
#include "CSPropModeAbsorb.h"
#include "CSModeFileParser.h"

using std::cerr;
using std::endl;

Operator_Ext_ModeAbsorb::Operator_Ext_ModeAbsorb(Operator* op) : Operator_Extension(op)
{
	Initialize();
}

Operator_Ext_ModeAbsorb::~Operator_Ext_ModeAbsorb()
{
	for (int n = 0; n < 2; ++n)
	{
		Delete2DArray<double>(m_E_OverlapW[n], m_numLines_E);
		Delete2DArray<double>(m_E_SubtractW[n], m_numLines_E);
		Delete2DArray<double>(m_H_OverlapW[n], m_numLines_H);
		Delete2DArray<double>(m_H_SubtractW[n], m_numLines_H);
	}
}

Operator_Ext_ModeAbsorb::Operator_Ext_ModeAbsorb(Operator* op, Operator_Ext_ModeAbsorb* op_ext) : Operator_Extension(op, op_ext)
{
	Initialize();
}

Operator_Extension* Operator_Ext_ModeAbsorb::Clone(Operator* op)
{
	return new Operator_Ext_ModeAbsorb(op, this);
}

void Operator_Ext_ModeAbsorb::Initialize()
{
	m_ny = m_nyP = m_nyPP = -1;

	for (int i = 0; i < 3; ++i)
	{
		m_sheetX0[i] = 0;
		m_sheetX1[i] = 0;
	}

	m_normalSignPositive = true;
	m_ZWave = 376.73;
	m_dirSign = 1.0;

	m_numLines_E[0] = m_numLines_E[1] = 0;
	m_numLines_H[0] = m_numLines_H[1] = 0;

	for (int n = 0; n < 2; ++n)
	{
		m_E_OverlapW[n] = NULL;
		m_E_SubtractW[n] = NULL;
		m_H_OverlapW[n] = NULL;
		m_H_SubtractW[n] = NULL;
	}
}

bool Operator_Ext_ModeAbsorb::SetInitParams(CSPrimitives* prim, CSPropModeAbsorb* prop)
{
	CSPrimBox* cSheet = dynamic_cast<CSPrimBox*>(prim);
	if (!cSheet)
	{
		cerr << "Operator_Ext_ModeAbsorb::SetInitParams(): Warning: primitive is not a box, skipping."
			 << " ID: " << prim->GetID() << " @ Property: " << prop->GetName() << endl;
		return false;
	}

	int Snap_Dimension =
		m_Op->SnapBox2Mesh(
			cSheet->GetStartCoord()->GetCoords(m_Op->m_MeshType),
			cSheet->GetStopCoord()->GetCoords(m_Op->m_MeshType),
			m_sheetX0,
			m_sheetX1,
			false,
			true);

	if (Snap_Dimension <= 0)
	{
		if (Snap_Dimension >= -1)
			cerr << "Operator_Ext_ModeAbsorb::SetInitParams(): Warning: snapping failed! Dimension is: " << Snap_Dimension << " skipping."
				 << " ID: " << prim->GetID() << " @ Property: " << prop->GetName() << endl;
		return false;
	}

	// Find the normal direction (collapsed dimension)
	unsigned int sheetCheck = 0;
	m_ny = -1;
	for (int dimIdx = 0; dimIdx < 3; dimIdx++)
	{
		unsigned int Ncells = m_sheetX1[dimIdx] - m_sheetX0[dimIdx] + 1;
		sheetCheck += (Ncells == 1);
		m_ny = (Ncells == 1) ? dimIdx : m_ny;
	}

	if (sheetCheck != 1)
	{
		cerr << "Operator_Ext_ModeAbsorb::SetInitParams(): Warning: primitive is not a sheet! Skipping."
			 << " ID: " << prim->GetID() << " @ Property: " << prop->GetName() << endl;
		return false;
	}

	m_normalSignPositive = prop->GetNormalSignPositive();
	m_ZWave = prop->GetWaveImpedance();
	m_EModeFileName = prop->GetEModeFileName();
	m_HModeFileName = prop->GetHModeFileName();

	prim->SetPrimitiveUsed(true);

	return true;
}

bool Operator_Ext_ModeAbsorb::BuildExtension()
{
	m_nyP  = (m_ny + 1) % 3;
	m_nyPP = (m_ny + 2) % 3;

	// Directional sign for forward/backward wave decomposition.
	m_dirSign = m_normalSignPositive ? 1.0 : -1.0;

	unsigned int Ncells[3];
	for (int i = 0; i < 3; ++i)
		Ncells[i] = m_sheetX1[i] - m_sheetX0[i] + 1;

	// --- E-field mode on primary mesh (dualMesh=false) ---
	m_numLines_E[0] = Ncells[m_nyP];
	m_numLines_E[1] = Ncells[m_nyPP];

	for (int n = 0; n < 2; ++n)
	{
		m_E_OverlapW[n]  = Create2DArray<double>(m_numLines_E);
		m_E_SubtractW[n] = Create2DArray<double>(m_numLines_E);
	}

	{
		CSModeFileParser modeFile;
		modeFile.ParseFile(m_EModeFileName);

		// First pass: parse mode values and compute normalization
		double** modeTmp[2];
		for (int n = 0; n < 2; ++n)
			modeTmp[n] = Create2DArray<double>(m_numLines_E);

		unsigned int pos[3] = {0, 0, 0};
		double discLine[3] = {0, 0, 0};
		pos[m_ny] = m_sheetX0[m_ny];
		discLine[m_ny] = m_Op->GetDiscLine(m_ny, pos[m_ny], false);
		double norm = 0;

		for (unsigned int posP = 0; posP < m_numLines_E[0]; ++posP)
		{
			pos[m_nyP] = m_sheetX0[m_nyP] + posP;
			discLine[m_nyP] = m_Op->GetDiscLine(m_nyP, pos[m_nyP], false);
			for (unsigned int posPP = 0; posPP < m_numLines_E[1]; ++posPP)
			{
				pos[m_nyPP] = m_sheetX0[m_nyPP] + posPP;
				discLine[m_nyPP] = m_Op->GetDiscLine(m_nyPP, pos[m_nyPP], false);

				double locCoord[3] = {
					discLine[0] - m_Op->GetDiscLine(0, m_sheetX0[0], false),
					discLine[1] - m_Op->GetDiscLine(1, m_sheetX0[1], false),
					discLine[2] - m_Op->GetDiscLine(2, m_sheetX0[2], false)
				};

				modeFile.LinInterp2(locCoord[m_nyP], locCoord[m_nyPP],
					modeTmp[0][posP][posPP], modeTmp[1][posP][posPP]);

				double area = m_Op->GetNodeArea(m_ny, pos, false);
				for (int n = 0; n < 2; ++n)
					norm += modeTmp[n][posP][posPP] * modeTmp[n][posP][posPP] * area;
			}
		}

		// Normalize and compute overlap/subtract weights
		norm = sqrt(norm);
		if (norm > 0)
		{
			for (unsigned int posP = 0; posP < m_numLines_E[0]; ++posP)
			{
				pos[m_nyP] = m_sheetX0[m_nyP] + posP;
				for (unsigned int posPP = 0; posPP < m_numLines_E[1]; ++posPP)
				{
					pos[m_nyPP] = m_sheetX0[m_nyPP] + posPP;

					double area = m_Op->GetNodeArea(m_ny, pos, false);

					// Edge lengths for each transverse component
					double edgeLen_P  = m_Op->GetEdgeLength(m_nyP,  pos, false);
					double edgeLen_PP = m_Op->GetEdgeLength(m_nyPP, pos, false);

					double mP  = modeTmp[0][posP][posPP] / norm;
					double mPP = modeTmp[1][posP][posPP] / norm;

					// OverlapW: mode_norm * area / edgeLen
					// (multiply by Volt to get field*mode*area contribution to overlap)
					m_E_OverlapW[0][posP][posPP] = (edgeLen_P  > 0) ? mP  * area / edgeLen_P  : 0;
					m_E_OverlapW[1][posP][posPP] = (edgeLen_PP > 0) ? mPP * area / edgeLen_PP : 0;

					// SubtractW: * mode_norm * edgeLen
					// (multiply by overlap coeff to get voltage/current correction)
					m_E_SubtractW[0][posP][posPP] = mP  * edgeLen_P;
					m_E_SubtractW[1][posP][posPP] = mPP * edgeLen_PP;
				}
			}
		}

		for (int n = 0; n < 2; ++n)
			Delete2DArray<double>(modeTmp[n], m_numLines_E);
	}

	// --- H-field mode on dual mesh (dualMesh=true) ---
	m_numLines_H[0] = Ncells[m_nyP]  > 0 ? Ncells[m_nyP]  - 1 : 0;
	m_numLines_H[1] = Ncells[m_nyPP] > 0 ? Ncells[m_nyPP] - 1 : 0;

	for (int n = 0; n < 2; ++n)
	{
		m_H_OverlapW[n]  = Create2DArray<double>(m_numLines_H);
		m_H_SubtractW[n] = Create2DArray<double>(m_numLines_H);
	}

	{
		CSModeFileParser modeFile;
		modeFile.ParseFile(m_HModeFileName);

		double** modeTmp[2];
		for (int n = 0; n < 2; ++n)
			modeTmp[n] = Create2DArray<double>(m_numLines_H);

		unsigned int pos[3] = {0, 0, 0};
		double discLine[3] = {0, 0, 0};
		pos[m_ny] = m_sheetX0[m_ny];
		discLine[m_ny] = m_Op->GetDiscLine(m_ny, pos[m_ny], true);
		double norm = 0;

		for (unsigned int posP = 0; posP < m_numLines_H[0]; ++posP)
		{
			pos[m_nyP] = m_sheetX0[m_nyP] + posP;
			discLine[m_nyP] = m_Op->GetDiscLine(m_nyP, pos[m_nyP], true);
			for (unsigned int posPP = 0; posPP < m_numLines_H[1]; ++posPP)
			{
				pos[m_nyPP] = m_sheetX0[m_nyPP] + posPP;
				discLine[m_nyPP] = m_Op->GetDiscLine(m_nyPP, pos[m_nyPP], true);

				double locCoord[3] = {
					discLine[0] - m_Op->GetDiscLine(0, m_sheetX0[0], true),
					discLine[1] - m_Op->GetDiscLine(1, m_sheetX0[1], true),
					discLine[2] - m_Op->GetDiscLine(2, m_sheetX0[2], true)
				};

				modeFile.LinInterp2(locCoord[m_nyP], locCoord[m_nyPP],
					modeTmp[0][posP][posPP], modeTmp[1][posP][posPP]);

				double area = m_Op->GetNodeArea(m_ny, pos, true);
				for (int n = 0; n < 2; ++n)
					norm += modeTmp[n][posP][posPP] * modeTmp[n][posP][posPP] * area;
			}
		}

		norm = sqrt(norm);
		if (norm > 0)
		{
			for (unsigned int posP = 0; posP < m_numLines_H[0]; ++posP)
			{
				pos[m_nyP] = m_sheetX0[m_nyP] + posP;
				for (unsigned int posPP = 0; posPP < m_numLines_H[1]; ++posPP)
				{
					pos[m_nyPP] = m_sheetX0[m_nyPP] + posPP;

					double area = m_Op->GetNodeArea(m_ny, pos, true);

					double edgeLen_P  = m_Op->GetEdgeLength(m_nyP,  pos, true);
					double edgeLen_PP = m_Op->GetEdgeLength(m_nyPP, pos, true);

					double mP  = modeTmp[0][posP][posPP] / norm;
					double mPP = modeTmp[1][posP][posPP] / norm;

					m_H_OverlapW[0][posP][posPP] = (edgeLen_P  > 0) ? mP  * area / edgeLen_P  : 0;
					m_H_OverlapW[1][posP][posPP] = (edgeLen_PP > 0) ? mPP * area / edgeLen_PP : 0;

					m_H_SubtractW[0][posP][posPP] = mP  * edgeLen_P;
					m_H_SubtractW[1][posP][posPP] = mPP * edgeLen_PP;
				}
			}
		}

		for (int n = 0; n < 2; ++n)
			Delete2DArray<double>(modeTmp[n], m_numLines_H);
	}

	return true;
}

Engine_Extension* Operator_Ext_ModeAbsorb::CreateEngineExtention()
{
	Engine_Ext_ModeAbsorb* eng_ext = new Engine_Ext_ModeAbsorb(this);
	return eng_ext;
}

void Operator_Ext_ModeAbsorb::ShowStat(std::ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
	ostr << " E-mode lines: " << m_numLines_E[0] << " x " << m_numLines_E[1] << endl;
	ostr << " H-mode lines: " << m_numLines_H[0] << " x " << m_numLines_H[1] << endl;
	ostr << " Wave impedance: " << m_ZWave << " Ohm" << endl;
	ostr << " Direction sign: " << m_dirSign << endl;
}
