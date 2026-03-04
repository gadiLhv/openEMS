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

#ifndef OPERATOR_EXT_MODEABSORB_H
#define OPERATOR_EXT_MODEABSORB_H

#include "FDTD/operator.h"
#include "operator_extension.h"

class CSPropModeAbsorb;
class CSPrimitives;

class Operator_Ext_ModeAbsorb : public Operator_Extension
{
	friend class Engine_Ext_ModeAbsorb;
public:
	Operator_Ext_ModeAbsorb(Operator* op);
	~Operator_Ext_ModeAbsorb();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual std::string GetExtensionName() const
	{
		return std::string("Mode-matched waveguide absorber");
	}

	virtual void ShowStat(std::ostream &ostr) const;

	virtual bool SetInitParams(CSPrimitives* prim, CSPropModeAbsorb* prop);

protected:
	Operator_Ext_ModeAbsorb(Operator* op, Operator_Ext_ModeAbsorb* op_ext);
	void Initialize();

	int m_ny, m_nyP, m_nyPP;

	unsigned int m_sheetX0[3];
	unsigned int m_sheetX1[3];

	bool m_normalSignPositive;

	std::string m_EModeFileName;
	std::string m_HModeFileName;

	// E mode on primary mesh: pre-computed weights
	// OverlapW[n] = mode_norm * area / edgeLen (for computing overlap integral from voltages)
	// SubtractW[n] = mode_norm * edgeLen (for subtracting mode component from voltages)
	unsigned int m_numLines_E[2];
	double** m_E_OverlapW[2];
	double** m_E_SubtractW[2];

	// H mode on dual mesh: pre-computed weights
	unsigned int m_numLines_H[2];
	double** m_H_OverlapW[2];
	double** m_H_SubtractW[2];

	// Damping factor: fraction of mode projection to subtract per timestep.
	// Derived from local Courant number for stability: alpha = v*dt/delta.
	double m_alpha;
};

#endif // OPERATOR_EXT_MODEABSORB_H
