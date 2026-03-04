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

#ifndef ENGINE_EXT_MODEABSORB_H
#define ENGINE_EXT_MODEABSORB_H

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_ModeAbsorb;

class Engine_Ext_ModeAbsorb : public Engine_Extension
{
public:
	Engine_Ext_ModeAbsorb(Operator_Ext_ModeAbsorb* op_ext);
	virtual ~Engine_Ext_ModeAbsorb();

	// E-field absorption: compute overlap after voltage update, subtract in Apply2Voltages
	virtual void DoPostVoltageUpdates() {Engine_Ext_ModeAbsorb::DoPostVoltageUpdates(0);}
	virtual void DoPostVoltageUpdates(int threadID);
	virtual void Apply2Voltages() {Engine_Ext_ModeAbsorb::Apply2Voltages(0);}
	virtual void Apply2Voltages(int threadID);

	// H-field absorption: compute overlap after current update, subtract in Apply2Current
	virtual void DoPostCurrentUpdates() {Engine_Ext_ModeAbsorb::DoPostCurrentUpdates(0);}
	virtual void DoPostCurrentUpdates(int threadID);
	virtual void Apply2Current() {Engine_Ext_ModeAbsorb::Apply2Current(0);}
	virtual void Apply2Current(int threadID);

protected:
	Operator_Ext_ModeAbsorb* m_Op_MA;

	template <typename EngType>
	void DoPostVoltageUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2VoltagesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPostCurrentUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void Apply2CurrentImpl(EngType* eng, int threadID);

	int m_ny, m_nyP, m_nyPP;

	unsigned int m_posStart[3];

	// E mode weights (primary mesh)
	unsigned int m_numLines_E[2];
	double** m_E_OverlapW[2];
	double** m_E_SubtractW[2];

	// H mode weights (dual mesh)
	unsigned int m_numLines_H[2];
	double** m_H_OverlapW[2];
	double** m_H_SubtractW[2];

	// Overlap coefficients computed per-timestep
	double m_a_E;
	double m_a_H;
};

#endif // ENGINE_EXT_MODEABSORB_H
