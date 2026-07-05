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

#ifndef PROCESSMODEMATCH_H
#define PROCESSMODEMATCH_H

#include "processintegral.h"
#include "CSModeFileParser.h"

class CSFunctionParser;

//! Processing class to match a mode to a given analytic function and return the integral value.
/*!
  The analytric function has to be defined in drawing units.
  It will return the integral value and the mode purity as a secondary value.
  */
class ProcessModeMatch : public ProcessIntegral
{
public:
	ProcessModeMatch(Engine_Interface_Base* eng_if);
	virtual ~ProcessModeMatch();

	virtual std::string GetProcessingName() const;

	virtual std::string GetIntegralName(int row) const;

	virtual void InitProcess();
	virtual void Reset();

	//! Set the field type (0 electric field, 1 magnetic field)
	void SetFieldType(int type);
	//! Set the mode function in the given direction ny. For example: SetModeFunction(0,"cos(pi/1000*x)*sin(pi/500*y)");
	void SetModeFunction(int ny, std::string function);

	void		SetModeFileName(std::string fileName);
	std::string GetModeFileName() {return m_ModeFileName;};
	bool		GetFieldSourceIsFile() {return m_FieldSourceIsFile;};
	void		SetFieldSourceIsFile(bool isFile) {m_FieldSourceIsFile = isFile;};

	virtual int		GetNumberOfIntegrals() const {return 2;}
	virtual double*	CalcMultipleIntegrals();

	//! Access the precomputed mode distribution for field component 0 (nyP) or 1 (nyPP).
	//! Valid only after InitProcess() has been called.
	const double* const* GetModeDist(int component) const { return m_ModeDist[component]; }

	//! Fill \p lines with the number of grid lines in each transverse dimension.
	void GetNumLines(unsigned int lines[2]) const { lines[0] = m_numLines[0]; lines[1] = m_numLines[1]; }

	//! Fill \p s with the snapped start indices of the mode plane (valid after InitProcess).
	//! Consumers applying per-cell corrections MUST use these indices, not their own snap.
	void GetStartPos(unsigned int s[3]) const { for (int n=0;n<3;++n) s[n] = start[n]; }

	//! Set the local-frame origin (drawing units) used to look up mode-file coordinates.
	//! Must be the physical start corner of the defining primitive box (the same frame
	//! the excitation uses via shiftCoordsForModeFile). Without it the origin falls back
	//! to the snapped start mesh line, which on the dual mesh is offset by half a cell.
	void SetModeFileOrigin(const double origin[3])
	{ for (int n=0;n<3;++n) m_ModeFileOrigin[n] = origin[n]; m_ModeFileOriginSet = true; }

	//! Yee-consistent sampling mode (for modal absorbers).
	//! Fields are read RAW per Yee edge (no node interpolation) and the mode template
	//! is evaluated at each component's own edge-midpoint coordinate -- the same
	//! convention the mode-file excitation uses. This makes the measurement operator
	//! exactly consistent with the engine's per-edge correction: subtracting
	//! a*mode*edgeLength from the same edges is then measured back as exactly 'a'.
	//! Without this, the absorber's measure/subtract loop has a systematic mismatch
	//! that leaves unabsorbable residue at the plane (and can turn the loop unstable).
	void SetYeeConsistent(bool val) { m_YeeConsistent = val; }

protected:
	//normal direction of the mode plane
	int m_ny;

	int m_ModeFieldType;

	std::string m_ModeFunction[3];
	CSFunctionParser* m_ModeParser[2];

	std::string m_ModeFileName;
	bool		m_FieldSourceIsFile;

	double		m_ModeFileOrigin[3];
	bool		m_ModeFileOriginSet;
	bool		m_YeeConsistent;

	unsigned int m_numLines[2];
	double** m_ModeDist[2];
};

#endif // PROCESSMODEMATCH_H
