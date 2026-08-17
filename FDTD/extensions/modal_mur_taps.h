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

#ifndef MODAL_MUR_TAPS_H
#define MODAL_MUR_TAPS_H

#include <vector>

//! Tap design for the Dispersive Modal Mur one-way termination.
/*!
  The one-way condition needs no impedance. For a mode leaving the domain
  through the sheet plane, the modal amplitude there is simply the DELAYED
  amplitude one plane inside:

      a_k(w) = a_{k-1}(w) * exp(-j*beta(w)*dz)

  and the sheet plane's modal component is OVERWRITTEN with it.

  Two properties make this work where the impedance formulation could not:

    - |exp(-j*beta*dz)| = 1 in band, -> 1 at cutoff, and exp(-alpha*dz) < 1
      below cutoff. Bounded by 1 EVERYWHERE, so the filter is passive BY
      CONSTRUCTION: no drain cap, no stability wall, no passivity guard, and
      evanescent content is handled exactly (a real decay per cell) rather
      than being misread by a direction test that cannot work there.

    - Because it is an overwrite, information flows one way. Whatever sits
      behind the sheet can be driven but can never drive back, so the
      feedback loop that destabilised the injection form cannot form at all.

  beta comes from the EXACT DISCRETE LATTICE dispersion relation, not the
  continuum sqrt(k0^2 - kc^2). Using the continuum form leaves a mismatch
  that grows with cell size; the lattice form is exact for the grid the
  fields actually live on.
  */
namespace ModalMur
{
	//! Design the per-cell one-way delay FIR exp(-j*beta*dz).
	/*!
	  \param fc       modal cutoff [Hz]. May be zero or NEGATIVE; negative is
	                  taken as kc^2 = -(2*pi*fc/c0)^2 (TEM / quasi-TEM lines).
	  \param dz       spacing between the read plane and the sheet plane [m].
	  \param dt       FDTD timestep [s].
	  \param numTaps  FIR length (compile-time predef at the call site).
	  \param taps     out: numTaps filter coefficients, newest sample first.
	  \param maxAbsH  out: max |H(f)| of the TRUNCATED filter over [0, fNyq].
	                  Must stay at or just above 1; a value meaningfully above
	                  1 means the truncation has broken passivity and the run
	                  can grow without bound.
	  \param acausalFrac out: fraction of the untruncated kernel's energy that
	                  landed in the second half of the DFT period, i.e. at
	                  negative time. Small (~1e-4) for a well-posed target.
	  \return false if the parameters are unusable (non-positive dz/dt/numTaps).
	  */
	bool DesignDelayTaps(double fc, double dz, double dt, unsigned int numTaps,
	                     std::vector<double>& taps,
	                     double& maxAbsH, double& acausalFrac);
}

#endif // MODAL_MUR_TAPS_H
