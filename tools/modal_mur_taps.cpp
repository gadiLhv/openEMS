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

#include "modal_mur_taps.h"

#include <cmath>
#include <complex>

#ifndef __C0__
#define __C0__ 299792458.0
#endif

namespace
{
	//! DFT length used to sample the target and invert it back to time.
	/*!
	  The exact kernel decays only ALGEBRAICALLY (it is a Bessel tail ringing
	  at the cutoff -- the memory of the cutoff resonance), so it never truly
	  ends and the inverse transform always aliases a little. Sampling at
	  16x the tap count keeps that aliasing far below the truncation error
	  that the tail taper introduces anyway.
	  */
	const unsigned int MODAL_MUR_DFT_LEN = 8192;

	//! Fraction of the tap set covered by the raised-cosine tail taper.
	const double MODAL_MUR_TAPER_FRAC = 0.15;
}

bool ModalMur::DesignDelayTaps(double fc, double v, double dz, double dt,
                               unsigned int numTaps,
                               std::vector<double>& taps,
                               double& maxAbsH, double& acausalFrac)
{
	taps.clear();
	maxAbsH = 0.0;
	acausalFrac = 0.0;

	if ((dz <= 0.0) || (dt <= 0.0) || (v <= 0.0) || (numTaps == 0))
		return false;

	unsigned int Md = MODAL_MUR_DFT_LEN;
	while (Md < 4 * numTaps)
		Md *= 2;

	// Signed squared cutoff. A negative fc is deliberate and means kc^2 < 0:
	// TEM and quasi-TEM lines cut off at (or numerically just below) DC, and
	// this keeps the whole design real instead of forcing complex arithmetic.
	const double kc  = 2.0 * M_PI * fc / v;
	const double kc2 = (fc < 0.0) ? -(kc * kc) : (kc * kc);

	// ---- sample the target on the one-sided DFT grid ----------------------
	const unsigned int Nh = Md / 2;
	std::vector< std::complex<double> > D(Nh + 1);

	for (unsigned int k = 0; k <= Nh; ++k)
	{
		double w  = 2.0 * M_PI * (double)k / ((double)Md * dt);
		double sn = sin(0.5 * w * dt);

		// Exact discrete lattice dispersion relation, same relation the
		// S-parameter de-embedding uses.
		double argd = (2.0 / (v * dt)) * (2.0 / (v * dt)) * sn * sn - kc2;
		double sarg = 0.5 * dz * sqrt(fabs(argd));

		// THREE regimes, all with |D| <= 1. Note dz cancels out of every
		// exponent: (2/dz)*f(sarg)*dz = 2*f(sarg). dz survives only inside
		// sarg itself.
		if (argd < 0.0)
		{
			// (a) below the WAVEGUIDE cutoff: evanescent in z, real decay
			//     per cell. Never fires when kc^2 <= 0, i.e. for TEM.
			D[k] = std::complex<double>(exp(-2.0 * asinh(sarg)), 0.0);
		}
		else if (sarg <= 1.0)
		{
			// (b) propagating on the grid: unit modulus, pure delay
			double ph = 2.0 * asin(sarg < 1.0 ? sarg : 1.0);
			D[k] = std::polar(1.0, -ph);
		}
		else
		{
			// (c) above the GRID cutoff: sin(beta*dz/2) > 1 has no real
			//     solution, so the lattice cannot propagate this fast and
			//     beta*dz/2 = pi/2 - j*acosh(sarg), giving a decaying,
			//     sign-flipped response. Clamping sarg at 1 here instead
			//     (as a first cut did) forces unit modulus where the true
			//     response decays, which makes the target anticausal and
			//     pushes max|H| above 1 -- destroying the one property the
			//     whole scheme rests on.
			D[k] = std::complex<double>(-exp(-2.0 * acosh(sarg)), 0.0);
		}
	}

	// ---- inverse DFT back to a real kernel --------------------------------
	// Hermitian symmetry is implied rather than materialised, so the sum runs
	// over the one-sided grid with the interior bins counted twice. A cos/sin
	// table keyed on (k*n) mod Md keeps this an O(Md^2) table lookup instead
	// of O(Md^2) transcendental calls, which is the difference between
	// milliseconds and seconds at startup.
	std::vector<double> cosTab(Md), sinTab(Md);
	for (unsigned int m = 0; m < Md; ++m)
	{
		double ang = 2.0 * M_PI * (double)m / (double)Md;
		cosTab[m] = cos(ang);
		sinTab[m] = sin(ang);
	}

	std::vector<double> hFull(Md, 0.0);
	for (unsigned int n = 0; n < Md; ++n)
	{
		double acc = D[0].real();                       // k = 0, always real
		for (unsigned int k = 1; k < Nh; ++k)
		{
			unsigned int idx = (unsigned int)(((unsigned long long)k * n) % Md);
			acc += 2.0 * (D[k].real() * cosTab[idx] - D[k].imag() * sinTab[idx]);
		}
		// Nyquist bin, real, alternating sign
		acc += D[Nh].real() * ((n % 2) ? -1.0 : 1.0);
		hFull[n] = acc / (double)Md;
	}

	// Energy that landed at "negative time" (the upper half of the period).
	double eTot = 0.0, eTail = 0.0;
	for (unsigned int n = 0; n < Md; ++n)
	{
		eTot += hFull[n] * hFull[n];
		if (n >= Nh)
			eTail += hFull[n] * hFull[n];
	}
	acausalFrac = (eTot > 0.0) ? sqrt(eTail / eTot) : 0.0;

	// ---- truncate and taper -----------------------------------------------
	taps.assign(hFull.begin(), hFull.begin() + numTaps);

	unsigned int nTaper = (unsigned int)(MODAL_MUR_TAPER_FRAC * numTaps + 0.5);
	if (nTaper > numTaps)
		nTaper = numTaps;
	for (unsigned int m = 1; m <= nTaper; ++m)
		taps[numTaps - nTaper + m - 1] *= 0.5 * (1.0 + cos(M_PI * (double)m / (double)nTaper));

	// ---- verify the REALISED response -------------------------------------
	// What the truncated, tapered tap set actually does -- not what the target
	// asked for. Passivity is the property the scheme rests on, so it has to
	// be checked on the thing that will run, and it has to be checked here
	// rather than discovered as a diverging simulation an hour later.
	const unsigned int Nver = 4096;
	for (unsigned int k = 0; k <= Nver; ++k)
	{
		double wdt = M_PI * (double)k / (double)Nver;   // 0 .. pi  (0 .. fNyq)
		double re = 0.0, im = 0.0;
		for (unsigned int n = 0; n < numTaps; ++n)
		{
			double ang = wdt * (double)n;
			re += taps[n] * cos(ang);
			im -= taps[n] * sin(ang);
		}
		double mag = sqrt(re * re + im * im);
		if (mag > maxAbsH)
			maxAbsH = mag;
	}

	return true;
}
