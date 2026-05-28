// fwhIntegrate.H

// Problem dimensions
const scalar fourPi = 4.0 * constant::mathematical::pi;
const label  nTimes = timeValues.size();											// Number of sampled time steps
const label  nFaces = faceCentres.size();											// Number of faces on surface
const label  nObs   = observers.size();												// Number of observers

// At least three time steps are required to apply any 3-point derivative stencil
if (nTimes < 3)
{
	FatalErrorInFunction
		<< "At least 3 time steps are required for the quadratic derivative "
		<< "stencils. Collected: " << nTimes
		<< exit(FatalError);
}

// * * * * * * * * * * * * * *  Unit Face Normals  * * * * * * * * * * * * * //
//
// Normalise the area-weighted face normals Sf to obtain unit outward normals.
// The +VSMALL denominator guards against zero-area degenerate faces that
// would otherwise cause a divide-by-zero.

vectorField nHat(nFaces);															// Unit outward normals

forAll(nHat, fI)
{
	nHat[fI] = faceNormals[fI] / (mag(faceNormals[fI]) + VSMALL);					// Normalise; +VSMALL prevents div/0 on any zero-area faces
}

// * * * * * * * * *  Observer-Independent Surface Histories  * * * * * * * * //
//
// LvecHistory[tI][fI]  =  p*n_hat + rho*U*(U.n_hat)
//   This is the permeable-surface loading vector (units ~ Pa).
//   Gets values from the surface, no effect of observers here.
//
// UnHistory[tI][fI]    =  U . n_hat
//   Normal component of the surface velocity (m/s).
//   Enters the thickness (monopole) term directly.
//
// Pre-computing so that no extra effort in the iterative loop.

Info << nl << "Pre-computing surface loading and normal-velocity histories ..." << nl;

List<vectorField> LvecHistory(nTimes, vectorField(nFaces, vector::zero));			// [tI][fI] permeable loading vector at each time and face
List<scalarField> UnHistory  (nTimes, scalarField(nFaces, 0.0));					// [tI][fI] surface-normal velocity at each time and face

for (label tI = 0; tI < nTimes; ++tI)
{
	for (label fI = 0; fI < nFaces; ++fI)
	{
		const scalar p_t = pSurfHistory[tI][fI];									// Pressure at this face and time step
		const vector U_t = USurfHistory[tI][fI];									// Velocity vector at this face and time step
		const scalar Un  = U_t & nHat[fI];											// Normal component of velocity: U . n_hat (dot product)

		LvecHistory[tI][fI] = p_t * nHat[fI] + rho * U_t * Un;
		UnHistory  [tI][fI] = Un;													// Store Un separately for the thickness term
	}
}

// * * * * * * * * * *  Differentiate Histories  * * * * * * * * * * * * * * //
//
// dLvecHistory[tI][fI]  =  d(Lvec)/dt  at stored time tI
// dUnHistory  [tI][fI]  =  d(Un)/dt    at stored time tI
//
// The FW-H integrand involves time derivatives evaluated at the retarded time.
// These are pre-computed here (before the observer loop) using 3-point
// nonuniform Lagrange stencils so they can be interpolated cheaply later.
//
// Interior nodes : centred stencil
//
//   f'(t_i) =  -h2/(h1*(h1+h2)) * f_{i-1}
//             + (h2-h1)/(h1*h2)  * f_i
//             + h1/(h2*(h1+h2))  * f_{i+1}
//
// Left boundary (i=0) : forward one-sided, nodes t0, t1, t2:
//   (h1 = t1-t0,  h2 = t2-t0)
//
//   f'(t0) =  -(h1+h2)/(h1*h2) * f0
//             + h2/(h1*(h2-h1)) * f1
//             - h1/(h2*(h2-h1)) * f2
//
// Right boundary (i=N-1) : backward one-sided, nodes t_{N-2}, t_{N-1}, t_N:
//   (h1 = t_{N-1}-t_{N-2},  h2 = t_N-t_{N-2})
//
//   f'(t_N) =  (h2-h1)/(h1*h2)     * f_{N-2}
//             - h2/(h1*(h2-h1))     * f_{N-1}
//             + (2h2-h1)/(h2*(h2-h1)) * f_N
//
// All three formulas are exact for polynomials up to degree 2 and have been
// verified algebraically against f(t) = t^2 for both uniform and non-uniform
// time spacing.

Info << "Differentiating histories (nonuniform 3-point stencils) ..." << nl;

List<vectorField> dLvecHistory(nTimes, vectorField(nFaces, vector::zero));			// d(Lvec)/dt at each stored time, used in the loading term
List<scalarField> dUnHistory  (nTimes, scalarField(nFaces, 0.0));					// d(Un)/dt   at each stored time, used in the thickness term

// Helper lambda: applies one stencil row to every face simultaneously.
// iOut				output time index to fill.
// i0, i1, i2		indices into timeValues for the three stencil nodes.
// c0, c1, c2		corresponding Lagrange-derivative coefficients.
auto applyStencil = [&]
(
	label iOut,
	label i0, label i1, label i2,
	scalar c0, scalar c1, scalar c2
)
{
	for (label fI = 0; fI < nFaces; ++fI)
	{
		dUnHistory[iOut][fI] =														// Scalar derivative: weighted combination of three Un values
			c0 * UnHistory[i0][fI]
		  + c1 * UnHistory[i1][fI]
		  + c2 * UnHistory[i2][fI];

		dLvecHistory[iOut][fI] =													// Vector derivative: same weights applied component-wise to Lvec
			c0 * LvecHistory[i0][fI]
		  + c1 * LvecHistory[i1][fI]
		  + c2 * LvecHistory[i2][fI];
	}
};

// Left boundary stencil (i = 0): forward one-sided using nodes t0, t1, t2.
{
	const scalar h1 = timeValues[1] - timeValues[0];								// First local step
	const scalar h2 = timeValues[2] - timeValues[0];								// Full span

	const scalar c0 = -(h1 + h2) / (h1 * h2);										// Coefficient for f[0]
	const scalar c1 =  h2        / (h1 * (h2 - h1));								// Coefficient for f[1]
	const scalar c2 = -h1        / (h2 * (h2 - h1));								// Coefficient for f[2]

	applyStencil(0,  0, 1, 2,  c0, c1, c2);
}

// Interior centred stencil
for (label i = 1; i < nTimes - 1; ++i)
{
	const scalar h1 = timeValues[i]     - timeValues[i - 1];
	const scalar h2 = timeValues[i + 1] - timeValues[i];

	const scalar c0 = -h2          / (h1 * (h1 + h2));
	const scalar c1 =  (h2 - h1)   / (h1 * h2);
	const scalar c2 =  h1          / (h2 * (h1 + h2));

	applyStencil(i,  i-1, i, i+1,  c0, c1, c2);
}

// Right boundary stencil (i = N)
{
	const label  N  = nTimes - 1;
	const scalar h1 = timeValues[N - 1] - timeValues[N - 2];
	const scalar h2 = timeValues[N]     - timeValues[N - 2];

	const scalar c0 =  (h2 - h1)    / (h1 * h2);
	const scalar c1 = -h2            / (h1 * (h2 - h1));
	const scalar c2 =  (2*h2 - h1)  / (h2 * (h2 - h1));

	applyStencil(N,  N-2, N-1, N,  c0, c1, c2);
}

// * * * * * * * * * * * * *  Interpolation Helpers  * * * * * * * * * * * * //
//
// bracket() finds the interval [lo, hi] that brackets tau using binary search.

auto bracket = [&](scalar tau, label& lo, label& hi)
{
	if (tau <= timeValues[0])          { lo = 0;          hi = 1;          return; }
	if (tau >= timeValues[nTimes - 1]) { lo = nTimes - 2; hi = nTimes - 1; return; }
	lo = 0; hi = nTimes - 1;
	while (hi - lo > 1)
	{
		label mid = (lo + hi) / 2;
		(timeValues[mid] <= tau) ? (lo = mid) : (hi = mid);
	}
};

// * * * * * * *  Precomputed Coefficient Caches (high-order schemes) * * * * //

List<scalarField>   splMomentsScalar (nFaces);
List<List<vector>>  splMomentsVector (nFaces);
List<scalarField>   akimaSlopesScalar(nFaces);
List<List<vector>>  akimaSlopesVector(nFaces);
List<scalarField>   phSlopesScalar   (nFaces);
List<List<vector>>  phSlopesVector   (nFaces);

if (timeInterpolationScheme == "spline"
 ||  timeInterpolationScheme == "akimaCubic"
 ||  timeInterpolationScheme == "piecewiseHermite")
{
	Info << nl << "Precomputing interpolation coefficients (" << timeInterpolationScheme
	     << ") for " << nFaces << " faces ..." << nl;

	for (label fI = 0; fI < nFaces; ++fI)
	{
		if (timeInterpolationScheme == "spline")
		{
			#include "timeInterpolation/spline_setup.H"
		}
		else if (timeInterpolationScheme == "akimaCubic")
		{
			#include "timeInterpolation/akimaCubic_setup.H"
		}
		else // piecewiseHermite
		{
			#include "timeInterpolation/piecewiseHermite_setup.H"
		}
	}

	Info << "    Done." << nl;
}

// * * * * * * * * * * * * *  Generic Interpolation Lambda  * * * * * * * * * //
//
// A single generic (auto-parameter) lambda replaces the former interpScalar /
// interpVector pair.  Using 'const auto& hist' makes the lambda a template,
// which is required for 'if constexpr' inside spline.H / akimaCubic.H /
// piecewiseHermite.H to properly discard the branch whose moment-array type
// does not match the instantiation.  With a concrete (non-template) lambda the
// compiler instantiates both branches and fails on the type mismatch.
//
// Return type is deduced automatically:
//   hist = List<scalarField>  =>  returns scalar
//   hist = List<vectorField>  =>  returns vector

auto interp = [&]
(
    const auto& hist,
    scalar tau,
    label  fI
)
{
    if (timeInterpolationScheme == "nearest")
    {
        #include "timeInterpolation/nearest.H"
    }
    else if (timeInterpolationScheme == "linear")
    {
        #include "timeInterpolation/linear.H"
    }
    else if (timeInterpolationScheme == "hermiteCubic")
    {
        #include "timeInterpolation/hermiteCubic.H"
    }
    else if (timeInterpolationScheme == "spline")
    {
        #include "timeInterpolation/spline.H"
    }
    else if (timeInterpolationScheme == "akimaCubic")
    {
        #include "timeInterpolation/akimaCubic.H"
    }
    else if (timeInterpolationScheme == "piecewiseHermite")
    {
        #include "timeInterpolation/piecewiseHermite.H"
    }
    else
    {
        FatalErrorInFunction
            << "Unknown timeInterpolationScheme '" << timeInterpolationScheme
            << "' in fwhDict.sampling." << nl
            << "Valid options are: "
               "linear, nearest, hermiteCubic, spline, akimaCubic, piecewiseHermite"
            << exit(FatalError);
    }
    // Unreachable; provides a typed return for compiler deduction.
    return hist[0][fI];
};

// * * * * * * * * * * * * *  FW-H Surface Integration  * * * * * * * * * * * //

List<scalarField> pAcoustic(nObs, scalarField(nTimes, 0.0));

Info << nl << "Running FWH integration (Farassat 1A) ..." << nl;

forAll(observers, obsI)
{
	const vector& xObs = observers[obsI];
	Info << "  Observer " << obsI << " at " << xObs << nl;

	for (label tI = 0; tI < nTimes; ++tI)
	{
		const scalar t = timeValues[tI];
		scalar pSum = 0.0;

		for (label fI = 0; fI < nFaces; ++fI)
		{
			const vector r     = xObs - faceCentres[fI];
			const scalar rMag  = mag(r) + VSMALL;
			const vector rHat  = r / rMag;
			const scalar dS    = faceAreas[fI];

			const scalar tau = t - rMag / c0;

			// Use the single generic interp lambda for all four history arrays.
			const vector Lvec_tau  = interp(LvecHistory,  tau, fI);
			const vector dLvec_tau = interp(dLvecHistory, tau, fI);
			const scalar Un_tau    = interp(UnHistory,     tau, fI);
			const scalar dUn_tau   = interp(dUnHistory,    tau, fI);

			const scalar Lr     = Lvec_tau  & rHat;
			const scalar dLr_dt = dLvec_tau & rHat;

			const scalar pL =
				( dLr_dt / (c0 * rMag)
				+ Lr     / (rMag * rMag) )
				/ fourPi * dS;

			const scalar pT =
				rho * dUn_tau / rMag / fourPi * dS;

			pSum += pL + pT;
		}

		pAcoustic[obsI][tI] = pSum;
	}

	Info << "    Done.  p' range: ["
		 << min(pAcoustic[obsI]) << ", "
		 << max(pAcoustic[obsI]) << "] Pa" << nl;
}

Info << nl << "FWH integration complete." << nl;
