// hermiteCubic.H
// Non-uniform Catmull-Rom cubic Hermite interpolation.
//
// Tangents at the two endpoints of the bracketed interval [lo, hi] are
// estimated with the centred finite-difference rule applied to the nearest
// available neighbours.  At the boundary of the stored time window the
// one-sided (clamped) neighbour is used, which degrades gracefully to a
// linear extrapolation when both neighbours coincide with an endpoint.
//
// Expects in scope (provided by fwhIntegrate.H before the dispatch):
//   hist        const List<scalarField>& or const List<vectorField>&
//   tau         scalar   retarded time to evaluate
//   fI          label    face index
//   timeValues  DynamicList<scalar>
//   nTimes      label
//   bracket()   lambda

label lo, hi;
bracket(tau, lo, hi);

const scalar t_lo = timeValues[lo];
const scalar t_hi = timeValues[hi];
const scalar h    = t_hi - t_lo + VSMALL;   // interval width
const scalar u    = (tau - t_lo) / h;        // normalised position in [0,1]

// Clamped neighbour indices
const label im1 = max(lo - 1, label(0));
const label ip1 = min(hi + 1, nTimes - 1);

// Non-uniform Catmull-Rom tangent denominators
//   m_lo ~ df/dt at lo,  using the span from im1 to hi
//   m_hi ~ df/dt at hi,  using the span from lo  to ip1
const scalar dt_m = timeValues[hi]  - timeValues[im1] + VSMALL;
const scalar dt_p = timeValues[ip1] - timeValues[lo]  + VSMALL;

auto m_lo = (hist[hi][fI]  - hist[im1][fI]) / dt_m;
auto m_hi = (hist[ip1][fI] - hist[lo][fI])  / dt_p;

// Cubic Hermite basis polynomials H00, H10, H01, H11
const scalar u2 = u * u;
const scalar u3 = u2 * u;

return (  2.0*u3 - 3.0*u2 + 1.0) * hist[lo][fI]
     + (     u3 - 2.0*u2 + u   ) * (h * m_lo)
     + ( -2.0*u3 + 3.0*u2      ) * hist[hi][fI]
     + (     u3 -     u2       ) * (h * m_hi);
