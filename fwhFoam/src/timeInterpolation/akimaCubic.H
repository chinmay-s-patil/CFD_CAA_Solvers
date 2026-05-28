// akimaCubic.H
// Akima piecewise cubic interpolation — O(1) evaluation step.
//
// The Akima knot slopes m[i] must have been precomputed for this face by
// akimaCubic_setup.H before the integration loop.  This snippet evaluates
// the cubic Hermite formula on the bracketed interval.
//
// Expects in scope (provided by fwhIntegrate.H before the dispatch):
//   hist                 scalar or vector history   [tI][fI]
//   tau                  scalar   retarded time
//   fI                   label    face index
//   timeValues           DynamicList<scalar>
//   nTimes               label
//   bracket()            lambda
//   akimaSlopesScalar    List<scalarField>    precomputed Un slopes
//   akimaSlopesVector    List<List<vector>>   precomputed Lvec slopes

label lo, hi;
bracket(tau, lo, hi);

const scalar hAkima = timeValues[hi] - timeValues[lo] + VSMALL;
const scalar u      = (tau - timeValues[lo]) / hAkima;
const scalar u2     = u * u;
const scalar u3     = u2 * u;

// Cubic Hermite basis polynomials H00, H10, H01, H11
typedef std::remove_const_t<std::remove_reference_t<
    decltype(hist[0][fI])>> valType;

if constexpr (std::is_same_v<valType, scalar>)
{
    const scalarField& m = akimaSlopesScalar[fI];
    return   ( 2.0*u3 - 3.0*u2 + 1.0) * hist[lo][fI]
           + (     u3 - 2.0*u2 + u  ) * (hAkima * m[lo])
           + (-2.0*u3 + 3.0*u2      ) * hist[hi][fI]
           + (     u3 -     u2      ) * (hAkima * m[hi]);
}
else
{
    const List<vector>& m = akimaSlopesVector[fI];
    return   ( 2.0*u3 - 3.0*u2 + 1.0) * hist[lo][fI]
           + (     u3 - 2.0*u2 + u  ) * (hAkima * m[lo])
           + (-2.0*u3 + 3.0*u2      ) * hist[hi][fI]
           + (     u3 -     u2      ) * (hAkima * m[hi]);
}
