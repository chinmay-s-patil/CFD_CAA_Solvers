// spline.H
// Natural cubic spline interpolation
//
// The second-derivative moments M[i] must have been precomputed for this
// face by spline_setup.H before the integration loop.  This snippet just
// evaluates the cubic formula on the bracketed interval.
//
// Expects in scope (provided by fwhIntegrate.H before the dispatch):
//   hist                scalar or vector history   [tI][fI]
//   tau                 scalar   retarded time
//   fI                  label    face index
//   timeValues          DynamicList<scalar>
//   nTimes              label
//   bracket()           lambda
//   splMomentsScalar    List<scalarField>    precomputed Un moments
//   splMomentsVector    List<List<vector>>   precomputed Lvec moments

label lo, hi;
bracket(tau, lo, hi);

const scalar hSpl  = timeValues[hi] - timeValues[lo] + VSMALL;
const scalar a     = (timeValues[hi] - tau) / hSpl;
const scalar b     = (tau - timeValues[lo]) / hSpl;

// Select the correct moment array depending on whether hist is scalar or vector.
// The typedef valType is used to dispatch at compile time.
typedef std::remove_const_t<std::remove_reference_t<
    decltype(hist[0][fI])>> valType;

if constexpr (std::is_same_v<valType, scalar>)
{
    const scalarField& M = splMomentsScalar[fI];
    return   a * hist[lo][fI]
           + b * hist[hi][fI]
           + ((a*a*a - a) * M[lo] + (b*b*b - b) * M[hi]) * (hSpl * hSpl / 6.0);
}
else
{
    const List<vector>& M = splMomentsVector[fI];
    return   a * hist[lo][fI]
           + b * hist[hi][fI]
           + ((a*a*a - a) * M[lo] + (b*b*b - b) * M[hi]) * (hSpl * hSpl / 6.0);
}
