
label lo, hi;
bracket(tau, lo, hi);
const scalar alpha =
    (tau - timeValues[lo]) / (timeValues[hi] - timeValues[lo] + VSMALL);
if(alpha < 0.5)
{
    return hist[lo][fI];
}
else
{
    return hist[hi][fI];
}
