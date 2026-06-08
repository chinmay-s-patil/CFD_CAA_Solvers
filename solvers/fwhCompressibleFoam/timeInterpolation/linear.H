
label lo, hi;
bracket(tau, lo, hi);
const scalar alpha =
    (tau - timeValues[lo]) / (timeValues[hi] - timeValues[lo] + VSMALL);
return (1.0 - alpha)*hist[lo][fI] + alpha*hist[hi][fI];
