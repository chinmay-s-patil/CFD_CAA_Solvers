// dft.H
// Direct O(N^2) Discrete Fourier Transform (uniform sampling assumed)
//
// Frequency bin k represents f_k = k * df = k / (N * dt), where dt is the
// uniform time step.  The phase angle uses the sample index n, which is
// algebraically equivalent to using actual times only when spacing is uniform.
//
// For non-uniform time series (adaptive time-stepping, stride-skipped output)
// use nudft.H instead, which uses actual timeValues[n] in the phase angle.
//
// Normalisation convention:
//   |P(k)| = 2 * |X(k)| / N   for interior bins  (k = 1 .. Nhalf-2)
//   |P(0)| =     |X(0)| / N   (DC)
//   |P(Nhalf-1)| = |X(Nhalf-1)| / N   (Nyquist)
//
// Expects in scope (provided by writeOutputs.H before the dispatch):
//   sig    const scalarField&   acoustic pressure time history (length N)
//   N      label                number of samples
//   Nhalf  label                one-sided bin count  (N/2 + 1)
//   ftMag  scalarField          output: one-sided magnitude spectrum [Pa]

for (label k = 0; k < Nhalf; k++)
{
    scalar re = 0.0, im = 0.0;
    const scalar phi0 =
        -2.0 * constant::mathematical::pi * scalar(k) / scalar(N);

    for (label n = 0; n < N; n++)
    {
        re += sig[n] * Foam::cos(phi0 * scalar(n));
        im += sig[n] * Foam::sin(phi0 * scalar(n));
    }

    // Two-sided to one-sided: multiply interior bins by 2
    const scalar scale = (k == 0 || k == Nhalf - 1) ? 1.0 : 2.0;
    ftMag[k] = scale * Foam::sqrt(re*re + im*im) / scalar(N);
}
