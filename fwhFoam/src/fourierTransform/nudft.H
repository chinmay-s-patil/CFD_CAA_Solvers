// nudft.H
// Non-Uniform Discrete Fourier Transform (NUDFT), O(N²).
//
// Each frequency bin k represents physical frequency f_k = k * df = k / T,
// where T = t_{N-1} - t_0 is the total signal duration.  The phase angle
// for sample n is computed from the actual stored time t_n:
//
//   X(k) = sum_{n=0}^{N-1} x(t_n) * exp(-i * 2*pi * f_k * t_n)
//
// This is correct for any time sampling pattern — uniform, non-uniform,
// adaptive, or stride-skipped.  For uniformly-spaced data this reduces
// algebraically to the standard DFT (dft.H) when df = 1/(N*dt).
//
// Normalisation convention (same as dft.H and fft.H for easy comparison):
//   |P(k)| = 2 * |X(k)| / N   for interior bins  (k = 1 .. Nhalf-2)
//   |P(0)| =     |X(0)| / N   (DC)
//   |P(Nhalf-1)| = |X(Nhalf-1)| / N   (Nyquist)
//
// Expects in scope (provided by writeOutputs.H before the dispatch):
//   sig         const scalarField&   acoustic pressure time history (length N)
//   N           label                number of samples
//   Nhalf       label                one-sided bin count  (N/2 + 1)
//   T           scalar               total duration  t_{N-1} - t_0  [s]
//   timeValues  DynamicList<scalar>  actual sample times [s]
//   ftMag       scalarField          output: one-sided magnitude spectrum [Pa]

{
    const scalar twoPiOverT =
        2.0 * constant::mathematical::pi / max(T, VSMALL);

    for (label k = 0; k < Nhalf; k++)
    {
        scalar re = 0.0, im = 0.0;
        const scalar omega_k = scalar(k) * twoPiOverT;  // angular frequency of bin k

        for (label n = 0; n < N; n++)
        {
            // Phase angle uses the actual sample time — correct for any spacing
            const scalar phi = -omega_k * timeValues[n];
            re += sig[n] * Foam::cos(phi);
            im += sig[n] * Foam::sin(phi);
        }

        // Two-sided to one-sided: multiply interior bins by 2
        const scalar scale = (k == 0 || k == Nhalf - 1) ? 1.0 : 2.0;
        ftMag[k] = scale * Foam::sqrt(re*re + im*im) / scalar(N);
    }
}
