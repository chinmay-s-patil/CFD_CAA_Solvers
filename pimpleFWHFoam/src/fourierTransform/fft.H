// fft.H
// Radix-2 Cooley-Tukey iterative FFT.
//
// A self-contained O(N log N) real-to-complex FFT implemented entirely in
// standard C++17 / OpenFOAM primitives : no external library required.
//
// IMPORTANT: This algorithm assumes UNIFORM time spacing.  For non-uniform
// time series use fourierTransformType NUDFT (nudft.H) instead.  writeOutputs.H
// will emit a WarningInFunction before reaching this code if non-uniform
// spacing is detected.
//
// The algorithm:
//   1. Zero-pad the input signal to the next power of two, M = 2^ceil(log2(N)).
//   2. Bit-reverse permutation (in-place).
//   3. Butterfly passes from stage 1 up to log2(M).
//   4. Extract the one-sided (real-input) spectrum: bins 0 .. M/2, scaled
//      with the same 1/N amplitude convention as dft.H so that the two
//      results are directly comparable.
//
// Normalisation convention (same as dft.H):
//   |P(k)| = 2 * |X(k)| / N   for interior bins  (k = 1 .. Nhalf-2)
//   |P(0)| =     |X(0)| / N   (DC)
//   |P(Nhalf-1)| = |X(Nhalf-1)| / N   (Nyquist, if M is even)
//
// where X(k) is the raw DFT sum and N is the original (unpadded) length.
// The factor of 2 accounts for the energy folded from the negative-frequency
// mirror.
//
// Expects in scope (provided by writeOutputs.H before the dispatch):
//   sig    const scalarField&   acoustic pressure time history (length N)
//   N      label                number of samples
//   Nhalf  label                one-sided bin count  (N/2 + 1)
//   T      scalar               total duration t_{N-1} - t_0  [s]  (unused here;
//                               the frequency axis is managed by writeOutputs.H)
//   ftMag  scalarField          output: one-sided magnitude spectrum [Pa]
//
// Note: ftMag has Nhalf entries sized for the original N.  The FFT is
// computed on a zero-padded copy of length M >= N, but only Nhalf bins
// are written, preserving the same frequency axis as the DFT.

{
    // * * * * * * * * * * 1. Next power-of-two length * * * * * * * * * * //

    label M = 1;
    while (M < N) M <<= 1;      // M = 2^ceil(log2(N))

    // * * * 2. Copy input into complex work array, zero-pad to length M * //

    // We store interleaved real/imag as a List<scalar> of length 2*M.
    // re[k] = work[2*k], im[k] = work[2*k+1].

    List<scalar> work(2 * M, 0.0);
    for (label n = 0; n < N; n++)
    {
        work[2*n] = sig[n];     // real part from signal
        // work[2*n+1] = 0 (imaginary already zero)
    }

    // * * * * * * * * * * 3. Bit-reverse permutation  * * * * * * * * * * //

    {
        label j = 0;
        for (label i = 1; i < M; i++)
        {
            label bit = M >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j)
            {
                std::swap(work[2*i],   work[2*j]);
                std::swap(work[2*i+1], work[2*j+1]);
            }
        }
    }

    // * * * * * * * * * * 4. Cooley-Tukey butterfly passes  * * * * * * * //

    const scalar twoPi = 2.0 * constant::mathematical::pi;

    for (label len = 2; len <= M; len <<= 1)
    {
        const scalar theta = -twoPi / scalar(len);   // twiddle angle per stage
        const scalar wRe0  = Foam::cos(theta);
        const scalar wIm0  = Foam::sin(theta);

        for (label i = 0; i < M; i += len)
        {
            scalar wRe = 1.0, wIm = 0.0;            // twiddle factor W^0

            for (label k = 0; k < len / 2; k++)
            {
                // Even element
                const scalar uRe = work[2*(i + k)    ];
                const scalar uIm = work[2*(i + k) + 1];

                // Odd element * twiddle
                const scalar vRe = work[2*(i + k + len/2)    ] * wRe
                                 - work[2*(i + k + len/2) + 1] * wIm;
                const scalar vIm = work[2*(i + k + len/2)    ] * wIm
                                 + work[2*(i + k + len/2) + 1] * wRe;

                // Butterfly
                work[2*(i + k)          ] = uRe + vRe;
                work[2*(i + k) + 1      ] = uIm + vIm;
                work[2*(i + k + len/2)  ] = uRe - vRe;
                work[2*(i + k + len/2)+1] = uIm - vIm;

                // Advance twiddle: W^{k+1} = W^k * W^1
                const scalar wReNew = wRe * wRe0 - wIm * wIm0;
                wIm = wRe * wIm0 + wIm * wRe0;
                wRe = wReNew;
            }
        }
    }

    // * * * * * * * * * * 5. Extract one-sided magnitude spectrum * * * * //
    //
    // Nhalf entries, same frequency axis as the DFT (bin spacing = 1/(N*dt)).
    // We normalise by N (not M) to be amplitude-consistent with dft.H.

    for (label k = 0; k < Nhalf; k++)
    {
        const scalar re = work[2*k    ];
        const scalar im = work[2*k + 1];

        // Two-sided → one-sided: factor 2 for interior bins
        const scalar scale = (k == 0 || k == Nhalf - 1) ? 1.0 : 2.0;
        ftMag[k] = scale * Foam::sqrt(re*re + im*im) / scalar(N);
    }
}
