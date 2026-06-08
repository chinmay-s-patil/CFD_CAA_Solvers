// czt.H
// Chirp Z-Transform (CZT) based zoom FFT.
//
// Evaluates the Z-transform on a spiral arc in the z-plane that maps to a
// uniform frequency grid over [fStart, fEnd].  For zoom factors > ~4x this
// gives far better frequency resolution than a full-length FFT of the same
// signal at the cost of O((N+M) log(N+M)) via Bluestein's identity.
//
// Algorithm (Bluestein / Rabiner 1969):
//   Given N input samples x[n] sampled at rate fs, evaluate M frequency
//   lines uniformly spaced over [fStart, fEnd]:
//
//   Step 1.  Define:
//     W  = exp(-j 2pi df / fs)        df = (fEnd - fStart) / (M - 1)
//     A  = exp( j 2pi fStart / fs)    starting phasor
//
//   Step 2.  Apply Bluestein's identity  nk = -(k-n)^2/2 + n^2/2 + k^2/2:
//     X(k) = W^{k^2/2} * IFFT[ FFT[x[n] * A^{-n} * W^{n^2/2}]
//                               * FFT[W^{-m^2/2}] ]_k
//
//   The three inner transforms are computed on a zero-padded length
//   L = next power-of-two >= N + M - 1, making the overall complexity
//   O((N+M) log(N+M)).
//
// Normalisation:
//   Same 1/N amplitude convention as dft.H / fft.H / fftw.H so that the
//   zoom spectrum amplitude is directly comparable to the full-range FFT.
//   Interior bins carry the factor-of-2 one-sided correction.
//
// fs estimation:
//   Because time steps can be non-uniform the sample rate is estimated as
//   fs = (N-1) / (t_{N-1} - t_0).  A warning is emitted when the
//   coefficient of variation of dt exceeds 1% (matching the existing FT
//   uniformity check in writeOutputs.H).
//
// Expects in scope (provided by writeOutputs.H zoom block):
//   sig          const scalarField&        signal, length N
//   N            label
//   timeValues   DynamicList<scalar>       actual sample times [s]
//   zoomDict     const dictionary&         zoom{} sub-dict
//   zoomFTMag    scalarField               output: M-point magnitude [Pa]
//   zoomFreqs    scalarField               output: M-point frequency axis [Hz]
//   zoomM        label                     number of zoom lines (= nLines)
//   zoomFs       scalar                    estimated sample rate [Hz]
//   win          scalarField (nfft=N)      window coefficients (pre-built by caller)

{
    // -----------------------------------------------------------------------
    //  Frequency axis setup
    // -----------------------------------------------------------------------

    const scalar fStart = zoomDict.get<scalar>("fStart");
    const scalar fEnd   = zoomDict.get<scalar>("fEnd");
    const label  M      = zoomM;           // number of output lines

    const scalar df     = (fEnd - fStart) / max(scalar(M - 1), VSMALL);
    const scalar fs     = zoomFs;

    // Populate output frequency axis
    for (label k = 0; k < M; k++)
        zoomFreqs[k] = fStart + scalar(k) * df;

    // -----------------------------------------------------------------------
    //  Pre-process: apply window, mean removal, detrend
    //  (caller already built win[]; here we copy + window the signal)
    // -----------------------------------------------------------------------

    // Windowed, optionally mean-removed / detrended input
    List<scalar> xw(N, 0.0);
    {
        // mean removal
        scalar mu = 0.0;
        if (zoomDict.getOrDefault<Switch>("removeMean", true))
        {
            for (label n = 0; n < N; n++) mu += sig[n];
            mu /= scalar(N);
        }

        // linear detrend
        scalar aSlope = 0.0, bInt = 0.0;
        if (zoomDict.getOrDefault<Switch>("detrend", false))
        {
            const scalar t0 = timeValues[0];
            scalar sum_t  = 0, sum_y  = 0;
            scalar sum_tt = 0, sum_ty = 0;
            for (label n = 0; n < N; n++)
            {
                const scalar tn = timeValues[n] - t0;
                const scalar yn = sig[n] - mu;
                sum_t  += tn;  sum_y  += yn;
                sum_tt += tn*tn; sum_ty += tn*yn;
            }
            const scalar denom = N*sum_tt - sum_t*sum_t + VSMALL;
            aSlope = (N*sum_ty - sum_t*sum_y) / denom;
            bInt   = (sum_y - aSlope*sum_t) / N;
        }

        for (label n = 0; n < N; n++)
        {
            const scalar tn = timeValues[n] - timeValues[0];
            xw[n] = (sig[n] - mu - aSlope*tn - bInt) * win[n];
        }
    }

    // -----------------------------------------------------------------------
    //  Bluestein CZT
    // -----------------------------------------------------------------------

    // Angular increments (radians per sample)
    const scalar phi_W  = -2.0 * constant::mathematical::pi * df  / fs;   // W phase step
    const scalar phi_A  =  2.0 * constant::mathematical::pi * fStart / fs; // A start phase

    // Zero-padded FFT length: L >= N + M - 1, next power of two
    label L = 1;
    while (L < N + M - 1) L <<= 1;

    // Helper: in-place iterative Cooley-Tukey FFT (same logic as fft.H)
    // Works on interleaved real/imag List<scalar> of length 2*L.
    auto fftInPlace = [&](List<scalar>& w, bool inverse)
    {
        // Bit-reverse permutation
        {
            label j = 0;
            for (label i = 1; i < L; i++)
            {
                label bit = L >> 1;
                for (; j & bit; bit >>= 1) j ^= bit;
                j ^= bit;
                if (i < j)
                {
                    std::swap(w[2*i],   w[2*j]);
                    std::swap(w[2*i+1], w[2*j+1]);
                }
            }
        }
        // Butterfly passes
        const scalar sign = inverse ? 1.0 : -1.0;
        for (label len = 2; len <= L; len <<= 1)
        {
            const scalar theta = sign * 2.0*constant::mathematical::pi / scalar(len);
            const scalar wRe0  = Foam::cos(theta);
            const scalar wIm0  = Foam::sin(theta);
            for (label i = 0; i < L; i += len)
            {
                scalar wRe = 1.0, wIm = 0.0;
                for (label k = 0; k < len/2; k++)
                {
                    const scalar uRe = w[2*(i+k)];
                    const scalar uIm = w[2*(i+k)+1];
                    const scalar vRe = w[2*(i+k+len/2)]*wRe - w[2*(i+k+len/2)+1]*wIm;
                    const scalar vIm = w[2*(i+k+len/2)]*wIm + w[2*(i+k+len/2)+1]*wRe;
                    w[2*(i+k)]         = uRe + vRe;
                    w[2*(i+k)+1]       = uIm + vIm;
                    w[2*(i+k+len/2)]   = uRe - vRe;
                    w[2*(i+k+len/2)+1] = uIm - vIm;
                    const scalar wReN = wRe*wRe0 - wIm*wIm0;
                    wIm = wRe*wIm0 + wIm*wRe0;
                    wRe = wReN;
                }
            }
        }
        if (inverse)
            for (label i = 0; i < 2*L; i++) w[i] /= scalar(L);
    };

    // --- Build yn = x[n] * A^{-n} * W^{n^2/2}  (length N, zero-padded to L)
    List<scalar> yn(2*L, 0.0);
    for (label n = 0; n < N; n++)
    {
        const scalar phiA  = -phi_A * scalar(n);          // A^{-n}
        const scalar phiW  =  phi_W * scalar(n)*scalar(n) * 0.5;  // W^{n^2/2}
        const scalar phi   = phiA + phiW;
        yn[2*n]   = xw[n] * Foam::cos(phi);
        yn[2*n+1] = xw[n] * Foam::sin(phi);
    }

    // --- Build hm = W^{-m^2/2}  for m = 0..M-1, and m = L-(N-1)..L-1
    //     (the "chirp filter" in circular convolution notation)
    List<scalar> hm(2*L, 0.0);
    // Positive indices 0..M-1
    for (label m = 0; m < M; m++)
    {
        const scalar phi = -phi_W * scalar(m)*scalar(m) * 0.5;
        hm[2*m]   = Foam::cos(phi);
        hm[2*m+1] = Foam::sin(phi);
    }
    // Wrap-around indices L-(N-1)..L-1
    for (label n = 1; n < N; n++)
    {
        const label  idx = L - n;
        const scalar phi = -phi_W * scalar(n)*scalar(n) * 0.5;
        hm[2*idx]   = Foam::cos(phi);
        hm[2*idx+1] = Foam::sin(phi);
    }

    // --- FFT both sequences
    fftInPlace(yn, false);
    fftInPlace(hm, false);

    // --- Complex multiply in frequency domain
    List<scalar> conv(2*L, 0.0);
    for (label i = 0; i < L; i++)
    {
        conv[2*i]   = yn[2*i]*hm[2*i]   - yn[2*i+1]*hm[2*i+1];
        conv[2*i+1] = yn[2*i]*hm[2*i+1] + yn[2*i+1]*hm[2*i];
    }

    // --- IFFT
    fftInPlace(conv, true);

    // --- Extract M output bins: multiply by W^{k^2/2}
    for (label k = 0; k < M; k++)
    {
        const scalar phi  = phi_W * scalar(k)*scalar(k) * 0.5;
        const scalar cosP = Foam::cos(phi);
        const scalar sinP = Foam::sin(phi);
        const scalar re   = conv[2*k]*cosP - conv[2*k+1]*sinP;
        const scalar im   = conv[2*k]*sinP + conv[2*k+1]*cosP;

        // Amplitude: 1/N normalisation (one-sided factor already applied below)
        zoomFTMag[k] = Foam::sqrt(re*re + im*im) / scalar(N);
    }

    // One-sided correction: all zoom bins are interior (not DC or Nyquist of
    // the full signal), so multiply by 2 for amplitude consistency with fft.H
    for (label k = 0; k < M; k++)
        zoomFTMag[k] *= 2.0;
}