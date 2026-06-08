// heterodyne.H
// Zoom FFT via heterodyne frequency shift + FIR decimation + FFTW.
//
// Pipeline:
//   1. Frequency-shift the signal to baseband by multiplying by
//      exp(-j 2pi f_centre t_n).  For a real input signal the baseband
//      signal is complex; we keep both real and imaginary parts.
//   2. Apply a windowed-sinc FIR lowpass filter with cutoff BW/2,
//      where BW = fEnd - fStart.  This suppresses aliasing before decimation.
//   3. Decimate by factor D (auto or user-specified).
//   4. Run FFTW on the decimated complex baseband signal.
//   5. Shift the output frequency axis back to [fStart, fEnd].
//
// Decimation factor D:
//   auto: D = floor(fs / (oversampling * BW))   oversampling = 2.56 (default)
//   The oversampling factor ensures the FIR transition band fits between
//   BW/2 and fs_decimated/2 without aliasing.  2.56 is the standard choice
//   (matches Agilent / NI spectrum analyser convention).
//
// FIR design (windowed sinc):
//   h[n] = sinc(2 * fc * (n - M/2)) * win[n]
//   fc   = cutoff normalised to decimated fs (0.5 * BW / (fs/D))
//   win  = user-selected window (Hamming default, same options as FFT windows)
//   M    = lpfOrder (number of taps; must be even for linear phase)
//
// Normalisation:
//   Same 1/N amplitude convention as the other FT implementations so the
//   zoom spectrum magnitude is directly comparable.
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
//   win          scalarField (nfft=N)      signal window coefficients (pre-built)

{
    const scalar fStart = zoomDict.get<scalar>("fStart");
    const scalar fEnd   = zoomDict.get<scalar>("fEnd");
    const label  M      = zoomM;
    const scalar fs     = zoomFs;
    const scalar BW     = fEnd - fStart;
    const scalar fCentre = 0.5 * (fStart + fEnd);

    if (BW < VSMALL)
        FatalErrorInFunction
            << "zoom.fEnd must be greater than zoom.fStart."
            << exit(FatalError);

    // -----------------------------------------------------------------------
    //  Decimation factor
    // -----------------------------------------------------------------------

    const scalar oversampling = 2.56;
    label D = 0;
    {
        const word decStr =
            zoomDict.getOrDefault<word>("decimationFactor", "auto");
        if (decStr == "auto")
            D = max(label(1), label(Foam::floor(fs / (oversampling * BW))));
        else
            D = max(label(1), label(Foam::readScalar(IStringStream(decStr)())));
    }

    const scalar fsD = fs / scalar(D);    // decimated sample rate

    // -----------------------------------------------------------------------
    //  FIR lowpass filter design (windowed-sinc)
    // -----------------------------------------------------------------------

    const label lpfOrder = zoomDict.getOrDefault<label>("lpfOrder", 64);
    const label lpfTaps  = lpfOrder + 1;   // ensure odd tap count for symmetry
    const label lpfHalf  = lpfTaps / 2;

    // FIR window selection
    const word lpfWinType =
        zoomDict.getOrDefault<word>("lpfWindow", "Hamming");

    scalarField lpfWin(lpfTaps, 1.0);
    {
        const scalar pi = constant::mathematical::pi;
        const label  P  = lpfTaps - 1;

        if (lpfWinType == "Hamming")
            for (label n = 0; n < lpfTaps; n++)
                lpfWin[n] = 0.54 - 0.46*Foam::cos(2.0*pi*scalar(n)/scalar(P));
        else if (lpfWinType == "Hann" || lpfWinType == "Hanning")
            for (label n = 0; n < lpfTaps; n++)
                lpfWin[n] = 0.5*(1.0 - Foam::cos(2.0*pi*scalar(n)/scalar(P)));
        else if (lpfWinType == "Blackman")
            for (label n = 0; n < lpfTaps; n++)
            {
                const scalar phi = 2.0*pi*scalar(n)/scalar(P);
                lpfWin[n] = 0.42 - 0.5*Foam::cos(phi) + 0.08*Foam::cos(2.0*phi);
            }
        else if (lpfWinType == "rectangular")
        {
            // already 1.0
        }
        else
        {
            WarningInFunction
                << "Unknown zoom.lpfWindow '" << lpfWinType
                << "'. Using Hamming." << nl;
            for (label n = 0; n < lpfTaps; n++)
                lpfWin[n] = 0.54 - 0.46*Foam::cos(2.0*pi*scalar(n)/
                             scalar(P));
        }
    }

    // Normalised cutoff: fc = (BW/2) / fsD  — slightly less than 0.5 to
    // leave a guard band matched to the oversampling factor.
    const scalar fc = (BW / 2.0) / fsD;

    scalarField lpfH(lpfTaps, 0.0);
    {
        const scalar pi = constant::mathematical::pi;
        for (label n = 0; n < lpfTaps; n++)
        {
            const scalar m = scalar(n) - scalar(lpfHalf);
            if (Foam::mag(m) < VSMALL)
                lpfH[n] = 2.0 * fc * lpfWin[n];
            else
                lpfH[n] = Foam::sin(2.0*pi*fc*m) / (pi*m) * lpfWin[n];
        }
        // Normalise DC gain to 1
        scalar dcGain = 0.0;
        for (label n = 0; n < lpfTaps; n++) dcGain += lpfH[n];
        dcGain = max(dcGain, VSMALL);
        for (label n = 0; n < lpfTaps; n++) lpfH[n] /= dcGain;
    }

    // -----------------------------------------------------------------------
    //  Pre-process input signal (mean removal + detrend from zoom dict)
    // -----------------------------------------------------------------------

    List<scalar> xproc(N, 0.0);
    {
        scalar mu = 0.0;
        if (zoomDict.getOrDefault<Switch>("removeMean", true))
        {
            for (label n = 0; n < N; n++) mu += sig[n];
            mu /= scalar(N);
        }

        scalar aSlope = 0.0, bInt = 0.0;
        if (zoomDict.getOrDefault<Switch>("detrend", false))
        {
            const scalar t0 = timeValues[0];
            scalar sum_t = 0, sum_y = 0, sum_tt = 0, sum_ty = 0;
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
            xproc[n] = (sig[n] - mu - aSlope*tn - bInt) * win[n];
        }
    }

    // -----------------------------------------------------------------------
    //  Step 1: Heterodyne — shift to baseband
    //  xB[n] = xproc[n] * exp(-j 2pi fCentre t_n)
    //  Store as interleaved real/imag.
    // -----------------------------------------------------------------------

    List<scalar> xB(2*N, 0.0);
    for (label n = 0; n < N; n++)
    {
        const scalar phi = -2.0*constant::mathematical::pi*fCentre*timeValues[n];
        xB[2*n]   = xproc[n] * Foam::cos(phi);
        xB[2*n+1] = xproc[n] * Foam::sin(phi);
    }

    // -----------------------------------------------------------------------
    //  Step 2: FIR lowpass filter (linear convolution, direct form)
    //  Operate on complex baseband signal.
    // -----------------------------------------------------------------------

    List<scalar> xFilt(2*N, 0.0);
    for (label n = 0; n < N; n++)
    {
        scalar re = 0.0, im = 0.0;
        for (label k = 0; k < lpfTaps; k++)
        {
            const label idx = n - k;
            if (idx < 0 || idx >= N) continue;
            re += lpfH[k] * xB[2*idx];
            im += lpfH[k] * xB[2*idx+1];
        }
        xFilt[2*n]   = re;
        xFilt[2*n+1] = im;
    }

    // -----------------------------------------------------------------------
    //  Step 3: Decimate
    // -----------------------------------------------------------------------

    const label Nd = (N + D - 1) / D;    // decimated length (ceiling)
    List<scalar> xDec(2*Nd, 0.0);
    for (label n = 0; n < Nd; n++)
    {
        const label src = n * D;
        if (src < N)
        {
            xDec[2*n]   = xFilt[2*src];
            xDec[2*n+1] = xFilt[2*src+1];
        }
    }

    // -----------------------------------------------------------------------
    //  Step 4: FFTW on complex decimated signal
    //  Use fftw_plan_dft_1d (complex-to-complex) so both I and Q are used.
    // -----------------------------------------------------------------------

    // Zero-pad to next power of two for FFTW efficiency
    label Lp = 1;
    while (Lp < Nd) Lp <<= 1;

    fftw_complex* inC  = fftw_alloc_complex(static_cast<size_t>(Lp));
    fftw_complex* outC = fftw_alloc_complex(static_cast<size_t>(Lp));

    for (label n = 0; n < Lp; n++)
    {
        inC[n][0] = (n < Nd) ? xDec[2*n]   : 0.0;
        inC[n][1] = (n < Nd) ? xDec[2*n+1] : 0.0;
    }

    fftw_plan planC = fftw_plan_dft_1d
    (
        static_cast<int>(Lp),
        inC, outC,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );
    fftw_execute(planC);
    fftw_destroy_plan(planC);

    // -----------------------------------------------------------------------
    //  Step 5: Extract M lines from the complex spectrum, map to [fStart,fEnd]
    //
    //  The complex FFT of a baseband signal gives a two-sided spectrum
    //  centred at 0 Hz (= fCentre in the original domain).
    //  Positive frequencies [0, fsD/2]  map to [fCentre, fCentre + fsD/2]
    //  Negative frequencies [-fsD/2, 0] map to [fCentre - fsD/2, fCentre]
    //
    //  We linearly interpolate M lines uniformly spaced over [fStart, fEnd].
    // -----------------------------------------------------------------------

    const scalar dfZoom = (fEnd - fStart) / max(scalar(M - 1), VSMALL);
    const scalar dfDec  = fsD / scalar(Lp);     // bin spacing of decimated FFT

    for (label k = 0; k < M; k++)
    {
        const scalar fk      = fStart + scalar(k) * dfZoom;
        zoomFreqs[k]         = fk;

        // Frequency relative to centre, mapped to decimated FFT bin index.
        // Negative relative frequencies wrap to the upper half of the DFT.
        const scalar fRel    = fk - fCentre;
        scalar binF          = fRel / dfDec;
        if (binF < 0) binF  += scalar(Lp);

        // Linear interpolation between adjacent bins
        const label  binLo   = label(Foam::floor(binF)) % Lp;
        const label  binHi   = (binLo + 1) % Lp;
        const scalar alpha   = binF - Foam::floor(binF);

        const scalar re = (1.0-alpha)*outC[binLo][0] + alpha*outC[binHi][0];
        const scalar im = (1.0-alpha)*outC[binLo][1] + alpha*outC[binHi][1];

        // Normalise by N (original signal length), factor 2 for one-sided
        zoomFTMag[k] = 2.0 * Foam::sqrt(re*re + im*im) / scalar(N);
    }

    fftw_free(inC);
    fftw_free(outC);
}