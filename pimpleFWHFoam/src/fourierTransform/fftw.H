// fftw.H
// Windowed FFT via the FFTW3 library, with optional mean removal and
// linear detrend.  Supports Welch-style segment averaging when nfft < N.
//
// Pipeline applied to each segment before the DFT:
//   1. Extract segment of length nfft from the signal.
//   2. (Optional) subtract the segment mean                  [removeMean yes]
//   3. (Optional) subtract the linear least-squares trend    [detrend    yes]
//   4. Apply window function read from windows/<windowType>.H [windowType ...]
//   5. Call fftw_plan_dft_r2c_1d : the real-input, half-complex transform.
//   6. Accumulate power spectrum; convert to amplitude at the end.
//
// IMPORTANT: This algorithm assumes UNIFORM time spacing within each segment.
// For non-uniform time series use fourierTransformType NUDFT (nudft.H) instead.
// writeOutputs.H will emit a WarningInFunction before reaching this code if
// non-uniform spacing is detected.
//
// fwhDict.fourierTransform entries consumed here:
//   windowType    word     Hann | Hamming | Blackman | FlatTop |
//                          Bartlett | Tukey | rectangular
//   removeMean    Switch   yes / no   (default: yes)
//   detrend       Switch   yes / no   (default: no)
//   nfft          label | "auto"       segment length (auto => full signal)
//
// Normalisation: same 1/N amplitude convention as dft.H and fft.H so that
// all three produce directly comparable spectra.  When Welch averaging is
// used the reported amplitude is the RMS-averaged one-sided spectrum.
//
// Expects in scope (provided by writeOutputs.H before the dispatch):
//   sig               const scalarField&   signal (length N)
//   N                 label
//   Nhalf             label                N/2 + 1
//   T                 scalar               total duration [s]  (unused here;
//                                          the frequency axis is in writeOutputs.H)
//   timeValues        DynamicList<scalar>  actual sample times [s]
//   ftMag             scalarField          output magnitude spectrum [Pa]
//   fourierTransform  const dictionary&    sub-dict from fwhDict
//
// External dependency: FFTW3   (included via fwhFoam.C, linked with -lfftw3)

{
    // * * * * * * * * * Read FFTW-specific dict entries * * * * * * * * * //

    const word   windowType  = fourierTransform.getOrDefault<word>  ("windowType",  "rectangular");
    const Switch removeMean  = fourierTransform.getOrDefault<Switch> ("removeMean",  true);
    const Switch doDetrend   = fourierTransform.getOrDefault<Switch> ("detrend",     false);

    // nfft: segment length.  "auto" or absent => use full signal length.
    label nfft = N;
    if (fourierTransform.found("nfft"))
    {
        // Read as word first; if it is not "auto" parse as a label.
        const word nfftStr = fourierTransform.getOrDefault<word>("nfft", "auto");
        if (nfftStr != "auto")
        {
            nfft = Foam::readLabel(fourierTransform.lookup("nfft"));
            nfft = max(label(2), min(nfft, N));   // clamp to [2, N]
        }
    }

    const label nfftHalf = nfft / 2 + 1;    // one-sided FFTW output bins

    // * * * * * * * * * Build the window coefficient array  * * * * * * * //
    //
    // Each windows/*.H snippet sets win[n] for n in [0, nfft).
    // The snippet sees: nfft (label), n (loop variable), win (scalarField&).

    scalarField win(nfft, 1.0);     // default: rectangular (all ones)
    {
        const scalar pi = constant::mathematical::pi;
        (void)pi;   // consumed by window snippets; suppress warning until implemented


        if (windowType == "Hann" || windowType == "Hanning")
        {
            #include "windows/hann.H"
        }
        else if (windowType == "Hamming")
        {
            #include "windows/hamming.H"
        }
        else if (windowType == "Blackman")
        {
            #include "windows/blackman.H"
        }
        else if (windowType == "FlatTop")
        {
            #include "windows/flattop.H"
        }
        else if (windowType == "Bartlett")
        {
            #include "windows/bartlett.H"
        }
        else if (windowType == "Tukey")
        {
            #include "windows/tukey.H"
        }
        else if (windowType == "rectangular")
        {
            // win already initialised to 1.0
        }
        else
        {
            FatalErrorInFunction
                << "Unknown windowType '" << windowType
                << "' in fwhDict.fourierTransform." << nl
                << "Valid options: Hann, Hamming, Blackman, FlatTop, "
                   "Bartlett, Tukey, rectangular"
                << exit(FatalError);
        }
    }

    // Window normalisation factor (amplitude-preserving: divide by mean(win))
    scalar winNorm = 0.0;
    for (label n = 0; n < nfft; n++) winNorm += win[n];
    winNorm = max(winNorm, VSMALL);

    // * * * * * * * * * * Welch segment averaging setup * * * * * * * * * //
    //
    // 50 % overlap when nfft < N; only one segment when nfft == N.

    const label hop      = (nfft < N) ? nfft / 2 : N;
    const label nSeg     = (nfft < N) ? (N - nfft) / hop + 1 : 1;

    // Accumulate power spectrum across segments
    scalarField powerSum(nfftHalf, 0.0);

    // Allocate FFTW buffers once (reused across segments)
    double* in     = fftw_alloc_real   (static_cast<size_t>(nfft));
    fftw_complex* out = fftw_alloc_complex(static_cast<size_t>(nfftHalf));

    fftw_plan plan = fftw_plan_dft_r2c_1d
    (
        static_cast<int>(nfft),
        in, out,
        FFTW_ESTIMATE        // no measurement : single use in post-processing
    );

    for (label segI = 0; segI < nSeg; segI++)
    {
        const label start = segI * hop;

        // * * * * * * * * * * * Fill input buffer * * * * * * * * * * * * //

        for (label n = 0; n < nfft; n++)
        {
            in[n] = sig[min(start + n, N - 1)];
        }

        // * * * * * * * * Optional: remove segment mean * * * * * * * * * //

        if (removeMean)
        {
            double mu = 0.0;
            for (label n = 0; n < nfft; n++) mu += in[n];
            mu /= nfft;
            for (label n = 0; n < nfft; n++) in[n] -= mu;
        }

        // * * * * * * * * * * * Optional: linear detrend  * * * * * * * * //
        //
        // Fit y = a*n + b by least-squares and subtract.
        // For uniform samples the closed-form solution is:
        //   a = (sum(n*y) - N*mean_n*mean_y) / (sum(n^2) - N*mean_n^2)
        //   b = mean_y - a * mean_n

        if (doDetrend)
        {
            // Fit y = a*t + b by least-squares using actual segment times.
            // This is correct for non-uniform spacing; for uniform spacing it
            // is equivalent to the index-based formula.
            const double t0seg = timeValues[start];
            double sum_t = 0.0, sum_y = 0.0, sum_tt = 0.0, sum_ty = 0.0;
            for (label n = 0; n < nfft; n++)
            {
                const double tn = timeValues[min(start + n, N - 1)] - t0seg;
                sum_t  += tn;
                sum_y  += in[n];
                sum_tt += tn * tn;
                sum_ty += tn * in[n];
            }
            const double denom = nfft * sum_tt - sum_t * sum_t + VSMALL;
            const double a     = (nfft * sum_ty - sum_t * sum_y) / denom;
            const double b     = (sum_y - a * sum_t) / nfft;
            for (label n = 0; n < nfft; n++)
            {
                const double tn = timeValues[min(start + n, N - 1)] - t0seg;
                in[n] -= a * tn + b;
            }
        }

        // * * * * * * * * * * * * * Apply window  * * * * * * * * * * * * //

        for (label n = 0; n < nfft; n++) in[n] *= win[n];

        // * * * * * * * * * * * Execute FFTW plan * * * * * * * * * * * * //
        // * * * * * * * * * * * * * Execute FFTW plan * * * * * * * * * * * * //

        fftw_execute(plan);

        // * * * * * * * * Accumulate one-sided raw power  * * * * * * * * //
        //
        // No winNorm division here : correction is applied once in the
        // final stage loop so the single-segment, rectangular-window case
        // reproduces dft.H / fft.H exactly.

        for (label k = 0; k < nfftHalf; k++)
        {
            const double re = out[k][0];
            const double im = out[k][1];
            powerSum[k] += re*re + im*im;
        }
    }

    // * * * * * * * * * * * Destroy FFTW resources  * * * * * * * * * * * //

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    // * * * * * * * Convert averaged power to amplitude spectrum  * * * * //
    //
    // For a windowed FFT the amplitude-preserving normalisation is:
    //
    //   A(k) = scale * sqrt( mean(|Y(k)|^2) ) / winNorm
    //
    // where winNorm = sum(win) (coherent gain).  This makes the
    // single-segment rectangular case identical to dft.H / fft.H
    // (winNorm == nfft for rectangular, so the formula reduces to
    //  scale * |X(k)| / N as in those files).
    //
    // The result is mapped onto Nhalf bins (original N-based axis):
    // only min(nfftHalf, Nhalf) bins are filled; extras remain zero.

    const label kMax = min(nfftHalf, Nhalf);

    for (label k = 0; k < kMax; k++)
    {
        const scalar meanPower = powerSum[k] / scalar(nSeg);

        // Two-sided → one-sided factor
        const scalar scale = (k == 0 || k == nfftHalf - 1) ? 1.0 : 2.0;

        // Amplitude = scale * sqrt(meanPower) / winNorm
        ftMag[k] = scale * Foam::sqrt(max(meanPower, scalar(0.0))) / winNorm;
    }
}
