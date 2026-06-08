// writeZoomOutputs.H
// Zoom FFT post-processing block.
// Included at the end of the observer loop in writeOutputs.H when
// zoomPostProcess yes is set in fwhDict.
//
// Expects in scope:
//   sig, N, timeValues                  (from the existing FT block)
//   zoomPostProcess   Switch
//   fwhSubDict        const dictionary&
//   obsDir            fileName
//   obsI              label
//   surfaceName       word
//   observers         List<vector>
//   pRef              scalar

if (zoomPostProcess && N > 1)
{
    const dictionary& zoomDict = fwhSubDict.subDict("zoom");

    const word   zoomType     = zoomDict.get<word>("type");
    const label  zoomM        = zoomDict.get<label>("nLines");
    const scalar fStart       = zoomDict.get<scalar>("fStart");
    const scalar fEnd         = zoomDict.get<scalar>("fEnd");
    const Switch writeZoomFT  = zoomDict.getOrDefault<Switch>("writeZoomFT",  true);
    const Switch writeZoomSPL = zoomDict.getOrDefault<Switch>("writeZoomSPL", true);

    // ------------------------------------------------------------------
    //  Estimate sample rate from timeValues
    // ------------------------------------------------------------------

    const scalar T     = timeValues[N - 1] - timeValues[0];
    const scalar zoomFs = scalar(N - 1) / max(T, VSMALL);

    // Uniformity warning (same threshold as the main FT block)
    {
        const scalar dtMean = T / max(label(1), N - 1);
        scalar dtVar = 0.0;
        for (label i = 0; i < N - 1; ++i)
        {
            const scalar dev = (timeValues[i + 1] - timeValues[i]) - dtMean;
            dtVar += dev * dev;
        }
        const scalar dtCV =
            Foam::sqrt(dtVar / max(N - 1, 1)) / (dtMean + VSMALL);
        if (dtCV > 0.01)
        {
            WarningInFunction
                << "Zoom FFT for observer " << obsI
                << " (surface " << surfaceName << ")"
                << " has non-uniform time spacing (CV = "
                << 100.0*dtCV << " %).\n"
                << "The estimated fs = " << zoomFs
                << " Hz may be inaccurate.  Results are approximate." << nl;
        }
    }

    // Validate band
    if (fStart < 0 || fEnd <= fStart || fEnd > zoomFs * 0.5)
    {
        WarningInFunction
            << "Zoom band [" << fStart << ", " << fEnd << "] Hz is invalid "
            << "for observer " << obsI
            << " (fs = " << zoomFs << " Hz, Nyquist = " << 0.5*zoomFs
            << " Hz). Skipping zoom for this observer." << nl;
    }
    else
    {
        // ------------------------------------------------------------------
        //  Build signal window (re-uses window snippets from windows/ dir)
        // ------------------------------------------------------------------

        const word winType =
            zoomDict.getOrDefault<word>("windowType", "Hanning");

        scalarField win(N, 1.0);
        {
            const label nfft = N;       // window snippets use 'nfft'
            const scalar pi  = constant::mathematical::pi;
            (void)pi;

            if (winType == "Hann" || winType == "Hanning")
            {
                #include "windows/hann.H"
            }
            else if (winType == "Hamming")
            {
                #include "windows/hamming.H"
            }
            else if (winType == "Blackman")
            {
                #include "windows/blackman.H"
            }
            else if (winType == "FlatTop")
            {
                #include "windows/flattop.H"
            }
            else if (winType == "Bartlett")
            {
                #include "windows/bartlett.H"
            }
            else if (winType == "Tukey")
            {
                #include "windows/tukey.H"
            }
            else if (winType == "rectangular")
            {
                // already 1.0
            }
            else
            {
                WarningInFunction
                    << "Unknown zoom.windowType '" << winType
                    << "'. Using rectangular." << nl;
            }
        }

        // ------------------------------------------------------------------
        //  Allocate output arrays
        // ------------------------------------------------------------------

        scalarField zoomFTMag(zoomM, 0.0);
        scalarField zoomFreqs(zoomM, 0.0);

        // ------------------------------------------------------------------
        //  Dispatch to zoom method
        // ------------------------------------------------------------------

        if (zoomType == "CZT")
        {
            #include "fourierTransform/czt.H"
        }
        else if (zoomType == "heterodyne")
        {
            #include "fourierTransform/heterodyne.H"
        }
        else
        {
            FatalErrorInFunction
                << "Unknown zoom.type '" << zoomType << "'." << nl
                << "Valid options: CZT | heterodyne"
                << exit(FatalError);
        }

        // ------------------------------------------------------------------
        //  Write zoomFT.dat
        // ------------------------------------------------------------------

        if (writeZoomFT)
        {
            const fileName zftFile = obsDir / "zoomFT.dat";
            OFstream ozft(zftFile);

            ozft << "# fwhFoam: Zoom FFT magnitude spectrum" << nl
                 << "# Surface   " << surfaceName << nl
                 << "# Observer  " << obsI
                 << "  position  " << observers[obsI] << nl
                 << "# Method    " << zoomType << nl
                 << "# Band      [" << fStart << ", " << fEnd << "] Hz" << nl
                 << "# Lines     " << zoomM << nl
                 << "# Window    " << winType << nl
                 << "# " << nl
                 << "# frequency [Hz]    |P(f)| [Pa]" << nl;

            for (label k = 0; k < zoomM; k++)
                ozft << zoomFreqs[k] << token::TAB << zoomFTMag[k] << nl;

            Info << "    zoomFT.dat written  -> " << zftFile << nl;
        }

        // ------------------------------------------------------------------
        //  Write zoomSPL.dat
        // ------------------------------------------------------------------

        if (writeZoomSPL)
        {
            const fileName zsplFile = obsDir / "zoomSPL.dat";
            OFstream ozspl(zsplFile);

            ozspl << "# fwhFoam: Zoom SPL spectrum" << nl
                  << "# Surface   " << surfaceName << nl
                  << "# Observer  " << obsI
                  << "  position  " << observers[obsI] << nl
                  << "# Method    " << zoomType << nl
                  << "# Band      [" << fStart << ", " << fEnd << "] Hz" << nl
                  << "# Lines     " << zoomM << nl
                  << "# Window    " << winType << nl
                  << "# pRef      " << pRef << " Pa" << nl
                  << "# " << nl
                  << "# frequency [Hz]    SPL [dB]" << nl;

            for (label k = 0; k < zoomM; k++)
            {
                // All zoom bins are interior — convert amplitude to RMS
                const scalar p_rms =
                    max(zoomFTMag[k] / Foam::sqrt(2.0), scalar(1e-30));
                const scalar spl   =
                    20.0 * Foam::log10(p_rms / pRef);

                ozspl << zoomFreqs[k] << token::TAB << spl << nl;
            }

            Info << "    zoomSPL.dat written -> " << zsplFile << nl;
        }
    }
}
