// writeOutputs.H
// Writes acoustic pressure time history (and optionally FT / SPL) for each
// observer to postProcessing/fwhFoam/<surfaceName>/observer<N>/.
//
// Expects in scope (set by fwhFoam.C / fwhIntegrate.H):
//   surfaceName     word
//   rootDir         fileName           base postProcessing directory
//   timeValues      DynamicList<scalar>
//   pAcoustic       List<scalarField>  [observer][timeIndex]
//   observers       List<vector>
//   nObs            label
//   nTimes          label
//   writeSignal     Switch
//   writeFT         Switch
//   writeSPL        Switch
//   pRef            scalar

const fileName surfOutDir = rootDir / surfaceName;
mkDir(surfOutDir);

Info << nl
     << "======================================" << nl
     << "  fwhFoam output: " << surfaceName << nl
     << "======================================" << nl
     << "  Observers    : " << nObs   << nl
     << "  Time samples : " << nTimes << nl
     << "  Output dir   : " << surfOutDir << nl
     << "======================================" << nl;

forAll(observers, obsI)
{
    const fileName obsDir = surfOutDir / ("observer" + Foam::name(obsI));
    mkDir(obsDir);

    const scalarField& sig = pAcoustic[obsI];

    Info << nl << "  Observer " << obsI
         << "  position " << observers[obsI] << nl
         << "    samples : " << nTimes << nl
         << "    p' min  : " << min(sig) << " Pa" << nl
         << "    p' max  : " << max(sig) << " Pa" << nl;

    if (writeSignal)
    {
        const fileName sigFile = obsDir / "pTime.dat";
        OFstream os(sigFile);

        os  << "# fwhFoam: acoustic pressure time history" << nl
            << "# Surface  " << surfaceName << nl
            << "# Observer " << obsI
            << "  position " << observers[obsI] << nl
            << "# " << nl
            << "# time [s]          p' [Pa]" << nl;

        for (label tI = 0; tI < nTimes; tI++)
        {
            os  << timeValues[tI] << token::TAB
                << sig[tI]        << nl;
        }

        Info << "    pTime.dat written -> " << sigFile << nl;
    }

    if (writeFT && nTimes > 1)
    {
        const label N = nTimes;

        const scalar T     = timeValues[N - 1] - timeValues[0];
        const scalar df    = 1.0 / max(T, VSMALL);
        const label  Nhalf = N / 2 + 1;

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
            if (dtCV > 0.01
             && (   fourierTransformType == "DFT"
                 || fourierTransformType == "FFT"
                 || fourierTransformType == "FFTW"))
            {
                WarningInFunction
                    << "Time series for observer " << obsI
                    << " (surface " << surfaceName << ")"
                    << " has non-uniform spacing (CV = "
                    << 100.0*dtCV << " %).\n"
                    << fourierTransformType
                    << " assumes uniform sampling.  Switch to"
                    << " fourierTransformType NUDFT for correct results"
                    << " with non-uniform time series." << nl;
            }
        }

        scalarField ftMag(Nhalf, 0.0);

        if (fourierTransformType == "DFT")
        {
            #include "fourierTransform/dft.H"
        }
        else if (fourierTransformType == "NUDFT")
        {
            #include "fourierTransform/nudft.H"
        }
        else if (fourierTransformType == "FFT")
        {
            #include "fourierTransform/fft.H"
        }
        else if (fourierTransformType == "FFTW")
        {
            #include "fourierTransform/fftw.H"
        }
        else
        {
            FatalErrorInFunction
                << "Unknown fourierTransformType '" << fourierTransformType
                << "' in fwhDict." << nl
                << "Valid options are: DFT, NUDFT, FFT, FFTW"
                << exit(FatalError);
        }

        {
            const fileName ftFile = obsDir / "ft.dat";
            OFstream os(ftFile);

            os  << "# fwhFoam: FT magnitude spectrum" << nl
                << "# Surface  " << surfaceName << nl
                << "# Observer " << obsI
                << "  position " << observers[obsI] << nl
                << "# df = " << df << " Hz  (T = " << T << " s,  N = " << N << ")" << nl
                << "# " << nl
                << "# frequency [Hz]    |P(f)| [Pa]" << nl;

            for (label k = 0; k < Nhalf; k++)
            {
                os  << scalar(k) * df << token::TAB
                    << ftMag[k]       << nl;
            }

            Info << "    ft.dat written   -> " << ftFile << nl;
        }

        if (writeSPL)
        {
            const fileName splFile = obsDir / "spl.dat";
            OFstream os2(splFile);

            os2 << "# fwhFoam: Sound Pressure Level" << nl
                << "# Surface  " << surfaceName << nl
                << "# Observer " << obsI
                << "  position " << observers[obsI] << nl
                << "# pRef = " << pRef << " Pa" << nl
                << "# " << nl
                << "# frequency [Hz]    SPL [dB]" << nl;

            for (label k = 0; k < Nhalf; k++)
            {
                const scalar p_rms = max(
                    (k == 0) ? ftMag[k] : ftMag[k] / Foam::sqrt(2.0),
                    scalar(1e-30)
                );
                const scalar spl = 20.0 * Foam::log10(p_rms / pRef);

                os2 << scalar(k) * df << token::TAB
                    << spl             << nl;
            }

            Info << "    spl.dat written  -> " << splFile << nl;
        }
    }
    else if (writeFT)
    {
        WarningInFunction
            << "writeFT requested but fewer than 2 time samples available "
            << "for observer " << obsI << "  skipping FT/SPL." << nl;
    }
}

// * * * * * * * * * * * * * * OASPL Directivity * * * * * * * * * * * * * * //
//
// Computes OASPL = 20 log10(p_rms / pRef) at each observer from the time signal,
// then writes [theta_deg, OASPL_dB] to oaspl_directivity.dat.
//
// theta is the angle of each observer relative to oasplCentre, measured from
// oasplRefDirection in the plane perpendicular to their cross product (or just
// the atan2 of the projected components for 2-D cases).
//
// Requires in scope: writeOASPL, oasplCentre, oasplRefDirection, oasplFileName,
//                    observers, pAcoustic, nObs, nTimes, pRef, surfOutDir.

if (writeOASPL && nObs > 0 && nTimes > 0)
{
    const fileName oasplFile = surfOutDir / oasplFileName + ".dat";
    OFstream oas(oasplFile);

    oas << "# fwhFoam: OASPL directivity" << nl
        << "# Surface    " << surfaceName << nl
        << "# centre     " << oasplCentre << nl
        << "# refDir     " << oasplRefDirection << nl
        << "# " << nl
        << "# theta [deg]    OASPL [dB]    observer" << nl;

    // Build an orthonormal basis in the plane of the arc:
    //   e1 = normalised refDirection
    //   e2 = normalised (e1 × up) or any perpendicular
    // For a general 3-D arc the user places observers, so we compute the angle
    // as atan2( r·e2 , r·e1 ) where e2 is orthogonal to e1 in the arc plane.
    // We derive e2 as the component of the "most misaligned" axis perpendicular to e1.

    const vector e1 = oasplRefDirection / (mag(oasplRefDirection) + VSMALL);

    // Pick a helper axis not parallel to e1
    vector helper(0, 0, 1);
    if (mag(e1 ^ helper) < 0.1) helper = vector(0, 1, 0);
    const vector e2 = (helper - (helper & e1)*e1);
    const vector e2n = e2 / (mag(e2) + VSMALL);

    forAll(observers, obsI)
    {
        const vector r = observers[obsI] - oasplCentre;

        // Polar angle in the arc plane
        const scalar theta =
            Foam::atan2(r & e2n, r & e1) * 180.0 / constant::mathematical::pi;

        // p_rms from time signal
        const scalarField& sig = pAcoustic[obsI];
        scalar sumSq = 0.0;
        forAll(sig, tI) sumSq += sig[tI] * sig[tI];
        const scalar pRms = Foam::sqrt(sumSq / max(scalar(nTimes), SMALL));
        const scalar oasplVal =
            20.0 * Foam::log10(max(pRms, scalar(1e-30)) / pRef);

        oas << theta << token::TAB << oasplVal
            << token::TAB << "# obs " << obsI
            << "  pos " << observers[obsI] << nl;
    }

    Info << nl << "  OASPL directivity written -> " << oasplFile << nl;
}

Info << nl
     << "  Output complete: " << surfaceName << nl;
