if (startTime >= endTime)
{
    FatalErrorInFunction
        << "startTime (" << startTime << ") must be less than endTime ("
        << endTime << ")." << nl
        << "Check startTime and endTime in fwhDict."
        << exit(FatalError);
}

if (observers.size() < 1)
{
    FatalErrorInFunction
        << "No observers defined in fwhDict."
        << exit(FatalError);
}

if (samplingPtr->empty())
{
    FatalErrorInFunction
        << "No surfaces defined in fwhDict sampling sub-dict."
        << exit(FatalError);
}

if (interpolationScheme.empty())
{
    FatalErrorInFunction
        << "faceInterpolationScheme in fwhDict sampling sub-dict is empty."
        << exit(FatalError);
}

{
    const wordList validTimeSchemes
    {
        "linear", "nearest", "hermiteCubic",
        "spline", "akimaCubic", "piecewiseHermite"
    };
    if (!validTimeSchemes.found(timeInterpolationScheme))
    {
        FatalErrorInFunction
            << "Unknown timeInterpolationScheme '"
            << timeInterpolationScheme << "'." << nl
            << "Valid options are: " << validTimeSchemes
            << exit(FatalError);
    }
}

{
    const wordList validFTTypes{"DFT", "NUDFT", "FFT", "FFTW"};
    if (writeFT && !validFTTypes.found(fourierTransformType))
    {
        FatalErrorInFunction
            << "Unknown fourierTransform type '"
            << fourierTransformType << "'." << nl
            << "Valid options are: " << validFTTypes
            << exit(FatalError);
    }
}

// * * * * * * * * *  Validate volumetric field sub-dicts  * * * * * * * * * //

wordList sampledSurfaceNames(samplingPtr->size());
forAll(*samplingPtr, sI)
    sampledSurfaceNames[sI] = (*samplingPtr)[sI].name();

if (writeVolumePPrime)
{
    const dictionary& ppDict = fwhSubDict.subDict("pPrime");
    const word sName = ppDict.get<word>("surfaceName");
    if (!sampledSurfaceNames.found(sName))
    {
        FatalErrorInFunction
            << "pPrime.surfaceName '" << sName
            << "' not found in sampling.surfaces{}." << nl
            << "Available: " << sampledSurfaceNames
            << exit(FatalError);
    }
    if (ppDict.get<label>("writeInterval") < 1)
    {
        FatalErrorInFunction
            << "pPrime.writeInterval must be >= 1."
            << exit(FatalError);
    }
}

if (writeVolumeOASPL)
{
    const dictionary& od = fwhSubDict.subDict("volumeOASPL");
    const word sName = od.get<word>("surfaceName");
    if (!sampledSurfaceNames.found(sName))
    {
        FatalErrorInFunction
            << "volumeOASPL.surfaceName '" << sName
            << "' not found in sampling.surfaces{}." << nl
            << "Available: " << sampledSurfaceNames
            << exit(FatalError);
    }

    ITstream its = od.lookup("time");
    token tok;
    its >> tok;
    if (tok.isNumber())
    {
        const scalar tEnd = tok.scalarToken();
        if (tEnd <= startTime)
            FatalErrorInFunction
                << "volumeOASPL.time (" << tEnd
                << ") must be greater than startTime (" << startTime << ")."
                << exit(FatalError);
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        scalar v1, v2;
        its >> v1 >> v2;
        if (v1 >= v2)
            FatalErrorInFunction
                << "volumeOASPL.time range (" << v1 << " " << v2
                << "): minT must be less than maxT."
                << exit(FatalError);
    }
    else
    {
        FatalErrorInFunction
            << "volumeOASPL.time must be a scalar or a two-element list (minT maxT)."
            << exit(FatalError);
    }
}

if (writeOASPL && mag(oasplRefDirection) < SMALL)
{
    FatalErrorInFunction
        << "OASPL.refDirection is a zero vector. "
        << "Specify a non-zero direction for θ = 0°."
        << exit(FatalError);
}

// * * * * * * * * * * *  Validate zoom sub-dict  * * * * * * * * * * * * * //

if (zoomPostProcess)
{
    if (!fwhSubDict.found("zoom"))
    {
        FatalErrorInFunction
            << "zoomPostProcess is enabled but no zoom{} sub-dict found in fwhDict."
            << exit(FatalError);
    }

    const dictionary& zd = fwhSubDict.subDict("zoom");

    const wordList validZoomTypes{"CZT", "heterodyne"};
    const word zType = zd.get<word>("type");
    if (!validZoomTypes.found(zType))
    {
        FatalErrorInFunction
            << "Unknown zoom.type '" << zType << "'." << nl
            << "Valid options are: " << validZoomTypes
            << exit(FatalError);
    }

    const scalar fStart = zd.get<scalar>("fStart");
    const scalar fEnd   = zd.get<scalar>("fEnd");
    if (fStart < 0 || fEnd <= fStart)
    {
        FatalErrorInFunction
            << "zoom.fStart (" << fStart << ") must be >= 0 and "
            << "zoom.fEnd (" << fEnd << ") must be > fStart."
            << exit(FatalError);
    }

    const label nLines = zd.get<label>("nLines");
    if (nLines < 2)
    {
        FatalErrorInFunction
            << "zoom.nLines (" << nLines << ") must be >= 2."
            << exit(FatalError);
    }

    if (zd.found("windowType"))
    {
        const wordList validWin
        {
            "Hann", "Hanning", "Hamming", "Blackman",
            "FlatTop", "Bartlett", "Tukey", "rectangular"
        };
        const word wType = zd.get<word>("windowType");
        if (!validWin.found(wType))
        {
            FatalErrorInFunction
                << "Unknown zoom.windowType '" << wType << "'." << nl
                << "Valid options are: " << validWin
                << exit(FatalError);
        }
    }
}
