if (observers.size() < 1)
{
    FatalErrorInFunction
        << "No observers defined in fwhDict."
        << exit(FatalError);
}

if (fwhSamplingPtr->empty())
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

// Helper: collect all surface names from the sampling block
wordList sampledSurfaceNames(fwhSamplingPtr->size());
forAll(*fwhSamplingPtr, sI)
    sampledSurfaceNames[sI] = (*fwhSamplingPtr)[sI].name();

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

    // Validate the 'time' entry: must be scalar or 2-element list.
    // Use lookup() -> ITstream to avoid the deprecated
    // lookupEntry(word, bool, bool) overload (removed in OF >= 2306).
    ITstream its = od.lookup("time");
    token tok;
    its >> tok;
    if (tok.isNumber())
    {
        const scalar tEnd = tok.scalarToken();
        if (tEnd <= runTime.startTime().value())
            FatalErrorInFunction
                << "volumeOASPL.time (" << tEnd
                << ") must be greater than startTime."
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
