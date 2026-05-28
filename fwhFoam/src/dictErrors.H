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
