
Info << nl << nl << "Sampling Settings" << nl
     << "-------------------" << endl;
Info << "Surfaces             = " << fwhSamplingPtr->size() << endl;
forAll(*fwhSamplingPtr, surfI)
{
    Info << "  [" << surfI << "] " << (*fwhSamplingPtr)[surfI].name() << endl;
}
Info << "Interpolation Scheme = " << interpolationScheme << endl;
Info << "Surface Format       = " << surfaceFormat << endl;
Info << "Fields               = " << fwhFields << endl;
