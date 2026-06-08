
Info << nl << nl << "Sampling Settings" << nl
     << "-------------------" << endl;
Info << "Surfaces             = " << samplingPtr->size() << endl;
forAll(*samplingPtr, surfI)
{
    Info << "  [" << surfI << "] " << (*samplingPtr)[surfI].name() << endl;
}
Info << "Interpolation Scheme = " << interpolationScheme << endl;
Info << "Surface Format       = " << surfaceFormat << endl;
Info << "Fields               = " << fields << endl;
