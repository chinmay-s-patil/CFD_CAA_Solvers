// fwhPostProcess.H
// Runs after the PIMPLE time loop is finished.
// Integrates and writes FW-H outputs for each sampled surface.

Info << nl << "==================================================" << nl;
Info << "FW-H Post-Processing Phase" << nl;
Info << "==================================================" << nl;

for (label surfI = 0; surfI < fwhNSurfaces; ++surfI)
{
    sampledSurface& surf = (*fwhSamplingPtr)[surfI];
    const word surfaceName = surf.name();

    Info << nl
         << "----------------------------------------" << nl
         << "Processing Surface: " << surfaceName << nl
         << "----------------------------------------" << nl;

    // Set up references that fwhIntegrate.H expects
    const DynamicList<scalar>&      timeValues   = fwhTimeValues[surfI];
    const DynamicList<scalarField>& pSurfHistory = fwhPHistory[surfI];
    const DynamicList<vectorField>& USurfHistory = fwhUHistory[surfI];
    const pointField&               faceCentres  = fwhFaceCentres[surfI];
    const vectorField&              faceNormals  = fwhFaceNormals[surfI];
    const scalarField&              faceAreas    = fwhFaceAreas[surfI];

    if (timeValues.size() < 3)
    {
        WarningInFunction
            << "Not enough time steps sampled for surface " << surfaceName
            << " (need at least 3). Skipping." << nl;
        continue;
    }

    // Set up root dir for outputs
    const fileName rootDir = fwhRootDir;
    const fileName surfDir = rootDir / surfaceName;
    mkDir(surfDir);

    // Run the integration and outputs
    #include "fwhIntegrate.H"
    #include "writeVolumeFields.H"
    #include "writeOutputs.H"
}

Info << nl << "FW-H Post-Processing Complete." << nl;
