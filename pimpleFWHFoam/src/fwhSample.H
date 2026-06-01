// fwhSample.H
// Included inside the time loop after the PIMPLE correctors finish.
// Samples p and U onto each FW-H surface and appends to history lists.
// On first call grabs the surface geometry (assumed time-independent).

{
    const scalar t = runTime.value();

    for (label surfI = 0; surfI < fwhNSurfaces; ++surfI)
    {
        sampledSurface& surf = (*fwhSamplingPtr)[surfI];
        surf.update();

        // Grab geometry on first sample
        if (!fwhGeometryGrabbed[surfI])
        {
            fwhFaceCentres[surfI] = surf.Cf();
            fwhFaceNormals[surfI] = surf.Sf();
            fwhFaceAreas[surfI]   = surf.magSf();

            Pstream::combineReduce(fwhFaceCentres[surfI], ListOps::appendEqOp<point>());
            Pstream::combineReduce(fwhFaceNormals[surfI], ListOps::appendEqOp<vector>());
            Pstream::combineReduce(fwhFaceAreas[surfI],   ListOps::appendEqOp<scalar>());

            fwhGeometryGrabbed[surfI] = true;
        }

        // Sample p and U onto the surface
        autoPtr<interpolation<scalar>> pInterp =
            interpolation<scalar>::New(interpolationScheme, p);
        autoPtr<interpolation<vector>> UInterp =
            interpolation<vector>::New(interpolationScheme, U);

        scalarField pSurf(surf.sample(*pInterp));
        vectorField USurf(surf.sample(*UInterp));

        Pstream::combineReduce(pSurf, ListOps::appendEqOp<scalar>());
        Pstream::combineReduce(USurf, ListOps::appendEqOp<vector>());

        fwhTimeValues[surfI].append(t);
        fwhPHistory[surfI].append(pSurf);
        fwhUHistory[surfI].append(USurf);
    }

    Info << nl
         << "  FW-H sampled " << fwhNSurfaces << " surface(s) at t = " << t
         << ",  total steps stored = " << fwhTimeValues[0].size() << nl;
}
