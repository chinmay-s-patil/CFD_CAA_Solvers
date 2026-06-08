// writeVolumeFields.H
//
// Evaluates the FW-H surface integral at every mesh cell centroid to produce
// volumetric acoustic fields.  Included from fwhFoam.C after fwhIntegrate.H
// so the following are in scope from that file:
//
//   interp          : generic time-interpolation lambda (const auto& hist, tau, fI)
//   LvecHistory     : List<vectorField>  [tI][fI]  permeable loading vector
//   dLvecHistory    : List<vectorField>  [tI][fI]  d(Lvec)/dt
//   UnHistory       : List<scalarField>  [tI][fI]  surface-normal velocity
//   dUnHistory      : List<scalarField>  [tI][fI]  d(Un)/dt
//   faceCentres     : pointField         face centroids
//   faceAreas       : scalarField        face area magnitudes
//   nHat            : vectorField        unit face normals  (not used here directly)
//   nTimes, nFaces  : label
//   timeValues      : DynamicList<scalar>
//   surfaceName     : word               current surface name
//
// From createFields.H:
//   c0, rho, pRef, fourPi
//   writeVolumePPrime, writeVolumeOASPL
//   fwhSubDict      : dictionary&
//
// From fwhFoam.C loop:
//   mesh, runTime, rootDir

if (writeVolumePPrime || writeVolumeOASPL)
{

// * * * * *  Identify which features apply to this surface  * * * * * * * * //

const word ppSurfName =
    writeVolumePPrime
  ? fwhSubDict.subDict("pPrime").get<word>("surfaceName")
  : word::null;

const word oasplSurfName =
    writeVolumeOASPL
  ? fwhSubDict.subDict("volumeOASPL").get<word>("surfaceName")
  : word::null;

const bool doPP        = writeVolumePPrime && (ppSurfName    == surfaceName);
const bool doVolOASPL  = writeVolumeOASPL  && (oasplSurfName == surfaceName);

if (doPP || doVolOASPL)
{
    // * * *  Lambda: FWH integral at a single cell centroid for one time step
    //
    // Identical formula to fwhIntegrate.H observer loop.  The lambda captures
    // all history arrays and surface geometry from the enclosing scope.
    //
    auto fwhAtCell = [&](const vector& xCell, const scalar t) -> scalar
    {
        scalar pSum = 0.0;
        for (label fI = 0; fI < nFaces; ++fI)
        {
            const vector r    = xCell - faceCentres[fI];
            const scalar rMag = mag(r) + VSMALL;
            const vector rHat = r / rMag;
            const scalar dS   = faceAreas[fI];
            const scalar tau  = t - rMag / c0;

            const vector Lvec_tau  = interp(LvecHistory,  tau, fI);
            const vector dLvec_tau = interp(dLvecHistory, tau, fI);
            const scalar Un_tau    = interp(UnHistory,     tau, fI);
            const scalar dUn_tau   = interp(dUnHistory,    tau, fI);

            const scalar Lr  = Lvec_tau  & rHat;
            const scalar dLr = dLvec_tau & rHat;

            pSum +=
                (  (dLr / (c0 * rMag) + Lr / (rMag * rMag)) / fourPi
                 + rho * dUn_tau / rMag / fourPi
                ) * dS;
        }
        return pSum;
    };

    // * * * * * * Cell centroids (processor-local in parallel runs) * * * * //
    const pointField& CC    = mesh.C();
    const label       nCells = mesh.nCells();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // A) Volumetric p' field  :  one volScalarField per (written) time step
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (doPP)
    {
        const dictionary& ppDict   = fwhSubDict.subDict("pPrime");
        const word  fieldName      = ppDict.get<word>  ("fieldName");
        const label writeEvery     = ppDict.get<label> ("writeInterval");
        const Switch doField       = ppDict.get<Switch>("writeField");
        const Switch doVTK         = ppDict.get<Switch>("writeVTK");
        const Switch doRaw         = ppDict.get<Switch>("writeRaw");

        // postProcessing output dir for raw/VTK files
        const fileName ppDir = rootDir / surfaceName / "volumeFields";
        mkDir(ppDir);

        Info << nl << "Computing volumetric pPrime field for surface \""
             << surfaceName << "\" ..." << nl;

        label nWritten = 0;

        for (label tI = 0; tI < nTimes; ++tI)
        {
            if ((tI % writeEvery) != 0) continue;

            const scalar t = timeValues[tI];

            // Allocate field : processor-local, matches mesh.C()
            scalarField ppLocal(nCells, 0.0);

            forAll(ppLocal, cellI)
                ppLocal[cellI] = fwhAtCell(CC[cellI], t);

            // writeField: OF volScalarField to case time directory
            // The existing .foam file picks this up automatically alongside p, U.
            if (doField)
            {
                runTime.setTime(t, tI);

                volScalarField ppVol
                (
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(fieldName, dimPressure, 0.0)
                );

                // ppLocal is processor-local; forAll over internal cells only
                forAll(ppVol, cellI)
                    ppVol[cellI] = ppLocal[cellI];

                ppVol.write();
            }

            // writeVTK: legacy VTK file to postProcessing dir
            // Writes cell-centroid point cloud (POLYVERTEX) with scalar data.
            // Full mesh topology VTK requires foamToVTK; this provides a quick
            // point-cloud view useful for validating wave fronts in ParaView.
            if (doVTK && Pstream::master())
            {
                // Gather all cell centres and values to master
                pointField  allCC   (CC);
                scalarField allPP   (ppLocal);
                Pstream::combineReduce(allCC, ListOps::appendEqOp<point>());
                Pstream::combineReduce(allPP, ListOps::appendEqOp<scalar>());

                const fileName vtkPath =
                    ppDir / (fieldName + "_" + Foam::name(t) + ".vtk");
                OFstream vtk(vtkPath);

                const label N = allCC.size();
                vtk << "# vtk DataFile Version 3.0\n"
                    << "fwhFoam pPrime t=" << t << "\n"
                    << "ASCII\n"
                    << "DATASET POLYDATA\n"
                    << "POINTS " << N << " double\n";
                forAll(allCC, i)
                    vtk << allCC[i].x() << " "
                        << allCC[i].y() << " "
                        << allCC[i].z() << "\n";
                vtk << "VERTICES 1 " << N + 1 << "\n" << N;
                forAll(allCC, i) vtk << " " << i;
                vtk << "\n"
                    << "POINT_DATA " << N << "\n"
                    << "SCALARS " << fieldName << " double 1\n"
                    << "LOOKUP_TABLE default\n";
                forAll(allPP, i)
                    vtk << allPP[i] << "\n";
            }

            // writeRaw: flat binary blob + companion .info
            if (doRaw && Pstream::master())
            {
                scalarField allPP(ppLocal);
                Pstream::combineReduce(allPP, ListOps::appendEqOp<scalar>());

                const fileName rawPath =
                    ppDir / (fieldName + "_" + Foam::name(t) + ".raw");
                // Binary write via std::ofstream (OFstream is text-mode)
                {
                    std::ofstream raw(rawPath.c_str(),
                                      std::ios::binary | std::ios::trunc);
                    raw.write(
                        reinterpret_cast<const char*>(allPP.cdata()),
                        allPP.size() * sizeof(scalar)
                    );
                }
                // Companion text metadata
                OFstream info(rawPath + ".info");
                info << "time    " << t           << nl
                     << "nCells  " << allPP.size() << nl
                     << "scalar  " << sizeof(scalar) << " bytes" << nl
                     << "endian  little" << nl;
            }

            ++nWritten;
        }

        Info << "  Written " << nWritten << " pPrime snapshot(s)." << nl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // B) Volumetric OASPL field  :  RMS p' over time window -> single snapshot
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (doVolOASPL)
    {
        const dictionary& od  = fwhSubDict.subDict("volumeOASPL");
        const word  fieldName = od.get<word>  ("fieldName");
        const Switch doField  = od.get<Switch>("writeField");
        const Switch doVTK    = od.get<Switch>("writeVTK");

        // Parse time window: scalar → [startTime, tEnd]  or  (minT maxT)
        scalar tMin = timeValues[0];
        scalar tMax = timeValues[nTimes - 1];
        {
            ITstream its = od.lookup("time");
            token tok;
            its >> tok;
            if (tok.isNumber())
            {
                tMax = tok.scalarToken();
            }
            else if (tok.isPunctuation(token::BEGIN_LIST))
            {
                its >> tMin >> tMax;
            }
        }

        Info << nl << "Computing volumetric OASPL field for surface \""
             << surfaceName << "\"" << nl
             << "  time window: [" << tMin << ", " << tMax << "] s ..." << nl;

        // Accumulate sum(p'^2) and count over the window
        scalarField sumP2(nCells, 0.0);
        label nCount = 0;

        for (label tI = 0; tI < nTimes; ++tI)
        {
            const scalar t = timeValues[tI];
            if (t < tMin - SMALL || t > tMax + SMALL) continue;
            ++nCount;

            forAll(sumP2, cellI)
            {
                const scalar pp = fwhAtCell(CC[cellI], t);
                sumP2[cellI] += pp * pp;
            }
        }

        if (nCount == 0)
        {
            WarningInFunction
                << "volumeOASPL: no time steps found in window ["
                << tMin << ", " << tMax << "] for surface '"
                << surfaceName << "'. Skipping OASPL field write." << nl;
        }
        else
        {
            // OASPL(x) = 20 log10( p_rms(x) / pRef )
            scalarField oasplLocal(nCells, 0.0);
            forAll(oasplLocal, cellI)
            {
                const scalar pRms =
                    Foam::sqrt(sumP2[cellI] / scalar(nCount));
                oasplLocal[cellI] =
                    20.0 * Foam::log10(max(pRms, scalar(1e-30)) / pRef);
            }

            // writeField
            if (doField)
            {
                runTime.setTime(tMax, nTimes - 1);

                volScalarField oasplVol
                (
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(fieldName, dimless, 0.0)
                );

                forAll(oasplVol, cellI)
                    oasplVol[cellI] = oasplLocal[cellI];

                oasplVol.write();
                Info << "  OASPL field written to " << runTime.timeName()
                     << "/" << fieldName << nl;
            }

            // writeVTK: same point-cloud approach as pPrime
            if (doVTK && Pstream::master())
            {
                pointField  allCC    (CC);
                scalarField allOASPL (oasplLocal);
                Pstream::combineReduce(allCC,    ListOps::appendEqOp<point>());
                Pstream::combineReduce(allOASPL, ListOps::appendEqOp<scalar>());

                const fileName oasplDir = rootDir / surfaceName / "volumeFields";
                mkDir(oasplDir);

                const fileName vtkPath =
                    oasplDir / (fieldName + "_" + Foam::name(tMin)
                    + "_to_" + Foam::name(tMax) + ".vtk");
                OFstream vtk(vtkPath);

                const label N = allCC.size();
                vtk << "# vtk DataFile Version 3.0\n"
                    << "fwhFoam OASPL [" << tMin << "," << tMax << "]\n"
                    << "ASCII\n"
                    << "DATASET POLYDATA\n"
                    << "POINTS " << N << " double\n";
                forAll(allCC, i)
                    vtk << allCC[i].x() << " "
                        << allCC[i].y() << " "
                        << allCC[i].z() << "\n";
                vtk << "VERTICES 1 " << N + 1 << "\n" << N;
                forAll(allCC, i) vtk << " " << i;
                vtk << "\n"
                    << "POINT_DATA " << N << "\n"
                    << "SCALARS " << fieldName << "_dB double 1\n"
                    << "LOOKUP_TABLE default\n";
                forAll(allOASPL, i)
                    vtk << allOASPL[i] << "\n";

                Info << "  OASPL VTK written -> " << vtkPath << nl;
            }

            Info << "  volumeOASPL used " << nCount << " time step(s)." << nl;
        }
    }

} // end if (doPP || doVolOASPL)

} // end if (writeVolumePPrime || writeVolumeOASPL)
