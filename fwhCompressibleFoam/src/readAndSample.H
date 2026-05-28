// readAndSample.H
//
// For each surface in samplingPtr, populates:
//   timeValues      DynamicList<scalar>         simulation times used
//   pSurfHistory    DynamicList<scalarField>     pressure history  [t][face]
//   USurfHistory    DynamicList<vectorField>     velocity history  [t][face]
//   rhoSurfHistory  DynamicList<scalarField>     density history   [t][face]
//   faceCentres     pointField
//   faceNormals     vectorField
//   faceAreas       scalarField
//
// Cache files are written per surface under:
//   postProcessing/fwhCompressibleFoam/<surfaceName>/timeHistory_p.csv  (and Ux/Uy/Uz/rho)

#include <fstream>
#include <sstream>
#include <string>

DynamicList<scalar>      timeValues;
DynamicList<scalarField> pSurfHistory;
DynamicList<vectorField> USurfHistory;
DynamicList<scalarField> rhoSurfHistory;  // per-face density (compressible)

const word surfaceName = surf.name();

const fileName surfDir = rootDir / surfaceName;
mkDir(surfDir);

const fileName cacheP   = surfDir / "timeHistory_p.csv";
const fileName cacheUx  = surfDir / "timeHistory_Ux.csv";
const fileName cacheUy  = surfDir / "timeHistory_Uy.csv";
const fileName cacheUz  = surfDir / "timeHistory_Uz.csv";
const fileName cacheRho = surfDir / "timeHistory_rho.csv";

const bool cacheExists =
    isFile(cacheP) && isFile(cacheUx) && isFile(cacheUy)
    && isFile(cacheUz) && isFile(cacheRho);

const bool doSampling = bool(recomputeTimeHistory) || !cacheExists;

if (doSampling)
{
    if (cacheExists)
        Info << nl << "recomputeTimeHistory = yes -> overwriting cached files." << nl;
    else
        Info << nl << "Cache not found -> sampling from CFD time directories." << nl;

    label validTimeCount = 0;

    forAll(runTime.times(), timeI)
    {
        if (timeI == 0) continue;

        runTime.setTime(runTime.times()[timeI], timeI);
        const scalar t = runTime.value();

        if (t < startTime) continue;
        if (t > endTime)   break;

        if ((validTimeCount++ % timeStepStride) != 0) continue;

        Info << nl << "  Sampling t = " << runTime.timeName() << nl;

        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        surf.update();

        autoPtr<interpolation<scalar>> pInterp =
            interpolation<scalar>::New(interpolationScheme, p);
        autoPtr<interpolation<vector>> UInterp =
            interpolation<vector>::New(interpolationScheme, U);
        autoPtr<interpolation<scalar>> rhoInterp =
            interpolation<scalar>::New(interpolationScheme, rho);

        scalarField pSurf  (surf.sample(*pInterp));
        vectorField USurf  (surf.sample(*UInterp));
        scalarField rhoSurf(surf.sample(*rhoInterp));

        Pstream::combineReduce(pSurf,   ListOps::appendEqOp<scalar>());
        Pstream::combineReduce(USurf,   ListOps::appendEqOp<vector>());
        Pstream::combineReduce(rhoSurf, ListOps::appendEqOp<scalar>());

        Info << "  Surface faces = " << pSurf.size() << nl
             << "  p surf: min = " << min(pSurf)
             << "  max = " << max(pSurf) << nl;

        timeValues.append(t);
        pSurfHistory.append(pSurf);
        USurfHistory.append(USurf);
        rhoSurfHistory.append(rhoSurf);
    }

    Info << nl << "Finished sampling" << nl;

    if (writeTimeHistory && Pstream::master())
    {
        Info << nl << "Writing time history files ..." << nl;

        const label nT = timeValues.size();
        const label nF = nT > 0 ? pSurfHistory[0].size() : 0;

        auto writeCSV =
            [&](const fileName& path,
                const DynamicList<scalarField>& data,
                const word& fieldName)
        {
            OFstream os(path);
            os  << "# fwhFoam: " << fieldName << " time history" << nl
                << "# nTimes=" << nT << "  nFaces=" << nF << nl
                << "# time";
            for (label i = 0; i < nF; i++)
                os << token::TAB << "f" << i;
            os << nl;

            for (label k = 0; k < nT; k++)
            {
                os << timeValues[k];
                for (label i = 0; i < nF; i++)
                    os << token::TAB << data[k][i];
                os << nl;
            }
            Info << "  Wrote: " << path << nl;
        };

        DynamicList<scalarField> UxH(nT), UyH(nT), UzH(nT);
        for (label k = 0; k < nT; k++)
        {
            scalarField Ux(nF), Uy(nF), Uz(nF);
            forAll(USurfHistory[k], fI)
            {
                Ux[fI] = USurfHistory[k][fI].x();
                Uy[fI] = USurfHistory[k][fI].y();
                Uz[fI] = USurfHistory[k][fI].z();
            }
            UxH.append(Ux);
            UyH.append(Uy);
            UzH.append(Uz);
        }

        writeCSV(cacheP,   pSurfHistory,  "pressure");
        writeCSV(cacheUx,  UxH,           "Ux");
        writeCSV(cacheUy,  UyH,           "Uy");
        writeCSV(cacheUz,  UzH,           "Uz");
        writeCSV(cacheRho, rhoSurfHistory, "rho");
    }
    else if (!writeTimeHistory)
    {
        Info << nl
             << "writeTimeHistory = no -> skipping CSV write; "
             << "data held in memory only." << nl;
    }
}
else
{
    Info << nl
         << "recomputeTimeHistory = no and all cache files present." << nl
         << "Reading time history from:" << nl
         << "  " << cacheP  << nl
         << "  " << cacheUx << nl
         << "  " << cacheUy << nl
         << "  " << cacheUz << nl;

    auto readCSV =
        [&](const fileName& path,
            DynamicList<scalar>& times,
            bool fillTimes) -> DynamicList<scalarField>
    {
        std::ifstream ifs(path.c_str());
        if (!ifs.good())
            FatalErrorInFunction
                << "Cannot open cache file: " << path
                << exit(FatalError);

        DynamicList<scalarField> table;
        std::string line;
        bool headerSkipped = false;

        while (std::getline(ifs, line))
        {
            if (line.empty() || line[0] == '#') continue;
            if (!headerSkipped) { headerSkipped = true; continue; }

            std::istringstream ss(line);
            DynamicList<scalar> row;
            scalar val;
            while (ss >> val) row.append(val);
            if (row.size() < 2) continue;

            if (fillTimes) times.append(row[0]);

            const label nF = row.size() - 1;
            scalarField faceVals(nF);
            for (label i = 0; i < nF; i++)
                faceVals[i] = row[i + 1];

            table.append(faceVals);
        }
        return table;
    };

    DynamicList<scalar> sink;

    auto pTable   = readCSV(cacheP,   timeValues, true);
    auto UxTable  = readCSV(cacheUx,  sink,        false);
    auto UyTable  = readCSV(cacheUy,  sink,        false);
    auto UzTable  = readCSV(cacheUz,  sink,        false);
    auto rhoTable = readCSV(cacheRho, sink,        false);

    const label nT = timeValues.size();
    const label nF = nT > 0 ? pTable[0].size() : 0;

    if (UxTable.size() != nT || UyTable.size() != nT
    ||  UzTable.size() != nT || rhoTable.size() != nT)
    {
        FatalErrorInFunction
            << "Cache files have inconsistent row counts:" << nl
            << "  p   : " << nT              << nl
            << "  Ux  : " << UxTable.size()  << nl
            << "  Uy  : " << UyTable.size()  << nl
            << "  Uz  : " << UzTable.size()  << nl
            << "  rho : " << rhoTable.size() << nl
            << "Delete the cache files and re-run with recomputeTimeHistory yes."
            << exit(FatalError);
    }

    for (label k = 0; k < nT; k++)
    {
        pSurfHistory.append(pTable[k]);
        rhoSurfHistory.append(rhoTable[k]);

        vectorField USurf(nF);
        for (label fI = 0; fI < nF; fI++)
            USurf[fI] = vector(UxTable[k][fI], UyTable[k][fI], UzTable[k][fI]);

        USurfHistory.append(USurf);
    }

    Info << "  Loaded " << nT << " time steps, " << nF << " faces per step." << nl;
}

// Surface geometry (time-independent for a stationary FW-H surface)
{
    bool surfaceUpdated = false;
    forAll(runTime.times(), timeI)
    {
        if (timeI == 0) continue;
        runTime.setTime(runTime.times()[timeI], timeI);
        const scalar t = runTime.value();
        if (t >= startTime && t <= endTime)
        {
            surf.update();
            surfaceUpdated = true;
            break;
        }
    }

    if (!surfaceUpdated)
        FatalErrorInFunction
            << "No time steps found in sampling window ["
            << startTime << ", " << endTime << "]. "
            << "Check startTime/endTime in fwhDict."
            << exit(FatalError);
}

pointField  faceCentres(surf.Cf());
vectorField faceNormals(surf.Sf());
scalarField faceAreas  (surf.magSf());

Pstream::combineReduce(faceCentres, ListOps::appendEqOp<point>());
Pstream::combineReduce(faceNormals, ListOps::appendEqOp<vector>());
Pstream::combineReduce(faceAreas,   ListOps::appendEqOp<scalar>());

Info << nl
     << "Surface \"" << surfaceName << "\": "
     << timeValues.size() << " steps, "
     << faceCentres.size() << " faces" << nl;
