// printVolumeFieldSettings.H
//
// Prints pPrime volumetric field and OASPL directivity settings to Info.
// Included from fwhFoam.C after printSamplingSettings.H, before readAndSample.H.
// Requires in scope (set by createFields.H):
//   writeOASPL, oasplCentre, oasplRefDirection, oasplFileName
//   writeVolumePPrime, writeVolumeOASPL
//   fwhSubDict

Info << nl << "Volume Field / OASPL Settings" << nl
     << "------------------------------" << nl;
Info << "writeOASPL               = " << (writeOASPL       ? "yes" : "no") << nl;
Info << "writeVolumePPrime        = " << (writeVolumePPrime ? "yes" : "no") << nl;
Info << "writeVolumeOASPL         = " << (writeVolumeOASPL  ? "yes" : "no") << nl;

if (writeOASPL)
{
    Info << "  OASPL.centre         = " << oasplCentre       << nl;
    Info << "  OASPL.refDirection   = " << oasplRefDirection  << nl;
    Info << "  OASPL.fileName       = " << oasplFileName      << nl;
}

if (writeVolumePPrime)
{
    const dictionary& ppDict = fwhSubDict.subDict("pPrime");
    Info << "  pPrime.surfaceName   = " << ppDict.get<word>("surfaceName")    << nl;
    Info << "  pPrime.fieldName     = " << ppDict.get<word>("fieldName")      << nl;
    Info << "  pPrime.writeInterval = " << ppDict.get<label>("writeInterval") << nl;
    Info << "  pPrime.writeField    = " << ppDict.get<Switch>("writeField")   << nl;
    Info << "  pPrime.writeVTK      = " << ppDict.get<Switch>("writeVTK")     << nl;
    Info << "  pPrime.writeRaw      = " << ppDict.get<Switch>("writeRaw")     << nl;
}

if (writeVolumeOASPL)
{
    const dictionary& od = fwhSubDict.subDict("volumeOASPL");
    Info << "  volumeOASPL.surfaceName = " << od.get<word>("surfaceName") << nl;
    Info << "  volumeOASPL.fieldName   = " << od.get<word>("fieldName")   << nl;
    Info << "  volumeOASPL.writeField  = " << od.get<Switch>("writeField") << nl;
    Info << "  volumeOASPL.writeVTK    = " << od.get<Switch>("writeVTK")   << nl;

    // Print time window in a human-readable form
    ITstream its = od.lookup("time");
    token tok;
    its >> tok;
    if (tok.isNumber())
    {
        Info << "  volumeOASPL.time     = [startTime, " << tok.scalarToken() << " s]" << nl;
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        scalar v1, v2;
        its >> v1 >> v2;
        Info << "  volumeOASPL.time     = [" << v1 << " s, " << v2 << " s]" << nl;
    }
}

Info << nl;
