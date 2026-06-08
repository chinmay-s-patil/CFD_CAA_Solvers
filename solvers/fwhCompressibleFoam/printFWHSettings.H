Info<<nl<<nl<<"FWH Settings"<<nl
    <<"-------------------"<<endl;
Info<<"rhoRef                   = "<<rhoRef<<endl;
Info<<"c0                       = "<<c0<<endl;
Info<<"pRef                     = "<<pRef<<endl;
Info<<"Observers                = "<<observers<<endl;
Info<<"Time Period               = "<<startTime<<"s to "<< endTime << "s"<<endl;
Info<<"Time Step Stride         = "<<timeStepStride<<endl;
Info<<"Time Interp Scheme       = "<<timeInterpolationScheme<<endl;
Info<<nl;
Info<<"writeTimeHistory         = "<< (writeTimeHistory ? "yes" : "no")<<endl;
Info<<"recomputeTimeHistory     = "<< (recomputeTimeHistory ? "yes" : "no")<<endl;
Info<<nl;
Info<<"writeSignal              = "<< (writeSignal ? "yes" : "no")<<endl;
Info<<"writeSPL                 = "<< (writeSPL ? "yes" : "no")<<endl;
Info<<"writeFT                  = "<< (writeFT ? "yes" : "no")<<endl;
if (writeFT)
{
    Info<<"  FT type              = "<<fourierTransformType<<endl;
}
Info<<nl;
Info<<"writeOASPL               = "<< (writeOASPL ? "yes" : "no")<<endl;
Info<<"writeVolumePPrime        = "<< (writeVolumePPrime ? "yes" : "no")<<endl;
Info<<"writeVolumeOASPL         = "<< (writeVolumeOASPL ? "yes" : "no")<<endl;
Info<<nl;
Info<<"zoomPostProcess          = "<< (zoomPostProcess ? "yes" : "no")<<endl;
if (zoomPostProcess)
{
    const dictionary& zd = fwhSubDict.subDict("zoom");
    Info<<"  zoom.type            = "<< zd.get<word>("type") <<endl;
    Info<<"  zoom.fStart          = "<< zd.get<scalar>("fStart") <<" Hz"<<endl;
    Info<<"  zoom.fEnd            = "<< zd.get<scalar>("fEnd")   <<" Hz"<<endl;
    Info<<"  zoom.nLines          = "<< zd.get<label>("nLines")  <<endl;
    Info<<"  zoom.windowType      = "
        << zd.getOrDefault<word>("windowType", "Hanning") <<endl;
    Info<<"  zoom.removeMean      = "
        << (zd.getOrDefault<Switch>("removeMean", true) ? "yes" : "no") <<endl;
    Info<<"  zoom.detrend         = "
        << (zd.getOrDefault<Switch>("detrend", false) ? "yes" : "no") <<endl;
    Info<<"  zoom.writeZoomFT     = "
        << (zd.getOrDefault<Switch>("writeZoomFT",  true) ? "yes" : "no") <<endl;
    Info<<"  zoom.writeZoomSPL    = "
        << (zd.getOrDefault<Switch>("writeZoomSPL", true) ? "yes" : "no") <<endl;

    const word zType = zd.get<word>("type");
    if (zType == "heterodyne")
    {
        Info<<"  zoom.decimationFactor= "
            << zd.getOrDefault<word>("decimationFactor", "auto") <<endl;
        Info<<"  zoom.lpfOrder        = "
            << zd.getOrDefault<label>("lpfOrder", 64) <<endl;
        Info<<"  zoom.lpfWindow       = "
            << zd.getOrDefault<word>("lpfWindow", "Hamming") <<endl;
    }
}