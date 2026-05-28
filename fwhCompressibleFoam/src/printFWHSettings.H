
Info<<nl<<nl<<"FWH Settings"<<nl
    <<"-------------------"<<endl;
Info<<"rho                      = "<<rho<<endl;
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
