// Hann windowing

for (label n = 0; n < nfft; n++)
    win[n] = 0.54 - 0.46*Foam::cos(2.0*pi*scalar(n)/scalar(nfft - 1));
