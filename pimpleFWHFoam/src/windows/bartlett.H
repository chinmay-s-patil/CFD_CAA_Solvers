// Bartlett windowing

for (label n = 0; n < nfft; n++)
    win[n] = 1.0 - Foam::mag(2.0*scalar(n)/scalar(nfft - 1) - 1.0);
