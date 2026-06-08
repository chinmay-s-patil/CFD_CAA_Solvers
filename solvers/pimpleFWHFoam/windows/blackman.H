// Blackman windowing

for (label n = 0; n < nfft; n++)
{
    const scalar phi = 2.0*pi*scalar(n)/scalar(nfft - 1);
    win[n] = 0.42 - 0.5*Foam::cos(phi) + 0.08*Foam::cos(2.0*phi);
}
