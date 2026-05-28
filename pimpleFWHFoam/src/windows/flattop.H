// Flat-top windowing

const scalar a0 = 0.21557895, a1 = 0.41663158, a2 = 0.277263158;
const scalar a3 = 0.083578947, a4 = 0.006947368;
for (label n = 0; n < nfft; n++)
{
    const scalar phi = 2.0*pi*scalar(n)/scalar(nfft - 1);
    win[n] = a0 - a1*Foam::cos(phi)   + a2*Foam::cos(2.0*phi)
                - a3*Foam::cos(3.0*phi) + a4*Foam::cos(4.0*phi);
}
