// Tukey windowing

const scalar alpha = 0.5;   // taper fraction; expose via dict if needed
for (label n = 0; n < nfft; n++)
{
    const scalar x = scalar(n)/scalar(nfft - 1);
    if (x < alpha/2.0)
        win[n] = 0.5*(1.0 - Foam::cos(2.0*pi*x/alpha));
    else if (x <= 1.0 - alpha/2.0)
        win[n] = 1.0;
    else
        win[n] = 0.5*(1.0 - Foam::cos(2.0*pi*(1.0 - x)/alpha));
}
