// piecewiseHermite_setup.H
// Precomputation step for monotone piecewise cubic Hermite interpolation.
//
// Called ONCE per face (fI) before the observer / reception-time loop.
// Fills phSlopesScalar[fI] and phSlopesVector[fI] with the Fritsch-Carlson
// monotone knot slopes m[i] for both history arrays.
//
// Complexity: O(N) per face.
//
// Expects in scope:
//   fI                label
//   nTimes            label
//   timeValues        DynamicList<scalar>
//   LvecHistory       List<vectorField>   [tI][fI]
//   UnHistory         List<scalarField>   [tI][fI]
//   phSlopesScalar    List<scalarField>   [fI] output (sized nTimes)
//   phSlopesVector    List<List<vector>>  [fI] output (sized nTimes)

{
    const label N = nTimes - 1;

    scalarField hSpan(N);
    for (label i = 0; i < N; ++i)
        hSpan[i] = timeValues[i + 1] - timeValues[i] + VSMALL;

    // Helper: given scalar f[0..nTimes-1], return Fritsch-Carlson slopes m[].
    auto fcSlopes = [&](auto getFval) -> scalarField
    {
        scalarField delta(N);
        for (label i = 0; i < N; ++i)
            delta[i] = (getFval(i + 1) - getFval(i)) / hSpan[i];

        scalarField m(nTimes, 0.0);

        // Interior: non-uniform centred differences
        for (label i = 1; i < N; ++i)
        {
            const scalar dt_m = timeValues[i + 1] - timeValues[i - 1] + VSMALL;
            m[i] = (getFval(i + 1) - getFval(i - 1)) / dt_m;
        }

        // Boundary: one-sided
        if (nTimes >= 3)
        {
            const scalar dt_p = timeValues[2] - timeValues[0] + VSMALL;
            m[0] = (getFval(2) - getFval(0)) / dt_p;

            const scalar dt_m = timeValues[N] - timeValues[N - 2] + VSMALL;
            m[N] = (getFval(N) - getFval(N - 2)) / dt_m;
        }
        else
        {
            m[0] = delta[0];
            m[N] = delta[N - 1];
        }

        // Fritsch-Carlson monotonicity enforcement
        for (label i = 0; i < N; ++i)
        {
            if (Foam::mag(delta[i]) < VSMALL)
            {
                m[i] = m[i + 1] = 0.0;
                continue;
            }
            const scalar dMag  = Foam::mag(delta[i]) + VSMALL;
            const scalar alpha = Foam::mag(m[i])     / dMag;
            const scalar beta  = Foam::mag(m[i + 1]) / dMag;
            const scalar r2    = alpha * alpha + beta * beta;
            if (r2 > 9.0)
            {
                const scalar tau_fc = 3.0 / Foam::sqrt(r2);
                m[i]     *= tau_fc;
                m[i + 1] *= tau_fc;
            }
        }
        return m;
    };

// * * * * * * * * * * * * * *  Scalar slopes (Un)  * * * * * * * * * * * * * //
    phSlopesScalar[fI] = fcSlopes([&](label i) { return UnHistory[i][fI]; });

// * * * * * * * * * * * * *  Vector slopes (Lvec)  * * * * * * * * * * * * * //
    {
        List<vector>& mv = phSlopesVector[fI];
        mv.resize(nTimes, vector::zero);
        for (label comp = 0; comp < 3; ++comp)
        {
            scalarField mc = fcSlopes([&](label i) { return LvecHistory[i][fI][comp]; });
            for (label i = 0; i < nTimes; ++i)
                mv[i][comp] = mc[i];
        }
    }
}
