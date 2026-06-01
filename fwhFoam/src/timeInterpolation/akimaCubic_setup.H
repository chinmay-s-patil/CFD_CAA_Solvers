// akimaCubic_setup.H
// Precomputation step for Akima piecewise cubic interpolation.
//
// Called ONCE per face (fI) before the observer / reception-time loop.
// Fills akimaSlopesScalar[fI] and akimaSlopesVector[fI] with the Akima
// knot slopes m[i] for both history arrays.
//
// Complexity: O(N) per face.
//
// Expects in scope:
//   fI                   label
//   nTimes               label
//   timeValues           DynamicList<scalar>
//   LvecHistory          List<vectorField>   [tI][fI]
//   UnHistory            List<scalarField>   [tI][fI]
//   akimaSlopesScalar    List<scalarField>   [fI] output (sized nTimes)
//   akimaSlopesVector    List<List<vector>>  [fI] output (sized nTimes)

{
    const label N = nTimes - 1;   // number of sub-intervals

    // Helper: given a List of N+4 divided differences (offset by 2),
    // compute Akima slopes for nTimes knots.
    auto akimaSlopesFromDivDiff = [&](const scalarField& dArr) -> scalarField
    {
        scalarField m(nTimes, 0.0);
        for (label i = 0; i < nTimes; ++i)
        {
            const scalar dm2 = dArr[i];
            const scalar dm1 = dArr[i + 1];
            const scalar dp0 = dArr[i + 2];
            const scalar dp1 = dArr[i + 3];

            const scalar w1   = Foam::mag(dp1 - dp0);
            const scalar w2   = Foam::mag(dm1 - dm2);
            const scalar wSum = w1 + w2 + VSMALL;
            m[i] = (w1 * dm1 + w2 * dp0) / wSum;
        }
        return m;
    };

// * * * * * * * * * * * * * *  Scalar slopes (Un)  * * * * * * * * * * * * * //
    {
        scalarField dArr(N + 4, 0.0);
        for (label k = 0; k < N; ++k)
        {
            const scalar dt = timeValues[k + 1] - timeValues[k] + VSMALL;
            dArr[k + 2] = (UnHistory[k + 1][fI] - UnHistory[k][fI]) / dt;
        }
        dArr[1] = 2.0 * dArr[2] - dArr[3];
        dArr[0] = 2.0 * dArr[1] - dArr[2];
        dArr[N + 2] = 2.0 * dArr[N + 1] - dArr[N];
        dArr[N + 3] = 2.0 * dArr[N + 2] - dArr[N + 1];

        akimaSlopesScalar[fI] = akimaSlopesFromDivDiff(dArr);
    }

// * * * * * * * * * * * * *  Vector slopes (Lvec)  * * * * * * * * * * * * * //
    {
        List<vector>& mv = akimaSlopesVector[fI];
        mv.resize(nTimes, vector::zero);

        for (label comp = 0; comp < 3; ++comp)
        {
            scalarField dArr(N + 4, 0.0);
            for (label k = 0; k < N; ++k)
            {
                const scalar dt = timeValues[k + 1] - timeValues[k] + VSMALL;
                dArr[k + 2] = (LvecHistory[k + 1][fI][comp] - LvecHistory[k][fI][comp]) / dt;
            }
            dArr[1] = 2.0 * dArr[2] - dArr[3];
            dArr[0] = 2.0 * dArr[1] - dArr[2];
            dArr[N + 2] = 2.0 * dArr[N + 1] - dArr[N];
            dArr[N + 3] = 2.0 * dArr[N + 2] - dArr[N + 1];

            scalarField mc = akimaSlopesFromDivDiff(dArr);
            for (label i = 0; i < nTimes; ++i)
                mv[i][comp] = mc[i];
        }
    }
}
