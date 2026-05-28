// spline_setup.H
// Precomputation step for natural cubic spline interpolation.
//
// Called ONCE per face (fI) before the observer / reception-time loop.
// Fills splMomentsScalar[fI] and splMomentsVector[fI] with the N
// second-derivative moments M[i] for the two history arrays.
//
// Complexity: O(N) tridiagonal solve per face.
//
// Expects in scope:
//   fI                label           face being precomputed
//   nTimes            label
//   timeValues        DynamicList<scalar>
//   LvecHistory       List<vectorField>   [tI][fI]
//   UnHistory         List<scalarField>   [tI][fI]
//   splMomentsScalar  List<scalarField>   [fI] output (sized nTimes)
//   splMomentsVector  List<List<vector>>  [fI] output (sized nTimes)

{
    // Step widths
    scalarField hSpan(nTimes - 1);
    for (label i = 0; i < nTimes - 1; ++i)
        hSpan[i] = timeValues[i + 1] - timeValues[i] + VSMALL;

    // Helper: solve for moments given a RHS vector (Thomas sweep, natural BCs)
    // Returns M[0..nTimes-1].
    auto solveSplineMomentsScalar = [&](const scalarField& rhs) -> scalarField
    {
        scalarField diagMod(nTimes, 0.0);
        scalarField rhsMod (nTimes, 0.0);

        diagMod[0] = 1.0;
        rhsMod [0] = 0.0;

        for (label i = 1; i < nTimes - 1; ++i)
        {
            const scalar diag   = 2.0 * (hSpan[i - 1] + hSpan[i]);
            const scalar factor = hSpan[i - 1] / (diagMod[i - 1] + VSMALL);
            diagMod[i] = diag - factor * hSpan[i - 1];
            rhsMod [i] = rhs[i] - factor * rhsMod[i - 1];
        }
        diagMod[nTimes - 1] = 1.0;
        rhsMod [nTimes - 1] = 0.0;

        scalarField M(nTimes, 0.0);
        for (label i = nTimes - 2; i >= 1; --i)
            M[i] = (rhsMod[i] - hSpan[i] * M[i + 1]) / (diagMod[i] + VSMALL);
        return M;
    };

// * * * * * * * * * * * * * *  Scalar moments (Un)  * * * * * * * * * * * * * //
    {
        scalarField rhs(nTimes, 0.0);
        for (label i = 1; i < nTimes - 1; ++i)
        {
            rhs[i] =   6.0 * (UnHistory[i + 1][fI] - UnHistory[i    ][fI]) / hSpan[i]
                     - 6.0 * (UnHistory[i    ][fI] - UnHistory[i - 1][fI]) / hSpan[i - 1];
        }
        splMomentsScalar[fI] = solveSplineMomentsScalar(rhs);
    }

// * * * * * * * * * * * * *  Vector moments (Lvec)  * * * * * * * * * * * * * //
    //
    // The tridiagonal matrix is identical for all components, so we solve
    // three scalar systems (x, y, z) and repack into vector moments.
    {
        List<vector>& Mv = splMomentsVector[fI];
        Mv.resize(nTimes, vector::zero);

        for (label comp = 0; comp < 3; ++comp)
        {
            scalarField rhs(nTimes, 0.0);
            for (label i = 1; i < nTimes - 1; ++i)
            {
                rhs[i] =
                      6.0 * (LvecHistory[i + 1][fI][comp] - LvecHistory[i    ][fI][comp]) / hSpan[i]
                    - 6.0 * (LvecHistory[i    ][fI][comp] - LvecHistory[i - 1][fI][comp]) / hSpan[i - 1];
            }
            scalarField Mc = solveSplineMomentsScalar(rhs);
            for (label i = 0; i < nTimes; ++i)
                Mv[i][comp] = Mc[i];
        }
    }
}
