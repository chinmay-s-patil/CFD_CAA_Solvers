# fwhCompressibleFoam — changes from fwhFoam (incompressible)

This file documents every deliberate difference between `fwhFoam` (incompressible)
and `fwhCompressibleFoam`.  Files that were copied unchanged are not listed.

---

## Why a separate solver?

In incompressible flow the density `rho` is constant everywhere, so it can be
read once from `fwhDict` as a scalar and used directly in the FW-H formula.

In compressible flow `rho` varies in space and time, which changes two things:

1. The **loading vector** `Lvec = p*n + rho*U*Un` must use the local density at
   each surface face at each time step, not a global constant.

2. The **thickness (monopole) source** is the time derivative of the surface
   mass flux `rho*Un`, not just `rho_ref * dUn/dt`.  The full form is:

   ```
   p'_T = (1/4pi) * integral[ d(rho*Un)/dt / r ] dS
   ```

   For constant rho this reduces to `rho * dUn/dt / r`, which is what
   fwhFoam uses.  fwhCompressibleFoam computes the product `rho*Un` first,
   then differentiates it with the same 3-point stencils.

---

## File-by-file changes

### `Make/files`

| Before | After |
|--------|-------|
| `EXE = $(FOAM_USER_APPBIN)/fwhFoam` | `EXE = $(FOAM_USER_APPBIN)/fwhCompressibleFoam` |

---

### `fwhFoam.C`

- Application name in the header changed to `fwhCompressibleFoam`.
- `postProcessing/fwhFoam/` output directory changed to
  `postProcessing/fwhCompressibleFoam/`.

---

### `createFields.H`

| Before | After |
|--------|-------|
| `const scalar rho = fwhSubDict.get<scalar>("rho");` | `const scalar rhoRef = fwhSubDict.get<scalar>("rhoRef");` |

`rhoRef` is still read from `fwhDict` and is used for `pRef`-based OASPL
calculations.  The per-face density at each time step now comes from
`rhoSurfHistory` which is populated by `readAndSample.H`.

---

### `readAndSample.H`

**New field sampled:** `rho` (volScalarField) is read alongside `p` and `U`
at every time step and interpolated onto the FW-H surface.

**New history list:**
```cpp
DynamicList<scalarField> rhoSurfHistory;   // [tI][fI]  local density
```

**New cache file:** `timeHistory_rho.csv` written and read back in the same
way as the existing pressure/velocity CSVs.  The cache-existence check now
requires all five files to be present.

**Cache consistency check** now also verifies that the `rho` table has the
same number of rows as `p`.

---

### `fwhIntegrate.H`

**New history lists:**
```cpp
List<scalarField> rhoHistory    (nTimes, scalarField(nFaces, 0.0));
List<scalarField> rhoUnHistory  (nTimes, scalarField(nFaces, 0.0));
List<scalarField> dRhoUnHistory (nTimes, scalarField(nFaces, 0.0));
```

**Loading vector** — uses local density instead of the constant:
```cpp
// Before (incompressible)
LvecHistory[tI][fI] = p_t * nHat[fI] + rho * U_t * Un;

// After (compressible)
const scalar rho_t = rhoSurfHistory[tI][fI];
LvecHistory[tI][fI] = p_t * nHat[fI] + rho_t * U_t * Un;
```

**rho*Un product** is built before differentiation:
```cpp
rhoUnHistory[tI][fI] = rhoHistory[tI][fI] * UnHistory[tI][fI];
```

**`applyStencil` lambda** now also fills `dRhoUnHistory` with the same
Lagrange coefficients used for `dUnHistory` and `dLvecHistory`.

**Thickness term in the observer loop** — uses `dRhoUn` instead of `rho * dUn`:
```cpp
// Before (incompressible)
const scalar pT = rho * dUn_tau / rMag / fourPi * dS;

// After (compressible)
const scalar dRhoUn = interp(dRhoUnHistory, tau, fI);
const scalar pT = dRhoUn / rMag / fourPi * dS;
```

---

### `writeVolumeFields.H`

The `fwhAtCell` lambda (used for the volumetric p' and OASPL fields) mirrors
the same change as the observer loop above — it uses `dRhoUnHistory` for the
thickness term instead of the constant `rho * dUnHistory`.

---

## fwhDict changes

Replace `rho` with `rhoRef`:

```
// before
rho     1.225;

// after
rhoRef  1.225;    // ambient/reference density for OASPL reference scaling
```

The solver now also expects a `rho` field to exist in each CFD time directory
(standard output of rhoPimpleFoam, rhoSimpleFoam, etc.).

---

## What did NOT change

- The loading (dipole) term `pL` is identical in both solvers.  The compressible
  change only affects the thickness (monopole) source.
- All interpolation schemes, Fourier transform options, OASPL/pPrime field
  outputs, and dictionary structure are unchanged.
- The incompressible solver `fwhFoam` is unmodified.
