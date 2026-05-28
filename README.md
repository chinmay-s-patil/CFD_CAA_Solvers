# <p align="center">fwhFoam</p>

An OpenFOAM post-processing utility for far-field aeroacoustics using the
**Farassat 1A** formulation of the Ffowcs Williams–Hawkings (FW-H) acoustic
analogy.  It reads transient CFD surface data produced by `pimpleFoam` (or any
incompressible/compressible solver) and reconstructs the acoustic pressure
signal at arbitrary observer (microphone) locations — no modifications to the
flow solver are required.

---

## Physics

The acoustic pressure at observer **x** is split into thickness (monopole) and
loading (dipole) terms:

```
p'(x,t) = p'_T + p'_L
```

Both are evaluated on a stationary permeable control surface enclosing all noise
sources.  Sound is assumed to propagate through a quiescent medium at speed `c0`
(free-stream), making this approach ideal for low-Mach flows (M < 0.3: drones,
rotors, fans, airfoils).  The quadrupole (volume) term is neglected, which is
standard practice for incompressible CFD.

Retarded-time interpolation and nonuniform quadratic time-derivative stencils
(Farassat 1A) are used throughout, giving second-order accuracy even on
non-uniform time sequences.

---

## Workflow

```
1.  pimpleFoam          ← transient CFD (LES / DES / URANS)
         |
         └─ samples p, U on FW-H surface every timestep
                (via sampledSurfaces in fwhDict)
2.  fwhFoam             ← reads surface history, computes acoustics
         |
         ├─ postProcessing/fwhFoam/observer0/pTime.dat   (time signal)
         ├─ postProcessing/fwhFoam/observer0/fft.dat     (frequency spectrum)
         └─ postProcessing/fwhFoam/observer0/spl.dat     (SPL in dB)
```

---

## Build

```bash
cd src
wmake
```

Requires a sourced OpenFOAM v2412 (or compatible) environment.  The binary is
installed to `$FOAM_USER_APPBIN`.

---

## Configuration — `system/fwhDict`

```c
FWH
{
    rho             1.225;      // kg/m³  (free-stream density)
    c0              343;        // m/s    (speed of sound)
    pRef            2e-5;       // Pa     (reference pressure for SPL)

    startTime       0.1;        // skip initial transient (Can be set to 0, constant is ignored in code)
    endTime         1.0;
    timeStepStride  1;          // use every Nth timestep

    observers               // microphone positions (m)
    (
        (10  0  0)
        (20  0  0)
    );

    recomputeTimeHistory    yes;  // yes: re-sample; no: use cached CSVs
    writeTimeHistory        yes;  // write CSV cache to postProcessing/fwhFoam/

    sampling
    {
        surfaceName         FWHSurface;     // for info / printing only
        surfaceType         sphere;         // for info / printing only
        interpolationScheme cellPoint;
        surfaceFormat       vtk;
        fields              (p U);

        surfaces
        {
            FWHSurface
            {
                type    patch;
                patches (myFWHPatch);   // or use a sphere / sampledSurface
            }
        }
    }

    writeSignal yes;    // p'(t) vs t
    writeFFT    yes;    // |P(f)| vs f
    writeSPL    yes;    // SPL (dB) vs f
}
```

The `sampling` sub-dictionary is passed directly to OpenFOAM's `sampledSurfaces`
class, so any surface type supported there (sphere, patch, plane, …) works.

---

## Outputs

All results are written to `postProcessing/fwhFoam/observer<N>/`.

| File          | Content                                             |
| ------------- | --------------------------------------------------- |
| `pTime.dat` | Acoustic pressure `p'(t)`vs time                  |
| `fft.dat`   | One-sided magnitude spectrum `\|P(f)\|`vs frequency |
| `spl.dat`   | Sound Pressure Level (dB re 20 µPa) vs frequency   |

Time-history CSV cache (`timeHistory_p/Ux/Uy/Uz.csv`) is written alongside and
reused on subsequent runs when `recomputeTimeHistory no` is set — useful for
re-running the acoustic integration with different observers or `pRef` without
repeating expensive mesh I/O.

---

## Tutorials

| Case                             | Description                                                                                                            |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `tutorial/vortexShred`         | 2-D cylinder in crossflow (Re=100, LES). FW-H surface is the cylinder wall patch. Expected shedding peak at ≈ 0.1 Hz. |

Run a tutorial:

```bash
cd tutorial/vortexShred
./Allrun
fwhFoam
```

---

## Design Notes

* **No solver modification** — `fwhFoam` is a standalone post-processing
  utility.  The CFD solver writes surface data; `fwhFoam` reads it afterwards.
  This is the correct separation of concerns for acoustic analogies.
* **Parallel-safe** — surface fields are gathered across MPI ranks with
  `Pstream::combineReduce` before integration.
* **Stationary surface, stationary medium** — moving-surface (rotating rotor)
  and convective-medium corrections are not implemented.
* **FFT** — a direct O(N²) DFT is used, sufficient for typical post-processing
  sample counts.  Replace with FFTW for large datasets.

---

## References

Farassat, F. (2007). *Derivation of Formulations 1 and 1A of Farassat.*
NASA/TM-2007-214853.
