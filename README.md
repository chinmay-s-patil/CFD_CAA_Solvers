<!-- =========================================================== -->
> [!WARNING]
> **Work in Progress:  Not Yet Validated**
>
> The implementation is functionally complete and builds cleanly against OpenFOAM v2412, but the results
> have **not yet been validated** against reference data or an experimental benchmark.
> The validation process is still underway, please do contact me if any issues are encountered or
> any sugesstions you may have.
<!-- =========================================================== -->

---

# CFD–CAA Solvers:  FW-H Aeroacoustics for OpenFOAM

A collection of aeroacoustic solvers and post-processing utilities built
around the **Farassat 1A** formulation of the Ffowcs Williams–Hawkings
(FW-H) equation.  Given a transient CFD run, these tools reconstruct the
far-field acoustic pressure at arbitrary observer (microphone) locations
without modifying the flow solver itself.

All solvers operate on OpenFOAM v2412 and share the same `system/fwhDict`
configuration file described below.

---

## Background:  what the FW-H equation does

The Ffowcs Williams–Hawkings equation is an exact rearrangement of the
Navier–Stokes equations into the form of an inhomogeneous wave equation.
It decomposes the far-field acoustic pressure into three physically distinct
contributions, each tied to a different source mechanism:

```
p'(x, t) = p'_T  +  p'_L  +  p'_Q
                     
          thickness  loading  quadrupole
         (monopole) (dipole) (volume)
```

- **Thickness (monopole)**:  fluid displaced by the surface moving through
  the medium.  For a stationary surface this is the time rate of change of
  the mass flux through each surface element.
- **Loading (dipole)**:  unsteady surface forces (pressure + viscous stress)
  acting on the fluid.  This is almost always the dominant source for
  lifting surfaces, fans, and rotors.
- **Quadrupole (volume)**:  turbulent Reynolds stresses and Lighthill's
  stress tensor throughout the volume.  Significant only at high Mach
  numbers (M ≳ 0.3) and is expensive to compute; in practice it is
  neglected for low-Mach incompressible flows.

The Farassat 1A form evaluates the surface integrals in the *time domain*
using retarded-time interpolation, which is what all three solvers here
implement.

---

## Solvers at a Glance

Three executables are provided.  They share the physics kernel but differ
in how and when the CFD data is produced.

### `fwhFoam`:  post-processing utility (incompressible)

| Attribute | Detail |
|-----------|--------|
| Type | Standalone post-processor; run after the CFD solver finishes |
| Flow regime | Incompressible (constant `rho` from `fwhDict`) |
| Thickness term | Implemented:  monopole:  `d(rho·Un)/dt` with constant `rho` |
| Loading term | Implemented:  dipole:  `p·n̂ + ρ·U·Uₙ` |
| Quadrupole | Not implemented (valid for M < 0.3) |
| Viscous stress in loading | Implemented via the permeable-surface form of Farassat 1A |

`fwhFoam` reads the `p` and `U` surface history written by any OpenFOAM
flow solver, reconstructs the retarded-time surface integrals, and writes
the acoustic signals to `postProcessing/fwhFoam/`.  Because it is a
separate executable, no changes to `system/controlDict` are needed and the
acoustic post-processing can be re-run with different observer positions or
`pRef` values without repeating the CFD.

### `fwhCompressibleFoam`:  post-processing utility (compressible)

| Attribute | Detail |
|-----------|--------|
| Type | Standalone post-processor; run after the CFD solver finishes |
| Flow regime | Compressible:  `rho` varies in space and time |
| Thickness term | Implemented:  monopole:  `d(rho·Un)/dt` with *local* `rho` per face |
| Loading term | Implemented:  dipole:  `p·n̂ + ρ_local·U·Uₙ` |
| Quadrupole | Not Implemented |
| Viscous stress in loading | Implemented via the permeable-surface form of Farassat 1A |

The key difference from `fwhFoam` is that density is not treated as a
constant.  The solver additionally reads the `rho` field from the CFD time
directories, samples it onto the FW-H surface alongside `p` and `U`, and
uses the product `rho·Un` (differentiated together) for the thickness term.
This is the correct form for flows where density fluctuations are
non-negligible:  outputs of `rhoPimpleFoam`, `rhoSimpleFoam`, etc.

The `fwhDict` entry changes from `rho` to `rhoRef`; everything else is
identical to `fwhFoam`.

### `pimpleFWHFoam`:  coupled flow + acoustics solver

| Attribute | Detail |
|-----------|--------|
| Type | Full flow solver; replaces `pimpleFoam` in the case |
| Base solver | `pimpleFoam`:  incompressible, turbulent, moving-mesh PIMPLE |
| Flow regime | Incompressible (constant `rho`) |
| Thickness term | Implemented:  monopole |
| Loading term | Implemented: dipole |
| Quadrupole | Not Implemented |
| Viscous stress in loading | Implemented |
| Sampling | Inline:  surface data is collected every time step during the run |
| Post-processing | Runs automatically at the end of the time loop |

`pimpleFWHFoam` is `pimpleFoam` with the FW-H machinery wired directly
into the time loop.  The surface is sampled after every PIMPLE correction
via `fwhSample.H`, and the full acoustic integration is triggered once the
simulation ends via `fwhPostProcess.H`.  This avoids writing large volumes
of surface field data to disk while still capturing the complete time
history.  Results go to `postProcessing/pimpleFWHFoam/`.

---

## Build

Each solver is built independently.

```bash
cd fwhFoam/src && wmake
cd fwhCompressibleFoam/src && wmake
cd pimpleFWHFoam/src && wmake
```

A sourced OpenFOAM v2412 (or compatible v2x) environment is required.
Binaries are installed to `$FOAM_USER_APPBIN`.  FFTW3 must be available
(`libfftw3-dev` on Ubuntu).

---

## Configuration:  `system/fwhDict`

All solvers read `system/fwhDict`.  The full structure is shown below with
every available entry.

```c
FWH
{
    // Acoustic constants 
    rho             1.225;      // [kg/m³] fluid density (fwhFoam / pimpleFWHFoam)
    // rhoRef       1.225;      // use rhoRef instead for fwhCompressibleFoam
    c0              343.0;      // [m/s]   speed of sound in the far field
    pRef            2e-5;       // [Pa]    reference pressure for SPL (20 µPa for air)

    // Time window 
    startTime       0.1;        // skip initial transient [s]
    endTime         500.0;      // last time step to include [s]
    timeStepStride  1;          // use every N-th available time step

    // Retarded-time interpolation 
    // Options: linear | nearest | hermiteCubic | spline
    //          akimaCubic | piecewiseHermite
    timeInterpolationScheme  linear;

    // Observer (microphone) positions 
    observers
    (
        (0.5  0    0)
        (0.5  0.5  0)
        (10   0    0)
        (20   0    0)
    );

    // Output toggles 
    writeSignal             yes;    // acoustic pressure time signal: pTime.dat
    writeSPL                yes;    // Sound Pressure Level spectrum: spl.dat
    writeTimeHistory        no;     // cache CSV files to disk
    recomputeTimeHistory    yes;    // reuse cache if it already exists

    // OASPL directivity 
    writeOASPL              yes;

    OASPL
    {
        centre          (0 0 0);        // arc reference point
        refDirection    (1 0 0);        // θ = 0° direction
        fileName        oaspl_directivity;
    }

    // Volumetric p' field (p' at every mesh cell, per time step) 
    writeVolumePPrime       yes;

    pPrime
    {
        surfaceName     FWHSurface;     // must match a name in sampling.surfaces{}
        fieldName       pPrime;
        writeInterval   1;              // write every N collected time steps
        writeField      yes;            // OpenFOAM volScalarField (loadable in ParaView)
        writeVTK        no;             // VTK point cloud .vtk
        writeRaw        no;             // flat binary + .info metadata
    }

    // Volumetric OASPL field (RMS p': single snapshot) 
    writeVolumeOASPL        no;

    // volumeOASPL
    // {
    //    surfaceName     FWHSurface;
    //    time            (4.0 8.0);  // (minT maxT) or single scalar
    //    fieldName       OASPL;
    //    writeField      yes;
    // }

    // Fourier transform 
    fourierTransform
    {
        writeFT         yes;

        // Options: DFT | FFT | FFTW | NUDFT
        //  DFT  :  direct O(N²) transform, correct for any spacing
        //  FFT  :  Cooley–Tukey, fast, assumes uniform sampling
        //  FFTW :  FFT via FFTW3 library with windowing / detrending
        //  NUDFT:  non-uniform DFT, correct for variable time steps 
        //          (Particularly required for pimpleFWHFoam)
        type            FFTW;

        // Entries used only by FFTW:
        windowType      Hanning;        // Rectangular | Hann | Hamming | Bartlett | 
                                        // Blackman | FlatTop | Tukey
        detrend         yes;
        removeMean      yes;
        nfft            auto;           // or an explicit integer
        normalize       yes;
        normalizeWith   amplitude;      // amplitude | power
    }

    // FW-H sampling surface 
    // The sampling sub-dict is forwarded to OpenFOAM's sampledSurfaces,
    // so any surface type that sampledSurfaces supports works here.
    sampling
    {
        // Informational only (printed at startup)
        surfaceName             FWHSurface;
        surfaceType             cylinder;

        faceInterpolationScheme cellPoint;
        fields                  (p U);     // must include p and U unless renamed by external solver.

        surfaces                            // Multiple surfaces can be defined.
        {
            FWHSurface
            {
                type    patch;
                patches (cylinder);
            }
        }
    }
}
```

---

## Outputs

### Always written

These files are produced unconditionally for every observer and every
sampling surface.

```
postProcessing/
└ fwhFoam/                           or pimpleFWHFoam/ or fwhCompressibleFoam/
    └ <surfaceName>/
        └ observer<N>/
            ├ pTime.dat              acoustic pressure p'(t) vs time [s, Pa]
            ├ ft.dat                 one-sided magnitude spectrum |P(f)| [Hz, Pa]
            └ spl.dat                SPL in dB re pRef vs frequency [Hz, dB]
```

| File | Content |
|------|---------|
| `pTime.dat` | Raw acoustic pressure time signal, columns: `time [s]`, `p' [Pa]` |
| `ft.dat` | One-sided Fourier magnitude spectrum, columns: `f [Hz]`, `\|P(f)\| [Pa]` |
| `spl.dat` | Sound Pressure Level, columns: `f [Hz]`, `SPL [dB]` |

### Optional:  switched on in `fwhDict`

| Switch | File(s) written | Notes |
|--------|-----------------|-------|
| `writeOASPL yes` | `<surfaceName>/oaspl_directivity.dat` | `[θ°, OASPL dB]` per observer; useful for polar directivity plots |
| `writeTimeHistory yes` | `timeHistory_p/Ux/Uy/Uz.csv` (and `_rho.csv` for the compressible solver) | CSV cache of the sampled surface history; avoids re-reading the mesh on subsequent runs when `recomputeTimeHistory no` |
| `writeVolumePPrime yes` | `postProcessing/…/pPrime/<time>/` | `p'` as a `volScalarField` computed at every mesh cell; loads directly in ParaView via the `.foam` file |
| `writeVolumeOASPL yes` | `postProcessing/…/OASPL/` | RMS of `p'` over the specified time window written as a single `volScalarField` snapshot |

---

## References

Farassat, F. (2007). *Derivation of Formulations 1 and 1A of Farassat.*
NASA/TM-2007-214853.

Ffowcs Williams, J. E. & Hawkings, D. L. (1969). Sound generation by
turbulence and surfaces in arbitrary motion. *Philosophical Transactions of
the Royal Society of London A*, 264, 321–342.

---

# Contact Information

E-Mail: <patil.chinmay3031@gmail.com>

GitHub: <https://github.com/chinmay-s-patil>
