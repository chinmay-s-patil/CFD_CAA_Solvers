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

# CFD–CAA Solvers: FW-H Aeroacoustics for OpenFOAM

A suite of high-performance aeroacoustic solvers and post-processing tools for OpenFOAM (v2412+), implementing **Farassat Formulation 1A** of the **Ffowcs Williams–Hawkings (FW-H)** acoustic analogy for permeable and solid surfaces.

Given transient CFD flow data, these tools predict far-field acoustic signals, compute Sound Pressure Level (SPL) spectra, generate directivity maps, perform Zoom FFT analysis, and reconstruct full 3D volumetric acoustic pressure fields `p'(x, t)` and OASPL distributions across the mesh.

---

## Technical Overview & Governing Physics

The Ffowcs Williams–Hawkings (FW-H) equation is an exact rearrangement of the compressible Navier–Stokes equations into an inhomogeneous wave equation. For a stationary surface `S` (permeable or solid) defined by `f(x, t) = 0`, the acoustic pressure `p'(x, t)` at an observer position `x` and reception time `t` is decomposed into three source terms:

```
p'(x, t) = p'_T(x, t) + p'_L(x, t) + p'_Q(x, t)

Where:
  p'_T  =  Thickness (Monopole term: mass displacement)
  p'_L  =  Loading (Dipole term: surface pressure & momentum flux)
  p'_Q  =  Quadrupole (Volume term: turbulent Reynolds stresses)
```

```
                            ┌──────────────────────────────────────────────┐
                            │      Far-Field Acoustic Pressure p'(x,t)     │
                            └──────────────────────┬───────────────────────┘
                                                   │
                ┌──────────────────────────────────┼──────────────────────────────────┐
                │                                  │                                  │
   ┌────────────┴────────────┐        ┌────────────┴────────────┐        ┌────────────┴────────────┐
   │  Thickness (Monopole)   │        │     Loading (Dipole)    │        │   Quadrupole (Volume)   │
   │    p'_T(x,t) (Mass)     │        │   p'_L(x,t) (Force/Press) │        │  p'_Q(x,t) (Turbulence) │
   └─────────────────────────┘        └─────────────────────────┘        └─────────────────────────┘
```

### Farassat Formulation 1A Integrals

All solvers in this suite implement Farassat Formulation 1A evaluated in the time domain using retarded-time (`τ`) integration:

```
τ = t - r / c₀,    r = |x - y|,    r̂ = (x - y) / r
```

1. **Thickness (Monopole) Term `p'_T`**:
   - **Incompressible Flow (`fwhFoam` & `pimpleFWHFoam`)**: Density `ρ₀` is constant.
     ```
     p'_T(x, t) = (1 / 4π) · ∫ [ (ρ₀ · dUn/dt) / r ]_τ dS
     ```
   - **Compressible Flow (`fwhCompressibleFoam`)**: Density varies in space and time `ρ(y, τ)`.
     ```
     p'_T(x, t) = (1 / 4π) · ∫ [ (d(ρ·Un)/dt) / r ]_τ dS
     ```

2. **Loading (Dipole) Term `p'_L`**:
   Uses the permeable-surface formulation with local force/momentum vector `L = p·n̂ + ρ·U·Un`:
   ```
   p'_L(x, t) = (1 / 4π) · ∫ [ (dL/dt · r̂) / (c₀ · r) + (L · r̂) / r² ]_τ dS
   ```

3. **Quadrupole (Volume) Term `p'_Q`**:
   Neglected under the low-to-moderate Mach number assumption (`M < 0.3`), standard for permeable FW-H surfaces enclosing the non-linear vorticity region.

---

## Solvers Matrix

Three solver executables are included in the `solvers/` directory:

| Attribute | `fwhFoam` | `fwhCompressibleFoam` | `pimpleFWHFoam` |
| :--- | :--- | :--- | :--- |
| **Executable Type** | Standalone post-processor | Standalone post-processor | Coupled flow + acoustics solver |
| **Base CFD Solver** | Any OpenFOAM solver | Any OpenFOAM solver | `pimpleFoam` (PIMPLE / DyM) |
| **Flow Regime** | Incompressible | Compressible | Incompressible |
| **Density Handling** | Constant `ρ₀` (from `fwhDict`) | Time-varying per-face `ρ(x,t)` | Constant `ρ₀` (from `fwhDict`) |
| **Thickness Term** | `ρ₀ · dUn/dt` | `d(ρ·Un)/dt` | `ρ₀ · dUn/dt` |
| **Loading Term** | `p·n̂ + ρ₀·U·Un` | `p·n̂ + ρ(x,t)·U·Un` | `p·n̂ + ρ₀·U·Un` |
| **Data Sampling** | Reads written CFD time dirs | Reads written CFD time dirs (`p, U, rho`) | On-the-fly inside time loop |
| **Disk Storage** | High (surface fields on disk) | High (surface fields on disk) | Zero (samples in memory) |

---

## Solver Descriptions

### 1. `solvers/fwhFoam` — Incompressible Post-Processor
- Reads transient `p` and `U` fields from saved OpenFOAM CFD time directories.
- Evaluates surface integrals using constant density `ρ`.
- Writes acoustic time signals and spectra to `postProcessing/fwhFoam/`.
- Allows quick re-evaluation with modified observer arrays, reference pressure, or interpolation schemes without re-running CFD.

### 2. `solvers/fwhCompressibleFoam` — Compressible Post-Processor
- Designed for compressible flow runs (e.g. `rhoPimpleFoam`, `rhoSimpleFoam`, `dbnsFoam`).
- Reads `p`, `U`, and `rho` fields at every sampled time step.
- Solves the full compressible thickness formulation by differentiating `d(ρ·Un)/dt` at each face.
- Evaluates loading vectors using exact local surface density `ρ(y, τ)`.
- Caches `rho` time history alongside pressure and velocity (`timeHistory_rho.csv`).

### 3. `solvers/pimpleFWHFoam` — Coupled Flow & Aeroacoustics Solver
- Combines OpenFOAM's `pimpleFoam` (transient, incompressible, turbulent, moving-mesh solver) directly with the FW-H engine.
- Performs on-the-fly surface sampling after each PIMPLE corrector step (`fwhSample.H`).
- Triggers acoustic Farassat 1A integration automatically at simulation end (`fwhPostProcess.H`).
- Eliminates the need to save massive surface field datasets to disk.

---

## Advanced Numerical Features

### Retarded-Time Interpolation Engine
To evaluate source quantities at emission time `τ = t - r/c₀`, six high-order time interpolation schemes are available (`timeInterpolationScheme` in `fwhDict`):
- `linear`: Piecewise linear interpolation.
- `nearest`: Nearest neighbor interpolation.
- `hermiteCubic`: Standard cubic Hermite spline interpolation.
- `spline`: Natural cubic spline with precomputed 2nd-derivative moments per face.
- `akimaCubic`: Akima cubic interpolation, designed to suppress non-physical overshoots and oscillations near sharp flow transients.
- `piecewiseHermite`: PCHIP (Piecewise Cubic Hermite Interpolating Polynomial), shape-preserving monotonic interpolation.

### Non-Uniform Stencil Differentiation
Time derivatives (`dL/dt` and `dUn/dt`) are precomputed using 3-point non-uniform Lagrange polynomial stencils, exact for polynomials up to 2nd order:
- **Interior Nodes**: Centred 3-point stencil using `(t_{i-1}, t_i, t_{i+1})`.
- **Start Boundary (`t₀`)**: Forward 3-point one-sided stencil using `(t₀, t₁, t₂)`.
- **End Boundary (`t_N`)**: Backward 3-point one-sided stencil using `(t_{N-2}, t_{N-1}, t_N)`.

### Spectral & Fourier Analysis Suite
- `DFT`: Direct O(N²) Discrete Fourier Transform. Exact for non-uniform time step distributions.
- `FFT`: Cooley–Tukey Fast Fourier Transform (assumes uniform time sampling).
- `FFTW`: Link to the FFTW3 library with zero-padding (`nfft`), mean removal, linear detrending, and amplitude or power normalization.
- `NUDFT`: Non-Uniform Discrete Fourier Transform, essential for adaptive time-stepping simulations (`pimpleFWHFoam`).

### Zoom FFT Analysis (`zoomPostProcess yes`)
Provides zoomed-in frequency band analysis over `[fStart, fEnd]` with user-specified resolution lines (`nLines`):
- `CZT`: Chirp Z-Transform for high-resolution spectral zoom without padding.
- `heterodyne`: Heterodyne Zoom FFT via complex digital down-conversion, low-pass filtering, decimation, and standard FFT.

### Windowing Functions
Supported windowing types for spectral processing: `Rectangular`, `Hann` (or `Hanning`), `Hamming`, `Blackman`, `FlatTop`, `Bartlett`, `Tukey`.

### 3D Volumetric Acoustic Field Reconstruction
- **Volumetric p' Field (`writeVolumePPrime yes`)**: Maps the acoustic pressure `p'(x, t)` to every cell center across the entire 3D computational domain at specified intervals. Output as an OpenFOAM `volScalarField` in `postProcessing/.../pPrime/<time>/` (viewable in ParaView).
- **Volumetric OASPL Field (`writeVolumeOASPL yes`)**: Computes the RMS acoustic pressure over a specified time window `[minT, maxT]` and converts it to a 3D spatial field of OASPL (dB).

### Polar OASPL Directivity Analysis (`writeOASPL yes`)
Calculates OASPL directivity across observer points relative to a defined acoustic center (`centre`) and reference axis (`refDirection`), writing polar angular distributions.

### CSV Surface History Caching
When `writeTimeHistory yes` is set, sampled surface histories (`p`, `U_x`, `U_y`, `U_z`, and `rho`) are saved to CSV files in `postProcessing/.../<surfaceName>/`. Setting `recomputeTimeHistory no` skips mesh re-reading on subsequent runs.

---

## Directory Structure (Solvers Only)

```
solvers/
├── fwhFoam/                        # Incompressible post-processing solver
│   ├── fwhFoam.C                   # Main program entry
│   ├── createFields.H              # Dictionary reading & initialization
│   ├── dictErrors.H                # Configuration validation & error checking
│   ├── readAndSample.H             # Mesh & surface field sampling / CSV caching
│   ├── fwhIntegrate.H              # Farassat 1A integration kernel
│   ├── writeOutputs.H              # Output writer (pTime, ft, spl, OASPL)
│   ├── writeVolumeFields.H         # 3D volScalarField reconstruction (pPrime, OASPL)
│   ├── writeZoomOutputs.H          # Zoom FFT output writer (zoomFT, zoomSPL)
│   ├── printFWHSettings.H          # Console reporting helpers
│   ├── printSamplingSettings.H
│   ├── printVolumeFieldSettings.H
│   ├── Make/                       # wmake build specification
│   ├── fourierTransform/           # Spectral analysis header modules
│   │   ├── czt.H                   # Chirp Z-Transform algorithm
│   │   ├── dft.H                   # Discrete Fourier Transform
│   │   ├── fft.H                   # Cooley-Tukey FFT
│   │   ├── fftw.H                  # FFTW3 library interface
│   │   ├── heterodyne.H           # Quadrature heterodyne zoom FFT
│   │   └── nudft.H                 # Non-Uniform DFT
│   ├── timeInterpolation/          # Retarded-time interpolation modules
│   │   ├── linear.H
│   │   ├── nearest.H
│   │   ├── hermiteCubic.H
│   │   ├── spline.H / spline_setup.H
│   │   ├── akimaCubic.H / akimaCubic_setup.H
│   │   └── piecewiseHermite.H / piecewiseHermite_setup.H
│   ├── windows/                    # Spectral windowing definitions
│   │   ├── rectangular.H, hann.H, hamming.H, blackman.H, flattop.H, bartlett.H, tukey.H
│   └── explaination/               # Documentation & theory source
│
├── fwhCompressibleFoam/            # Compressible post-processing solver
│   ├── fwhFoam.C                   # Main program entry for compressible solver
│   ├── CHANGES_from_incompressible.md # Comprehensive differences log vs fwhFoam
│   ├── createFields.H              # Reads rhoRef instead of constant rho
│   ├── readAndSample.H             # Samples p, U, and rho; manages timeHistory_rho.csv
│   ├── fwhIntegrate.H              # Compressible d(rho*Un)/dt & local rho Lvec kernel
│   ├── writeVolumeFields.H         # Compressible 3D field reconstruction
│   └── ... (shares fourierTransform, timeInterpolation, windows modules)
│
└── pimpleFWHFoam/                  # Coupled flow + acoustics solver
    ├── pimpleFWHFoam.C             # PIMPLE loop integrated with FW-H sampling
    ├── createFWHFields.H           # FW-H initialization for coupled run
    ├── fwhSample.H                 # In-memory time-step surface sampling
    ├── fwhPostProcess.H            # End-of-run acoustic post-processing trigger
    └── ... (shares integration, spectral, and interpolation modules)
```

---

## Build Instructions

### Prerequisites
1. **OpenFOAM Environment**: OpenFOAM v2412 (or compatible OpenFOAM-v2x release) must be sourced.
2. **FFTW3 Library**: `libfftw3-dev` (Ubuntu/Debian) or equivalent development library installed.

### Compilation
Build each solver independently using OpenFOAM's `wmake` utility:

```bash
# Build Incompressible Post-Processor
cd solvers/fwhFoam
wmake

# Build Compressible Post-Processor
cd ../fwhCompressibleFoam
wmake

# Build Coupled PIMPLE Flow + Acoustics Solver
cd ../pimpleFWHFoam
wmake
```

Compiled binaries are installed directly to `$FOAM_USER_APPBIN`.

---

## Configuration Reference: `system/fwhDict`

All three solvers read the central configuration file located at `system/fwhDict`. Below is an annotated example illustrating all supported options:

```cpp
FWH
{
    // =========================================================================
    //  Acoustic Physical Constants
    // =========================================================================
    rho             1.225;      // Ambient fluid density [kg/m³] (fwhFoam & pimpleFWHFoam)
    // rhoRef       1.225;      // Ambient reference density [kg/m³] (fwhCompressibleFoam)
    c0              343.0;      // Speed of sound in ambient medium [m/s]
    pRef            2e-5;       // Reference pressure for SPL calculation [Pa] (20 µPa for air)

    // =========================================================================
    //  Time Window Controls
    // =========================================================================
    startTime       0.1;        // Start time for acoustic analysis [s] (skips transient)
    endTime         10.0;       // End time for acoustic analysis [s]
    timeStepStride  1;          // Stride ratio (use every N-th time step)

    // =========================================================================
    //  Retarded-Time Interpolation Scheme
    //  Options: linear | nearest | hermiteCubic | spline | akimaCubic | piecewiseHermite
    // =========================================================================
    timeInterpolationScheme  akimaCubic;

    // =========================================================================
    //  Observer Positions (Microphone Coordinates in meters)
    // =========================================================================
    observers
    (
        (1.0  0.0  0.0)
        (0.0  1.0  0.0)
        (5.0  0.0  0.0)
        (10.0 0.0  0.0)
    );

    // =========================================================================
    //  Output Switches & Cache Controls
    // =========================================================================
    writeSignal             yes;    // Write acoustic pressure time signal (pTime.dat)
    writeSPL                yes;    // Write Sound Pressure Level spectrum (spl.dat)
    writeTimeHistory        yes;   // Save sampled surface histories as CSV
    recomputeTimeHistory    no;    // Reuse cached CSV history if present

    // =========================================================================
    //  OASPL Directivity Map (Optional)
    // =========================================================================
    writeOASPL              yes;

    OASPL
    {
        centre          (0 0 0);          // Arc reference origin
        refDirection    (1 0 0);          // Axis corresponding to theta = 0 deg
        fileName        oaspl_directivity;// Output filename prefix
    }

    // =========================================================================
    //  3D Volumetric Transient Acoustic Pressure Field (Optional)
    // =========================================================================
    writeVolumePPrime       yes;

    pPrime
    {
        surfaceName     FWHSurface;       // Target surface defined in sampling.surfaces
        fieldName       pPrime;           // Output field name
        writeInterval   10;               // Write every N time steps
        writeField      yes;              // Output OpenFOAM volScalarField
        writeVTK        no;               // Output VTK point cloud
        writeRaw        no;               // Output flat binary raw format
    }

    // =========================================================================
    //  3D Volumetric OASPL Field Snapshot (Optional)
    // =========================================================================
    writeVolumeOASPL        no;

    // volumeOASPL
    // {
    //     surfaceName     FWHSurface;
    //     time            (1.0 5.0);        // Time interval (minT maxT) or scalar end time
    //     fieldName       OASPL;
    //     writeField      yes;
    // }

    // =========================================================================
    //  Fourier Transform Configuration
    // =========================================================================
    fourierTransform
    {
        writeFT         yes;

        // Transform algorithm: DFT | FFT | FFTW | NUDFT
        type            FFTW;

        // FFTW specific parameters:
        // windowType: Rectangular | Hann | Hamming | Bartlett | Blackman | FlatTop | Tukey
        windowType      Hann;
        detrend         yes;              // Remove linear trend from signal
        removeMean      yes;              // Subtract mean pressure
        nfft            auto;             // FFT size (auto or explicit power of 2)
        normalize       yes;              // Apply window energy normalization
        normalizeWith   amplitude;        // amplitude | power
    }

    // =========================================================================
    //  Zoom FFT Configuration (Optional)
    // =========================================================================
    zoomPostProcess     yes;

    zoom
    {
        type            CZT;              // CZT | heterodyne
        fStart          100.0;            // Band start frequency [Hz]
        fEnd            1000.0;           // Band end frequency [Hz]
        nLines          500;              // Number of spectral output lines
        windowType      Hann;             // Windowing function
        writeZoomFT     yes;              // Write zoomFT.dat
        writeZoomSPL    yes;              // Write zoomSPL.dat
    }

    // =========================================================================
    //  FW-H Sampling Surface Definitions
    // =========================================================================
    sampling
    {
        surfaceName             FWHSurface;
        surfaceType             patch;
        faceInterpolationScheme cellPoint;
        fields                  (p U);

        surfaces
        {
            FWHSurface
            {
                type    patch;
                patches (cylinderSurface);
            }
        }
    }
}
```

---

## Output Data Structure

All output files are written under the case `postProcessing/` directory:

```
postProcessing/
└── <solverName>/                      # fwhFoam, fwhCompressibleFoam, or pimpleFWHFoam
    └── <surfaceName>/
        ├── oaspl_directivity.dat      # [theta_deg, OASPL_dB] per observer (if writeOASPL yes)
        ├── timeHistory_p.csv          # Sampled surface pressure cache (if writeTimeHistory yes)
        ├── timeHistory_Ux.csv         # Sampled surface velocity X cache
        ├── timeHistory_Uy.csv         # Sampled surface velocity Y cache
        ├── timeHistory_Uz.csv         # Sampled surface velocity Z cache
        ├── timeHistory_rho.csv        # Sampled surface density cache (fwhCompressibleFoam)
        │
        ├── observer0/                 # Directory per observer coordinate
        │   ├── pTime.dat              # Acoustic pressure signal: [time (s), p' (Pa)]
        │   ├── ft.dat                 # One-sided Fourier magnitude: [f (Hz), |P(f)| (Pa)]
        │   ├── spl.dat                # Sound Pressure Level: [f (Hz), SPL (dB re pRef)]
        │   ├── zoomFT.dat             # Zoom FFT magnitude: [f (Hz), |P(f)| (Pa)]
        │   └── zoomSPL.dat            # Zoom SPL spectrum: [f (Hz), SPL (dB)]
        │
        ├── observer1/
        │   └── ...
        │
        ├── pPrime/                     # 3D volumetric p' field (if writeVolumePPrime yes)
        │   ├── <time_0>/
        │   │   └── pPrime             # OpenFOAM volScalarField loadable in ParaView
        │   └── <time_1>/
        │       └── pPrime
        │
        └── OASPL/                     # 3D volumetric OASPL snapshot (if writeVolumeOASPL yes)
            └── <time>/
                └── OASPL              # OpenFOAM volScalarField loadable in ParaView
```

---

## Mathematical References

1. **Farassat, F.** (2007). *Derivation of Formulations 1 and 1A of Farassat*. NASA/TM-2007-214853.
2. **Ffowcs Williams, J. E., & Hawkings, D. L.** (1969). *Sound generation by turbulence and surfaces in arbitrary motion*. Philosophical Transactions of the Royal Society of London A, 264(1151), 321–342.
3. **Akima, H.** (1970). *A new method of interpolation and smooth curve fitting based on local procedures*. Journal of the ACM, 17(4), 589–602.

---

## Author & Contact

- **Author**: Chinmay Patil
- **GitHub**: [chinmay-s-patil](https://github.com/chinmay-s-patil)
- **E-Mail**: [patil.chinmay3031@gmail.com](mailto:patil.chinmay3031@gmail.com)
