# Building exstruct

`exstruct` extracts high-dissipation (or any thresholded) structures from a 3D
scalar field stored in an HDF5 file. It is a single Fortran source file,
[exstruct.f90](exstruct.f90), that links against the **HDF5 Fortran** library.

There are two ways to build it: the **CMake** build (recommended, portable) and
the legacy **Makefile** (hard-coded for `nvfortran`).

## Prerequisites

- A Fortran compiler — `gfortran` or `nvfortran` are both supported.
- **HDF5 with Fortran bindings**, built with the *same* compiler you use here.
  The HDF5 Fortran modules (`hdf5.mod`) are compiler-specific and cannot be read
  by a different compiler, so a mismatch is the most common build failure.
- CMake ≥ 3.18 (for the CMake build; required for the `Fortran_PREPROCESS`
  property, since the source uses `#ifdef`).
- Optionally, an OpenMP-capable compiler for threaded builds.

## CMake build (recommended)

```sh
cmake -B build -DHDF5_ROOT=/path/to/hdf5
cmake --build build -j
```

The executable is produced as `build/exstruct.exe`.

### Common options

| Option | Effect | Default |
| --- | --- | --- |
| `-DHDF5_ROOT=<prefix>` | Point CMake at your HDF5 install | (autodetect) |
| `-DCMAKE_Fortran_COMPILER=nvfortran` (or `gfortran`) | Choose the compiler | (autodetect) |
| `-DEXSTRUCT_OPENMP=ON` | Enable OpenMP threading | `OFF` |
| `-DCMAKE_BUILD_TYPE=Release\|Debug\|RelWithDebInfo\|MinSizeRel` | Build type | `Release` |

> **Important:** the HDF5 Fortran modules must have been built with the same
> Fortran compiler used here. If CMake picks up an HDF5 built with a different
> compiler, set both `-DCMAKE_Fortran_COMPILER=...` and `-DHDF5_ROOT=...` so they
> agree.

Example, threaded build with `nvfortran`:

```sh
cmake -B build \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DHDF5_ROOT=/home/giorgos/libs/hdf5 \
  -DEXSTRUCT_OPENMP=ON
cmake --build build -j
```

## Legacy Makefile build

The [Makefile](Makefile) is hard-coded for `nvfortran` and a fixed HDF5 prefix.
Edit `HDF5ROOT` (and `ARCH`/`FC` if needed) to match your system, then:

```sh
make            # produces ./exstruct.exe
make clean      # removes *.o and *.mod
```

This path is kept for convenience on the original development machine; prefer
the CMake build elsewhere.

## Running

```sh
exstruct.exe -f FILE -n FIELD [-s SDEVS] [-v VOLUME]
```

| Flag | Meaning | Default |
| --- | --- | --- |
| `-f`, `--file FILE` | Path to the input HDF5 file (required) | — |
| `-n`, `--field FIELD` | Name of the 3D dataset to read (required) | — |
| `-s`, `--sdevs N` | Threshold at `mean + N*stddev` (positive integer) | `3` |
| `-v`, `--volume N` | Output only structures with more than N points | `0` |
| `-h`, `--help` | Show help and exit | — |

Options also accept the `--key=value` form.

Example, using the bundled sample [in.h5](in.h5) (which contains a `/diss`
dataset):

```sh
./build/exstruct.exe -f in.h5 -n diss -s 3 -v 0
```

### Outputs

Files are written to the current working directory:

- `struct.bin` — unformatted binary: the structure count, then per structure the
  point count and its `(i,j,k)` integer coordinates.
- `gs.dat` — text file, one line per structure giving its point count.
- `out.<n>.vtk` — one ASCII VTK polydata file per structure whose size exceeds
  `--volume`.

All of these (plus `*.mod`, `*.exe`, the `build/` directory, and the sample HDF5
files) are listed in [.gitignore](.gitignore) and are not tracked.

## Troubleshooting

- **`Cannot open module file hdf5.mod` / unreadable module:** the HDF5 Fortran
  modules were built with a different compiler than the one in use. Rebuild HDF5
  with the same compiler, or point `-DHDF5_ROOT` / `-DCMAKE_Fortran_COMPILER` at
  a matching install.
- **HDF5 not found by CMake:** pass `-DHDF5_ROOT=<prefix>` explicitly.
- **`error: the field dataset is not 3-dimensional`** (at runtime): the dataset
  named by `-n` exists but is not a rank-3 array — check the field name with
  `h5dump -n FILE`.
