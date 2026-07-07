# Added Derived Quantities

This document describes the derived-quantity operators added to the ADIOS2 derived
variables toolkit (`source/adios2/toolkit/derived/`).

Derived variables are expressions evaluated over one or more existing variables at
write time. An expression such as `"GRADIENT(temperature)"` registers a new variable
whose data is computed from the operands via the operators below.

## Summary of additions

Four new operators were added, plus a correctness fix to the existing `POW` operator.

| Operator | Aliases | Operands | Output |
|----------|---------|----------|--------|
| `GRADIENT` | `GRAD` | 1 scalar 3D field (+ optional axis constant) | Vector field (`d0,d1,d2,3`), or scalar partial derivative |
| `MEAN` | — | 1 field (any dims) | Single scalar value (per writer block) |
| `VARIANCE` | `VAR` | 1 field (any dims) | Single scalar value (per writer block), floating point |
| `SPECTRUM` | `FFT` | 3 fields (`ux, uy, uz`), 3D, power-of-two dims | 1D radially-binned energy spectrum `E(k)` |

## Operators

### `GRADIENT(f)` / `GRAD(f)` / `GRADIENT(f, axis)`

Computes the spatial gradient of a 3D scalar field using **2nd-order central
differences in index space** (grid spacing = 1), with one-sided (clamped)
differences at the domain edges. This matches the discretisation convention already
used by the existing `CURL` operator.

- `GRADIENT(f)` — returns a vector field with a trailing component dimension of size
  3. Layout is `[i, j, k, component]` with component 0/1/2 = `d/d(dim0)`, `d/d(dim1)`,
  `d/d(dim2)` (component contiguous / last).
- `GRADIENT(f, axis)` — returns a scalar field containing the single partial
  derivative `d(f)/d(dim=axis)`, where `axis ∈ {0, 1, 2}`.

Only implemented for 3D arrays.

### `MEAN(f)`

Arithmetic mean (reduction) of a field to a single scalar value.

> **Note:** Under MPI decomposition this is a **per-writer-block** mean. A true
> global mean requires a size-weighted combine of the per-block values on the reader
> side. Accumulation is done in `double` precision regardless of input type.

### `VARIANCE(f)` / `VAR(f)`

Population variance (reduction) of a field to a single scalar value:

```
VAR(f) = (1/N) * sum_i (f_i - mean)^2
```

computed with a numerically stable **two-pass** algorithm (mean first, then sum of
squared deviations) accumulated in `double` precision. Divides by `N` (population
variance), matching the reduction convention of `MEAN`.

- Output is floating point — `double` for all inputs, `long double` for `long double`
  input — regardless of the input type, so integer fields yield a fractional variance.
- Like `MEAN`, this is a **per-writer-block** statistic. Under MPI decomposition a true
  global variance requires a size-weighted combine of the per-block `(mean, variance,
  count)` triples on the reader side (parallel / pooled-variance formula).

Intended use: a cheap statistical **trigger signal**. For example, monitoring
`VARIANCE(V)` of a Gray-Scott reaction-diffusion field detects the transition from the
perturbation-dominated phase to the reactive-pattern phase as a sharp rise in variance.

### `SPECTRUM(ux, uy, uz)` / `FFT(ux, uy, uz)`

Radially-binned kinetic-energy spectrum `E(k)` of a velocity field, computed via a
global 3D FFT:

```
E(k) = 0.5 * sum_{round|kvec|=k} (|ux^|^2 + |uy^|^2 + |uz^|^2) / N^2
```

normalised so that `sum_k E(k) = <0.5|u|^2>` (Parseval's theorem).

- Requires **power-of-two dimensions** on all three axes.
- Like `MEAN`, this is a **per-writer-block global operation**.
- Output is always `double` (energy), regardless of the input velocity type.
- Output is a 1D array of length `round(sqrt(3) * (nmax/2)) + 1` wavenumber shells,
  where `nmax = max(d0, d1, d2)`.

The FFT is an in-house iterative radix-2 Cooley–Tukey implementation (forward
transform, numpy sign convention), applied along each axis in turn.

## Bug fix: `POW` exponent

The `POW` operator previously truncated its exponent to an integer
(`static_cast<size_t>(std::stoull(...))`). It now parses the exponent as a
floating-point value (`std::stod`), so non-integer exponents such as `POW(x, 1.5)`
and `POW(x, 0.5)` are honoured instead of being truncated.

## Files changed

- `Expression.h` — new `ExpressionOperator` enum values: `OP_GRAD`, `OP_MEAN`, `OP_SPECTRUM`.
- `Expression.cpp` — operator properties, string→operator mappings (incl. aliases
  `GRAD`, `FFT`), and operator→function bindings.
- `Function.h` — declarations for the new `*Func` and `*DimsFunc` functions.
- `Function.cpp` — implementations:
  - `ApplyGradient` / `GradientFunc` / `GradDimsFunc`
  - `MeanFunc` / `MeanDimsFunc`
  - `fft1d` / `fft3d` / `spectrumBins` / `ApplySpectrum` / `SpectrumFunc` / `SpectrumDimsFunc`
  - `PowFunc` exponent fix

## Type support

All operators dispatch over the standard ADIOS2 primitive types via
`ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG`. `SPECTRUM` requires floating-point
input (`FloatTypeFunc`); `GRADIENT` and `MEAN` preserve the input type (`SameTypeFunc`).
