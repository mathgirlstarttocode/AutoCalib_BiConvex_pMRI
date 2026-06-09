# AutoCalib_BiConvex_pMRI

Code accompanying:

Ni, Y., Strohmer, T. Auto-Calibration and Biconvex Compressive Sensing with Applications to Parallel MRI. J Fourier Anal Appl 32, 32 (2026). https://doi.org/10.1007/s00041-025-10223-1

## Repository structure

- `src/`: MATLAB helper functions and operators.
- `lsense.m`: main lifted l1,2 reconstruction routine.
- `exp/`: experiment scripts and notebooks.
- `data/`: small/preprocessed data files used by the experiments.
- `output/`: generated reconstruction outputs. Large outputs are not tracked.
- `external/`: third-party utilities included with attribution.

## Dependencies

MATLAB scripts require MATLAB with the Image Processing Toolbox and Signal Processing Toolbox. Some scripts use `parfor`, which requires the Parallel Computing Toolbox.

The synthetic phantom generator requires the external `MRIPhantomv0-8` package if regenerating phantom data from scratch.

The BART/NLINV baseline notebook requires the BART toolbox.

Python notebooks require NumPy, SciPy, Matplotlib, CVXPY, PyWavelets, and optionally MOSEK.

## Notes

The small synthetic transition-curve notebooks use direct CVXPY solves. The larger image reconstruction scripts use the custom lifted l1,2 FISTA-style solver in `lsense.m`.

