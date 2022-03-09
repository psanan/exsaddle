# exSaddle

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4040667.svg)](https://doi.org/10.5281/zenodo.4040667)

By Dave A. May and Patrick Sanan.

PETSc codes to test preconditioners for Stokes and mixed-form linear elasticity saddle-point systems.

Requires [PETSc](https://mcs.anl.gov/petsc) 3.12 or later.

This includes a family of Q2-Q1 (Taylor-Hood) codes for Stokes and mixed-form linear elasticity systems, suitable for tests with both monolithic and segregated preconditioners constructed by composing Approximate Block Factorization (ABF) and/or multigrid methods.

You will also find modifications of PETSc KSP tutorial examples ex23 and ex42.

Also see [PCILUPACK](github.com/psanan/pcilupack) which is a plug-in
version of the ILUPACK-based preconditioners here, useful for trying
out the preconditioner with an existing PETSc-based code.

A [paper](https://se.copernicus.org/articles/11/2031/2020/) and its [included supplement](https://se.copernicus.org/articles/11/2031/2020/se-11-2031-2020-supplement.pdf) include further information, experiments, and instructions regarding this code.

If appropriate, please cite the paper and/or this repository.

```
@article{SananMayBollhoeferSchenk2020,
  title        = {Pragmatic solvers for {3D} {S}tokes and elasticity problems with heterogeneous coefficients: {E}valuating modern incomplete {LDL$^T$} preconditioners},
  author       = {Sanan, Patrick and May, Dave A. and Bollh\"{o}fer, Matthias and Schenk, Olaf},
  year         = 2020,
  journal      = {Solid Earth},
  volume       = 11,
  number       = 6,
  pages        = {2031--2045},
  doi          = {10.5194/se-11-2031-2020},
  url          = {https://se.copernicus.org/articles/11/2031/2020/}
}

@software{SananMay2020,
  author       = {Sanan, Patrick and May, Dave A.},
  title        = {exSaddle},
  month        = sep,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {0.1.1},
  doi          = {10.5281/zenodo.4040667},
  url          = {https://doi.org/10.5281/zenodo.4040667}
}
```
