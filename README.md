# exSaddle
By Dave A. May and Patrick Sanan.

PETSc codes to test preconditioners for Stokes and Lame' saddle-point systems.

Requires [PETSc](https://mcs.anl.gov/petsc) 3.12 or later.

This includes a family of Q2-Q1 (Taylor-Hood) codes for Stokes and mixed-form linear elasticity systems, suitable for tests with both monolithic and segregated preconditioners constructed by composing Approximate Block Factorization (ABF) and/or multigrid methods.

You will also find modifications of PETSc KSP tutorial examples ex23 and ex42.

Also see [PCILUPACK](github.com/psanan/pcilupack) which is a plug-in
version of the ILUPACK-based preconditioners here, useful for trying
out the preconditioner with an existing PETSc-based code.

A [paper](https://se.copernicus.org/preprints/se-2020-79/>) and its [included supplement](https://se.copernicus.org/preprints/se-2020-79/se-2020-79-supplement.pdf) include further information, experiments, and instructions regarding this code.

If appropriate, please cite this repository (using the DOI which will be added above) and/or the paper.

```
@article{SananMayBollhoeferSchenk2020,
  title   = {{P}ragmatic Solvers for {3D} {S}tokes and Elasticity Problems with Heterogeneous Coefficients: {E}valuating Modern Incomplete {$LDL^T$} Preconditioners},
  author  = {Patrick Sanan and Dave A. May and Matthias Bollh\"{o}fer and Olaf Schenk},
  journal = {Solid Earth preprint},
  url     = {https://se.copernicus.org/preprints/se-2020-79/},
  year    = {2020},
}
```
(Note that this information will be updated [here](github.com/psanan/exsaddle) once a final version of the paper is available online.)
