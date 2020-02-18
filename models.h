#ifndef MODELS_H_
#define MODELS_H_

#include "exSaddle.h"
#include "femixedspace.h"

#include "petscdm.h"

#if defined(LAME)
#define DEFAULT_MODEL 6
#else
#define DEFAULT_MODEL 2
#endif

PetscErrorCode ISCreate_BCList(DM dm,PetscBool global,IS *_is,PetscReal **_vals);
#if defined(LAME)
PetscErrorCode Lame_EvaluateCoefficients(PetscReal coor[],PetscReal *mu,PetscReal *lambda,PetscReal Fu[],PetscReal Fp[]);
#else
PetscErrorCode Stokes_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[]);
#endif
PetscErrorCode ComputeReferenceSolution(DM dm_saddle,Vec *Xref);

#endif
