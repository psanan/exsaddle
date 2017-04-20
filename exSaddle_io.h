#ifndef EXSADDLE_IO_H_
#define EXSADDLE_IO_H_

#include "petscdm.h"
#include "petscksp.h"

#include "exSaddle.h"

PetscErrorCode SaddleReportSolutionDiagnostics(DM dm_saddle,Vec X);
PetscErrorCode DumpOperator(Mat A,const char* name);
PetscErrorCode DumpPreconditioner(KSP ksp,const char* name);
PetscErrorCode DumpSolution(Vec X, const char* name);
PetscErrorCode DumpPreconditionedOperator(KSP ksp,const char* name);
PetscErrorCode ViewFields(DM dm_saddle,Vec X,const char* tag);

#endif
