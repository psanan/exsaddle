#ifndef FEMIXED_SPACE_H_
#define FEMIXED_SPACE_H_

#include "petscdm.h"

#include "exSaddle.h"

#if defined(LAME)
typedef struct {
  PetscScalar mu,lambda;
  PetscScalar Fu[NSD];
  PetscScalar Fp;
} LameCoefficient;
#define EXSADDLE_NCOEFF NSD+3
#else
typedef struct {
  PetscScalar eta;
  PetscScalar Fu[NSD];
  PetscScalar Fp;
} StokesCoefficient;
#define EXSADDLE_NCOEFF NSD+2
#endif

typedef struct {
  PetscInt n_qpoints;
  PetscReal *qp_coor;
  PetscReal *qp_weight;
} VQuadraturePoint;

struct _p_FEMixedSpace {
  PetscInt          n_u_elements_domain; /* Global velocity elements */
  PetscInt          n_p_elements_domain; /* Global pressure elements */
  PetscInt          n_u_elements;        /* Number of local velocity elements */
  PetscInt          n_p_elements;        /* Number of local pressure elements (warning: assumed to be equal to n_u_elements in some places)*/
  PetscInt          *u_el_nd_map;        /* Q2 element -> node map (local indexing) */
  PetscInt          *p_el_nd_map;        /* Q1 element -> node map (local indexing) */
#if NSD ==2 
  PetscInt          lmx,lmy;
#elif NSD == 3
  PetscInt          lmx,lmy,lmz;
#endif
  VQuadraturePoint  *vol_quadrature;
#if defined(LAME)
  LameCoefficient   *coeff_cell; /* indexed as qp_wise_coeff[eidx*n_qpoints + qidx] */
  LameCoefficient   *coeff_qp;   /* indexed as cell_wise_coeff[eidx] */
#else
  StokesCoefficient *coeff_cell; /* indexed as qp_wise_coeff[eidx*n_qpoints + qidx] */
  StokesCoefficient *coeff_qp;   /* indexed as cell_wise_coeff[eidx] */
#endif
  IS                u_is_global,u_is_local;
  PetscReal         *u_bc_global,*u_bc_local;
  PetscInt          l_nbc;
  PetscInt          *bc_local_idx_g; 
  DM                dm_saddle;
  PetscReal         *e_centroid,*e_dx;
};

typedef struct _p_FEMixedSpace *FEMixedSpace;

PetscErrorCode FEMixedSpaceCreate(FEMixedSpace *space);
PetscErrorCode DMCreate_SaddleQ2Q1(MPI_Comm comm,PetscInt mx,PetscInt my,PetscInt mz,DM *dm_stokes,FEMixedSpace fespace);
PetscErrorCode DMDASetUniformCoordinates_Saddle(DM dm_saddle,PetscReal x0,PetscReal x1,PetscReal y0,PetscReal y1,PetscReal z0,PetscReal z1);
PetscErrorCode FEMixedSpaceQuadratureCreate(FEMixedSpace space,PetscBool alloc_cell_coeff,PetscBool alloc_qp_coeff);
PetscErrorCode FEMixedSpaceDefineQPwiseProperties(FEMixedSpace space,DM dmv);
PetscErrorCode FEMixedSpaceDefineQPwiseProperties_Q1Projection(PetscInt nlevels,FEMixedSpace _space[],DM _dmscalar[]);
PetscErrorCode FEMixedSpaceDestroy(FEMixedSpace *space);
PetscErrorCode DMCreateDomainDecomposition_DMDAFEQ2Q1(DM dm,PetscInt *len,char ***namelist,IS **innerislist,IS **outerislist,DM **dmlist);
PetscErrorCode MatAssemble_Saddle_NULL(FEMixedSpace space,DM dm_saddle,Mat A);
PetscErrorCode MatAssemble_Saddle(FEMixedSpace space,DM dm_saddle,Mat A,Vec rhs_diri);
PetscErrorCode VecAssemble_F1_qp(FEMixedSpace space,Vec F);
PetscErrorCode VecAssemble_F2_qp(FEMixedSpace space,Vec F);
PetscErrorCode ImposeDirichletValuesIS(Vec x,IS is,PetscScalar vals[]);
PetscErrorCode ZeroDirichletValuesIS(Vec x,IS is);
PetscErrorCode MatAssemble_Schur(FEMixedSpace space,DM dmp,Mat S);
PetscErrorCode FEMixedSpaceBCISCreate(FEMixedSpace space,DM dmsaddle);

#endif
