#include "models.h"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

/*****************************************************************************/
static PetscErrorCode SolCx_ISCreate_BCList(DM dmv,PetscBool global,IS *_is,PetscReal **_vals)
{
  PetscErrorCode   ierr;
  PetscInt         i,j,si,sj,sk,ni,nj,nk,M,N,P;
#if NSD == 3
  PetscInt         k,length;
#endif
  PetscInt         *idx;
  PetscInt         L,cnt,dof_to_constrain;
  IS               is;
  PetscBool        is_free_slip = PETSC_FALSE;
  PetscReal        *vals;
  static PetscBool been_here = PETSC_FALSE;

  PetscFunctionBeginUser;

  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"Boundary Conditions: SolCx\n");
    been_here = PETSC_TRUE;
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-freesliphack",&is_free_slip,NULL);CHKERRQ(ierr);

  ierr = DMDAGetInfo(dmv,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  if (global) {
    ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  } else {
    ierr = DMDAGetGhostCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  }

#if NSD == 2
  PetscMalloc(sizeof(PetscReal)*(ni+nj+ni+nj),&vals);
  PetscMemzero(vals,sizeof(PetscReal)*(ni+nj+ni+nj));
  
  L = 0;
  if (si == 0) { L += nj * 1; }
  if (sj == 0) { L += ni * 1; }
  
  if (si+ni == M) { L += nj * 1; }
  if (is_free_slip) {
    if (sj+nj == N) { L += ni * 1; }
  }
  PetscMalloc(sizeof(PetscInt)*L,&idx);
  
  cnt = 0;
  if (si == 0) {
    i = 0; dof_to_constrain = 0;
    for (j=0; j<nj; ++j) {
      idx[cnt] = NSD*(i + j*ni )+dof_to_constrain; ++cnt;
    }
  }
  if (sj == 0) {
    j = 0; dof_to_constrain = 1;
    for (i=0; i<ni; ++i) {
      idx[cnt] = NSD*(i + j*ni )+dof_to_constrain; ++cnt;
    }
  }
  
  if (si+ni == M) {
    i = ni-1; dof_to_constrain = 0;
    for (j=0; j<nj; ++j) {
      idx[cnt] = NSD*(i + j*ni )+dof_to_constrain; ++cnt;
    }
  }
  if (is_free_slip) {
    if (sj+nj == N) {
      j = nj-1; dof_to_constrain = 1;
      for (i=0; i<ni; ++i) {
        idx[cnt] = NSD*(i + j*ni )+dof_to_constrain; ++cnt;
      }
    }
  }
#elif NSD == 3 
  length = 2 * (ni*nj + ni*nk + nj*nk); /* 6 faces */
  PetscMalloc(sizeof(PetscReal)*length,&vals);
  PetscMemzero(vals,sizeof(PetscReal)*length);
  
  L = 0;
  if (si == 0) { L += nj * nk; } 
  if (sj == 0) { L += ni * nk; }
  if (sk == 0) { L += ni * nj; }
  
  if (si+ni == M) { L += nj * nk; }
  if (is_free_slip) {
    if (sj+nj == N) { L += ni * nk; }
  }
  if (sk+nk == P) { L += ni * nj; }
  PetscMalloc(sizeof(PetscInt)*L,&idx);
  

  cnt = 0;
  if (si == 0) {
    i = 0; dof_to_constrain = 0;
    for (j=0; j<nj; ++j) {
      for (k=0; k<nk; ++k) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj)+dof_to_constrain; ++cnt; 
      }
    }
  }
  if (sj == 0) {
    j = 0; dof_to_constrain = 1;
    for (i=0; i<ni; ++i) {
      for (k=0; k<nk; ++k) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj)+dof_to_constrain; ++cnt;
      }
    }
  }
  if (sk == 0) {
    k = 0; dof_to_constrain = 2;
    for (i=0; i<ni; ++i) {
      for (j=0; j<nj; ++j) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj)+dof_to_constrain; ++cnt;
      }
    }
  }

  if (si+ni == M) {
    i = ni-1; dof_to_constrain = 0;
    for (j=0; j<nj; ++j) {
      for (k=0; k<nk; ++k) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj)+dof_to_constrain; ++cnt;
      }
    }
  }
  if (is_free_slip) {
    if (sj+nj == N) {
      j = nj-1; dof_to_constrain = 1;
      for (i=0; i<ni; ++i) {
        for (k=0; k<nk; ++k) {
          idx[cnt] = NSD*(i + j*ni + k*ni*nj)+dof_to_constrain; ++cnt;
        }
      }
    }
  }
  if (sk+nk == P) {
    k = nk-1; dof_to_constrain = 2;
    for (i=0; i<ni; ++i) {
      for (j=0; j<nj; ++j) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj)+dof_to_constrain; ++cnt;
      }
    }
  }
#endif

  ierr = ISCreateGeneral(PETSC_COMM_SELF,L,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  PetscFree(idx);
  *_is   = is;
  *_vals = vals;
  
  PetscFunctionReturn(0);
}

#if defined(LAME) || NSD==3
/*****************************************************************************/
static PetscErrorCode FixedBase_ISCreate_BCList(DM dmv,PetscBool global,IS *_is,PetscReal **_vals)
{
  PetscErrorCode   ierr;
  PetscInt         i,j,si,sj,sk,ni,nj,nk,M,N,P;
#if NSD == 3
  PetscInt         k;
#endif
  PetscInt         *idx;
  PetscInt         L,cnt,dof_to_constrain;
  IS               is;
  PetscReal        *vals;
  static PetscBool been_here = PETSC_FALSE;
 
  PetscFunctionBeginUser;

  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"Boundary Conditions: FixedBase\n");
    been_here = PETSC_TRUE;
  }

  ierr = DMDAGetInfo(dmv,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  if (global) {
    ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  } else {
    ierr = DMDAGetGhostCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  }

#if NSD == 2
  PetscMalloc(sizeof(PetscReal)*(NSD*ni),&vals);
  PetscMemzero(vals,sizeof(PetscReal)*(NSD*ni));
 
  /* Constrain only the bottom face */
  L = (sj == 0) ? NSD * ni : 0;
  PetscMalloc(sizeof(PetscInt)*L,&idx);
  
  cnt = 0;
  if (sj == 0) {
    for(dof_to_constrain = 0; dof_to_constrain<NSD; ++dof_to_constrain){
      for (i=0; i<ni; ++i) {
        j = 0;
        idx[cnt] = NSD*(i + j*ni )+dof_to_constrain; ++cnt;
      }
    }
  }

#elif NSD == 3 
  PetscMalloc(sizeof(PetscReal)*(NSD*ni*nk),&vals);
  PetscMemzero(vals,sizeof(PetscReal)*(NSD*ni*nk));
 
  /* Constrain only the bottom face */
  L = (sj == 0) ? NSD * ni * nk : 0;
  PetscMalloc(sizeof(PetscInt)*L,&idx);
  
  cnt = 0;
  if (sj == 0) {
    j = 0;
    for(dof_to_constrain = 0; dof_to_constrain<NSD; ++dof_to_constrain){
      for (i=0; i<ni; ++i) {
        for (k=0; k<nk; ++k) {
          idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; ++cnt;
        }
      }
    }
  }
#endif

  ierr = ISCreateGeneral(PETSC_COMM_SELF,L,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  PetscFree(idx);
  *_is   = is;
  *_vals = vals;
  
  PetscFunctionReturn(0);
}
#endif

#if defined(LAME)
/*****************************************************************************/
static PetscErrorCode Compression_ISCreate_BCList(DM dmv,PetscBool global,IS *_is,PetscReal **_vals)
{
  PetscErrorCode   ierr;
  PetscInt         i,j,si,sj,sk,ni,nj,nk,M,N,P;
#if NSD == 3
  PetscInt         k;
#endif
  PetscInt         *idx;
  PetscInt         L,cnt,dof_to_constrain;
  IS               is;
  PetscReal        *vals;
  static PetscBool been_here = PETSC_FALSE;
  const PetscReal  displacement = 0.1;
 
  PetscFunctionBeginUser;

  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"Boundary Conditions: Compression\n");
    been_here = PETSC_TRUE;
  }

  ierr = DMDAGetInfo(dmv,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  if (global) {
    ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  } else {
    ierr = DMDAGetGhostCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  }

  /* Constrain the left and right faces */
  L = 0;
  if (si == 0)    L += NSD * nj * nk;
  if (si+ni == N) L += NSD * nj * nk;
  PetscMalloc(sizeof(PetscInt)*L,&idx);
  PetscMalloc(      sizeof(PetscScalar)*L,&vals);
  PetscMemzero(vals,sizeof(PetscScalar)*L);
  cnt = 0;
#if NSD == 2
  if (si == 0) {
    i = 0;
    for(dof_to_constrain = 0; dof_to_constrain<NSD; ++dof_to_constrain){
      for (j=0; j<nj; ++j) {
        idx[cnt] = NSD*(i + j*ni)+dof_to_constrain; 
        if (dof_to_constrain == 0 ){
          vals[cnt] = displacement;
        }
        ++cnt;
      }
    }
  }
  if (si+ni == N) {
    i = ni-1;
    for(dof_to_constrain = 0; dof_to_constrain<NSD; ++dof_to_constrain){
      for (j=0; j<nj; ++j) {
        idx[cnt] = NSD*(i + j*ni)+dof_to_constrain; 
        if (dof_to_constrain == 0){
          vals[cnt] = - displacement;
        }
        ++cnt;
      }
    }
  }
#elif NSD == 3 
  if (si == 0) {
    i = 0;
    for(dof_to_constrain = 0; dof_to_constrain<NSD; ++dof_to_constrain){
      for (j=0; j<nj; ++j) {
        for (k=0; k<nk; ++k) {
          idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
          if (dof_to_constrain == 0 ){
            vals[cnt] = displacement;
          }
          ++cnt;
        }
      }
    }
  }
  if (si+ni == N) {
    i = ni-1;
    for(dof_to_constrain = 0; dof_to_constrain<NSD; ++dof_to_constrain){
      for (j=0; j<nj; ++j) {
        for (k=0; k<nk; ++k) {
          idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
          if (dof_to_constrain == 0){
            vals[cnt] = - displacement;
          }
          ++cnt;
        }
      }
    }
  }
#endif

  ierr = ISCreateGeneral(PETSC_COMM_SELF,L,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  PetscFree(idx);
  *_is   = is;
  *_vals = vals;
  
  PetscFunctionReturn(0);
}
#endif

#if defined(LAME) && NSD == 3
/*****************************************************************************/
static PetscErrorCode Compression2_ISCreate_BCList(DM dmv,PetscBool global,IS *_is,PetscReal **_vals)
{
  PetscErrorCode   ierr;
  PetscInt         i,j,k,si,sj,sk,ni,nj,nk,M,N,P;
  PetscInt         *idx;
  PetscInt         L,cnt;
  IS               is;
  PetscReal        *vals;
  static PetscBool been_here = PETSC_FALSE;
  const PetscReal  displacement = 0.1;
 
  PetscFunctionBeginUser;

  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"Boundary Conditions: Compression2\n");
    been_here = PETSC_TRUE;
  }

  ierr = DMDAGetInfo(dmv,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  if (global) {
    ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  } else {
    ierr = DMDAGetGhostCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  }

  /* Free slip on all faces except the top, which is unconstrained. Apply
     opposing displacements on the left and right faces */
  L = 0;
  if (si == 0)    L += nj * nk; /* left face */
  if (si+ni == N) L += nj * nk; /* right face */
  if (sj == 0)    L += ni * nk; /* bottom face */
  if (sk == 0)    L += nj * nk; /* back face */
  if (sk+nk == P) L += nj * nk; /* front face */
  PetscMalloc(sizeof(PetscInt)*L,&idx);
  PetscMalloc(      sizeof(PetscScalar)*L,&vals);
  PetscMemzero(vals,sizeof(PetscScalar)*L);
  cnt = 0;
  /* Displace on left face */
  if (si == 0) {
    const PetscInt dof_to_constrain = 0;

    i = 0;
    for (j=0; j<nj; ++j) {
      for (k=0; k<nk; ++k) {
        idx[cnt]  = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
        vals[cnt] = displacement;
        ++cnt;
      }
    }
  }
  /* Displace (in the opposite direction) on the right face */
  if (si+ni == N) {
    const PetscInt dof_to_constrain = 0;

    i = ni-1;
    for (j=0; j<nj; ++j) {
      for (k=0; k<nk; ++k) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
        vals[cnt] = - displacement;
        ++cnt;
      }
    }
  }

  /* Free slip on the bottom face */
  if (sj == 0) {
    const PetscInt dof_to_constrain = 1;

    j = 0;
    for (i=0; i<ni; ++i) {
      for (k=0; k<nk; ++k) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
        vals[cnt] = 0.0;
        ++cnt;
      }
    }
  }

  /* Free slip on back face */
  if (sk == 0) {
    const PetscInt dof_to_constrain = 2;

    k = 0;
    for (i=0; i<ni; ++i) {
      for (j=0; j<nj; ++j) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
        vals[cnt] = 0.0;
        ++cnt;
      }
    }
  }

  /* Free slip on front face */
  if (sk+nk == P) {
    const PetscInt dof_to_constrain = 2;

    k = nk-1;
    for (i=0; i<ni; ++i) {
      for (j=0; j<nj; ++j) {
        idx[cnt] = NSD*(i + j*ni + k*ni*nj )+dof_to_constrain; 
        vals[cnt] = 0.0;
        ++cnt;
      }
    }
  }


  ierr = ISCreateGeneral(PETSC_COMM_SELF,L,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  PetscFree(idx);
  *_is   = is;
  *_vals = vals;
  
  PetscFunctionReturn(0);
}
#endif

#ifndef LAME
#if NSD==2
/**********************************************************/
/* Helper functions for MMS1 */
static PetscScalar MMS1_SOLX(PetscReal x, PetscReal y){ return 20*x*y*y*y; }
static PetscScalar MMS1_SOLY(PetscReal x, PetscReal y){ return 5*(x*x*x*x-y*y*y*y); }
static PetscScalar MMS1_SOLP(PetscReal x, PetscReal y){ return 60*x*x*y-20*y*y*y; }

static PetscErrorCode StokesMMS1_ISCreate_BCList(DM dmv,PetscBool global,IS *_is,PetscReal **_vals)
{
  PetscErrorCode   ierr;
  PetscInt         i,j,d,si,sj,sk,ni,nj,nk,M,N,P,*idx,L,cnt;
  IS               is;
  PetscReal        *vals;
  static PetscBool been_here = PETSC_FALSE;
  Vec              coord_l;
  PetscScalar      *LA_coord_l;
  PetscReal        c[NSD];
 
  PetscFunctionBeginUser;

  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"Boundary Conditions: StokesMMS1\n");
    been_here = PETSC_TRUE;
  }

  ierr = DMDAGetInfo(dmv,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  if (global) {
    ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  } else {
    ierr = DMDAGetGhostCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  }
  ierr = DMGetCoordinatesLocal(dmv,&coord_l);CHKERRQ(ierr);
  ierr = VecGetArray(coord_l,&LA_coord_l);CHKERRQ(ierr);
 
  /* Constrain all 4 faces in both velocity components*/
  L = 0;
  if (si    == 0) L += NSD * nj;
  if (si+ni == N) L += NSD * nj;
  if (sj    == 0) L += NSD * ni;
  if (sj+nj == M) L += NSD * ni;

  PetscMalloc(sizeof(PetscInt)*L,&idx);
  PetscMalloc(sizeof(PetscReal)*L,&vals);
  PetscMemzero(vals,sizeof(PetscReal)*L);

  cnt = 0;
  if (si == 0) {
    i = 0;
    for (j=0; j<nj; ++j) {
      for (d=0; d<NSD; ++d){
        c[d] = LA_coord_l[NSD*(i+j*ni)+d];
      }
      for(d = 0; d<NSD; ++d){
        idx[cnt] = NSD*(i + j*ni)+d;
        vals[cnt] = d == 0 ? MMS1_SOLX(c[0],c[1]) : MMS1_SOLY(c[0],c[1]);
        // DEBUG
#if 0
        {
          PetscMPIInt rank;
          ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] debug. i=%d,j=%d,dof=%d,c=[%g,%g],val=%g \n",rank,i,j,d,c[0],c[1],vals[cnt]);CHKERRQ(ierr);
        }
#endif
        ++cnt;
      }
    }
  }

  if (si+ni == N) {
    i = ni-1;
    for (j=0; j<nj; ++j) {
      for (d=0; d<NSD; ++d){
        c[d] = LA_coord_l[NSD*(i+j*ni)+d];
      }
      for(d = 0; d<NSD; ++d){
        idx[cnt] = NSD*(i + j*ni)+d; 
        vals[cnt] = d == 0 ? MMS1_SOLX(c[0],c[1]) : MMS1_SOLY(c[0],c[1]);
        // DEBUG
#if 0 
        {
          PetscMPIInt rank;
          ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] debug. i=%d,j=%d,dof=%d,c=[%g,%g],val=%g \n",rank,i,j,d,c[0],c[1],vals[cnt]);CHKERRQ(ierr);
        }
#endif
        ++cnt;
      }
    }
  }

  if (sj == 0) {
    j = 0;
    for (i=0; i<ni; ++i) {
      for (d=0; d<NSD; ++d){
        c[d] = LA_coord_l[NSD*(i+j*ni)+d];
      }
      for(d = 0; d<NSD; ++d){
        idx[cnt] = NSD*(i + j*ni)+d; 
        vals[cnt] = d == 0 ? MMS1_SOLX(c[0],c[1]) : MMS1_SOLY(c[0],c[1]);

        // DEBUG
#if 0
        {
          PetscMPIInt rank;
          ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] debug. i=%d,j=%d,dof=%d,c=[%g,%g],val=%g \n",rank,i,j,d,c[0],c[1],vals[cnt]);CHKERRQ(ierr);
        }
#endif
        ++cnt;
      }
    }
  }

  if (sj+nj == M) {
    j = nj-1;
    for (i=0; i<ni; ++i) {
      for (d=0; d<NSD; ++d){
        c[d] = LA_coord_l[NSD*(i+j*ni)+d];
      }
      for(d = 0; d<NSD; ++d){
        idx[cnt] = NSD*(i + j*ni)+d;
        vals[cnt] = d == 0 ? MMS1_SOLX(c[0],c[1]) : MMS1_SOLY(c[0],c[1]);

        // DEBUG
#if 0
        {
          PetscMPIInt rank;
          ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] debug. i=%d,j=%d,dof=%d,c=[%g,%g],val=%g \n",rank,i,j,d,c[0],c[1],vals[cnt]);CHKERRQ(ierr);
        }
#endif
        ++cnt;
      }
    }
  }

  ierr = VecRestoreArray(coord_l,&LA_coord_l);CHKERRQ(ierr);

  ierr = ISCreateGeneral(PETSC_COMM_SELF,L,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  PetscFree(idx);
  *_is   = is;
  *_vals = vals;
  
  PetscFunctionReturn(0);
}
#endif
#endif

/*****************************************************************************/
/* Note: using integer values to define models is bad design. It should be done
   with an enumerated type */
PetscErrorCode ISCreate_BCList(DM dmv,PetscBool global,IS *_is,PetscReal **_vals)
{
  PetscErrorCode ierr;
  PetscInt       model = DEFAULT_MODEL;
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-model",&model,NULL);CHKERRQ(ierr);

  switch (model) {
#if defined(LAME)
    case 8:
      ierr = FixedBase_ISCreate_BCList(dmv,global,_is,_vals);CHKERRQ(ierr);
      break;
    case 9:
    case 10:
      ierr = Compression_ISCreate_BCList(dmv,global,_is,_vals);CHKERRQ(ierr);
      break;
#endif
#if NSD==3
    case 11:
      ierr = FixedBase_ISCreate_BCList(dmv,global,_is,_vals);CHKERRQ(ierr);
      break;
#endif
#if defined(LAME) && NSD==3
    case 12:
      ierr = Compression2_ISCreate_BCList(dmv,global,_is,_vals);CHKERRQ(ierr);
      break;
#endif
#if !defined(LAME) && NSD==2
    case 101:
      ierr = StokesMMS1_ISCreate_BCList(dmv,global,_is,_vals);CHKERRQ(ierr);
      break;
#endif
    default:
      ierr = SolCx_ISCreate_BCList(dmv,global,_is,_vals);CHKERRQ(ierr);
      break;
  }

  PetscFunctionReturn(0);
}

#if defined(LAME)
/*
Lame' (Elasticity)
*/
/*****************************************************************************/
PetscErrorCode LameOneSinker_EvaluateCoefficients(PetscReal coor[],PetscReal *mu,PetscReal *lambda, PetscReal Fu[],PetscReal Fp[])
{
  PetscErrorCode   ierr;
  static PetscReal opts_mu0,opts_mu1,opts_lambda0,opts_lambda1,opts_rad;
  static PetscBool been_here = PETSC_FALSE;
  PetscReal        mu_qp,lambda_qp,rho_qp;
  PetscBool        inside = PETSC_FALSE;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: LameOneSinker\n");
    opts_mu0     = 1.0;
    opts_mu1     = 1.0;
    opts_lambda0 = 1.0;
    opts_lambda1 = 2.0;
    opts_rad     = 0.25;
    ierr = PetscOptionsGetReal(NULL,NULL,"-mu0",&opts_mu0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-mu1",&opts_mu1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-lambda0",&opts_lambda0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-lambda1",&opts_lambda1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_r",&opts_rad,0);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: mu0 %1.4e\n",opts_mu0);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: mu1 %1.4e\n",opts_mu1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: lambda0 %1.4e\n",opts_lambda0);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: lambda1 %1.4e\n",opts_lambda1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: rad %1.4e\n",opts_rad);CHKERRQ(ierr);
    been_here = PETSC_TRUE;
  }
  
  mu_qp = opts_mu0;
  lambda_qp = opts_lambda0;
  rho_qp = 1.0;
  
  {
#if NSD == 2
    const PetscReal sep2 = (coor[0] - 0.5)*(coor[0] - 0.5) + (coor[1] - 0.5)*(coor[1] - 0.5);
#elif NSD == 3 
    const PetscReal sep2 = (coor[0] - 0.5)*(coor[0] - 0.5) + (coor[1] - 0.5)*(coor[1] - 0.5) + (coor[2] - 0.5)*(coor[2] - 0.5);
#endif
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;
  }
  
  if (inside) {
    rho_qp = 2.0; /* 2x density inside the inclusion */
    mu_qp = opts_mu1;
    lambda_qp = opts_lambda1;
  }
  
  if (lambda) {
    *lambda = lambda_qp;
  }
  if (mu) {
    *mu = mu_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
  PetscFunctionReturn(0);
}

/******************************************************************************/
PetscErrorCode LameXSinker_EvaluateCoefficients(PetscReal coor[],PetscReal *mu,PetscReal *lambda,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_mu0,opts_mu1,opts_lambda0,opts_lambda1,opts_rad;
  static PetscBool been_here = PETSC_FALSE;
  static PetscInt  opts_numinc;
  const PetscReal posx[8] = { 0.27 , 0.6  , 0.7  , 0.2 , 0.85  , 0.4 , 0.16 , 0.55 };
  const PetscReal posy[8] = { 0.63 , 0.83 , 0.33 , 0.2 , 0.65  , 0.3 , 0.84 , 0.54 };
#if NSD == 3
  const PetscReal posz[8] = { 0.50 , 0.40 , 0.30 , 0.70 , 0.65  , 0.4 , 0.8  , 0.50 };
#endif
  PetscReal      mu_qp,lambda_qp,rho_qp;
  PetscBool      inside = PETSC_FALSE;
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: LameXSinker\n");
    opts_mu0     = 1.0;
    opts_mu1     = 1.0;
    opts_lambda0 = 1.0;
    opts_lambda1 = 1.0;
    opts_rad     = 0.05;
    opts_numinc  = 3;
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-mu0",&opts_mu0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-mu1",&opts_mu1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-lambda0",&opts_lambda0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-lambda1",&opts_lambda1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_r",&opts_rad,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-sinker_n",&opts_numinc,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: mu0 %1.4e\n",opts_mu0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: mu1 %1.4e\n",opts_mu1);
    PetscPrintf(PETSC_COMM_WORLD,"  params: lambda0 %1.4e\n",opts_lambda0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: lambda1 %1.4e\n",opts_lambda1);
    PetscPrintf(PETSC_COMM_WORLD,"  params: num sinkers %d\n",opts_numinc);
    PetscPrintf(PETSC_COMM_WORLD,"  params: sinker radius %1.4e\n",opts_rad);
    been_here = PETSC_TRUE;
  }

  if (opts_numinc > 8) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many sinkers");
  }
  if (opts_rad > 0.05) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Sinker Radius too big"); /* perhaps to cautious */
  }
  
  mu_qp     = opts_mu0;
  lambda_qp = opts_lambda0;
  rho_qp    = 1.0;

  for(i=0; i<opts_numinc; ++i){
#if NSD == 2
    const PetscReal d2 = (coor[0] - posx[i])*(coor[0] - posx[i]) + (coor[1] - posy[i])*(coor[1] - posy[i]);
#elif NSD == 3
    const PetscReal d2 = (coor[0] - posx[i])*(coor[0] - posx[i]) + (coor[1] - posy[i])*(coor[1] - posy[i]) + (coor[2] - posz[i])*(coor[2] - posz[i]);
#endif
    if  (d2 < opts_rad * opts_rad) {
      inside = PETSC_TRUE;
      break;
    }
  } 
  
  if (inside) {
    mu_qp     = opts_mu1;
    lambda_qp = opts_lambda1;
    rho_qp = 1.1;
  }
  
  if (mu) {
    *mu = mu_qp;
  }
  if (lambda) {
    *lambda = lambda_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}

/********************************************************************/
PetscErrorCode LameHomogeneous_EvaluateCoefficients(PetscReal coor[],PetscReal *mu,PetscReal *lambda, PetscReal Fu[],PetscReal Fp[])
{
  PetscErrorCode   ierr;
  static PetscReal opts_mu0,opts_lambda0;
  static PetscBool been_here = PETSC_FALSE;
  PetscReal        mu_qp,lambda_qp,rho_qp;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: LameHomogeneous\n");
    opts_mu0     = 1.0;
    opts_lambda0 = 1.0;
    ierr = PetscOptionsGetReal(NULL,NULL,"-mu0",&opts_mu0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-lambda0",&opts_lambda0,0);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: mu0 %1.4e\n",opts_mu0);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: lambda0 %1.4e\n",opts_lambda0);CHKERRQ(ierr);
    been_here = PETSC_TRUE;
  }
  
  mu_qp = opts_mu0;
  lambda_qp = opts_lambda0;
  rho_qp = 1.0;
  
  if (lambda) {
    *lambda = lambda_qp;
  }
  if (mu) {
    *mu = mu_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode Lame_EvaluateCoefficients(PetscReal coor[],PetscReal *mu,PetscReal *lambda,PetscReal Fu[],PetscReal Fp[])
{
  PetscErrorCode ierr;
  PetscInt model = DEFAULT_MODEL;
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-model",&model,NULL);CHKERRQ(ierr);
  
  switch (model) {
    case 2:
      ierr = LameXSinker_EvaluateCoefficients(coor,mu,lambda,Fu,Fp);CHKERRQ(ierr);
      break;
    case 6:
    case 8: 
    case 10: 
    case 12:
      ierr = LameOneSinker_EvaluateCoefficients(coor,mu,lambda,Fu,Fp);CHKERRQ(ierr);
      break;
    case 9: 
      ierr = LameHomogeneous_EvaluateCoefficients(coor,mu,lambda,Fu,Fp);CHKERRQ(ierr);
      break;
    default :
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Elasticity Model %d not implemented",model);
      break;
  }

  PetscFunctionReturn(0);
}
#endif

#ifndef LAME
/*
Stokes
*/
/*****************************************************************************/
PetscErrorCode StokesSolCx_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_eta0,opts_eta1,opts_xc;
  static PetscInt  opts_nz;
  static PetscBool been_here = PETSC_FALSE;
  PetscReal eta_qp;
  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: StokesSolCx\n");
    opts_eta0 = 1.0;
    opts_eta1 = 1.0;
    opts_xc   = 0.5;
    opts_nz   = 1;
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_xc",&opts_xc,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-solcx_nz",&opts_nz,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);
    PetscPrintf(PETSC_COMM_WORLD,"  params: xc   %1.4e\n",opts_xc);
    PetscPrintf(PETSC_COMM_WORLD,"  params: nz   %D\n",opts_nz);
    
    been_here = PETSC_TRUE;
  }
  
  eta_qp = opts_eta0;
  if (coor[0] > opts_xc) eta_qp = opts_eta1;
  if (eta) {
    *eta = eta_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = PetscSinReal(opts_nz*PETSC_PI*coor[1])*PetscCosReal(1.0*PETSC_PI*coor[0]);
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode StokesThreeSinker_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_eta0,opts_eta1,opts_rad;
  static PetscBool been_here = PETSC_FALSE;
  PetscReal        eta_qp,rho_qp;
  PetscBool        inside = PETSC_FALSE;
  PetscErrorCode   ierr;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: StokesThreeSinker\n");
    opts_eta0 = 1.0;
    opts_eta1 = 1.0;
    opts_rad = 0.1;
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_r",&opts_rad,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);
    PetscPrintf(PETSC_COMM_WORLD,"  params: rad  %1.4e\n",opts_rad);
    been_here = PETSC_TRUE;
  }
  
  eta_qp = opts_eta0;
  rho_qp = 1.0;

  {
    PetscReal sep2;
#if NSD == 2
    sep2 = (coor[0] - 0.27)*(coor[0] - 0.27) + (coor[1] - 0.63)*(coor[1] - 0.63);
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;

    sep2 = (coor[0] - 0.6)*(coor[0] - 0.6) + (coor[1] - 0.83)*(coor[1] - 0.83);
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;

    sep2 = (coor[0] - 0.7)*(coor[0] - 0.7) + (coor[1] - 0.33)*(coor[1] - 0.33);
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;
#elif NSD == 3
    sep2 = (coor[0] - 0.27)*(coor[0] - 0.27) + (coor[1] - 0.63)*(coor[1] - 0.63) + (coor[2] - 0.5)*(coor[2] - 0.5);
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;

    sep2 = (coor[0] - 0.6)*(coor[0] - 0.6) + (coor[1] - 0.83)*(coor[1] - 0.83) + (coor[2] - 0.5)*(coor[2] - 0.5);
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;

    sep2 = (coor[0] - 0.7)*(coor[0] - 0.7) + (coor[1] - 0.33)*(coor[1] - 0.33) + (coor[2] - 0.5)*(coor[2] - 0.5);
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;
#endif
  }
  
  if (inside) {
    eta_qp = opts_eta1;
    rho_qp = 1.1;
  }
  
  if (eta) {
    *eta = eta_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
  PetscFunctionReturn(0);
}
/*****************************************************************************/
PetscErrorCode StokesXSinker_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_eta0,opts_eta1,opts_rad;
  static PetscBool been_here = PETSC_FALSE;
  static PetscInt  opts_numinc;
  const PetscReal posx[8] = { 0.27 , 0.6  , 0.7  , 0.2 , 0.85  , 0.4 , 0.16 , 0.55 };
  const PetscReal posy[8] = { 0.63 , 0.83 , 0.33 , 0.2 , 0.65  , 0.3 , 0.84 , 0.54 };
#if NSD == 3
  const PetscReal posz[8] = { 0.50 , 0.40 , 0.30 , 0.70 , 0.65  , 0.4 , 0.8  , 0.50 };
#endif
  PetscReal eta_qp,rho_qp;
  PetscBool inside = PETSC_FALSE;
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: StokesXSinker\n");
    opts_eta0   = 1.0;
    opts_eta1   = 1.0;
    opts_rad    = 0.05;
    opts_numinc = 3;
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_r",&opts_rad,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-sinker_n",&opts_numinc,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);
    PetscPrintf(PETSC_COMM_WORLD,"  params: num sinkers %d\n",opts_numinc);
    PetscPrintf(PETSC_COMM_WORLD,"  params: sinker radius %1.4e\n",opts_rad);
    been_here = PETSC_TRUE;
  }

  if (opts_numinc > 8) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many sinkers");
  }
  if (opts_rad > 0.05) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Sinker Radius too big"); /* perhaps to cautious */
  }
  
  eta_qp = opts_eta0;
  rho_qp = 1.0;

  for(i=0; i<opts_numinc; ++i){
#if NSD == 2
    const PetscReal d2 = (coor[0] - posx[i])*(coor[0] - posx[i]) + (coor[1] - posy[i])*(coor[1] - posy[i]);
#elif NSD == 3
    const PetscReal d2 = (coor[0] - posx[i])*(coor[0] - posx[i]) + (coor[1] - posy[i])*(coor[1] - posy[i]) + (coor[2] - posz[i])*(coor[2] - posz[i]);
#endif
    if  (d2 < opts_rad * opts_rad) {
      inside = PETSC_TRUE;
      break;
    }
  } 
  
  if (inside) {
    eta_qp = opts_eta1;
    rho_qp = 1.1;
  }
  
  if (eta) {
    *eta = eta_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode StokesOneSinker_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_eta0,opts_eta1,opts_rad,opts_x,opts_y;
#if NSD==3
  static PetscReal opts_z;
#endif 
  static PetscBool been_here = PETSC_FALSE;
  PetscReal        eta_qp,rho_qp;
  PetscBool        inside = PETSC_FALSE;
	PetscErrorCode   ierr;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: StokesOneSinker\n");
    opts_eta0 = 1.0;
    opts_eta1 = 1.0;
    opts_rad = 0.25;
    opts_x    = 0.5;
    opts_y    = 0.5;
#if NSD==3
    opts_z    = 0.5;
#endif
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_r",&opts_rad,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_x",&opts_x,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_y",&opts_y,0);CHKERRQ(ierr);
#if NSD==3
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_z",&opts_z,0);CHKERRQ(ierr);
#endif
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: x %1.4e\n",opts_x);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: y %1.4e\n",opts_y);CHKERRQ(ierr);
#if NSD==3
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: z %1.4e\n",opts_z);CHKERRQ(ierr);
#endif
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  params: rad %1.4e\n",opts_rad);CHKERRQ(ierr);
    been_here = PETSC_TRUE;
  }
  
  eta_qp = opts_eta0;
  rho_qp = 1.0;
  
  {
    PetscReal sep2;
#if NSD == 2
    sep2 = (coor[0] - opts_x)*(coor[0] - opts_x) + (coor[1] - opts_y)*(coor[1] - opts_y);
#elif NSD == 3
    sep2 = (coor[0] - opts_x)*(coor[0] - opts_x) + (coor[1] - opts_y)*(coor[1] - opts_y) + (coor[2] - opts_z)*(coor[2] - opts_z);
#endif
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;
  }
  
  if (inside) {
    eta_qp = opts_eta1;
    rho_qp = 1.1;
  }
  
  if (eta) {
    *eta = eta_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}

/*****************************************************************************/
#if NSD == 3
/*
   (3d only)
 Generates location of spherical inclusions in the domain [0,Lx]x[0,Ly]x[0,Lz] assuming the
 inclusions cannot overlap, and they must be some min distance from the boundary edges.
 */
PetscErrorCode GenerateInclusionOrigins(PetscInt ninclusions,PetscReal rmax,PetscReal Lx,PetscReal Ly,PetscReal Lz, PetscReal min_sep_wall,PetscReal min_sep_region,PetscReal **_pos)
{
	PetscReal   *pos;
	PetscInt    p,found=0,overlap,attempt,loops=0;
	PetscInt    max_attempts = 50000;
 
  PetscOptionsGetInt(NULL,NULL,"-max_attempts",&max_attempts,NULL);
  PetscPrintf(PETSC_COMM_WORLD,"# GenerateInclusionOrigins:\n");
  PetscPrintf(PETSC_COMM_WORLD,"#   nregions       %D\n",ninclusions);
  PetscPrintf(PETSC_COMM_WORLD,"#   radius         %1.4e\n",rmax);
  PetscPrintf(PETSC_COMM_WORLD,"#   Lx,Ly,Lz       %1.4e,%1.4e,%1.4e\n",Lx,Ly,Lz);
  PetscPrintf(PETSC_COMM_WORLD,"#   min_sep        %1.4e (in terms of region radii)\n",min_sep_region/rmax);
  PetscPrintf(PETSC_COMM_WORLD,"#   min_wall sep   %1.4e (in terms of region radii)\n",min_sep_wall/rmax);
	
	PetscMalloc1(NSD*ninclusions,&pos);
	
  if (ninclusions >= max_attempts) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Error: max attempts is < ninclusions. Increase -max_attempts");
  }
  
	srand(0);
	loops = 0;
  
START_INCLUSION:
	/* if we didn't locate ninclusions in max_attempts, we throw ALL previous generated inclusions away and start again */
  
	++loops;
	found   = 0;
	attempt = 0;
  
	while (found < ninclusions) {
		PetscReal dx,dy,dz,range[NSD],cp[NSD],xp,yp,zp;
    
		if (attempt == max_attempts) { goto START_INCLUSION; }
		
    /* generate random centoid within [0,Lx]x[0,Ly]x[0,Lz] */
		xp = (PetscReal)( rand()/( (double)RAND_MAX ) );
		yp = (PetscReal)( rand()/( (double)RAND_MAX ) );
		zp = (PetscReal)( rand()/( (double)RAND_MAX ) );
    
		xp = xp * Lx;
		yp = yp * Ly;
		zp = zp * Lz;
		
    ++attempt;
    
    /* check of sphere is further from wall than min_sep_wall * radius */
		dx = min_sep_wall*rmax;
		range[0] = xp - dx;
		if (range[0] < 0.0) { continue; }
		range[0] = xp + dx;
		if (range[0] > Lx) { continue; }
		
		dy = min_sep_wall*rmax;
		range[1] = yp - dy;
		if (range[1] < 0.0) { continue; }
		range[1] = yp + dy;
		if (range[1] > Ly) { continue; }
    
		dz = min_sep_wall*rmax;
		range[2] = zp - dz;
		if (range[2] < 0.0) { continue; }
		range[2] = zp + dz;
		if (range[2] > Lz) { continue; }

		/* check if spheres overlap with each other */
		cp[0] = xp;
		cp[1] = yp;
		cp[2] = zp;
		overlap = 0;
		for (p=0; p<found; p++) {
			PetscReal sep;
			
			sep = PetscSqrtReal(  (pos[NSD*p+0]-cp[0])*(pos[NSD*p+0]-cp[0])
                          + (pos[NSD*p+1]-cp[1])*(pos[NSD*p+1]-cp[1])
                          + (pos[NSD*p+2]-cp[2])*(pos[NSD*p+2]-cp[2]) );
      
      /* check centroids are min_sep_region * radius distance from each other */
			if (sep < 2.0*rmax + min_sep_region*rmax) {
				overlap = 1;
				break;
			}
		}
		if (overlap == 1) { continue; }
		
		pos[NSD*found+0] = xp;
		pos[NSD*found+1] = yp;
		pos[NSD*found+2] = zp;
		
    ++found;
	}
  
	PetscPrintf(PETSC_COMM_WORLD,"# GenerateInclusionOrigins: performed %D loops: made %D attempts and correctly defined %D inclusions\n",loops,attempt,ninclusions);
	
	*_pos   = pos;
	PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode SinkerPtatin_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
	PetscErrorCode   ierr;
  static PetscReal opts_eta0,opts_eta1,opts_rad,*centroidpos;
  static PetscBool been_here = PETSC_FALSE;
  static PetscInt  opts_numinc;
  PetscReal        eta_qp,rho_qp;
  PetscBool        inside;
  PetscInt         k;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: SinkerPtatin\n");
    opts_eta0 = 1.0;
    opts_eta1 = 1.1;
    opts_rad = 0.05;
    opts_numinc = 3;
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-sinker_r",&opts_rad,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-sinker_n",&opts_numinc,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);
    
    ierr = GenerateInclusionOrigins(opts_numinc,opts_rad,1.0,1.0,1.0,1.5,1.5,&centroidpos);CHKERRQ(ierr);
    
    been_here = PETSC_TRUE;
  }
  
  eta_qp = opts_eta0;
  rho_qp = 1.0;
  
  inside = PETSC_FALSE;
  for (k=0; k<opts_numinc; ++k) {
    PetscReal *pos_k;
    PetscReal sep2;
    
    pos_k = &centroidpos[NSD*k];
    
#if NSD == 2
    sep2 = (coor[0] - pos_k[0])*(coor[0] - pos_k[0]) + (coor[1] - pos_k[1])*(coor[1] - pos_k[1]);
#elif NSD == 3
    sep2 = (coor[0] - pos_k[0])*(coor[0] - pos_k[0]) + (coor[1] - pos_k[1])*(coor[1] - pos_k[1]) + (coor[2] - pos_k[2])*(coor[2] - pos_k[2]);
#endif
    if (sep2 < opts_rad*opts_rad) inside = PETSC_TRUE;
    
    if (inside) break;
  }
  
  if (inside) {
    eta_qp = opts_eta1;
    rho_qp = 1.1;
  }
  
  if (eta) {
    *eta = eta_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = -rho_qp;
#if NSD == 3
    Fu[2] = 0.0;
#endif
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}
#endif

/*****************************************************************************/
PetscErrorCode StokesSolCx3d_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_eta0,opts_eta1,opts_xc;
  static PetscInt  opts_nz,opts_nz2; 
  static PetscBool been_here = PETSC_FALSE;
  PetscReal        eta_qp;
	PetscErrorCode   ierr;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: StokesSolCx3d\n");
    opts_eta0 = 1.0;
    opts_eta1 = 1.0;
    opts_xc   = 0.5;
    opts_nz   = 1;
    opts_nz2  = 1;
    
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-solcx_xc",&opts_xc,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-solcx_nz",&opts_nz,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);
    PetscPrintf(PETSC_COMM_WORLD,"  params: xc   %1.4e\n",opts_xc);
    PetscPrintf(PETSC_COMM_WORLD,"  params: nz   %D\n",opts_nz);
    PetscPrintf(PETSC_COMM_WORLD,"  params: nz2  %D\n",opts_nz2);
    
    been_here = PETSC_TRUE;
  }
  
  eta_qp = opts_eta0;
  if (coor[0] > opts_xc) eta_qp = opts_eta1;
  if (eta) {
    *eta = eta_qp;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = PetscSinReal(opts_nz*PETSC_PI*coor[1])*PetscCosReal(1.0*PETSC_PI*coor[0])*PetscSinReal(opts_nz2*PETSC_PI*coor[2]);
    Fu[2] = 0.0;
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}

#if NSD==2
/*****************************************************************************/
static PetscErrorCode StokesMMS1_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscBool been_here = PETSC_FALSE;
  
  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: StokesMMS1\n");
    been_here = PETSC_TRUE;
  }

  /* This is a solution, taken from Elman/Silvester/Wathen,
     of converging 2D flow with stream function 
     phi(x,y) = 5xy^4 - x^5
     
     u^x(x,y) = 20xy^3
     u^y(x,y) = 5 (x^4 - y^4)
       p(x,y) = 60 x^2y - 20 y^3 + constant
     
     It satisfies 

       \nabla^2 u + \grad p = 0,
                    \div  p = 0

     in the domain, hence here eta = 1 and the forcing
     is identically zero. Note that the Dirichlet boundary
     conditions are not zero. See the corresponding function
     which defines those.
     
     Note that this problem has a constant-pressure nullspace. 
     When a representative needs to be chosen, we will opt for the
     one which has zero mean of the discretized values (project
     out the component in the constant-pressure direction). In general
     this doesn't correspond to any mesh-independent choice of the 
     constant above.
     */

  if (eta) {
    *eta = 1.0;
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = 0.0;
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}
#endif

#if NSD==3
/*****************************************************************************/
/*
This is to experiment with a case in which it's important to have a preconditioner/smoother with strong coupling in the x direction.
   */
static PetscErrorCode StokesPseudoIce_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
  static PetscReal opts_eta0,opts_eta1;
  static PetscBool been_here = PETSC_FALSE;
  PetscReal        size_x;
	PetscErrorCode   ierr;

  PetscFunctionBeginUser;
  if (!been_here) {
    PetscPrintf(PETSC_COMM_WORLD,"ModelType: PseudoIce\n");
    opts_eta0 = 1.0;
    opts_eta1 = 10000.0;

    ierr = PetscOptionsGetReal(NULL,NULL,"-eta0",&opts_eta0,0);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-eta1",&opts_eta1,0);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta0 %1.4e\n",opts_eta0);
    PetscPrintf(PETSC_COMM_WORLD,"  params: eta1 %1.4e\n",opts_eta1);

    been_here = PETSC_TRUE;
  }

  size_x = 1.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-size_x",&size_x,NULL);CHKERRQ(ierr);/* This is a hack */
  if (eta) { 
    PetscReal xrel = coor[0]/size_x;
    *eta = xrel*opts_eta0 + (1-xrel)*opts_eta1; /* simple linear interp. This doesn't actually hit the boundary values */
  }
  if (Fu) {
    Fu[0] = 0.0;
    Fu[1] = 0.0; 
    Fu[2] = 1.0; /* Uniform forcing in z direction */
  }
  if (Fp) {
    Fp[0] = 0.0;
  }
	PetscFunctionReturn(0);
}
#endif
/*****************************************************************************/
PetscErrorCode Stokes_EvaluateCoefficients(PetscReal coor[],PetscReal *eta,PetscReal Fu[],PetscReal Fp[])
{
	PetscErrorCode ierr;
  PetscInt model = DEFAULT_MODEL;
  
  PetscOptionsGetInt(NULL,NULL,"-model",&model,NULL);
  
  switch (model) {
    case 0:
      ierr = StokesSolCx_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
    case 1:
      ierr = StokesThreeSinker_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
    case 2:
      ierr = StokesXSinker_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
#if NSD == 3
    case 5:
      ierr = StokesSolCx3d_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
#endif
    case 6:
      ierr = StokesOneSinker_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
#if NSD == 3
    case 7:
      ierr = SinkerPtatin_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
    case 11:
      ierr =  StokesPseudoIce_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
#endif
#if NSD == 2
    case 101:
      ierr = StokesMMS1_EvaluateCoefficients(coor,eta,Fu,Fp);CHKERRQ(ierr);
      break;
#endif
    default :
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Stokes Model %d not implemented",model);
      break;
  }
  PetscFunctionReturn(0);
}

#endif

/********************************************************************/
#ifndef LAME
#if NSD==2


PetscErrorCode StokesMMS1_ComputeReferenceSolution(DM dm_saddle,Vec *Xref)
{
  PetscErrorCode ierr;
  DM             dmv,dmp;
  Vec            Xv,Xp;
  PetscReal      c[NSD];
  PetscInt       i,j,d;

  PetscFunctionBeginUser;

  /* Create the vector */
  ierr = DMCreateGlobalVector(dm_saddle,Xref);CHKERRQ(ierr);

  /* Get Velocity and Pressure DMs */
  ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);

  /* Get Velocity and Pressure Vectors */
  ierr = DMCompositeGetAccess(dm_saddle,*Xref,&Xv,&Xp);CHKERRQ(ierr);

  /* Fill in Velocity values */
  {
    PetscInt         si,sj,sk,ni,nj,nk,M,N,P;
    const PetscReal  *LA_coord_l;
    PetscScalar      ***arr;
    Vec              coord_l; 

    ierr = DMGetCoordinatesLocal(dmv,&coord_l);CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmv,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dmv,Xv,&arr);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coord_l,&LA_coord_l);CHKERRQ(ierr);
    for(i=si;i<si+ni;++i){
      for(j=sj;j<sj+nj;++j){
        for (d=0; d<NSD; ++d){
          c[d] = LA_coord_l[NSD*(i+j*ni)+d];
        }
        for(d=0; d<NSD; ++d){
          PetscScalar val = d == 0 ? MMS1_SOLX(c[0],c[1]) : MMS1_SOLY(c[0],c[1]);
          arr[j][i][d] = val;
          // DEBUG
#if 0
          {
            PetscMPIInt rank;
            ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] debugsolv. i=%d,j=%d,d=%d,c=[%g,%g],val=%g \n",rank,i,j,d,c[0],c[1],val);CHKERRQ(ierr);
          }
#endif
        }
      }
    }
    ierr = VecRestoreArrayRead(coord_l,&LA_coord_l);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dmv,Xv,&arr);CHKERRQ(ierr);
  }

  /* Repeat for Pressure */
  {
    PetscInt         si,sj,sk,ni,nj,nk,M,N,P;
    const PetscReal  *LA_coord_l;
    Vec              coord_l; 
    PetscScalar      **arr;

    ierr = DMGetCoordinatesLocal(dmp,&coord_l);CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmp,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    ierr = DMDAGetCorners(dmp,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dmp,Xp,&arr);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coord_l,&LA_coord_l);CHKERRQ(ierr);
    for(i=si;i<si+ni;++i){
      for(j=sj;j<sj+nj;++j){
        for (d=0; d<NSD; ++d){
          c[d] = LA_coord_l[NSD*(i+j*ni)+d];
        }
        PetscScalar val = MMS1_SOLP(c[0],c[1]);
        arr[j][i] = val;
        // DEBUG
#if 0
        {
          PetscMPIInt rank;
          ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] debugsolp. i=%d,j=%d,d=%d,c=[%g,%g],val=%g \n",rank,i,j,d,c[0],c[1],val);CHKERRQ(ierr);
        }
#endif
      }
    }
    ierr = VecRestoreArrayRead(coord_l,&LA_coord_l);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dmp,Xp,&arr);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccess(dm_saddle,*Xref,&Xv,&Xp);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#endif
#endif

/********************************************************************/
/* Compute Reference solution (for MMS) tests */
PetscErrorCode ComputeReferenceSolution(DM dm_saddle,Vec *Xref)
{
  PetscFunctionBeginUser;
	PetscErrorCode ierr;
  PetscInt       model = DEFAULT_MODEL;
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-model",&model,NULL);CHKERRQ(ierr);
  
  switch (model) {
#ifndef LAME
#if NSD==2
    case 101:
      ierr = StokesMMS1_ComputeReferenceSolution(dm_saddle,Xref);CHKERRQ(ierr);
      break;
#endif
#endif
    default:    
      *Xref = NULL;
  }
  PetscFunctionReturn(0);

}
