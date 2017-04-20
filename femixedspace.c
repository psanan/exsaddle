#include "femixedspace.h"
#include "models.h"

#include <petscdmda.h>
#include <petscdmcomposite.h>
#include <petsc/private/dmimpl.h> 

/*****************************************************************************/
PetscErrorCode SaddlePreallocation_MPI(DM dmv,DM dmp,PetscInt nnz[],PetscInt onnz[],PetscInt max_nnz_per_row_diag, PetscInt max_nnz_per_row_offdiag)
{
  PetscErrorCode ierr;
  PetscInt       idx,ii,jj,i,j,si,sj,ni,nj,nnz_per_row=0,nnz_per_row_diag,nnz_per_row_offdiag,u_ni,u_nj,Ni,Nj;
#if NSD == 2
  PetscBool      vi,vj;
#elif NSD == 3
  PetscInt       kk,k,sk,nk,u_nk,Nk,nmod;
#endif
  PetscBool      sub_domain_boundary;

#if NSD == 2
  /* fill velocity rows */
  ierr = DMDAGetInfo(dmv,NULL,&Ni,&Nj,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dmv,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
  u_ni = ni;
  u_nj = nj;
  for (j=sj; j<sj+nj; ++j) {
    for (i=si; i<si+ni; ++i) {
      ii = i - si;
      jj = j - sj;
      idx = ii + jj * ni;
      
      vi = (i%2 == 0);
      vj = (j%2 == 0);
      if (vi && vj) { /* vertex */
        nnz_per_row = 2*(5*5) + 9;
      } else if (vi || vj) { /* edge */
        nnz_per_row = 2*(5*3) + 6;
      } else { /* center */
        nnz_per_row = 2*(3*3) + 4;
      }
      nnz_per_row_diag = nnz_per_row > max_nnz_per_row_diag ? max_nnz_per_row_diag : nnz_per_row;
      
      nnz[2*idx+0] = nnz_per_row_diag;
      nnz[2*idx+1] = nnz_per_row_diag;

      sub_domain_boundary = PETSC_FALSE;
      if (ii < 3) { sub_domain_boundary = PETSC_TRUE; }
      if (jj < 3) { sub_domain_boundary = PETSC_TRUE; }
      if (ii >= ni-3) { sub_domain_boundary = PETSC_TRUE; }
      if (jj >= nj-3) { sub_domain_boundary = PETSC_TRUE; }
     
      if (sub_domain_boundary) {
        nnz_per_row_offdiag = nnz_per_row > max_nnz_per_row_offdiag ? max_nnz_per_row_offdiag : nnz_per_row;
        onnz[2*idx+0] = nnz_per_row_offdiag;
        onnz[2*idx+1] = nnz_per_row_offdiag;
      }

    }
  }
  
  /* fill pressure rows */
  ierr = DMDAGetInfo(dmp,NULL,&Ni,&Nj,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dmp,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
  for (j=sj; j<sj+nj; ++j) {
    for (i=si; i<si+ni; ++i) {
      ii = i - si;
      jj = j - sj;
      idx = ii + jj * ni   +  u_ni*u_nj*2;
      
      nnz_per_row = 2*(5*5) + 9;
      nnz_per_row_diag = nnz_per_row > max_nnz_per_row_diag ? max_nnz_per_row_diag : nnz_per_row;
      nnz[idx] = nnz_per_row_diag;

      sub_domain_boundary = PETSC_FALSE;
      
      if (ii < 2) { sub_domain_boundary = PETSC_TRUE; }
      if (jj < 2) { sub_domain_boundary = PETSC_TRUE; }
      if (ii >= ni-2) { sub_domain_boundary = PETSC_TRUE; }
      if (jj >= nj-2) { sub_domain_boundary = PETSC_TRUE; }
     
      if (sub_domain_boundary) {
        nnz_per_row_offdiag = nnz_per_row > max_nnz_per_row_offdiag ? max_nnz_per_row_offdiag : nnz_per_row;
        onnz[idx] = nnz_per_row_offdiag;
      }

    }
  }
#elif NSD == 3 
  /* fill velocity rows */
  ierr = DMDAGetInfo(dmv,NULL,&Ni,&Nj,&Nk,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  u_ni = ni;
  u_nj = nj;
  u_nk = nk;
  for (k=sk; k<sk+nk; ++k) {
    for (j=sj; j<sj+nj; ++j) {
      for (i=si; i<si+ni; ++i) {
        ii = i - si;
        jj = j - sj;
        kk = k - sk;
        idx = ii + jj * ni + kk * ni * nj; 

        nmod = (i%2) + (j%2) + (k%2); /* How many dims align with the Q1 lattice */
        switch (nmod) {
          case 0 : /* vertex */
            nnz_per_row = 3*(5*5*5) + (3*3*3);
            break;
          case 1 : /* edge */
            nnz_per_row = 3*(5*5*3) + (3*3*2); 
            break;
          case 2 : /* face center */
            nnz_per_row = 3*(5*3*3) + (3*2*2); 
            break;
          case 3 : /* cell center */
            nnz_per_row = 3*(3*3*3) + (2*2*2);
            break;
        }
        nnz_per_row_diag = nnz_per_row > max_nnz_per_row_diag ? max_nnz_per_row_diag : nnz_per_row;

        nnz[3*idx+0] = nnz_per_row_diag;
        nnz[3*idx+1] = nnz_per_row_diag;
        nnz[3*idx+2] = nnz_per_row_diag;

        sub_domain_boundary = PETSC_FALSE;
        if (ii < 3) { sub_domain_boundary = PETSC_TRUE; }
        if (jj < 3) { sub_domain_boundary = PETSC_TRUE; }
        if (kk < 3) { sub_domain_boundary = PETSC_TRUE; }
        if (ii >= ni-3) { sub_domain_boundary = PETSC_TRUE; }
        if (jj >= nj-3) { sub_domain_boundary = PETSC_TRUE; }
        if (kk >= nk-3) { sub_domain_boundary = PETSC_TRUE; }

        if (sub_domain_boundary) {
          nnz_per_row_offdiag = nnz_per_row > max_nnz_per_row_offdiag ? max_nnz_per_row_offdiag : nnz_per_row;
          onnz[3*idx+0] = nnz_per_row_offdiag;
          onnz[3*idx+1] = nnz_per_row_offdiag;
          onnz[3*idx+2] = nnz_per_row_offdiag;
        }

      }
    }
  }
  
  /* fill pressure rows */
  ierr = DMDAGetInfo(dmp,NULL,&Ni,&Nj,&Nk,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dmp,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  for (k=sk; k<sk+nk; ++k) {
    for (j=sj; j<sj+nj; ++j) {
      for (i=si; i<si+ni; ++i) {
        ii = i - si;
        jj = j - sj;
        kk = k - sk;
        idx = ii + jj * ni + kk * ni * nj   +  u_ni*u_nj*u_nk*3;

        nnz_per_row = 3*(5*5*5) + (3*3*3);

        nnz_per_row_diag = nnz_per_row > max_nnz_per_row_diag ? max_nnz_per_row_diag : nnz_per_row;
        nnz[idx] = nnz_per_row_diag;

        sub_domain_boundary = PETSC_FALSE;

        if (ii < 2) { sub_domain_boundary = PETSC_TRUE; }
        if (jj < 2) { sub_domain_boundary = PETSC_TRUE; }
        if (kk < 2) { sub_domain_boundary = PETSC_TRUE; }
        if (ii >= ni-2) { sub_domain_boundary = PETSC_TRUE; }
        if (jj >= nj-2) { sub_domain_boundary = PETSC_TRUE; }
        if (kk >= nk-2) { sub_domain_boundary = PETSC_TRUE; }

        if (sub_domain_boundary) {
          nnz_per_row_offdiag = nnz_per_row > max_nnz_per_row_offdiag ? max_nnz_per_row_offdiag : nnz_per_row;
          onnz[idx] = nnz_per_row_offdiag;
        }

      }
    }
  }
#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
static PetscErrorCode SaddlePreallocation_SEQ(DM dmv,DM dmp,PetscInt nnz[],PetscInt max_nnz_per_row)
{
  PetscErrorCode ierr;
  PetscInt idx,ii,jj,i,j,si,sj,ni,nj,nnz_per_row=0,u_ni,u_nj;
#if NSD == 2
  PetscBool vi,vj;
#elif NSD == 3
  PetscInt kk,k,sk,nk,u_nk,nmod;
#endif
  
  /* fill velocity rows */
#if NSD == 2
  ierr = DMDAGetCorners(dmv,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
  u_ni = ni;
  u_nj = nj;
  for (j=sj; j<sj+nj; ++j) {
    for (i=si; i<si+ni; ++i) {
      ii = i - si;
      jj = j - sj;
      idx = ii + jj * ni;
      
      vi = vj = PETSC_FALSE;
      if (i%2 == 0) { vi = PETSC_TRUE; } 
      if (j%2 == 0) { vj = PETSC_TRUE; }
      if (vi && vj) { /* vertex */
        nnz_per_row = 2*(5*5) + 9;
      } else if (vi || vj) { /* edge */
        nnz_per_row = 2*(5*3) + 6;
      } else { /* center */
        nnz_per_row = 2*(3*3) + 4;
      }
      if (nnz_per_row > max_nnz_per_row ) nnz_per_row = max_nnz_per_row;

      nnz[2*idx+0] = nnz_per_row;
      nnz[2*idx+1] = nnz_per_row;
    }
  }
#elif NSD == 3
  ierr = DMDAGetCorners(dmv,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  u_ni = ni;
  u_nj = nj;
  u_nk = nk;
  for (k=sk; k<sk+nk; ++k) {
    for (j=sj; j<sj+nj; ++j) {
      for (i=si; i<si+ni; ++i) {
        ii = i - si;
        jj = j - sj;
        kk = k - sk;
        idx = ii + jj * ni + kk * ni * nj;
        nmod = (i%2) + (j%2) + (k%2); /* How many dims align with the Q1 lattice */
        switch (nmod) {
          case 0: /* vertex */
            nnz_per_row = 3*(5*5*5) + (3*3*3);
            break;
          case 1 : /* edge */
            nnz_per_row = 3*(5*5*3) + (3*3*2); 
            break;
          case 2 : /* face center */
            nnz_per_row = 3*(5*3*3) + (3*2*2); 
            break;
          case 3 : /* cell center */
            nnz_per_row = 3*(3*3*3) + (2*2*2);
            break;
        }
        if (nnz_per_row > max_nnz_per_row ) nnz_per_row = max_nnz_per_row;
        nnz[3*idx+0] = nnz_per_row;
        nnz[3*idx+1] = nnz_per_row;
        nnz[3*idx+2] = nnz_per_row;
      }
    }
  }
#endif

  /* fill pressure rows */
#if NSD == 2
  ierr = DMDAGetCorners(dmp,&si,&sj,NULL,&ni,&nj,NULL);CHKERRQ(ierr);
  for (j=sj; j<sj+nj; ++j) {
    for (i=si; i<si+ni; ++i) {
      ii = i - si;
      jj = j - sj;
      idx = ii + jj * ni   +  u_ni*u_nj*2;
      
      nnz_per_row = 2*(5*5) + 9;
      if (nnz_per_row > max_nnz_per_row ) nnz_per_row = max_nnz_per_row;
      nnz[idx] = nnz_per_row;
    }
  }
#elif NSD == 3
  ierr = DMDAGetCorners(dmp,&si,&sj,&sk,&ni,&nj,&nk);CHKERRQ(ierr);
  for (k=sk; k<sk+nk; ++k) {
    for (j=sj; j<sj+nj; ++j) {
      for (i=si; i<si+ni; ++i) {
        ii = i - si;
        jj = j - sj;
        kk = k - sk;
        idx = ii + jj * ni + kk * ni * nj   +  u_ni*u_nj*u_nk*3; 

        nnz_per_row = 3*(5*5*5) + (3*3*3);
        if (nnz_per_row > max_nnz_per_row ) nnz_per_row = max_nnz_per_row;
        nnz[idx] = nnz_per_row;
      }
    }
  }
#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
/* Note: we are lazy here and maintain separate 2d and 3d versions of this
         routing. In general this should be avoided */
#if NSD == 2
static PetscErrorCode DMDAPatchCreateGlobalIS2d(DM dmv,DM dmp,MatStencil lowerQ2,MatStencil upperQ2,MatStencil lowerQ1,MatStencil upperQ1,IS *is)
{
  PetscErrorCode ierr;
  const PetscInt *oai,*oaj; /* ownership array */
  const PetscInt *oaiQ1,*oajQ1; /* ownership array */
  PetscInt       *ri,*rj; /* ranges */
  PetscInt       *riQ1,*rjQ1; /* ranges */
  PetscInt       *offset;
  PetscInt       *offsetQ1;
  PetscInt       pi,pj,pM,pN;
  PetscInt       li,lj;
  PetscInt       i,j,cnt,sum,length,nudofs;
  PetscInt       lidx,idx;
  PetscInt       *indices;
  
  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(dmv,NULL,NULL,NULL,NULL,&pM,&pN,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(dmv,&oai,&oaj,NULL);CHKERRQ(ierr); /* Not multiplied by dof/node = 2 */
  ierr = DMDAGetOwnershipRanges(dmp,&oaiQ1,&oajQ1,NULL);CHKERRQ(ierr);
  
  PetscMalloc1(pM+1,&ri);
  PetscMalloc1(pN+1,&rj);
  PetscMalloc1(pM*pN,&offset);

  PetscMalloc1(pM+1,&riQ1);
  PetscMalloc1(pN+1,&rjQ1);
  PetscMalloc1(pM*pN,&offsetQ1);

  /* Compute the total number of local dof in the elements owned by this rank */
  length = 2*(upperQ2.i-lowerQ2.i)*(upperQ2.j-lowerQ2.j) + (upperQ1.i-lowerQ1.i)*(upperQ1.j-lowerQ1.j);
  PetscMalloc1(length,&indices);
 
  /* Accumulate the number of dof (for each component) in each direction for the Q2 elements */
  sum = 0; for (pi=0; pi<pM; ++pi) { ri[pi] = sum; sum += oai[pi]; } ri[pi] = sum;
  sum = 0; for (pj=0; pj<pN; ++pj) { rj[pj] = sum; sum += oaj[pj]; } rj[pj] = sum;

  /* Accumulate the number of dof in each direction for the Q1 elements */
  sum = 0; for (pi=0; pi<pM; ++pi) { riQ1[pi] = sum; sum += oaiQ1[pi]; } riQ1[pi] = sum;
  sum = 0; for (pj=0; pj<pN; ++pj) { rjQ1[pj] = sum; sum += oajQ1[pj]; } rjQ1[pj] = sum;

  /* Compute an offset for each rank that gives the number of velocity nodes (not dof) 
     on preceding ranks */
  sum = 0;
  cnt = 0;
  { 
    /* This defines the ordering of offset */
    for (pj=0; pj<pN; ++pj) {
      for (pi=0; pi<pM; ++pi) {
        offset[cnt] = sum;
        sum += oai[pi] * oaj[pj];
        ++cnt;
      }
    }
  }

  /* Compute an offset for each rank that gives the number of pressure nodes on preceding ranks */
  sum = 0;
  cnt = 0;
  {
    for (pj=0; pj<pN; ++pj) {
      for (pi=0; pi<pM; ++pi) {
        offsetQ1[cnt] = sum;
        sum += oaiQ1[pi] * oajQ1[pj];
        ++cnt;
      }
    }
  }

  /* Compute the global indices for the velocity nodes on this subdomain */ 
  cnt = 0;
  {
    
    for (j=lowerQ2.j; j<upperQ2.j; ++j) {

      /* Compute this point's processor index in the j direction 
         (might now be this proc because of overlap) 
        */
      for (pj=0; pj<pN; ++pj) { if ((j >= rj[pj]) && (j < rj[pj+1])) { break; } }
      
      for (i=lowerQ2.i; i<upperQ2.i; ++i) {
       /* Compute this point's processor index in the i direction 
         (might now be this proc because of overlap) to
         determine which rank pi,pj,pk contains point i,j,k 
        */
        for (pi=0; pi<pM; ++pi) { if ((i >= ri[pi]) && (i < ri[pi+1])) { break; } }
        
        /* shift by the start index of pi*/
        li = i - ri[pi];
        lj = j - rj[pj];
        
        lidx = li + lj * oai[pi]; /* velocity point index (doesn't take 2 dof into account) */
        idx = 2*lidx + 2*offset[pi + pj * pM] + offsetQ1[pi + pj * pM]; 
        
        indices[2*cnt+0] = idx + 0;
        indices[2*cnt+1] = idx + 1;
        
        ++cnt;
      }
    }
  }
  nudofs = 2 * cnt;

  cnt = nudofs;
  {
    
    for (j=lowerQ1.j; j<upperQ1.j; ++j) {
      for (pj=0; pj<pN; ++pj) { if ((j >= rjQ1[pj]) && (j < rjQ1[pj+1])) { break; } }
      
      for (i=lowerQ1.i; i<upperQ1.i; ++i) {
        /* determine which rank pi,pj,pk contains point i,j,k */
        for (pi=0; pi<pM; ++pi) { if ((i >= riQ1[pi]) && (i < riQ1[pi+1])) { break; } }
        
        /* shift by the start index of pi */
        li = i - riQ1[pi];
        lj = j - rjQ1[pj];
        
        lidx = li + lj * oaiQ1[pi];
        
        /* shift by TOTAL uv-dofs p-dofs + however many uv dofs there are on this sub-domain */
        idx = lidx + 2*offset[pi + pj * pM] + offsetQ1[pi + pj * pM] + 2*oai[pi]*oaj[pj];
        
        indices[cnt] = idx;
        
        ++cnt;
      }
    }
  }
  if (cnt != length) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Counts do not match");

  
  ierr = ISCreateGeneral(PETSC_COMM_SELF,length,indices,PETSC_COPY_VALUES,is);CHKERRQ(ierr);
  
  PetscFree(indices);
  PetscFree(offset);
  PetscFree(rj);
  PetscFree(ri);
  PetscFree(offsetQ1);
  PetscFree(rjQ1);
  PetscFree(riQ1);
  
  PetscFunctionReturn(0);
}
#elif NSD == 3
PetscErrorCode DMDAPatchCreateGlobalIS3d(DM dmv,DM dmp,MatStencil lowerQ2,MatStencil upperQ2,MatStencil lowerQ1,MatStencil upperQ1,IS *is)
{
  PetscErrorCode ierr;
  const PetscInt *oai,*oaj,*oak; /* ownership array */
  const PetscInt *oaiQ1,*oajQ1,*oakQ1; /* ownership array */
  PetscInt       *ri,*rj,*rk; /* ranges */
  PetscInt       *riQ1,*rjQ1,*rkQ1; /* ranges */
  PetscInt       *offset;
  PetscInt       *offsetQ1;
  PetscInt       pi,pj,pk,pM,pN,pP;
  PetscInt       li,lj,lk;
  PetscInt       i,j,k,cnt,sum,length,nudofs;
  PetscInt       lidx,idx;
  PetscInt       *indices;
  
  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(dmv,NULL,NULL,NULL,NULL,&pM,&pN,&pP,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(dmv,&oai,&oaj,&oak);CHKERRQ(ierr); /* Not multiplied by dof/node = 3 */
  ierr = DMDAGetOwnershipRanges(dmp,&oaiQ1,&oajQ1,&oakQ1);CHKERRQ(ierr);

  PetscMalloc1(pM+1,&ri);
  PetscMalloc1(pN+1,&rj);
  PetscMalloc1(pP+1,&rk);
  PetscMalloc1(pM*pN*pP,&offset);

  PetscMalloc1(pM+1,&riQ1);
  PetscMalloc1(pN+1,&rjQ1);
  PetscMalloc1(pP+1,&rkQ1);
  PetscMalloc1(pM*pN*pP,&offsetQ1);

  /* Compute the total number of local dof in the elements owned by this rank */
  length = 3*(upperQ2.i-lowerQ2.i)*(upperQ2.j-lowerQ2.j)*(upperQ2.k-lowerQ2.k) + (upperQ1.i-lowerQ1.i)*(upperQ1.j-lowerQ1.j)*(upperQ1.k-lowerQ1.k);
  PetscMalloc1(length,&indices);
 
  /* Accumulate the number of dof (for each component) in each direction for the Q2 elements */
  sum = 0; for (pi=0; pi<pM; ++pi) { ri[pi] = sum; sum += oai[pi]; } ri[pi] = sum;
  sum = 0; for (pj=0; pj<pN; ++pj) { rj[pj] = sum; sum += oaj[pj]; } rj[pj] = sum;
  sum = 0; for (pk=0; pk<pP; ++pk) { rk[pk] = sum; sum += oak[pk]; } rk[pk] = sum;

  /* Accumulate the number of dof in each direction for the Q1 elements */
  sum = 0; for (pi=0; pi<pM; ++pi) { riQ1[pi] = sum; sum += oaiQ1[pi]; } riQ1[pi] = sum;
  sum = 0; for (pj=0; pj<pN; ++pj) { rjQ1[pj] = sum; sum += oajQ1[pj]; } rjQ1[pj] = sum;
  sum = 0; for (pk=0; pk<pP; ++pk) { rkQ1[pk] = sum; sum += oakQ1[pk]; } rkQ1[pk] = sum;

  /* Compute an offset for each rank that gives the number of velocity nodes (not dof) 
     on preceding ranks */
  sum = 0;
  cnt = 0;
  { 
    /* This defines the ordering of offset */
    for (pk=0; pk<pP; ++pk) {
      for (pj=0; pj<pN; ++pj) {
        for (pi=0; pi<pM; ++pi) {
          offset[cnt] = sum;
          sum += oai[pi] * oaj[pj] * oak[pk];
          ++cnt;
        }
      }
    }
  }

  /* Compute an offset for each rank that gives the number of pressure nodes on preceding ranks */
  sum = 0;
  cnt = 0;
  {
    for (pk=0; pk<pP; ++pk) {
      for (pj=0; pj<pN; ++pj) {
        for (pi=0; pi<pM; ++pi) {
          offsetQ1[cnt] = sum;
          sum += oaiQ1[pi] * oajQ1[pj] * oakQ1[pk];
          ++cnt;
        }
      }
    }
  }

  /* Compute the global indices for the velocity nodes on this subdomain */ 
  cnt = 0;
  {
    for (k=lowerQ2.k; k<upperQ2.k; ++k) {
      for (pk=0; pk<pP; ++pk) { if ((k >= rk[pk]) && (k < rk[pk+1])) { break; } }

    for (j=lowerQ2.j; j<upperQ2.j; ++j) {

      /* Compute this point's processor index in the j direction 
         (might now be this proc because of overlap) 
        */
      for (pj=0; pj<pN; ++pj) { if ((j >= rj[pj]) && (j < rj[pj+1])) { break; } }
      
      for (i=lowerQ2.i; i<upperQ2.i; ++i) {
       /* Compute this point's processor index in the i direction 
         (might now be this proc because of overlap) to
         determine which rank pi,pj,pk contains point i,j,k 
        */
        for (pi=0; pi<pM; ++pi) { if ((i >= ri[pi]) && (i < ri[pi+1])) { break; } }
        
        /* shift by the start index of pi*/
        li = i - ri[pi];
        lj = j - rj[pj];
        lk = k - rk[pk];
        
        lidx = li + lj * oai[pi] + lk * oai[pi] * oaj[pj] ; /* velocity point index (doesn't take 2 dof into account) */
        idx = 3*lidx + 3*offset[pi + pj * pM + pk * pM * pN] + offsetQ1[pi + pj * pM + pk * pM * pN]; 
        
        indices[3*cnt+0] = idx + 0;
        indices[3*cnt+1] = idx + 1;
        indices[3*cnt+2] = idx + 2;
        
        ++cnt;
      }
    }

    }
  }
  nudofs = 3 * cnt;

  cnt = nudofs;
  {
    for (k=lowerQ1.k; k<upperQ1.k; ++k) {
      for (pk=0; pk<pP; ++pk) { if ((k >= rkQ1[pk]) && (k < rkQ1[pk+1])) { break; } }

      for (j=lowerQ1.j; j<upperQ1.j; ++j) {
        for (pj=0; pj<pN; ++pj) { if ((j >= rjQ1[pj]) && (j < rjQ1[pj+1])) { break; } }

        for (i=lowerQ1.i; i<upperQ1.i; ++i) {
          /* determine which rank pi,pj,pk contains point i,j,k */
          for (pi=0; pi<pM; ++pi) { if ((i >= riQ1[pi]) && (i < riQ1[pi+1])) { break; } }

          /* shift by the start index of pi */
          li = i - riQ1[pi];
          lj = j - rjQ1[pj];
          lk = k - rkQ1[pk];

          lidx = li + lj * oaiQ1[pi] + lk * oaiQ1[pi] * oajQ1[pj];

          /* shift by TOTAL uv-dofs p-dofs + however many uv dofs there are on this sub-domain */
          idx = lidx + 3*offset[pi + pj * pM + pk * pM * pN] + offsetQ1[pi + pj * pM + pk * pM * pN] + 3*oai[pi]*oaj[pj]*oak[pk];

          indices[cnt] = idx;

          ++cnt;
        }
      }
    }
  }
  if (cnt != length) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Counts do not match");

  ierr = ISCreateGeneral(PETSC_COMM_SELF,length,indices,PETSC_COPY_VALUES,is);CHKERRQ(ierr);
  
  PetscFree(indices);
  PetscFree(offset);
  PetscFree(rk);
  PetscFree(rj);
  PetscFree(ri);
  PetscFree(offsetQ1);
  PetscFree(rkQ1);
  PetscFree(rjQ1);
  PetscFree(riQ1);
  
  PetscFunctionReturn(0);
}
#endif

/*****************************************************************************/
static PetscErrorCode DMCreateMatrix_SaddleAIJ(DM dm,Mat *A)
{
  PetscErrorCode ierr;
  Mat            B;
  PetscInt       Ni_u,Nj_u,Nk_u,dof_u,ni_u,nj_u,nk_u,Ni_p,Nj_p,Nk_p,dof_p,ni_p,nj_p,nk_p,m,n,M,N;
  DM             dmv,dmp;
  PetscMPIInt    comm_size;
  PetscInt       *nnz,*onnz;

  PetscBool sbaijhack = PETSC_FALSE;
  
  ierr = DMCompositeGetEntries(dm,&dmv,&dmp);CHKERRQ(ierr);

  ierr = DMDAGetInfo(dmv,NULL,&Ni_u,&Nj_u,&Nk_u,NULL,NULL,NULL,&dof_u,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dmv,NULL,NULL,NULL,&ni_u,&nj_u,&nk_u);CHKERRQ(ierr);

  ierr = DMDAGetInfo(dmp,NULL,&Ni_p,&Nj_p,&Nk_p,NULL,NULL,NULL,&dof_p,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(dmp,NULL,NULL,NULL,&ni_p,&nj_p,&nk_p);CHKERRQ(ierr);
  
#if NSD == 2
  m = n = (ni_u * nj_u)*dof_u + (ni_p * nj_p)*dof_p;
  M = N = (Ni_u * Nj_u)*dof_u + (Ni_p * Nj_p)*dof_p;
#elif NSD == 3
  m = n = (ni_u * nj_u * nk_u)*dof_u + (ni_p * nj_p * nk_p)*dof_p;
  M = N = (Ni_u * Nj_u * Nk_u)*dof_u + (Ni_p * Nj_p * Nk_p)*dof_p;
#endif
  
  ierr = MatCreate(PetscObjectComm((PetscObject)dm),&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,m,n,M,N);CHKERRQ(ierr);
  ierr = MatSetType(B,MATAIJ);CHKERRQ(ierr);

  /* Some hacks to change the matrix type. If these are useful, they should be promoted to 
     documented options! */
  {
    PetscBool viennaclhack2 = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,NULL,"-viennaclhack2",&viennaclhack2,NULL);CHKERRQ(ierr);
    if (viennaclhack2){
      ierr = MatSetType(B,MATAIJVIENNACL);CHKERRQ(ierr);
    }
  }
  ierr = PetscOptionsGetBool(NULL,NULL,"-sbaijhack",&sbaijhack,NULL);CHKERRQ(ierr);
  if (sbaijhack){
    ierr = MatSetType(B,MATSBAIJ);CHKERRQ(ierr);
  }
  
  PetscMalloc1(m,&nnz);  PetscMemzero(nnz,sizeof(PetscInt)*m);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&comm_size);CHKERRQ(ierr);
  if (comm_size == 1) {
    PetscInt max_nnz_per_row = m;
    ierr = SaddlePreallocation_SEQ(dmv,dmp,nnz,max_nnz_per_row);CHKERRQ(ierr);
    if(sbaijhack){
      const PetscInt bs=1;
      ierr = MatSeqSBAIJSetPreallocation(B,bs,0,nnz);CHKERRQ(ierr); /* Note that nnz is too big since we only need upper tri */
    }else{
      ierr = MatSeqAIJSetPreallocation(B,0,nnz);CHKERRQ(ierr);
    }
  } else {
    PetscMalloc1(m,&onnz); PetscMemzero(onnz,sizeof(PetscInt)*m);
    PetscInt max_nnz_per_row_diag = m, max_nnz_per_row_offdiag = M-m;
    ierr = SaddlePreallocation_MPI(dmv,dmp,nnz,onnz,max_nnz_per_row_diag,max_nnz_per_row_offdiag);CHKERRQ(ierr);
    if (sbaijhack) {
      const PetscInt bs=1;
      ierr = MatMPISBAIJSetPreallocation(B,bs,0,nnz,0,onnz);CHKERRQ(ierr); /* Note that nnz is too big since we only need upper tri */
    } else {
      ierr = MatMPIAIJSetPreallocation(B,0,nnz,0,onnz);CHKERRQ(ierr);
    }
    PetscFree(onnz);
  }

  PetscFree(nnz);
  *A = B;
  PetscFunctionReturn(0);
}

/*****************************************************************************/
static PetscErrorCode DMDAGetProcessorsIJK(DM dm,PetscInt *pI,PetscInt *pJ,PetscInt *pK)
{
  PetscErrorCode ierr;
  PetscInt dim,m,n,p;
  PetscMPIInt size,rank;
  MPI_Comm comm;
  
  PetscFunctionBeginUser;
  
  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  ierr = DMDAGetInfo(dm,&dim,0,0,0,&m,&n,&p,0,0,0,0,0,0);CHKERRQ(ierr);
  if (dim == 1) {
    *pI = rank;
  }
  if (dim == 2) {
    *pJ = rank/m;
    *pI = rank - (*pJ)*m;
  }
  if (dim == 3) {
    PetscInt rIJ;
    
    *pK = rank/(m*n);
    rIJ = rank - (*pK)*m*n;
    *pJ = rIJ/m;
    *pI = rIJ - (*pJ)*m;
  }
  
  PetscFunctionReturn(0);
}

/*****************************************************************************/
static PetscErrorCode DMDAFEGetCornerElementQ2_3D(DM da,PetscInt *sei,PetscInt *sej,PetscInt *sek)
{
  PetscInt       i,j,k;
  PetscInt       si,sj,sk,m,n,p,width;
  PetscInt       sig,sjg,skg,mg,ng,pg;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&width,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  if (width != 2) SETERRQ(PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Stencil width must be 2 for Q2");
  ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&m,&n,&p);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&sig,&sjg,&skg,&mg,&ng,&pg);CHKERRQ(ierr);
  
  /* i-direction */
  if (sei) {
    for (i=si; i<si+m; ++i) {
      if (i%2 == 0 && i == si && i != 0) { continue; } /* reject first ghost if its's even */
      if (i%2 == 0) { *sei = i; break; }
    }
  }
  /* j-direction */
  if (sej) {
    for (j=sj; j<sj+n; ++j) {
      if (j%2 == 0 && j == sj && j != 0) { continue; } /* reject first ghost if its's even */
      if (j%2 == 0) { *sej = j; break; }
    }
  }
  /* k-direction */
  if (sek) {
    for (k=sk; k<sk+p; ++k) {
      if (k%2 == 0 && k == sk && k != 0) { continue; } /* reject first ghost if its's even */
      if (k%2 == 0) { *sek = k; break; }
    }
  }
  PetscFunctionReturn(0);
}

/*****************************************************************************/
static PetscErrorCode DMDAFEPatchCreateGlobalIS_Q2Q1(DM dm,IS *is)
{
  PetscErrorCode ierr;
  PetscInt       MQ2,NQ2,MQ1,NQ1,si=0,sj=0,mx,my;
#if NSD == 3
  PetscInt       PQ2,PQ1,sk=0,mz;
#endif
  MatStencil     lowerQ2,upperQ2,lowerQ1,upperQ1;
  PetscInt       overlap = 0;
  const char     *prefix;
  DM             dmv,dmp;
  FEMixedSpace   space;
  
  PetscFunctionBeginUser;
  ierr = DMGetOptionsPrefix(dm,&prefix);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,prefix,"-dmdafe_overlap",&overlap,NULL);CHKERRQ(ierr);
  
  ierr = DMGetApplicationContext(dm,(void**)&space);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(dm,&dmv,&dmp);CHKERRQ(ierr);
#if NSD == 2
  ierr = DMDAGetInfo(dmv,NULL,&MQ2,&NQ2,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dmp,NULL,&MQ1,&NQ1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAFEGetCornerElementQ2_3D(dmv,&si,&sj,NULL);CHKERRQ(ierr);
#elif NSD == 3
  ierr = DMDAGetInfo(dmv,NULL,&MQ2,&NQ2,&PQ2,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dmp,NULL,&MQ1,&NQ1,&PQ1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAFEGetCornerElementQ2_3D(dmv,&si,&sj,&sk);CHKERRQ(ierr);
#endif
  mx = space->lmx;
  my = space->lmy;
#if NSD == 3
  mz = space->lmz;
#endif
  
  lowerQ2.i = si - 2 * overlap; if (lowerQ2.i < 0) lowerQ2.i = 0;
  lowerQ2.j = sj - 2 * overlap; if (lowerQ2.j < 0) lowerQ2.j = 0;
#if NSD == 3
  lowerQ2.k = sk - 2 * overlap; if (lowerQ2.k < 0) lowerQ2.k = 0;
#endif
  
  upperQ2.i = si + 2 * (mx+overlap)+1; if (upperQ2.i > MQ2) upperQ2.i = MQ2;
  upperQ2.j = sj + 2 * (my+overlap)+1; if (upperQ2.j > NQ2) upperQ2.j = NQ2;
#if NSD == 3
  upperQ2.k = sk + 2 * (mz+overlap)+1; if (upperQ2.k > PQ2) upperQ2.k = PQ2;
#endif

  lowerQ1.i = si/2 - 1 * overlap; if (lowerQ1.i < 0) lowerQ1.i = 0;
  lowerQ1.j = sj/2 - 1 * overlap; if (lowerQ1.j < 0) lowerQ1.j = 0;
#if NSD == 3
  lowerQ1.k = sk/2 - 1 * overlap; if (lowerQ1.k < 0) lowerQ1.k = 0;
#endif
  
  upperQ1.i = si/2 + 1 * (mx+overlap)+1; if (upperQ1.i > MQ1) upperQ1.i = MQ1;
  upperQ1.j = sj/2 + 1 * (my+overlap)+1; if (upperQ1.j > NQ1) upperQ1.j = NQ1;
#if NSD == 3
  upperQ1.k = sk/2 + 1 * (mz+overlap)+1; if (upperQ1.k > PQ1) upperQ1.k = PQ1;
#endif
  
  ierr = PetscInfo1(dm,"Q2 sub-domain: overlap %D\n",overlap);CHKERRQ(ierr);
  ierr = PetscInfo2(dm,"Q2 sub-domain: i range [%D,%D]\n",lowerQ2.i,upperQ2.i);CHKERRQ(ierr);
  ierr = PetscInfo2(dm,"Q2 sub-domain: j range [%D,%D]\n",lowerQ2.j,upperQ2.j);CHKERRQ(ierr);
#if NSD == 3
  ierr = PetscInfo2(dm,"Q2 sub-domain: k range [%D,%D]\n",lowerQ2.k,upperQ2.k);CHKERRQ(ierr);
#endif
  ierr = PetscInfo2(dm,"Q1 sub-domain: i range [%D,%D]\n",lowerQ1.i,upperQ1.i);CHKERRQ(ierr);
  ierr = PetscInfo2(dm,"Q1 sub-domain: j range [%D,%D]\n",lowerQ1.j,upperQ1.j);CHKERRQ(ierr);
#if NSD == 3
  ierr = PetscInfo2(dm,"Q1 sub-domain: k range [%D,%D]\n",lowerQ1.k,upperQ1.k);CHKERRQ(ierr);
#endif
#if NSD == 2
  ierr = DMDAPatchCreateGlobalIS2d(dmv,dmp,lowerQ2,upperQ2,lowerQ1,upperQ1,is);CHKERRQ(ierr);
#elif NSD == 3
  ierr = DMDAPatchCreateGlobalIS3d(dmv,dmp,lowerQ2,upperQ2,lowerQ1,upperQ1,is);CHKERRQ(ierr);
#endif

  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode DMCreateDomainDecomposition_DMDAFEQ2Q1(DM dm,PetscInt *len,char ***namelist,IS **innerislist,IS **outerislist,DM **dmlist)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  *len         = 1;
  *namelist    = NULL;
  *outerislist = NULL;
  *dmlist      = NULL;
  ierr = PetscMalloc1(1,innerislist);CHKERRQ(ierr);
  ierr = DMDAFEPatchCreateGlobalIS_Q2Q1(dm,&(*innerislist)[0]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode FEMixedSpaceCreate(FEMixedSpace *space)
{
  FEMixedSpace fes;
  
  PetscFunctionBeginUser;
  PetscMalloc(sizeof(struct _p_FEMixedSpace),&fes);
  PetscMemzero(fes,sizeof(struct _p_FEMixedSpace));
  *space = fes;
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode _DMCreate_SaddleQ1_BuildElementNodeMap(DM dm,PetscInt s_el[],PetscInt e_el[],PetscInt **_u_e_map)
{
  PetscErrorCode ierr;
  PetscInt       *u_e_map;
  PetscInt       s_l[3],m_l[3];
  PetscInt       s_g[3],m_g[3];
  PetscInt       i,j,el_i,el_j,eidx;
#if NSD == 3
  PetscInt       k,el_k;
#endif

  ierr = DMDAGetGhostCorners(dm,&s_l[0],&s_l[1],&s_l[2],&m_l[0],&m_l[1],&m_l[2]);CHKERRQ(ierr); /* local space */
  ierr = DMDAGetCorners(dm,&s_g[0],&s_g[1],&s_g[2],&m_g[0],&m_g[1],&m_g[2]);CHKERRQ(ierr); /* global space */
  
  /* number of Q1 elements */
  el_i = (e_el[0] - s_el[0]);
  el_j = (e_el[1] - s_el[1]);
#if NSD == 3
  el_k = (e_el[2] - s_el[2]);
#endif

#if NSD == 2
  PetscMalloc(sizeof(PetscInt)*4*el_i*el_j,&u_e_map);

  eidx = 0;
  for (j=0; j<el_j; ++j) { 
    for (i=0; i<el_i; ++i) {
      PetscInt s0[NSD];
      
      s0[1] = j + s_el[1] - s_l[1]; /* This gives the local index */ 
      s0[0] = i + s_el[0] - s_l[0]; /* This gives the local index */

      u_e_map[4*eidx + 0] = (s0[0]+0) + (s0[1]+0)*m_l[0]; /* m_l[0] is the local number of points in the x direction (including ghosts, including unused ones)*/
      u_e_map[4*eidx + 1] = (s0[0]+1) + (s0[1]+0)*m_l[0];
      u_e_map[4*eidx + 2] = (s0[0]+0) + (s0[1]+1)*m_l[0];
      u_e_map[4*eidx + 3] = (s0[0]+1) + (s0[1]+1)*m_l[0];
#ifdef DEBUGGING_OUTPUT
      {
        PetscMPIInt rank;
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%D] Q1  p[%d-%d,] -> %d %d %d %d (local indices) \n",rank,rank,eidx,u_e_map[4*eidx + 0],u_e_map[4*eidx + 1],u_e_map[4*eidx + 2],u_e_map[4*eidx + 3]);CHKERRQ(ierr);
        ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);
        MPI_Barrier(PETSC_COMM_WORLD);
      }
#endif
      eidx++;
    }
  }
#elif NSD == 3 
  PetscMalloc(sizeof(PetscInt)*8*el_i*el_j*el_k,&u_e_map);

  eidx = 0;
  for (k=0; k<el_k; ++k) { 
    for (j=0; j<el_j; ++j) { 
      for (i=0; i<el_i; ++i) {
        PetscInt s0[NSD];

        s0[2] = k + s_el[2] - s_l[2]; /* This gives the local index */ 
        s0[1] = j + s_el[1] - s_l[1]; /* This gives the local index */ 
        s0[0] = i + s_el[0] - s_l[0]; /* This gives the local index */

        /* m_l[0] is the local number of points in the x direction (including ghosts, including unused ones)*/
        u_e_map[8*eidx + 0] = (s0[0]+0) + (s0[1]+0)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1]; 
        u_e_map[8*eidx + 1] = (s0[0]+1) + (s0[1]+0)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[8*eidx + 2] = (s0[0]+0) + (s0[1]+1)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[8*eidx + 3] = (s0[0]+1) + (s0[1]+1)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[8*eidx + 4] = (s0[0]+0) + (s0[1]+0)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1]; 
        u_e_map[8*eidx + 5] = (s0[0]+1) + (s0[1]+0)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[8*eidx + 6] = (s0[0]+0) + (s0[1]+1)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[8*eidx + 7] = (s0[0]+1) + (s0[1]+1)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];

        eidx++;
      }
    }
  }
#endif
  *_u_e_map = u_e_map;
  
  PetscFunctionReturn(0);
}
/*****************************************************************************/
PetscErrorCode _DMCreate_SaddleQ2_BuildElementNodeMap(DM dm,PetscInt s_el[],PetscInt e_el[],PetscInt **_u_e_map)
{
  PetscErrorCode ierr;
  PetscInt *u_e_map;
  PetscInt s_l[3],m_l[3];
  PetscInt s_g[3],m_g[3];
  PetscInt i,j,el_i,el_j,eidx;
#if NSD == 3
  PetscInt k,el_k;
#endif

  ierr = DMDAGetGhostCorners(dm,&s_l[0],&s_l[1],&s_l[2],&m_l[0],&m_l[1],&m_l[2]);CHKERRQ(ierr); /* local space */
  ierr = DMDAGetCorners(dm,&s_g[0],&s_g[1],&s_g[2],&m_g[0],&m_g[1],&m_g[2]);CHKERRQ(ierr); /* global space */

  /* number of Q2 elements */
  el_i = (e_el[0] - s_el[0])/2; /* note that e_el and s_el refer to global node numberings */
  el_j = (e_el[1] - s_el[1])/2;
#if NSD == 3
  el_k = (e_el[2] - s_el[2])/2;
#endif

#if NSD == 2
  PetscMalloc(sizeof(PetscInt)*9*el_i*el_j,&u_e_map);
  
  eidx = 0;
  for (j=0; j<el_j; ++j) {
    for (i=0; i<el_i; ++i) {
      PetscInt s0[NSD];
      
      s0[1] = 2*j + s_el[1] - s_l[1]; /* The local index */ 
      s0[0] = 2*i + s_el[0] - s_l[0]; /* The local index */

      u_e_map[9*eidx + 0] = (s0[0]+0) + (s0[1]+0)*m_l[0]; /* m_l[0] is the number of local x points, including ghosts (including unused ghosts) */
      u_e_map[9*eidx + 1] = (s0[0]+1) + (s0[1]+0)*m_l[0];
      u_e_map[9*eidx + 2] = (s0[0]+2) + (s0[1]+0)*m_l[0];

      u_e_map[9*eidx + 3] = (s0[0]+0) + (s0[1]+1)*m_l[0];
      u_e_map[9*eidx + 4] = (s0[0]+1) + (s0[1]+1)*m_l[0];
      u_e_map[9*eidx + 5] = (s0[0]+2) + (s0[1]+1)*m_l[0];

      u_e_map[9*eidx + 6] = (s0[0]+0) + (s0[1]+2)*m_l[0];
      u_e_map[9*eidx + 7] = (s0[0]+1) + (s0[1]+2)*m_l[0];
      u_e_map[9*eidx + 8] = (s0[0]+2) + (s0[1]+2)*m_l[0];

      eidx++;
    }
  }
#elif NSD == 3
  PetscMalloc(sizeof(PetscInt)*Q2_BASIS*el_i*el_j*el_k,&u_e_map);

  eidx = 0;
  for (k=0; k<el_k; ++k) {
    for (j=0; j<el_j; ++j) {
      for (i=0; i<el_i; ++i) {
        PetscInt s0[NSD];
 
        s0[0] = 2*i + s_el[0] - s_l[0]; /* The local index */
        s0[1] = 2*j + s_el[1] - s_l[1]; /* The local index */ 
        s0[2] = 2*k + s_el[2] - s_l[2]; /* The local index */
        /* Note the factor of 2 - this is because we are using macro elements */

        /* m_l[0] is the number of local x points, including ghosts (including unused ghosts) */
        u_e_map[27*eidx + 0 ] = (s0[0]+0) + (s0[1]+0)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 1 ] = (s0[0]+1) + (s0[1]+0)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 2 ] = (s0[0]+2) + (s0[1]+0)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 3 ] = (s0[0]+0) + (s0[1]+1)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 4 ] = (s0[0]+1) + (s0[1]+1)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 5 ] = (s0[0]+2) + (s0[1]+1)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 6 ] = (s0[0]+0) + (s0[1]+2)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 7 ] = (s0[0]+1) + (s0[1]+2)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 8 ] = (s0[0]+2) + (s0[1]+2)*m_l[0] + (s0[2]+0)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 9 ] = (s0[0]+0) + (s0[1]+0)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 10] = (s0[0]+1) + (s0[1]+0)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 11] = (s0[0]+2) + (s0[1]+0)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 12] = (s0[0]+0) + (s0[1]+1)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 13] = (s0[0]+1) + (s0[1]+1)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 14] = (s0[0]+2) + (s0[1]+1)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 15] = (s0[0]+0) + (s0[1]+2)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 16] = (s0[0]+1) + (s0[1]+2)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 17] = (s0[0]+2) + (s0[1]+2)*m_l[0] + (s0[2]+1)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 18] = (s0[0]+0) + (s0[1]+0)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 19] = (s0[0]+1) + (s0[1]+0)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 20] = (s0[0]+2) + (s0[1]+0)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 21] = (s0[0]+0) + (s0[1]+1)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 22] = (s0[0]+1) + (s0[1]+1)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 23] = (s0[0]+2) + (s0[1]+1)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];

        u_e_map[27*eidx + 24] = (s0[0]+0) + (s0[1]+2)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 25] = (s0[0]+1) + (s0[1]+2)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];
        u_e_map[27*eidx + 26] = (s0[0]+2) + (s0[1]+2)*m_l[0] + (s0[2]+2)*m_l[0]*m_l[1];

        eidx++;
      }
    }
  }
#endif
  *_u_e_map = u_e_map;
  
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode _DMCreate_SaddleQ1_BuildElementLayout(DM dm,PetscInt s_el[],PetscInt e_el[])
{
  PetscErrorCode ierr;
  PetscInt       s_l[3],m_l[3];
  PetscInt       s_g[3],m_g[3];
  PetscInt       sizes[3],i,dim;
  PetscMPIInt    rank;
  MPI_Comm       comm;
  
  PetscFunctionBeginUser;

  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  ierr = DMDAGetGhostCorners(dm,&s_l[0],&s_l[1],&s_l[2],&m_l[0],&m_l[1],&m_l[2]);CHKERRQ(ierr); /* local space */
  ierr = DMDAGetCorners(dm,&s_g[0],&s_g[1],&s_g[2],&m_g[0],&m_g[1],&m_g[2]);CHKERRQ(ierr); /* global space */
  ierr = DMDAGetInfo(dm,&dim,&sizes[0],&sizes[1],&sizes[2],NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
 
  /* For all but the last rank in each dim, the number of elements is the same as the number of global points.
     For the last rank in each dim, it is one less 
     Note that this means that are not using the ghost points on the lower end of either range - we start
     the elements at the same place as the global indices. Except for the last rank in each dim,
     we do rely on ghost points at the high end of the range.
     */
  for (i=0;i<dim;++i){
    s_el[i] = s_g[i];
    e_el[i] = s_g[i] + m_g[i] == sizes[i] ? s_g[i] + m_g[i] - 1 : s_g[i] + m_g[i];
  }

  PetscFunctionReturn(0);
}

/*****************************************************************************/
static PetscErrorCode _DMCreate_SaddleQ2_BuildElementLayout(DM dm,PetscInt s_el[],PetscInt e_el[])
{
  PetscErrorCode ierr;
  PetscInt       d,i,dim;
  PetscInt       s_l[3],m_l[3];
  PetscInt       s_g[3],m_g[3];
  PetscMPIInt    rank;
  MPI_Comm       comm;
  
  PetscFunctionBeginUser;

  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  
  ierr = DMDAGetGhostCorners(dm,&s_l[0],&s_l[1],&s_l[2],&m_l[0],&m_l[1],&m_l[2]);CHKERRQ(ierr); /* local space */
  ierr = DMDAGetCorners(dm,&s_g[0],&s_g[1],&s_g[2],&m_g[0],&m_g[1],&m_g[2]);CHKERRQ(ierr); /* global space */
  ierr = DMDAGetInfo(dm,&dim,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
 
  for (i=0; i<dim; ++i){
    s_el[i]  = -1;
    e_el[i]  = -1;
  }

  /* We start the elements either at the same place as the global numbering starts (thus ignoring both ghost points
     at the low end of the range, or 1 before (thus using 1 ghost point on either end of the range) */

  /* check global space index first, then check local space */
  for (d=0; d<dim; ++d) {
    if (s_g[d]%2 == 0) {
      s_el[d] = s_g[d];
    } else {
      if (s_g[d] - 1 >= s_l[d]) {
        s_el[d] = s_g[d] - 1;
      } else {
        SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"[start dim %d]: Cannot generate consistent macro element",d);
      }
    }

    if ((s_g[d] + m_g[d])%2 == 0) {
      e_el[d] = s_g[d] + m_g[d];
    } else {
      if ((s_g[d] + m_g[d] - 1) < (s_l[d] + m_l[d])) {
        e_el[d] = s_g[d] + m_g[d] - 1;
      } else {
        SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"[end dim %d]: Cannot generate consistent macro element",d);
      }
    }
    if ((e_el[d] - s_el[d])%2 != 0) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"[end dim %d]: Cannot generate consistent macro element: non-divisible by 2",d);
    }
#ifdef DEBUGGING_OUTPUT 
    PetscSynchronizedPrintf(comm,"Q2 [%D] dim %d [start-end] %D - %D lrange [%D - %D] : nmacro_el %D\n",rank,d,
                          s_el[d],e_el[d], s_l[d],s_l[d]+m_l[d]-1,(e_el[d] - s_el[d])/2);
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
#endif
  }

  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode DMCreate_SaddleQ2Q1(MPI_Comm comm,PetscInt mx,PetscInt my,PetscInt mz,DM *dm_saddle,FEMixedSpace fespace)
{
  PetscErrorCode ierr;
  DM             dm_stk,dm_vel,dm_p;
  PetscInt       stencil_width;
  PetscInt       s_el[3],e_el[3],s_el_p[3],e_el_p[3];
  PetscInt       np[3],*nel_i,*nel_j;
#if NSD == 3
  PetscInt       *nel_k;
#endif
  
  PetscFunctionBeginUser;
  
  /* standard Q1 mesh */
  stencil_width = 2;

  /* Its annoying this setup won't let me run 1 Q2 element per rank if comm.size > 1 */
#if NSD == 2
  ierr = DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                      2*mx+1,2*my+1,       PETSC_DECIDE,PETSC_DECIDE,             U_DOFS,stencil_width,NULL,NULL,     &dm_vel);CHKERRQ(ierr);
#elif NSD == 3
  ierr = DMDACreate3d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                      2*mx+1,2*my+1,2*mz+1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,U_DOFS,stencil_width,NULL,NULL,NULL,&dm_vel);CHKERRQ(ierr);
#endif
  ierr = DMDASetUniformCoordinates(dm_vel,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
  ierr = DMDAGetInfo(dm_vel,0,0,0,0,&np[0],&np[1],&np[2],0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = _DMCreate_SaddleQ2_BuildElementLayout(dm_vel,s_el,e_el);CHKERRQ(ierr);
  
  /* generate a list of number of Q2 elements in i,j,k on each proc - list is duplicated on all procs */
  PetscMalloc(sizeof(PetscInt)*np[0],&nel_i); PetscMemzero(nel_i,sizeof(PetscInt)*np[0]);
  PetscMalloc(sizeof(PetscInt)*np[1],&nel_j); PetscMemzero(nel_j,sizeof(PetscInt)*np[1]);
#if NSD == 3
  PetscMalloc(sizeof(PetscInt)*np[2],&nel_k); PetscMemzero(nel_k,sizeof(PetscInt)*np[2]);
#endif
  {
    PetscInt pI=0,pJ=0,pK=0;
    PetscInt el_i,el_j;
#if NSD == 3
    PetscInt el_k;
#endif
    PetscInt *recv,np_max;
    PetscInt i;
    
    np_max = np[0]; 
    for (i=1; i<NSD; ++i) {
      if (np[i] > np_max) np_max = np[i];
    }
    PetscMalloc(sizeof(PetscInt)*np_max,&recv);
    
    el_i = (e_el[0] - s_el[0])/2;
    el_j = (e_el[1] - s_el[1])/2;
#if NSD == 3
    el_k = (e_el[2] - s_el[2])/2;
#endif
    ierr = DMDAGetProcessorsIJK(dm_vel,&pI,&pJ,&pK);CHKERRQ(ierr);
    
    nel_i[pI] = el_i;
    nel_j[pJ] = el_j;
#if NSD == 3
    nel_k[pK] = el_k;
#endif
    
    PetscMemzero(recv,sizeof(PetscInt)*np_max);
    ierr = MPI_Allreduce(nel_i,recv,np[0],MPIU_INT,MPI_MAX,comm);CHKERRQ(ierr);
    PetscMemcpy(nel_i,recv,sizeof(PetscInt)*np[0]);
    
    PetscMemzero(recv,sizeof(PetscInt)*np_max);
    ierr = MPI_Allreduce(nel_j,recv,np[1],MPIU_INT,MPI_MAX,comm);CHKERRQ(ierr);
    PetscMemcpy(nel_j,recv,sizeof(PetscInt)*np[1]);
    
#if NSD == 3
    PetscMemzero(recv,sizeof(PetscInt)*np_max);
    ierr = MPI_Allreduce(nel_k,recv,np[2],MPIU_INT,MPI_MAX,comm);CHKERRQ(ierr);
    PetscMemcpy(nel_k,recv,sizeof(PetscInt)*np[2]);
#endif
    
    PetscFree(recv);
  }

  /* The number of P dofs in each direction on each node. This is almost identical to the number
     of Q2 elements per rank computed above, except that the last rank in each dimension gets
     an extra point */
  {
    PetscInt *npoints_i,*npoints_j,i,j;
#if NSD == 3
    PetscInt *npoints_k,k;    
#endif

    PetscMalloc(sizeof(PetscInt)*np[0],&npoints_i);
    PetscMalloc(sizeof(PetscInt)*np[1],&npoints_j);
#if NSD == 3
    PetscMalloc(sizeof(PetscInt)*np[2],&npoints_k);
#endif

    for(i=0;i<np[0];++i){
      npoints_i[i] = i == np[0]-1 ? nel_i[i] + 1: nel_i[i];
    }
    for(j=0;j<np[1];++j){
      npoints_j[j] = j == np[1]-1 ? nel_j[j] + 1: nel_j[j];
    }
#if NSD == 3
    for(k=0;k<np[2];++k){
      npoints_k[k] = k == np[2]-1 ? nel_k[k] + 1: nel_k[k];
    }
#endif

    stencil_width = 1;
#if NSD == 2
    ierr = DMDACreate2d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
        mx+1,my+1,np[0],np[1],P_DOFS,stencil_width,npoints_i,npoints_j,&dm_p);CHKERRQ(ierr);
#elif NSD == 3 
    ierr = DMDACreate3d(comm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
        mx+1,my+1,mz+1,np[0],np[1],np[2],P_DOFS,stencil_width,npoints_i,npoints_j,npoints_k,&dm_p);CHKERRQ(ierr);
#endif

    PetscFree(npoints_i);
    PetscFree(npoints_j);
#if NSD == 3
    PetscFree(npoints_k);
#endif
  }
  
  ierr = _DMCreate_SaddleQ1_BuildElementLayout(dm_p,s_el_p,e_el_p);CHKERRQ(ierr);

  /* setup the fespace */
  {
    PetscInt *u_e_map,*p_e_map;
    PetscInt el_i,el_j;
#if NSD == 3
    PetscInt el_k;
#endif
    
#if NSD == 2
    fespace->n_u_elements_domain = mx * my;
    fespace->n_p_elements_domain = mx * my;
#elif NSD == 3
    fespace->n_u_elements_domain = mx * my * mz;
    fespace->n_p_elements_domain = mx * my * mz;
#endif
    
    ierr = _DMCreate_SaddleQ2_BuildElementNodeMap(dm_vel,s_el,e_el,&u_e_map);CHKERRQ(ierr);
    ierr = _DMCreate_SaddleQ1_BuildElementNodeMap(dm_p,s_el_p,e_el_p,&p_e_map);CHKERRQ(ierr);

    fespace->u_el_nd_map = u_e_map;
    fespace->p_el_nd_map = p_e_map;

    el_i = (e_el[0] - s_el[0])/2;
    el_j = (e_el[1] - s_el[1])/2;
#if NSD == 3
    el_k = (e_el[2] - s_el[2])/2;
#endif
#if NSD == 2
    fespace->n_u_elements = el_i * el_j;
#elif NSD == 3
    fespace->n_u_elements = el_i * el_j  * el_k;
#endif
    fespace->lmx = el_i;
    fespace->lmy = el_j;
#if NSD == 3
    fespace->lmz = el_k;
#endif
    
    el_i = (e_el_p[0] - s_el_p[0]);
    el_j = (e_el_p[1] - s_el_p[1]);
#if NSD == 3
    el_k = (e_el_p[2] - s_el_p[2]);
#endif
#if NSD == 2
    fespace->n_p_elements = el_i * el_j;
#elif NSD == 3
    fespace->n_p_elements = el_i * el_j * el_k;
#endif
  }

	ierr = DMCompositeCreate(comm,&dm_stk);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(dm_stk,dm_vel);CHKERRQ(ierr);
	ierr = DMCompositeAddDM(dm_stk,dm_p);CHKERRQ(ierr);
  dm_stk->ops->creatematrix = DMCreateMatrix_SaddleAIJ;
	
	ierr = DMDASetFieldName(dm_vel,0,"ux");CHKERRQ(ierr);
	ierr = DMDASetFieldName(dm_vel,1,"uy");CHKERRQ(ierr);
#if NSD == 3
	ierr = DMDASetFieldName(dm_vel,2,"uz");CHKERRQ(ierr);
#endif
  ierr = DMDASetFieldName(dm_p,0,"p");CHKERRQ(ierr);
  
  ierr = DMDestroy(&dm_vel);CHKERRQ(ierr);
  ierr = DMDestroy(&dm_p);CHKERRQ(ierr);
  
  {
    Vec X;
    IS  *is_saddle_field;
    
    ierr = DMGetGlobalVector(dm_stk,&X);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(dm_stk,&X);CHKERRQ(ierr);
    ierr = DMCompositeGetGlobalISs(dm_stk,&is_saddle_field);CHKERRQ(ierr);
    
    ierr = ISDestroy(&is_saddle_field[0]);CHKERRQ(ierr);
    ierr = ISDestroy(&is_saddle_field[1]);CHKERRQ(ierr);
    ierr = PetscFree(is_saddle_field);CHKERRQ(ierr);
  }
  
  PetscFree(nel_i);
  PetscFree(nel_j);
#if NSD == 3
  PetscFree(nel_k);
#endif

  fespace->dm_saddle = dm_stk;
  *dm_saddle = dm_stk;
  
  PetscFunctionReturn(0);
}

/****************************************************************/
PetscErrorCode DMDASetUniformCoordinates_Saddle(DM dm_saddle,PetscReal x0,PetscReal x1,PetscReal y0,PetscReal y1,PetscReal z0,PetscReal z1)
{
  DM             dmv,dmp;
  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;
  ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(dmv,x0,x1,y0,y1,z0,z1);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(dmp,x0,x1,y0,y1,z0,z1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode FEMixedSpaceQuadratureCreate(FEMixedSpace space,PetscBool alloc_cell_coeff,PetscBool alloc_qp_coeff)
{
  PetscInt nqp;
  
  PetscFunctionBeginUser;
  PetscMalloc(sizeof(VQuadraturePoint),&space->vol_quadrature);
  PetscMemzero(space->vol_quadrature,sizeof(VQuadraturePoint));
  
  nqp = MAX_QUADP;
  space->vol_quadrature->n_qpoints = nqp;
  PetscMalloc(sizeof(PetscReal)*nqp*NSD,&space->vol_quadrature->qp_coor);
  PetscMalloc(sizeof(PetscReal)*nqp,  &space->vol_quadrature->qp_weight);
  {
    PetscReal xi1d[] = { -0.774596669241483, 0.0, 0.774596669241483 };
    PetscReal wt1d[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
    PetscInt i,j;
#if NSD == 3
    PetscInt k;
#endif
    PetscInt qidx;
    qidx = 0;
#if NSD == 3
    for (k=0; k<3; ++k) {
#endif
      for (j=0; j<3; ++j) {
        for (i=0; i<3; ++i) {
          space->vol_quadrature->qp_coor[NSD*qidx  ] = xi1d[i];
          space->vol_quadrature->qp_coor[NSD*qidx+1] = xi1d[j];
#if NSD == 3
          space->vol_quadrature->qp_coor[NSD*qidx+2] = xi1d[k];
#endif
#if NSD == 2
          space->vol_quadrature->qp_weight[qidx]     = wt1d[i] * wt1d[j];
#elif NSD == 3 
          space->vol_quadrature->qp_weight[qidx]     = wt1d[i] * wt1d[j] * wt1d[k];
#endif
          ++qidx;
        }
      }
#if NSD == 3
    }
#endif
  }
  
  if (alloc_cell_coeff) {
#ifdef LAME
    PetscMalloc(sizeof(LameCoefficient)*space->n_u_elements,&space->coeff_cell);
    PetscMemzero(space->coeff_cell,sizeof(LameCoefficient)*space->n_u_elements);
#else
    PetscMalloc(sizeof(StokesCoefficient)*space->n_u_elements,&space->coeff_cell);
    PetscMemzero(space->coeff_cell,sizeof(StokesCoefficient)*space->n_u_elements);
#endif
  } else {
    space->coeff_cell = NULL;
  }
  
  if (alloc_qp_coeff) {
    nqp = space->vol_quadrature->n_qpoints;
#ifdef LAME
    PetscMalloc(sizeof(LameCoefficient)*nqp*space->n_u_elements,&space->coeff_qp);
    PetscMemzero(space->coeff_qp,sizeof(LameCoefficient)*nqp*space->n_u_elements);
#else
    PetscMalloc(sizeof(StokesCoefficient)*nqp*space->n_u_elements,&space->coeff_qp);
    PetscMemzero(space->coeff_qp,sizeof(StokesCoefficient)*nqp*space->n_u_elements);
#endif
  } else {
    space->coeff_qp = NULL;
  }
  
  PetscFunctionReturn(0);
}


/*****************************************************************************/
PetscErrorCode FEMixedSpaceBCISCreate(FEMixedSpace space,DM dm_saddle)
{
  IS is_global,  is_local;
  PetscScalar    *u_bc_global,*u_bc_local;
  DM             dmv,dmp;
  PetscErrorCode ierr;
  
  ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
  
  ierr = ISCreate_BCList(dmv,PETSC_TRUE,&is_global,&u_bc_global);CHKERRQ(ierr);
  ierr = ISCreate_BCList(dmv,PETSC_FALSE,&is_local,&u_bc_local);CHKERRQ(ierr);

  space->u_is_global = is_global;
  space->u_is_local  = is_local;
  
  space->u_bc_global = u_bc_global;
  space->u_bc_local  = u_bc_local;
  
  {
    ISLocalToGlobalMapping *ltogms;
    const                  PetscInt *g_idx_u;
    PetscInt               i,k,len_is_local;
    const PetscInt         *is_local_idx;
    PetscInt               *bc_local_idx_g;
    
    ierr = DMCompositeGetISLocalToGlobalMappings(space->dm_saddle,&ltogms);CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltogms[0],&g_idx_u);CHKERRQ(ierr);
    ierr = ISGetSize(is_local,&len_is_local);CHKERRQ(ierr);
    PetscMalloc1(len_is_local,&bc_local_idx_g);
    
    ierr = ISGetIndices(is_local,&is_local_idx);CHKERRQ(ierr);
    for (k=0; k<len_is_local; ++k) {
      bc_local_idx_g[k] = g_idx_u[is_local_idx[k]];
    }
    ierr = ISRestoreIndices(is_local,&is_local_idx);CHKERRQ(ierr);
    
    ierr = ISLocalToGlobalMappingRestoreIndices(ltogms[0],&g_idx_u);CHKERRQ(ierr);
    for (i=0; i<2; ++i){
      ierr = ISLocalToGlobalMappingDestroy(&ltogms[i]);CHKERRQ(ierr);
    }
    PetscFree(ltogms);
    
    space->l_nbc = len_is_local;
    space->bc_local_idx_g = bc_local_idx_g;
  }
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode EvaluateBasis_Q1(PetscReal _xi[],PetscReal Ni[])
{
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];
#if NSD == 3
	PetscReal zeta = _xi[2];
#endif
#if NSD == 2
  Ni[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
  Ni[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
  Ni[2] = 0.25 * (1.0 - xi) * (1.0 + eta);
  Ni[3] = 0.25 * (1.0 + xi) * (1.0 + eta);
#elif NSD == 3 
	Ni[0] = 0.125 * ( 1.0 - xi ) * ( 1.0 - eta ) * ( 1.0 - zeta );
	Ni[1] = 0.125 * ( 1.0 + xi ) * ( 1.0 - eta ) * ( 1.0 - zeta );
	Ni[2] = 0.125 * ( 1.0 - xi ) * ( 1.0 + eta ) * ( 1.0 - zeta );
	Ni[3] = 0.125 * ( 1.0 + xi ) * ( 1.0 + eta ) * ( 1.0 - zeta );
	Ni[4] = 0.125 * ( 1.0 - xi ) * ( 1.0 - eta ) * ( 1.0 + zeta );
	Ni[5] = 0.125 * ( 1.0 + xi ) * ( 1.0 - eta ) * ( 1.0 + zeta );
	Ni[6] = 0.125 * ( 1.0 - xi ) * ( 1.0 + eta ) * ( 1.0 + zeta );
	Ni[7] = 0.125 * ( 1.0 + xi ) * ( 1.0 + eta ) * ( 1.0 + zeta );
#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode EvaluateBasis_Q2(PetscReal _xi[],PetscReal Ni[])
{
#if NSD == 2
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];
  
  Ni[0] = 0.5*eta*(eta-1.0)   * 0.5*xi*(xi-1.0);   /*0-0*/
  Ni[1] = 0.5*eta*(eta-1.0)   * (1.0+xi)*(1.0-xi); /*0-1*/
  Ni[2] = 0.5*eta*(eta-1.0)   * 0.5*(1.0+xi)*xi;   /*0-2*/
  
  Ni[3] = (1.0+eta)*(1.0-eta) * 0.5*xi*(xi-1.0);   /*1-0*/
  Ni[4] = (1.0+eta)*(1.0-eta) * (1.0+xi)*(1.0-xi); /*1-1*/
  Ni[5] = (1.0+eta)*(1.0-eta) * 0.5*(1.0+xi)*xi;   /*1-2*/
  
  Ni[6] = 0.5*(1.0+eta)*eta   * 0.5*xi*(xi-1.0);   /*2-0*/
  Ni[7] = 0.5*(1.0+eta)*eta   * (1.0+xi)*(1.0-xi); /*2-1*/
  Ni[8] = 0.5*(1.0+eta)*eta   * 0.5*(1.0+xi)*xi;   /*2-2*/
  PetscFunctionReturn(0);
#elif NSD == 3
	PetscInt  i,j,k,d,cnt;
	PetscReal basis_NI[3][3];
	
	for(d=0; d<3; ++d) {
		double xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); /* 0.5 * ( xi^2 - xi ) */
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); /* 1 - xi^2            */
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; /* 0.5 * ( xi^2 + xi ) */
	}
	
	cnt = 0;
	for(k=0; k<3; ++k) {
		for(j=0; j<3; ++j) {
			for(i=0; i<3; ++i) {
				Ni[cnt] = basis_NI[0][i] * basis_NI[1][j] * basis_NI[2][k];
				++cnt;
			}
		}
	}
#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode EvaluateBasisTransformation(PetscInt nbasis,
                                          PetscReal GNxi[],PetscReal GNeta[],PetscReal GNzeta[],
                                          PetscReal coords[],PetscReal *det_J)
{
  PetscInt  k;
#if NSD == 2
  PetscReal J[][2] = { {0.0,0.0} , {0.0,0.0} };
#elif NSD == 3 
  PetscReal J[][3] = { {0.0,0.0,0.0} , {0.0,0.0,0.0} , {0.0,0.0,0.0} };
#endif

  PetscFunctionBeginUser;
#if NSD == 2
  for (k=0; k<nbasis; ++k) {
    PetscInt  idx = NSD*k;
    PetscReal xc = coords[  idx];
    PetscReal yc = coords[++idx];
    
    /* J_ij = d(x^j) / d(xi^i) ~ \sum_I GN^i[I] * x^j[I] */
    J[0][0] += GNxi[k] * xc;
    J[0][1] += GNxi[k] * yc;
    
    J[1][0] += GNeta[k] * xc;
    J[1][1] += GNeta[k] * yc;
  }
  *det_J = J[0][0]*J[1][1] - J[0][1]*J[1][0];
#elif NSD == 3 
  for (k=0; k<nbasis; ++k) {
    PetscInt  idx = NSD*k;
    PetscReal xc = coords[idx  ];
    PetscReal yc = coords[idx+1];
    PetscReal zc = coords[idx+2];
    
    /* J_ij = d(x^j) / d(xi^i) ~ \sum_I GN^i[I] * x^j[I] */
    J[0][0] += GNxi[k] * xc;
    J[0][1] += GNxi[k] * yc;
    J[0][2] += GNxi[k] * zc;
    
    J[1][0] += GNeta[k] * xc;
    J[1][1] += GNeta[k] * yc;
    J[1][2] += GNeta[k] * zc;

    J[2][0] += GNzeta[k] * xc;
    J[2][1] += GNzeta[k] * yc;
    J[2][2] += GNzeta[k] * zc;
  }

  *det_J = 
		  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
		- J[0][1]*(J[1][0]*J[2][2] + J[1][2]*J[2][0]) 
		+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode EvaluateBasisDerivGlobal(PetscInt nbasis,
                                       PetscReal GNxi[],PetscReal GNeta[],PetscReal GNzeta[],
                                       PetscReal GNx[],PetscReal GNy[],PetscReal GNz[],
                                       PetscReal coords[],PetscReal *det_J)
{
  PetscInt  k;
#if NSD == 2
	PetscReal detJ;
  PetscReal iJ[2][2];
  PetscReal J[][2] = { {0.0,0.0} , {0.0,0.0} };
#elif NSD == 3 
  PetscReal t4, t6, t8, t10, t12, t14, t17;
  PetscReal J[3][3],iJ[3][3];
#endif
  
  PetscFunctionBeginUser;
#if NSD == 2
  for (k=0; k<nbasis; ++k) {
    PetscInt  idx = NSD*k;
    PetscReal xc = coords[  idx];
    PetscReal yc = coords[++idx];
    
    /* J_ij = d(x^j) / d(xi^i) ~ \sum_I GN^i[I] * x^j[I] */
    J[0][0] += GNxi[k] * xc;
    J[0][1] += GNxi[k] * yc;
    
    J[1][0] += GNeta[k] * xc;
    J[1][1] += GNeta[k] * yc;
  }
  detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

  iJ[0][0] =  J[1][1]/detJ;
  iJ[0][1] = -J[0][1]/detJ;
  iJ[1][0] = -J[1][0]/detJ;
  iJ[1][1] =  J[0][0]/detJ;
  
  *det_J = detJ;
  for (k=0; k<nbasis; ++k) {
    GNx[k] = GNxi[k]*iJ[0][0] + GNeta[k]*iJ[0][1];
    GNy[k] = GNxi[k]*iJ[1][0] + GNeta[k]*iJ[1][1];
  }
#elif NSD == 3 

  if (!GNxi) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNxi cannot be NULL");
  if (!GNeta) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNeta cannot be NULL");
  if (!GNzeta) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNzeta cannot be NULL");

  if (!GNx) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNx cannot be NULL");
  if (!GNy) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNy cannot be NULL");
  if (!GNz) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNz cannot be NULL");

  J[0][0] = J[0][1] = J[0][2] = 0.0;
  J[1][0] = J[1][1] = J[1][2] = 0.0;
  J[2][0] = J[2][1] = J[2][2] = 0.0;

  for (k=0; k<nbasis; ++k) {
    PetscReal xc = coords[NSD*k+0];
    PetscReal yc = coords[NSD*k+1];
    PetscReal zc = coords[NSD*k+2];

    J[0][0] += GNxi[k] * xc;
    J[0][1] += GNxi[k] * yc;
    J[0][2] += GNxi[k] * zc;

    J[1][0] += GNeta[k] * xc;
    J[1][1] += GNeta[k] * yc;
    J[1][2] += GNeta[k] * zc;

    J[2][0] += GNzeta[k] * xc;
    J[2][1] += GNzeta[k] * yc;
    J[2][2] += GNzeta[k] * zc;
  }
  /* flops = [NPE] * 18 */

  t4  = J[2][0] * J[0][1];
  t6  = J[2][0] * J[0][2];
  t8  = J[1][0] * J[0][1];
  t10 = J[1][0] * J[0][2];
  t12 = J[0][0] * J[1][1];
  t14 = J[0][0] * J[1][2]; /* 6 */
  t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  /* 12 */

  iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  /* 4 */
  iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; /* 5 */
  iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  /* 4 */
  iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;/* 6 */
  iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                /* 4 */
  iJ[1][2] = -(-t10 + t14) * t17;                            /* 4 */
  iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; /* 5 */
  iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               /* 5 */
  iJ[2][2] = (-t8 + t12) * t17;                              /* 3 */
  /* flops = [NQP] * 58 */
  /* shape function derivatives */
  for (k=0; k<nbasis; ++k) {
    GNx[k] = iJ[0][0]*GNxi[k] + iJ[0][1]*GNeta[k] + iJ[0][2]*GNzeta[k];

    GNy[k] = iJ[1][0]*GNxi[k] + iJ[1][1]*GNeta[k] + iJ[1][2]*GNzeta[k];

    GNz[k] = iJ[2][0]*GNxi[k] + iJ[2][1]*GNeta[k] + iJ[2][2]*GNzeta[k];
  }

  *det_J =
		  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
		- J[0][1]*(J[1][0]*J[2][2] + J[1][2]*J[2][0]) 
		+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode EvaluateBasisDerivLocal_Q1(PetscReal _xi[],PetscReal GNxi[],PetscReal GNeta[],PetscReal GNzeta[])
{
#if NSD == 2
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];
  
  /* xi deriv */
  GNxi[0] = -0.25 * (1.0 - eta);
  GNxi[1] =  0.25 * (1.0 - eta);
  GNxi[2] = -0.25 * (1.0 + eta);
  GNxi[3] =  0.25 * (1.0 + eta);
  
  /* eta deriv */
  GNeta[0] = -0.25 * (1.0 - xi);
  GNeta[1] = -0.25 * (1.0 + xi);
  GNeta[2] =  0.25 * (1.0 - xi);
  GNeta[3] =  0.25 * (1.0 + xi);
#elif NSD == 3 
	PetscReal xi   = _xi[0];
	PetscReal eta  = _xi[1];
	PetscReal zeta = _xi[2];
 
  PetscFunctionBeginUser;
  /* xi deriv */
	GNxi[0] = - 0.125 * ( 1.0 - eta ) * ( 1.0 - zeta );
	GNxi[1] =   0.125 * ( 1.0 - eta ) * ( 1.0 - zeta );
	GNxi[2] = - 0.125 * ( 1.0 + eta ) * ( 1.0 - zeta );
	GNxi[3] =   0.125 * ( 1.0 + eta ) * ( 1.0 - zeta );
	
	GNxi[4] = - 0.125 * ( 1.0 - eta ) * ( 1.0 + zeta );
	GNxi[5] =   0.125 * ( 1.0 - eta ) * ( 1.0 + zeta );
	GNxi[6] = - 0.125 * ( 1.0 + eta ) * ( 1.0 + zeta );
	GNxi[7] =   0.125 * ( 1.0 + eta ) * ( 1.0 + zeta );

  /* eta deriv */
	GNeta[0] = - 0.125 * ( 1.0 - xi ) * ( 1.0 - zeta );
	GNeta[1] = - 0.125 * ( 1.0 + xi ) * ( 1.0 - zeta );
	GNeta[2] =   0.125 * ( 1.0 - xi ) * ( 1.0 - zeta );
	GNeta[3] =   0.125 * ( 1.0 + xi ) * ( 1.0 - zeta );
	
	GNeta[4] = - 0.125 * ( 1.0 - xi ) * ( 1.0 + zeta );
	GNeta[5] = - 0.125 * ( 1.0 + xi ) * ( 1.0 + zeta );
	GNeta[6] =   0.125 * ( 1.0 - xi ) * ( 1.0 + zeta );
	GNeta[7] =   0.125 * ( 1.0 + xi ) * ( 1.0 + zeta );

  /* zeta deriv */
	GNzeta[0] = -0.125 * ( 1.0 - xi ) * ( 1.0 - eta );
	GNzeta[1] = -0.125 * ( 1.0 + xi ) * ( 1.0 - eta );
	GNzeta[2] = -0.125 * ( 1.0 - xi ) * ( 1.0 + eta );
	GNzeta[3] = -0.125 * ( 1.0 + xi ) * ( 1.0 + eta );
	
	GNzeta[4] = 0.125 * ( 1.0 - xi ) * ( 1.0 - eta );
	GNzeta[5] = 0.125 * ( 1.0 + xi ) * ( 1.0 - eta );
	GNzeta[6] = 0.125 * ( 1.0 - xi ) * ( 1.0 + eta );
	GNzeta[7] = 0.125 * ( 1.0 + xi ) * ( 1.0 + eta );
#endif
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode EvaluateBasisDerivLocal_Q2(PetscReal _xi[],PetscReal GNxi[],PetscReal GNeta[],PetscReal GNzeta[])
{
#if NSD == 2
  PetscReal   xi   = _xi[0];
  PetscReal   eta  = _xi[1];
#elif NSD == 3
	PetscScalar basis_NI[3][3];
	PetscScalar basis_GNI[3][3];
	PetscInt i,j,k,d,cnt;
#endif

  PetscFunctionBeginUser;
#if NSD == 2
  
  /* xi deriv */
  GNxi[0] = 0.5*eta*(eta-1.0)   * (xi-0.5);               /*0-0*/
  GNxi[1] = 0.5*eta*(eta-1.0)   * ( - 2.0 * xi );         /*0-1*/
  GNxi[2] = 0.5*eta*(eta-1.0)   * 0.5*( 1.0 + 2.0 * xi ); /*0-2*/
  
  GNxi[3] = (1.0+eta)*(1.0-eta) * (xi-0.5);               /*1-0*/
  GNxi[4] = (1.0+eta)*(1.0-eta) * ( - 2.0 * xi );         /*1-1*/
  GNxi[5] = (1.0+eta)*(1.0-eta) * 0.5*( 1.0 + 2.0 * xi ); /*1-2*/
  
  GNxi[6] = 0.5*(1.0+eta)*eta   * (xi-0.5);               /*2-0*/
  GNxi[7] = 0.5*(1.0+eta)*eta   * ( - 2.0 * xi );         /*2-1*/
  GNxi[8] = 0.5*(1.0+eta)*eta   * 0.5*( 1.0 + 2.0 * xi ); /*2-2*/
  
  /* eta deriv */
  GNeta[0] = (eta - 0.5) * 0.5*xi*(xi-1.0);   /*0-0*/
  GNeta[1] = (eta - 0.5) * (1.0+xi)*(1.0-xi); /*0-1*/
  GNeta[2] = (eta - 0.5) * 0.5*(1.0+xi)*xi;   /*0-2*/
  
  GNeta[3] = (-2.0*eta) * 0.5*xi*(xi-1.0);   /*1-0*/
  GNeta[4] = (-2.0*eta) * (1.0+xi)*(1.0-xi); /*1-1*/
  GNeta[5] = (-2.0*eta) * 0.5*(1.0+xi)*xi;   /*1-2*/
  
  GNeta[6] = 0.5*(1.0 + 2.0*eta) * 0.5*xi*(xi-1.0); /*2-0*/
  GNeta[7] = 0.5*(1.0 + 2.0*eta) * (1.0+xi)*(1.0-xi); /*2-1*/
  GNeta[8] = 0.5*(1.0 + 2.0*eta) * 0.5*(1.0+xi)*xi; /*2-2*/
#elif NSD == 3 

  if (!GNxi) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNxi cannot be NULL");
  if (!GNeta) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNeta cannot be NULL");
  if (!GNzeta) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"GNzeta cannot be NULL");
	for(d=0; d<3; ++d) {
		PetscScalar xi = _xi[d];
		
		basis_NI[d][0] = 0.5 * xi * (xi-1.0); /* 0.5 * ( xi^2 - xi ) */
		basis_NI[d][1] = (1.0+xi) * (1.0-xi); /* 1 - xi^2            */
		basis_NI[d][2] = 0.5 * (1.0+xi) * xi; /* 0.5 * ( xi^2 + xi ) */
		
		basis_GNI[d][0] = 0.5 * ( 2.0*xi - 1.0 );
		basis_GNI[d][1] = - 2.0*xi;
		basis_GNI[d][2] = 0.5 * ( 2.0*xi + 1.0 );
	}
	
	cnt = 0;
	for(k=0; k<3; ++k) {
		for(j=0; j<3; ++j) {
			for(i=0; i<3; ++i) {
				GNxi[cnt]   = basis_GNI[0][i]  *  basis_NI[1][j]  *  basis_NI[2][k];
				GNeta[cnt]  = basis_NI[0][i]   *  basis_GNI[1][j] *  basis_NI[2][k];
				GNzeta[cnt] = basis_NI[0][i]   *  basis_NI[1][j]  *  basis_GNI[2][k];
				++cnt;
			}
		}
	}
#endif
  PetscFunctionReturn(0);
}
/*****************************************************************************/
PetscErrorCode FEMixedSpaceDefineQPwiseProperties(FEMixedSpace space,DM dmv)
{
  PetscErrorCode ierr;
  Vec            coord_l;
  PetscScalar    *LA_coord_l;
  PetscReal      el_coor[NSD*U_BASIS];
  PetscInt       e,i,d;
  PetscScalar    Fu[U_DOFS],Fp;
#ifdef LAME
  PetscScalar    mu_qp,lambda_qp;
#else
  PetscScalar    eta_qp;
#endif
  PetscInt       q,nqp;
  PetscReal      *xi_qp,Ni[MAX_QUADP][U_BASIS];
  
  PetscFunctionBeginUser;
  nqp   = space->vol_quadrature->n_qpoints;
  xi_qp = space->vol_quadrature->qp_coor;
  
  for (q=0; q<nqp; ++q) {
    ierr = EvaluateBasis_Q2(&xi_qp[NSD*q],Ni[q]);CHKERRQ(ierr);
  }

  ierr = DMGetCoordinatesLocal(dmv,&coord_l);CHKERRQ(ierr);
#ifdef  DEBUGGING_OUTPUT 
  {
    PetscInt VecSizeLocal;
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = VecGetLocalSize(coord_l,&VecSizeLocal);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] coordinate vector size %d\n",rank,VecSizeLocal);CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);
  }
#endif
  ierr = VecGetArray(coord_l,&LA_coord_l);CHKERRQ(ierr);

  for (e=0; e<space->n_u_elements; ++e) {
    PetscInt *eidx = &space->u_el_nd_map[U_BASIS*e];
    
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = eidx[i];
      for (d=0; d<NSD; ++d) el_coor[NSD*i+d] = LA_coord_l[NSD*lnidx+d];
    }
    
    for (q=0; q<nqp; ++q) {
      PetscReal xqp[] = { 0.0, 0.0, 0.0 };
      
      for (i=0; i<U_BASIS; ++i) {
        for (d=0; d<NSD; ++d) xqp[d] += Ni[q][i] * el_coor[NSD*i+d];
#ifdef DEBUGGING_OUTPUT
  {
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] (%d) element %d (coords %f %f %f), qp %d, xqp=[%f,%f,%f]\n",rank,i,e,el_coor[NSD*i+0],el_coor[NSD*i+1],-999.9,q,xqp[0],xqp[1],xqp[2]);CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);CHKERRQ(ierr);
  }
#endif
      }

#ifdef LAME 
      ierr = Lame_EvaluateCoefficients  (xqp,&mu_qp,&lambda_qp,Fu,&Fp);CHKERRQ(ierr);
      space->coeff_qp[nqp*e+q].mu     = mu_qp;
      space->coeff_qp[nqp*e+q].lambda = lambda_qp;
#else
      ierr = Stokes_EvaluateCoefficients(xqp,&eta_qp,          Fu,&Fp);CHKERRQ(ierr);
      space->coeff_qp[nqp*e+q].eta   = eta_qp;
#endif
      for (d=0; d<U_DOFS; ++d) space->coeff_qp[nqp*e+q].Fu[d] = Fu[d];
      space->coeff_qp[nqp*e+q].Fp    = Fp;
    }
  }
  
  ierr = VecRestoreArray(coord_l,&LA_coord_l);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*****************************************************************************/
/* Note: in the old versions, there is also an unused cell-averaged procedure, if needed */
PetscErrorCode FEMixedSpaceDefineQPwiseProperties_Q1Projection(PetscInt nlevels,FEMixedSpace _space[],DM _dmscalar[])
{
  PetscErrorCode ierr;
  PetscInt       e,i,j,k,d;
#ifdef LAME
  PetscScalar    mu_qp,lambda_qp;
#else
  PetscScalar    eta_qp;
#endif
  PetscScalar    Fu_qp[NSD],Fp_qp;
  PetscInt       q,nqp;
  Vec            coeffs[EXSADDLE_NCOEFF],coeffs_l[EXSADDLE_NCOEFF],scale,coeffsLv[MG_DEPTH][EXSADDLE_NCOEFF];
  PetscScalar    el_coeffs[EXSADDLE_NCOEFF][Q1_BASIS];
  PetscScalar    el_scale[Q1_BASIS];
  DM             dmscalar;
  PetscReal      *xi_qp,Ni[MAX_QUADP][Q1_BASIS];
  FEMixedSpace   space;
  PetscScalar    *LA_values[EXSADDLE_NCOEFF];
  PetscBool      view_coeffs=PETSC_FALSE;
  
  PetscFunctionBeginUser;
  
  ierr = PetscOptionsGetBool(NULL,NULL,"-view_coeffs",&view_coeffs,NULL);CHKERRQ(ierr);

  dmscalar = _dmscalar[nlevels-1];
  space    = _space[nlevels-1];

  nqp   = space->vol_quadrature->n_qpoints;
  xi_qp = space->vol_quadrature->qp_coor;
  for (q=0; q<nqp; ++q) {
    ierr = EvaluateBasis_Q1(&xi_qp[NSD*q],Ni[q]);CHKERRQ(ierr);
  }

  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = DMCreateGlobalVector(dmscalar,&coeffs[i]);CHKERRQ(ierr);
  }
  ierr = DMCreateGlobalVector(dmscalar,&scale);CHKERRQ(ierr);
  
  /* fine level */
  for (e=0; e<space->n_u_elements; ++e) {
    PetscInt *eidx = &space->p_el_nd_map[P_BASIS*e];
    
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      PetscMemzero(el_coeffs[i],sizeof(PetscScalar)*Q1_BASIS);
    }
    PetscMemzero(el_scale,sizeof(PetscScalar)*Q1_BASIS);

    for (q=0; q<nqp; ++q) {
      for (i=0; i<Q1_BASIS; ++i) {
#ifdef LAME
        el_coeffs[0][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].mu;
        el_coeffs[4][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].lambda;
#if NSD == 3
        el_coeffs[5][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].Fu[2];
#endif 
#else
        el_coeffs[0][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].eta;
#if NSD == 3
        el_coeffs[4][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].Fu[2];
#endif 
#endif
        el_coeffs[1][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].Fu[0];
        el_coeffs[2][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].Fu[1];
        el_coeffs[3][i] += Ni[q][i] * space->coeff_qp[nqp*e+q].Fp;
        el_scale[i]     += Ni[q][i];
      }
    }
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = VecSetValuesLocal(coeffs[i],Q1_BASIS,eidx,el_coeffs[i],ADD_VALUES);CHKERRQ(ierr);
    }
    ierr = VecSetValuesLocal(scale,Q1_BASIS,eidx,el_scale,ADD_VALUES);CHKERRQ(ierr);
  }
  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = VecAssemblyBegin(coeffs[i]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(coeffs[i]);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(scale);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(scale);CHKERRQ(ierr);
  
  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = VecPointwiseDivide(coeffs[i],coeffs[i],scale);CHKERRQ(ierr);
  }

  /* Destroy scale vector that we are finished with */
  ierr = VecDestroy(&scale);CHKERRQ(ierr);

  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = DMCreateLocalVector(dmscalar,&coeffs_l[i]);CHKERRQ(ierr);
  }
  
  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = DMGlobalToLocalBegin(dmscalar,coeffs[i],INSERT_VALUES,coeffs_l[i]);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dmscalar,coeffs[i],INSERT_VALUES,coeffs_l[i]);CHKERRQ(ierr);
  }
  
  /* interpolate to quad points */
  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = VecGetArray(coeffs_l[i],&LA_values[i]);CHKERRQ(ierr);
  }
  for (e=0; e<space->n_u_elements; ++e) { /* Iterate over local elements */
    PetscInt *eidx = &space->p_el_nd_map[P_BASIS*e]; 
  
    /* get dofs */
    for (i=0; i<Q1_BASIS; ++i) {
      PetscInt lnidx = eidx[i];
      
      for (j=0; j<EXSADDLE_NCOEFF; ++j) {
          el_coeffs[j][i] = LA_values[j][lnidx];
      }
    }
    
    for (q=0; q<nqp; ++q) {
#ifdef LAME
      mu_qp = 0.0;
      lambda_qp = 0.0;
#else
      eta_qp = 0.0;
#endif
      for (d=0; d<U_DOFS; ++d) Fu_qp[d] = 0.0;
      Fp_qp = 0.0;
      for (i=0; i<Q1_BASIS; ++i) {
#ifdef LAME
        mu_qp      += Ni[q][i] * el_coeffs[0][i];
        lambda_qp  += Ni[q][i] * el_coeffs[4][i];
#if NSD == 3
        Fu_qp[2]   += Ni[q][i] * el_coeffs[5][i];
#endif 
#else
        eta_qp   += Ni[q][i] * el_coeffs[0][i];
#if NSD == 3
        Fu_qp[2]   += Ni[q][i] * el_coeffs[4][i];
#endif 
#endif
        Fu_qp[0] += Ni[q][i] * el_coeffs[1][i];
        Fu_qp[1] += Ni[q][i] * el_coeffs[2][i];
        Fp_qp    += Ni[q][i] * el_coeffs[3][i];
      }
#ifdef LAME
      space->coeff_qp[nqp*e+q].mu     = mu_qp;
      space->coeff_qp[nqp*e+q].lambda = lambda_qp;
#else
      space->coeff_qp[nqp*e+q].eta    = eta_qp;
#endif
      for (d=0; d<U_DOFS; ++d) space->coeff_qp[nqp*e+q].Fu[d] = Fu_qp[d];
      space->coeff_qp[nqp*e+q].Fp = Fp_qp;
    }
  }

  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = VecRestoreArray(coeffs_l[i],&LA_values[i]);CHKERRQ(ierr);
  }

  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    ierr = VecDestroy(&coeffs_l[i]);CHKERRQ(ierr);
  }
  if (view_coeffs) {
    /* View (global) vectors */
    char name[PETSC_MAX_PATH_LEN];
    PetscViewer viewer;

    PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"coeff_Lv%D.vts",nlevels-1);
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
#ifdef LAME
    ierr = PetscObjectSetName((PetscObject)coeffs[0],"mu");CHKERRQ(ierr);
    ierr = VecView(coeffs[0],viewer);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)coeffs[4],"lambda");CHKERRQ(ierr);
    ierr = VecView(coeffs[4],viewer);CHKERRQ(ierr);
#if NSD == 3
    ierr = PetscObjectSetName((PetscObject)coeffs[5],"Fu2");CHKERRQ(ierr);
    ierr = VecView(coeffs[5],viewer);CHKERRQ(ierr);
#endif
#else
    ierr = PetscObjectSetName((PetscObject)coeffs[0],"eta");CHKERRQ(ierr);
    ierr = VecView(coeffs[0],viewer);CHKERRQ(ierr);
#if NSD == 3
    ierr = PetscObjectSetName((PetscObject)coeffs[4],"Fu2");CHKERRQ(ierr);
    ierr = VecView(coeffs[4],viewer);CHKERRQ(ierr);
#endif
#endif
    ierr = PetscObjectSetName((PetscObject)coeffs[1],"Fu0");CHKERRQ(ierr);
    ierr = VecView(coeffs[1],viewer);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)coeffs[2],"Fu1");CHKERRQ(ierr);
    ierr = VecView(coeffs[2],viewer);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)coeffs[3],"Fp");CHKERRQ(ierr);
    ierr = VecView(coeffs[3],viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Set the global coeffs vectors for the finest level, to be interpolated */
  for (i=0; i<EXSADDLE_NCOEFF; ++i) {
    coeffsLv[nlevels-1][i] = coeffs[i];
  }
  for (k=0; k<nlevels-1; ++k) {
   
    /* Note that we could create and destroy these as we go, as
       we only use two at a time, but we haven't reached that level of 
       optimization yet */
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = DMCreateGlobalVector(_dmscalar[k],&coeffsLv[k][i]);CHKERRQ(ierr);
    }
  }
  
  for (k=nlevels-2; k>=0; k--) {
    Mat interpolate;
    Vec interp_scale;
    Vec coeffs_c[EXSADDLE_NCOEFF],coeffs_c_l[EXSADDLE_NCOEFF];

    ierr = DMCreateInterpolation(_dmscalar[k],_dmscalar[k+1],&interpolate,NULL);CHKERRQ(ierr);
    ierr = DMCreateInterpolationScale(_dmscalar[k],_dmscalar[k+1],interpolate,&interp_scale);CHKERRQ(ierr);

    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = MatRestrict(interpolate,coeffsLv[k+1][i],coeffsLv[k][i]);CHKERRQ(ierr);
      ierr = VecPointwiseMult(coeffsLv[k][i],coeffsLv[k][i],interp_scale);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&interpolate);CHKERRQ(ierr);
    ierr = VecDestroy(&interp_scale);CHKERRQ(ierr);
    
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      coeffs_c[i] = coeffsLv[k][i];
    }

    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = DMCreateLocalVector(_dmscalar[k],&coeffs_c_l[i]);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(_dmscalar[k],coeffs_c[i],INSERT_VALUES,coeffs_c_l[i]);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(_dmscalar[k],coeffs_c[i],INSERT_VALUES,coeffs_c_l[i]);CHKERRQ(ierr);
    }
    
    /* interpolate to quad points */
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = VecGetArray(coeffs_c_l[i],&LA_values[i]);CHKERRQ(ierr);
    }
    for (e=0; e<_space[k]->n_u_elements; ++e) {
      PetscInt *eidx = &_space[k]->p_el_nd_map[P_BASIS*e];
      
      /* get dofs */
      for (i=0; i<Q1_BASIS; ++i) {
        PetscInt lnidx = eidx[i];
        
        for (j=0; j<EXSADDLE_NCOEFF; ++j) {
            el_coeffs[j][i] = LA_values[j][lnidx];
        }
      }
      
      for (q=0; q<nqp; ++q) {
#ifdef LAME
        mu_qp = 0.0;
        lambda_qp = 0.0;
#else
        eta_qp = 0.0;
#endif
        for (d=0; d<U_DOFS; ++d) Fu_qp[d] = 0.0;
        Fp_qp = 0.0;
        for (i=0; i<Q1_BASIS; ++i) {
#ifdef LAME
          mu_qp     += Ni[q][i] * el_coeffs[0][i];
          lambda_qp += Ni[q][i] * el_coeffs[4][i];
#if NSD == 3
          Fu_qp[2] += Ni[q][i] * el_coeffs[5][i];
#endif
#else
          eta_qp   += Ni[q][i] * el_coeffs[0][i];
#if NSD == 3
          Fu_qp[2] += Ni[q][i] * el_coeffs[4][i];
#endif
#endif
          Fu_qp[0] += Ni[q][i] * el_coeffs[1][i];
          Fu_qp[1] += Ni[q][i] * el_coeffs[2][i];
          Fp_qp    += Ni[q][i] * el_coeffs[3][i];
        }
#ifdef LAME
        _space[k]->coeff_qp[nqp*e+q].mu     = mu_qp;
        _space[k]->coeff_qp[nqp*e+q].lambda = lambda_qp;
#else
        _space[k]->coeff_qp[nqp*e+q].eta    = eta_qp;
#endif
        for (d=0; d<U_DOFS; ++d) _space[k]->coeff_qp[nqp*e+q].Fu[d] = Fu_qp[d];
        _space[k]->coeff_qp[nqp*e+q].Fp = Fp_qp;
      }
    }
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = VecRestoreArray(coeffs_c_l[i],&LA_values[i]);CHKERRQ(ierr);
    }

    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = VecDestroy(&coeffs_c_l[i]);CHKERRQ(ierr);
    }
   
    if (view_coeffs) {
      char name[PETSC_MAX_PATH_LEN];
      PetscViewer viewer;

      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"coeff_Lv%D.vts",k);
      ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
#ifdef LAME
      ierr = PetscObjectSetName((PetscObject)coeffs_c[0],"mu");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[0],viewer);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)coeffs_c[4],"lambda");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[4],viewer);CHKERRQ(ierr);
#if NSD == 3
      ierr = PetscObjectSetName((PetscObject)coeffs_c[5],"Fu2");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[5],viewer);CHKERRQ(ierr);
#endif
#else
      ierr = PetscObjectSetName((PetscObject)coeffs_c[0],"eta");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[0],viewer);CHKERRQ(ierr);
#if NSD == 3
      ierr = PetscObjectSetName((PetscObject)coeffs_c[4],"Fu2");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[4],viewer);CHKERRQ(ierr);
#endif
#endif
      ierr = PetscObjectSetName((PetscObject)coeffs_c[1],"Fu0");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[1],viewer);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)coeffs_c[2],"Fu1");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[2],viewer);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)coeffs_c[3],"Fp");CHKERRQ(ierr);
      ierr = VecView(coeffs_c[3],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
  }
  
  /* Note that in theory we could create and destroy these as we go, as we only use two at a time,
     but that is more optimization than is warranted currently */
  for (k=0; k<nlevels; ++k) {
    for (i=0; i<EXSADDLE_NCOEFF; ++i) {
      ierr = VecDestroy(&coeffsLv[k][i]);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode FEMixedSpaceDestroy(FEMixedSpace *space)
{
  FEMixedSpace   fes;
  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;

  if (!space) PetscFunctionReturn(0);

  fes = *space;
  
  PetscFree(fes->vol_quadrature->qp_coor);
  PetscFree(fes->vol_quadrature->qp_weight);
  PetscFree(fes->vol_quadrature);
  
  if (fes->bc_local_idx_g) PetscFree(fes->bc_local_idx_g); 
  ierr = ISDestroy(&fes->u_is_global);CHKERRQ(ierr);
  ierr = ISDestroy(&fes->u_is_local);CHKERRQ(ierr);
  
  if (fes->u_bc_global) PetscFree(fes->u_bc_global);
  if (fes->u_bc_local) PetscFree(fes->u_bc_local);
  
  if (fes->coeff_cell) PetscFree(fes->coeff_cell);
  if (fes->coeff_qp) PetscFree(fes->coeff_qp);

  if (fes->u_el_nd_map) PetscFree(fes->u_el_nd_map);
  if (fes->p_el_nd_map) PetscFree(fes->p_el_nd_map);
  
  if (fes->e_centroid) PetscFree(fes->e_centroid);
  if (fes->e_dx) PetscFree(fes->e_dx);
  
  PetscFree(fes);
  
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode MatAssemble_Saddle_NULL(FEMixedSpace space,DM dm_saddle,Mat A)
{
  PetscErrorCode          ierr;
  DM                      dmv,dmp;
  PetscInt                e,i,d;
  PetscScalar             el_A11[U_DOFS*U_DOFS*U_BASIS*U_BASIS];
  PetscInt                el_u_idx[U_DOFS*U_BASIS];
  PetscInt                el_p_idx[P_DOFS*P_BASIS];
  ISLocalToGlobalMapping *ltogms;
  const PetscInt         *g_idx_u;
  const PetscInt         *g_idx_p;
  
  PetscFunctionBeginUser;

  PetscMemzero(el_A11,sizeof(PetscScalar)*NSD*NSD*U_BASIS*U_BASIS); /* zeros, so this one block (the biggest) is fine */

  ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
  
  /* Get mappings of local dof to global dof */
  ierr = DMCompositeGetISLocalToGlobalMappings(dm_saddle,&ltogms);CHKERRQ(ierr); 
  ierr = ISLocalToGlobalMappingGetIndices(ltogms[0],&g_idx_u);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetIndices(ltogms[1],&g_idx_p);CHKERRQ(ierr);
 
  for (e=0; e<space->n_u_elements; ++e) {
    PetscInt *e_u_idx = &space->u_el_nd_map[U_BASIS*e];
    PetscInt *e_p_idx = &space->p_el_nd_map[P_BASIS*e];
    
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = e_u_idx[i];
      for (d=0; d<U_DOFS; ++d) el_u_idx[U_DOFS*i+d] = U_DOFS*lnidx + d;
    }

    for (i=0; i<P_BASIS; ++i) {
      PetscInt lnidx = e_p_idx[i];
      
      el_p_idx[i] = lnidx;
    }

    /* Convert local to global indices */
    PetscInt  el_u_idx_g[U_DOFS*U_BASIS];
    PetscInt  el_p_idx_g[P_DOFS*P_BASIS];
    for (i=0;i<U_DOFS*U_BASIS;++i) el_u_idx_g[i] = g_idx_u[el_u_idx[i]];    
    for (i=0;i<P_DOFS*P_BASIS;++i) el_p_idx_g[i] = g_idx_p[el_p_idx[i]]; /* el_p_idx \in [0,...,npress_dof-1] */

    ierr = MatSetValues(A,U_DOFS*U_BASIS,el_u_idx_g,U_DOFS*U_BASIS,el_u_idx_g,el_A11,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(A,U_DOFS*U_BASIS,el_u_idx_g,P_DOFS*P_BASIS,el_p_idx_g,el_A11,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(A,P_DOFS*P_BASIS,el_p_idx_g,U_DOFS*U_BASIS,el_u_idx_g,el_A11,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(A,P_DOFS*P_BASIS,el_p_idx_g,P_DOFS*P_BASIS,el_p_idx_g,el_A11,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = ISLocalToGlobalMappingRestoreIndices(ltogms[0],&g_idx_u);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingRestoreIndices(ltogms[1],&g_idx_p);CHKERRQ(ierr);
  for (i=0;i<2;++i){
    ierr = ISLocalToGlobalMappingDestroy(&ltogms[i]);CHKERRQ(ierr);
  }
  PetscFree(ltogms);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode MatAssemble_Saddle(FEMixedSpace space,DM dm_saddle,Mat A,Vec rhs_diri)
{
  PetscErrorCode         ierr;
  DM                     dmv,dmp;
  PetscInt               e,i,j,k,d;
  PetscInt               el_u_idx[U_DOFS*U_BASIS];
  PetscInt               el_p_idx[P_DOFS*P_BASIS];
  PetscInt               q,nqp;
  PetscReal              *w_qp,*xi_qp;
#ifdef LAME
  PetscScalar            mu_c;
#else
  PetscScalar            eta_c;
#endif
  Vec                    coor_u_l;
  PetscScalar            *LA_coor_u;
  PetscReal              el_coor_u[U_DOFS*U_BASIS];
  PetscReal              Nu[MAX_QUADP][U_BASIS],GNuxi[MAX_QUADP][U_BASIS],GNueta[MAX_QUADP][U_BASIS],GNux[U_BASIS],GNuy[U_BASIS],Np[MAX_QUADP][P_BASIS];
#if NSD == 3
  PetscReal              GNuzeta[MAX_QUADP][U_BASIS],GNuz[U_BASIS];
#endif
#ifdef LAME
  Vec                    coor_p_l;
  PetscScalar            *LA_coor_p;
  PetscReal              el_coor_p[NSD*P_BASIS];
  PetscReal              GNpxi[MAX_QUADP][P_BASIS],GNpeta[MAX_QUADP][P_BASIS];
#if NSD == 3
  PetscReal              GNpzeta[MAX_QUADP][P_BASIS];
#endif
#endif
  PetscScalar            el_A11[U_DOFS*U_DOFS*U_BASIS*U_BASIS];
  PetscScalar            el_A12[U_DOFS*P_DOFS*U_BASIS*P_BASIS];
  PetscScalar            el_A21[P_DOFS*U_DOFS*P_BASIS*U_BASIS];
#ifdef LAME
  PetscScalar            el_A22[P_DOFS*P_DOFS*P_BASIS*P_BASIS];
#endif
  const PetscInt         *g_idx_u;
  const PetscInt         *g_idx_p;
  ISLocalToGlobalMapping *ltogms;
  
  PetscFunctionBeginUser;
  
  ierr = MatZeroEntries(A);CHKERRQ(ierr);
  ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);

  nqp   = space->vol_quadrature->n_qpoints;
  w_qp  = space->vol_quadrature->qp_weight;
  xi_qp = space->vol_quadrature->qp_coor;
 
  for (q=0; q<nqp; ++q) {
    ierr = EvaluateBasis_Q2(&xi_qp[NSD*q],Nu[q]);CHKERRQ(ierr);
#if NSD == 2
    ierr = EvaluateBasisDerivLocal_Q2(&xi_qp[NSD*q],GNuxi[q],GNueta[q],0         );CHKERRQ(ierr);
#elif NSD == 3
    ierr = EvaluateBasisDerivLocal_Q2(&xi_qp[NSD*q],GNuxi[q],GNueta[q],GNuzeta[q]);CHKERRQ(ierr);
#endif
    ierr = EvaluateBasis_Q1(&xi_qp[NSD*q],Np[q]);CHKERRQ(ierr);
#ifdef LAME
#if NSD == 2
    ierr = EvaluateBasisDerivLocal_Q1(&xi_qp[NSD*q],GNpxi[q],GNpeta[q],0         );CHKERRQ(ierr);
#elif NSD == 3
    ierr = EvaluateBasisDerivLocal_Q1(&xi_qp[NSD*q],GNpxi[q],GNpeta[q],GNpzeta[q]);CHKERRQ(ierr);
#endif
#endif
  }

  /* Get mapping from local to global indices */
  ierr = DMCompositeGetISLocalToGlobalMappings(dm_saddle,&ltogms);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetIndices(ltogms[0],&g_idx_u);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetIndices(ltogms[1],&g_idx_p);CHKERRQ(ierr);
  
  ierr = DMGetCoordinatesLocal(dmv,&coor_u_l);CHKERRQ(ierr);
  ierr = VecGetArray(coor_u_l,&LA_coor_u);CHKERRQ(ierr);

#ifdef LAME
  ierr = DMGetCoordinatesLocal(dmp,&coor_p_l);CHKERRQ(ierr);
  ierr = VecGetArray(coor_p_l,&LA_coor_p);CHKERRQ(ierr);
#endif
  
  for (e=0; e<space->n_u_elements; ++e) {
    PetscInt *e_u_idx = &space->u_el_nd_map[U_BASIS*e];
    PetscInt *e_p_idx = &space->p_el_nd_map[P_BASIS*e];
    PetscInt  el_u_idx_g[U_DOFS*U_BASIS];
    PetscInt  el_p_idx_g[P_DOFS*P_BASIS];
    
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = e_u_idx[i];
      for (d=0; d<U_DOFS; ++d) el_u_idx[U_DOFS*i+d] = U_DOFS*lnidx + d;
    }
    
    for (i=0; i<P_BASIS; ++i) {
      PetscInt lnidx = e_p_idx[i];
      el_p_idx[i] = lnidx;
    }
    
    /* get element coordinates */
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = e_u_idx[i];
      for (d=0; d<NSD; ++d) el_coor_u[NSD*i+d] = LA_coor_u[NSD*lnidx+d];
    }
#ifdef LAME
    for (i=0; i<P_BASIS; ++i) {
      PetscInt lnidx = e_p_idx[i];
      for (d=0; d<NSD; ++d) el_coor_p[NSD*i+d] = LA_coor_p[NSD*lnidx+d];
    }
#endif
    
    PetscMemzero(el_A11,sizeof(PetscScalar)*U_DOFS*U_DOFS*U_BASIS*U_BASIS);
    PetscMemzero(el_A12,sizeof(PetscScalar)*U_DOFS*P_DOFS*U_BASIS*P_BASIS);
    PetscMemzero(el_A21,sizeof(PetscScalar)*P_DOFS*U_DOFS*P_BASIS*U_BASIS);
#ifdef LAME
    PetscMemzero(el_A22,sizeof(PetscScalar)*P_DOFS*P_DOFS*P_BASIS*P_BASIS);
#endif

    /* Assemble A11 */
    /* Note that we assemble the full operator, even though it is 
       symmetric. This could be changed if we used a symmetric 
       matrix type */
    for (q=0; q<nqp; ++q) {
      PetscScalar fac;
      PetscReal   detJ;
#if NSD == 2
      PetscReal B[3][U_DOFS*U_BASIS],D[3];
#elif NSD == 3 
      PetscReal B[6][U_DOFS*U_BASIS],D[6];
#endif
      
      /* evaluate derivatives */
#if NSD == 2
      ierr = EvaluateBasisDerivGlobal(U_BASIS,GNuxi[q],GNueta[q],0,         GNux,GNuy,0,   el_coor_u,&detJ);CHKERRQ(ierr);
      for (i=0; i<U_BASIS; ++i) {
        B[0][U_DOFS*i] = GNux[i];  B[0][U_DOFS*i+1] = 0.0;
        B[1][U_DOFS*i] = 0.0;      B[1][U_DOFS*i+1] = GNuy[i];
        B[2][U_DOFS*i] = GNuy[i];  B[2][U_DOFS*i+1] = GNux[i];
      }
#elif NSD == 3 
      ierr = EvaluateBasisDerivGlobal(U_BASIS,GNuxi[q],GNueta[q],GNuzeta[q],GNux,GNuy,GNuz,el_coor_u,&detJ);CHKERRQ(ierr);
      for (i=0; i<U_BASIS; ++i) {
        B[0][U_DOFS*i] = GNux[i];  B[0][U_DOFS*i+1] = 0.0;     B[0][U_DOFS*i+2] = 0.0;
        B[1][U_DOFS*i] = 0.0;      B[1][U_DOFS*i+1] = GNuy[i]; B[1][U_DOFS*i+2] = 0.0;
        B[2][U_DOFS*i] = 0.0;      B[2][U_DOFS*i+1] = 0.0;     B[2][U_DOFS*i+2] = GNuz[i];

        B[3][U_DOFS*i] = GNuy[i];  B[3][U_DOFS*i+1] = GNux[i]; B[3][U_DOFS*i+2] = 0.0;
        B[4][U_DOFS*i] = GNuz[i];  B[4][U_DOFS*i+1] = 0.0;     B[4][U_DOFS*i+2] = GNux[i];
        B[5][U_DOFS*i] = 0.0;      B[5][U_DOFS*i+1] = GNuz[i]; B[5][U_DOFS*i+2] = GNuy[i];
      }
#endif
     
#ifdef LAME
      mu_c  = space->coeff_qp[nqp*e+q].mu;
      fac = mu_c  * w_qp[q] * detJ;
#else
      eta_c = space->coeff_qp[nqp*e+q].eta;
      fac = eta_c * w_qp[q] * detJ;
#endif


#if NSD == 2
      D[0] = 2.0 * fac;
      D[1] = 2.0 * fac;
      D[2] = 1.0 * fac;
      
      /* k_ij = Bt_ik D_kk B_kj = B_ik D_kk B_kj */
      for (i=0; i<U_DOFS*U_BASIS; ++i) {
        for (j=0; j<U_DOFS*U_BASIS; ++j) {
          for (k=0; k<3; ++k) {
            el_A11[i*(U_DOFS*U_BASIS)+j] += B[k][i] * D[k] * B[k][j];
          }
        }
      }
#elif NSD == 3 
      D[0] = 2.0 * fac;
      D[1] = 2.0 * fac;
      D[2] = 2.0 * fac;

      D[3] = 1.0 * fac;
      D[4] = 1.0 * fac;
      D[5] = 1.0 * fac;
      
      /* k_ij = Bt_ik D_kk B_kj = B_ik D_kk B_kj */
      for (i=0; i<U_DOFS*U_BASIS; ++i) {
        for (j=0; j<U_DOFS*U_BASIS; ++j) {
          for (k=0; k<6; ++k) {
            el_A11[i*(U_DOFS*U_BASIS)+j] += B[k][i] * D[k] * B[k][j];
          }
        }
      }
#endif
    }

    /* Assemble A12 */
  for (q=0; q<nqp; ++q) {
    PetscReal fac,detJ;

#if NSD == 2
    ierr = EvaluateBasisDerivGlobal(U_BASIS,GNuxi[q],GNueta[q],0,         GNux,GNuy,0,   el_coor_u,&detJ);CHKERRQ(ierr);
#elif NSD == 3 
    ierr = EvaluateBasisDerivGlobal(U_BASIS,GNuxi[q],GNueta[q],GNuzeta[q],GNux,GNuy,GNuz,el_coor_u,&detJ);CHKERRQ(ierr);
#endif

    fac = w_qp[q] * detJ;
    for (i=0; i<U_BASIS; ++i) {
      for (j=0; j<P_BASIS; ++j) {
        el_A12[(U_DOFS*i+0)*(P_BASIS)+j] -= GNux[i] * Np[q][j] * fac;
        el_A12[(U_DOFS*i+1)*(P_BASIS)+j] -= GNuy[i] * Np[q][j] * fac;
#if NSD == 3
        el_A12[(U_DOFS*i+2)*(P_BASIS)+j] -= GNuz[i] * Np[q][j] * fac;
#endif
      }
    }

    for (i=0; i<U_BASIS; ++i) {
      for (j=0; j<P_BASIS; ++j) {
        for(d=0; d<U_DOFS; ++d){
          el_A21[j*(U_DOFS*U_BASIS)+(U_DOFS*i+d)] = el_A12[(U_DOFS*i+d)*(P_BASIS)+j];
        }
      }
    }
  }
#ifdef LAME
    /* Assemble A22 */
    for (q=0; q<nqp; ++q) {
      PetscReal fac,detJ,lambda_c;
#if NSD == 2
      ierr = EvaluateBasisTransformation(P_BASIS,GNpxi[q],GNpeta[q],0,el_coor_p,&detJ);CHKERRQ(ierr);
#elif NSD == 3 
      ierr = EvaluateBasisTransformation(P_BASIS,GNpxi[q],GNpeta[q],GNpzeta[q],el_coor_p,&detJ);CHKERRQ(ierr);
#endif
      lambda_c = space->coeff_qp[nqp*e+q].lambda;
      fac = w_qp[q] * detJ / lambda_c;

      for (i=0; i<P_BASIS; ++i) {
        for (j=0; j<P_BASIS; ++j) {
          el_A22[i*(P_BASIS)+j] -= Np[q][i] * Np[q][j] * fac;
        }
      }
    }
#endif

    /* Convert to global indices and insert */
    for (i=0;i<U_DOFS*U_BASIS;++i) el_u_idx_g[i] = g_idx_u[el_u_idx[i]];
    for (i=0;i<P_DOFS*P_BASIS;++i) el_p_idx_g[i] = g_idx_p[el_p_idx[i]]; /* el_p_idx \in [0,...,npress_dof-1] */
    ierr = MatSetValues(A, U_DOFS*U_BASIS,el_u_idx_g, U_DOFS*U_BASIS,el_u_idx_g, el_A11,ADD_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(A, U_DOFS*U_BASIS,el_u_idx_g, P_DOFS*P_BASIS,el_p_idx_g, el_A12,ADD_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(A, P_DOFS*P_BASIS,el_p_idx_g, U_DOFS*U_BASIS,el_u_idx_g, el_A21,ADD_VALUES);CHKERRQ(ierr);
#ifdef LAME
    ierr = MatSetValues(A, P_DOFS*P_BASIS,el_p_idx_g, P_DOFS*P_BASIS,el_p_idx_g, el_A22,ADD_VALUES);CHKERRQ(ierr);
#endif
  }
  ierr = VecRestoreArray(coor_u_l,&LA_coor_u);CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = ISLocalToGlobalMappingRestoreIndices(ltogms[0],&g_idx_u);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingRestoreIndices(ltogms[1],&g_idx_p);CHKERRQ(ierr);
  for (i=0;i<2;++i){
    ierr = ISLocalToGlobalMappingDestroy(&ltogms[i]);CHKERRQ(ierr);
  }
  PetscFree(ltogms);

  if(rhs_diri){
    Vec temp;
    ierr = VecDuplicate(rhs_diri,&temp);CHKERRQ(ierr);
    ierr = VecZeroEntries(temp);CHKERRQ(ierr);
    ierr = ImposeDirichletValuesIS(temp,space->u_is_global,space->u_bc_global);CHKERRQ(ierr);
    ierr = MatMult(A,temp,rhs_diri);CHKERRQ(ierr); /* Note A has not had any rows or columns zeroed yet */
    ierr = VecDestroy(&temp);CHKERRQ(ierr);
    ierr = VecScale(rhs_diri,-1.0);CHKERRQ(ierr);
    ierr = ZeroDirichletValuesIS(rhs_diri,space->u_is_global);CHKERRQ(ierr);
  }

  ierr = MatZeroRowsColumns(A,space->l_nbc,space->bc_local_idx_g,1.0,NULL,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode VecAssemble_F1_qp(FEMixedSpace space,Vec F)
{
  PetscErrorCode ierr;
  DM             dmv,dmp;
  Vec            coor_u_l;
  PetscScalar    *LA_coor_u;
  PetscReal      el_coor_u[NSD*U_BASIS];
  PetscInt       e,i,d;
  PetscInt       q,nqp;
  PetscReal      *w_qp,*xi_qp;
  PetscReal      detJ,Ni[MAX_QUADP][U_BASIS],GNxi[MAX_QUADP][U_BASIS],GNeta[MAX_QUADP][U_BASIS],GNzeta[MAX_QUADP][U_BASIS];
  PetscScalar    el_f[NSD*U_BASIS];
  PetscInt       el_u_idx[NSD*U_BASIS];
  
  PetscFunctionBeginUser;

  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  nqp   = space->vol_quadrature->n_qpoints;
  w_qp  = space->vol_quadrature->qp_weight;
  xi_qp = space->vol_quadrature->qp_coor;
  
  ierr = DMCompositeGetEntries(space->dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dmv,&coor_u_l);CHKERRQ(ierr);
  ierr = VecGetArray(coor_u_l,&LA_coor_u);CHKERRQ(ierr);
  
  for (q=0; q<nqp; ++q) {
    ierr = EvaluateBasis_Q2(&xi_qp[NSD*q],Ni[q]);CHKERRQ(ierr);
    ierr = EvaluateBasisDerivLocal_Q2(&xi_qp[NSD*q],GNxi[q],GNeta[q],GNzeta[q]);CHKERRQ(ierr);
  }
  
  for (e=0; e<space->n_u_elements; ++e) {
    PetscInt *e_u_idx = &space->u_el_nd_map[U_BASIS*e];
    
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = e_u_idx[i];
      for (d=0; d<U_DOFS; ++d) el_u_idx[U_DOFS*i+d] = NSD*lnidx+d;
    }
    
    /* get element coordinates */
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = e_u_idx[i];
      for (d=0; d<NSD; ++d) el_coor_u[NSD*i+d] = LA_coor_u[NSD*lnidx+d];
    }
    
    PetscMemzero(el_f,sizeof(PetscScalar)*NSD*U_BASIS);
    for (q=0; q<nqp; ++q) {
      PetscScalar fac,Fu[U_DOFS];
      
      /* element coefficient */
      for (d=0; d<U_DOFS; ++d) Fu[d] = space->coeff_qp[nqp*e + q].Fu[d];
      
      /* evaluate derivatives */
      ierr = EvaluateBasisTransformation(U_BASIS,GNxi[q],GNeta[q],GNzeta[q],el_coor_u,&detJ);CHKERRQ(ierr);
      
      fac = w_qp[q] * detJ;
      
      for (i=0; i<U_BASIS; ++i) {
        for (d=0; d<U_DOFS; ++d) el_f[U_DOFS*i+d] += Ni[q][i] * Fu[d] * fac;
      }
    }
    ierr = VecSetValuesLocal(F,U_DOFS*U_BASIS,el_u_idx,el_f,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(coor_u_l,&LA_coor_u);CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode VecAssemble_F2_qp(FEMixedSpace space,Vec F)
{
  PetscErrorCode ierr;
  DM             dmv,dmp;
  Vec            coor_u_l;
  PetscScalar    *LA_coor_u;
  PetscReal      el_coor_u[NSD*U_BASIS];
  PetscInt       e,i,d;
  PetscInt       q,nqp;
  PetscScalar    *w_qp;
  PetscReal      *xi_qp;
  PetscReal      detJ,Ni_u[MAX_QUADP][U_BASIS],GNxi[MAX_QUADP][U_BASIS],GNeta[MAX_QUADP][U_BASIS],GNzeta[MAX_QUADP][U_BASIS],Ni_p[MAX_QUADP][P_BASIS];
  PetscScalar    el_f[P_BASIS];
  
  PetscFunctionBeginUser;

  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  nqp   = space->vol_quadrature->n_qpoints;
  w_qp  = space->vol_quadrature->qp_weight;
  xi_qp = space->vol_quadrature->qp_coor;
  
  ierr = DMCompositeGetEntries(space->dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dmv,&coor_u_l);CHKERRQ(ierr);
  ierr = VecGetArray(coor_u_l,&LA_coor_u);CHKERRQ(ierr);
  
  for (q=0; q<nqp; ++q) {
    ierr = EvaluateBasis_Q2(&xi_qp[NSD*q],Ni_u[q]);CHKERRQ(ierr);
    ierr = EvaluateBasisDerivLocal_Q2(&xi_qp[NSD*q],GNxi[q],GNeta[q],GNzeta[q]);CHKERRQ(ierr);
    ierr = EvaluateBasis_Q1(&xi_qp[NSD*q],Ni_p[q]);CHKERRQ(ierr);
  }
  
  for (e=0; e<space->n_u_elements; ++e) {
    PetscInt *e_u_idx = &space->u_el_nd_map[U_BASIS*e];
    PetscInt *e_p_idx = &space->p_el_nd_map[P_BASIS*e];
    
    /* get element coordinates */
    for (i=0; i<U_BASIS; ++i) {
      PetscInt lnidx = e_u_idx[i];
      for (d=0; d<NSD; ++d) el_coor_u[NSD*i+d] = LA_coor_u[NSD*lnidx+d];
    }
    
    PetscMemzero(el_f,sizeof(PetscReal)*P_BASIS);
    for (q=0; q<nqp; ++q) {
      PetscScalar fac,Fp;
      
      /* element coefficient */
      Fp = space->coeff_qp[nqp*e + q].Fp;
      
      /* evaluate derivatives */
      ierr = EvaluateBasisTransformation(U_BASIS,GNxi[q],GNeta[q],GNzeta[q],el_coor_u,&detJ);CHKERRQ(ierr);
      
      fac = w_qp[q] * detJ;
      
      for (i=0; i<P_BASIS; ++i) {
        el_f[i] += Ni_p[q][i] * Fp * fac;
      }
    }
    ierr = VecSetValuesLocal(F,P_BASIS,e_p_idx,el_f,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(coor_u_l,&LA_coor_u);CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode ZeroDirichletValuesIS(Vec x,IS is)
{
  PetscErrorCode ierr;
  PetscScalar    *LA_x;
  PetscInt       i,n,N;
  const PetscInt *idx;
  
  PetscFunctionBeginUser;
  ierr = VecGetLocalSize(x,&N);CHKERRQ(ierr);
  ierr = ISGetSize(is,&n);CHKERRQ(ierr);
  ierr = ISGetIndices(is,&idx);CHKERRQ(ierr);
  ierr = VecGetArray(x,&LA_x);CHKERRQ(ierr);
  for (i=0; i<n; ++i) {
    if (idx[i] >= N) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"BCIS not compatible with vector: idx[] > vec.row.n");
    }
    LA_x[idx[i]] = 0.0;
  }
  ierr = VecRestoreArray(x,&LA_x);CHKERRQ(ierr);
  ierr = ISRestoreIndices(is,&idx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode ImposeDirichletValuesIS(Vec x,IS is,PetscScalar vals[])
{
  PetscErrorCode ierr;
  PetscScalar    *LA_x;
  PetscInt       i,n,N;
  const PetscInt *idx;
  
  PetscFunctionBeginUser;
  ierr = VecGetLocalSize(x,&N);CHKERRQ(ierr);
  ierr = ISGetSize(is,&n);CHKERRQ(ierr);
  ierr = ISGetIndices(is,&idx);CHKERRQ(ierr);
  ierr = VecGetArray(x,&LA_x);CHKERRQ(ierr);
  for (i=0; i<n; ++i) {
    if (idx[i] >= N) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"BCIS not compatible with vector: idx[] > vec.row.n");
    }
    LA_x[idx[i]] = vals[i];
  }
  ierr = VecRestoreArray(x,&LA_x);CHKERRQ(ierr);
  ierr = ISRestoreIndices(is,&idx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode MatAssemble_Schur(FEMixedSpace space,DM dmp,Mat S)
{
  PetscErrorCode ierr;
  PetscInt               e,i,j,d;
  PetscInt               el_p_idx[P_DOFS*P_BASIS];
  PetscInt               q,nqp;
  PetscScalar            *w_qp;
  PetscReal              *xi_qp;
  Vec                    coor_p_l;
  PetscScalar            *LA_coor_p;
  PetscReal              el_coor_p[NSD*P_BASIS];
  PetscReal              Np[MAX_QUADP][P_BASIS];
  PetscReal              GNpxi[MAX_QUADP][P_BASIS],GNpeta[MAX_QUADP][P_BASIS];
#if NSD == 3
  PetscReal              GNpzeta[MAX_QUADP][P_BASIS];
#endif
  PetscScalar            el_A22[P_DOFS*P_DOFS*P_BASIS*P_BASIS];
  const PetscInt         *g_idx_p;
  ISLocalToGlobalMapping ltogm;
  
  PetscFunctionBeginUser;
  
  ierr = MatZeroEntries(S);CHKERRQ(ierr);
  
  nqp   = space->vol_quadrature->n_qpoints;
  w_qp  = space->vol_quadrature->qp_weight;
  xi_qp = space->vol_quadrature->qp_coor;
  
  for (q=0; q<nqp; ++q) {
    ierr = EvaluateBasis_Q1(&xi_qp[NSD*q],Np[q]);CHKERRQ(ierr);
#if NSD == 2
    ierr = EvaluateBasisDerivLocal_Q1(&xi_qp[NSD*q],GNpxi[q],GNpeta[q],0         );CHKERRQ(ierr);
#elif NSD == 3
    ierr = EvaluateBasisDerivLocal_Q1(&xi_qp[NSD*q],GNpxi[q],GNpeta[q],GNpzeta[q]);CHKERRQ(ierr);
#endif
  }
  
  ierr = DMGetLocalToGlobalMapping(dmp,&ltogm);CHKERRQ(ierr);
  ierr = ISLocalToGlobalMappingGetIndices(ltogm,&g_idx_p);CHKERRQ(ierr);
  
  ierr = DMGetCoordinatesLocal(dmp,&coor_p_l);CHKERRQ(ierr);
  ierr = VecGetArray(coor_p_l,&LA_coor_p);CHKERRQ(ierr);

  for (e=0; e<space->n_p_elements; ++e) {
    PetscInt *e_p_idx = &space->p_el_nd_map[P_BASIS*e];
    PetscInt  el_p_idx_g[P_DOFS*P_BASIS];
    
    for (i=0; i<P_BASIS; ++i) {
      PetscInt lnidx = e_p_idx[i];
      el_p_idx[i] = lnidx;
    }

    /* get element coordinates */
    for (i=0; i<P_BASIS; ++i) {
      PetscInt lnidx = e_p_idx[i];
      for (d=0; d<NSD; ++d) el_coor_p[NSD*i+d] = LA_coor_p[NSD*lnidx+d];
    }
    
    PetscMemzero(el_A22,sizeof(PetscScalar)*P_DOFS*P_DOFS*P_BASIS*P_BASIS);

    /* Assemble A22 

    Stokes requires we evaluate
     
     S = - \int (1/\eta) p q \, dV
     
     and Lame' uses the same plus the non-zero 2,2 block 

     S = - \int (1/\eta + 1/\mu) p q \, dV

      Note: in some other examples (ex42 and ex43), instead the weighting is by the inverse of the average viscosity over the element
     
    */
    
    for (q=0; q<nqp; ++q) {
      PetscReal   detJ;
      PetscScalar fac;
#ifdef LAME
      PetscScalar coeff_inv = 1.0/space->coeff_qp[nqp*e+q].lambda + 1.0/space->coeff_qp[nqp*e+q].mu;
#else
      PetscScalar coeff_inv = 1.0/space->coeff_qp[nqp*e+q].eta;
#endif
#if NSD == 2
      ierr = EvaluateBasisTransformation(P_BASIS,GNpxi[q],GNpeta[q],0,el_coor_p,&detJ);CHKERRQ(ierr);
#elif NSD == 3 
      ierr = EvaluateBasisTransformation(P_BASIS,GNpxi[q],GNpeta[q],GNpzeta[q],el_coor_p,&detJ);CHKERRQ(ierr);
#endif
      fac = w_qp[q] * detJ;
      
      for (i=0; i<P_BASIS; ++i) {
        for (j=0; j<P_BASIS; ++j) {
          el_A22[i*P_BASIS+j] -= coeff_inv * Np[q][i] * Np[q][j] * fac;
        }
      }
    }
    
    for (i=0;i<P_DOFS*P_BASIS;++i) {
      el_p_idx_g[i] = g_idx_p[el_p_idx[i]];
    }
    
    ierr = MatSetValues(S, P_DOFS*P_BASIS,el_p_idx_g, P_DOFS*P_BASIS,el_p_idx_g, el_A22,ADD_VALUES);CHKERRQ(ierr);
    
  }
  ierr = VecRestoreArray(coor_p_l,&LA_coor_p);CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(S,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(S,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = ISLocalToGlobalMappingRestoreIndices(ltogm,&g_idx_p);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}
