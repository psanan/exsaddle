#include "exSaddle_io.h"

#include <petscdmcomposite.h>
#include <petscdmda.h>

/*****************************************************************************/
PetscErrorCode SaddleReportSolutionDiagnostics(DM dm_saddle,Vec X)
{
  PetscErrorCode ierr;
  Vec            Xu,Xp;
  PetscReal      nrm[3];
  PetscReal      range[NSD];
 
  PetscFunctionBeginUser;
  ierr = DMCompositeGetAccess(dm_saddle,X,&Xu,&Xp);CHKERRQ(ierr);
  VecStrideNormAll(Xu,NORM_1,nrm);CHKERRQ(ierr);
#if NSD == 2
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v|_1   %+1.6e , %+1.6e \n",nrm[0],nrm[1]);CHKERRQ(ierr);
#elif NSD == 3 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v,w|_1   %+1.6e , %+1.6e , %+1.6e\n",nrm[0],nrm[1],nrm[2]);CHKERRQ(ierr);
#endif
  ierr = VecStrideNormAll(Xu,NORM_2,nrm);CHKERRQ(ierr);
#if NSD == 2
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v|_2   %+1.6e , %+1.6e \n",nrm[0],nrm[1]);CHKERRQ(ierr);
#elif NSD == 3 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v,w|_2   %+1.6e , %+1.6e , %+1.6e\n",nrm[0],nrm[1],nrm[2]);CHKERRQ(ierr);
#endif
  ierr = VecStrideNormAll(Xu,NORM_INFINITY,nrm);CHKERRQ(ierr);
#if NSD == 2
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v|_inf %+1.6e , %+1.6e \n",nrm[0],nrm[1]);CHKERRQ(ierr);
#elif NSD == 3 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v,w|_inf %+1.6e , %+1.6e , %+1.6e\n",nrm[0],nrm[1],nrm[2]);CHKERRQ(ierr);
#endif
  ierr = VecStrideMinAll(Xu,NULL,range);CHKERRQ(ierr);
#if NSD == 2
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v|_min %+1.6e , %+1.6e \n",range[0],range[1]);CHKERRQ(ierr);
#elif NSD == 3 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v,w|_min %+1.6e , %+1.6e , %+1.6e\n",range[0],range[1],range[2]);CHKERRQ(ierr);
#endif
  ierr = VecStrideMaxAll(Xu,NULL,range);CHKERRQ(ierr);
#if NSD == 2
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v|_max %+1.6e , %+1.6e \n",range[0],range[1]);CHKERRQ(ierr);
#elif NSD == 3 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u,v,w|_max %+1.6e , %+1.6e , %+1.6e\n",range[0],range[1],range[2]);CHKERRQ(ierr);
#endif
  ierr = VecNorm(Xp,NORM_1,&nrm[0]);CHKERRQ(ierr);
  ierr = VecNorm(Xp,NORM_2,&nrm[1]);CHKERRQ(ierr);
  ierr = VecNorm(Xp,NORM_INFINITY,&nrm[2]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|p|_1          %+1.6e\n",nrm[0]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|p|_2          %+1.6e\n",nrm[1]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|p|_inf        %+1.6e\n",nrm[2]);CHKERRQ(ierr);
  ierr = VecMin(Xp,NULL,&range[0]);CHKERRQ(ierr);
  ierr = VecMax(Xp,NULL,&range[1]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|p|_min        %+1.6e\n",range[0]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|p|_max        %+1.6e\n",range[1]);CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(dm_saddle,X,&Xu,&Xp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode DumpOperator(Mat A,const char* name)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Dumping operator to %s. This could be very slow!\n",name);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = MatView(A,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished dumping operator to %s.\n",name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode DumpSolution(Vec X,const char* name)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Dumping solution vector to %s.\n",name);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(X,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished dumping vector to %s.\n",name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode DumpPreconditioner(KSP ksp,const char* name)
{
  PetscErrorCode ierr;
  PC             pc;
  PetscViewer    viewer;
  Mat            B;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Dumping explicit PC to %s. This could be very slow!\n",name);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCComputeExplicitOperator(pc,&B);CHKERRQ(ierr); /* Slow */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = MatView(B,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished dumping explicit PC to %s.\n",name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*****************************************************************************/
PetscErrorCode DumpPreconditionedOperator(KSP ksp,const char* name)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  Mat            M;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Dumping preconditioned operator to %s. This could be very slow!\n",name);CHKERRQ(ierr);
  ierr = KSPComputeExplicitOperator(ksp,&M);CHKERRQ(ierr); /* Slow */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = MatView(M,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished dumping preconditioned operator to %s.\n",name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/******************************************************************************/
PetscErrorCode ViewFields(DM dm_saddle,Vec X,const char* tag)
{
  PetscErrorCode ierr;
  PetscViewer  viewer;
  Vec          Xu,Xp,scalar;
  DM           dmn,dmp,dmv;
  char         name[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
  ierr = DMDAGetReducedDMDA(dmv,1,&dmn);CHKERRQ(ierr); /* nodal dmda */
  ierr = DMDASetFieldName(dmn,0,"");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dmn,&scalar);CHKERRQ(ierr);

  ierr = DMCompositeGetAccess(dm_saddle,X,&Xu,&Xp);CHKERRQ(ierr);
#if NSD == 2 
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%suv.vts",tag);CHKERRQ(ierr);
#elif NSD == 3
  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%suvw.vts",tag);CHKERRQ(ierr);
#endif
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)scalar,"u");CHKERRQ(ierr);
  ierr = VecStrideGather(Xu,0,scalar,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecView(scalar,viewer);CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject)scalar,"v");CHKERRQ(ierr);
  ierr = VecStrideGather(Xu,1,scalar,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecView(scalar,viewer);CHKERRQ(ierr);

#if NSD == 3
  ierr = PetscObjectSetName((PetscObject)scalar,"w");CHKERRQ(ierr);
  ierr = VecStrideGather(Xu,2,scalar,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecView(scalar,viewer);CHKERRQ(ierr);
#endif

  ierr = VecDestroy(&scalar);CHKERRQ(ierr);
  ierr = DMDestroy(&dmn);CHKERRQ(ierr); 
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%sp.vts",tag);CHKERRQ(ierr);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Xp,"p");CHKERRQ(ierr);
  ierr = VecView(Xp,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(dm_saddle,X,&Xu,&Xp);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
