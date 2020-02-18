static char help[] = "Solves the "
#ifdef LAME
"Lame' system "
#else
"incompressible, variable viscosity Stokes equations " 
#endif
"in " 
#if NSD==2
"2d "
#elif NSD == 3
"3d "
#endif
"using a Q2-Q1 (Taylor-Hood) mixed finite element formulation. \n"
"-mx : number elements in x-direction \n"
"-my : number elements in y-direction \n"
"-size_x : box dimension in x-direction (note that coefficient structure is not scaled)\n"
"-size_y : box dimension in y-direction (note that coefficient structure is not scaled)\n"
#if NSD == 3
"-size_z : boz dimension in z-direction (note that coefficient structure is not scaled)\n"
#endif
"-model 0 : solcx \n"
"-model 1 : 3 sinkers \n"
"-model 2 : X<=8 sinkers, controlled by -sinker_n \n"
#if NSD == 3
"-model 5: ad-doc 3d solcx-type test\n"
#endif
"-model 6 : single sinker \n"
#if NSD == 3
"-model 7 : lots of sinkers, controlled by -sinker_n \n"
#endif
#ifdef LAME
"-model 8 : one inclusion, fixed base, free elsewhere \n"
"-model 9 : no inclusions, compression \n"
"-model 10: single inclusion, compression \n"
#endif
#if !defined(LAME) && NSD == 3
"-model 11: fixed base, viscosity varying across x\n"
#endif
#if defined(LAME) && NSD == 3
"-model 12: single inclusion, compression with free slip on all faces except top\n"
#endif
#ifdef LAME
"-mu0 : Lame''s second parameter (shear modulus) of background \n"
"-mu1 : Lame''s second parameter (shear modulus) of inclusion\n"
"-lambda0 : Lame''s first parameter of background\n"
"-lambda1 : Lame''s first parameter of inclusion \n"
#else
"-eta0 : viscosity of background \n"
"-eta1 : viscosity of inclusion \n"
#endif
"-sinker_r : radius of the sinker \n"
#if NSD==2
"-sinker_x,-sinker_y : position of sinker (-model 6 only) \n"
#elif NSD==3
"-sinker_x,-sinker_y,-sinker_z : position of sinker (-model 6 only) \n"
#endif
"-nlevels <n> : number of multigrid levels \n"
"-mg : set up for use with top level PCMG\n"
"-fs : set up for top level PCFIELDSPLIT\n"
"-fs_coarse : set up a fieldsplit PC on the coarse MG solver. Don't forget to set -saddle_mg_coarse_ksp_convergence_test default if needbe\n"
"-set_ksp_dm : set the dm for the top-level KSP (needed for 1-level ASM)\n"
"-diagnostics : report information about u,p solution \n"
"-view_fields : dump paraview files for u,v,p \n"
"-view_coeffs : dump coefficient and forcing fields as .vts\n"
"-dump_solution : save PETSc binary file of computed solution vector\n"
"-dump_operator : save PETSc binary file(s) of the operator(s)\n"
"-dump_smoother : save PETSc binary file(s) of smoother(s) [MEANINGLESS if the smoother is nonlinear]\n"
"-dump_preconditioner : save PETSc binary file(s) of explicit PC. SLOW\n"
"-dump_preconditioned_operator : save PETSc binary file of the preconditioned operator. SLOW\n"
"-dump_scaled_mass_matrix : dump the scaled mass matrix, if it is computed (with -fs)\n"
"-refinefactor : (badly named) how much, in each direction, the grid is coarsened to form the MG hierarchy. \n"
"                (for the -fs option, don't use this, use the (badly named) -da_refine_x and -da_refine_y with MG on the viscous block) \n"
"-twosolves : perform a second solve, with top level monitoring off, in a new logging stage.\n"
"-check_solution : compute error wrt a reference solution, if the model provides one\n"
"-constant_pressure_nullspace : define a constant-pressure nullspace on the top-level operator\n"
"-options_file_yaml <options.yml> : load a YAML options file (useful for highly nested solves)\n"
"\n";

/* Note: there are several older versions of this code in old_exSaddle/,
         which include some debugging output which may be of use if 
         doing further development */

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>
#include <petsc/private/dmimpl.h> /* To override ops->createdomaindecomposition */

#ifdef EXSADDLE_WITH_PCILUPACK
#include "pcilupack.h"
#endif
#ifdef EXSADDLE_WITH_PCILDL
#include "pcildl.h"
#endif

#include "exSaddle.h"
#include "femixedspace.h"
#include "exSaddle_io.h"
#include "models.h"

PetscErrorCode SaddleSolve_Q2Q1();
static PetscErrorCode ExtraSolves(KSP ksp, Vec F, Vec X);

/*****************************************************************************/
int main(int argc,char **args)
{
  PetscErrorCode ierr;
  
  ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
#ifdef EXSADDLE_WITH_PCILUPACK
  ierr = PCRegister(PCILUPACK,     PCCreate_ILUPACK     );CHKERRQ(ierr);
#endif
#ifdef EXSADDLE_WITH_PCILDL
  ierr = PCRegister(PCILDL,        PCCreate_ILDL        );CHKERRQ(ierr);
#endif
  ierr = SaddleSolve_Q2Q1();CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}

/*****************************************************************************/
/* A main solve routine that sets up the saddle point system, sets up the solver,
   and produces output */
PetscErrorCode SaddleSolve_Q2Q1()
{
  PetscErrorCode ierr;
  PetscInt       k;
  DM             dm_saddle,dm_saddle_levels[MG_DEPTH],dmp,dmv;
  FEMixedSpace   fespace,fespace_levels[MG_DEPTH];
  Vec            X,F,Xref = NULL;
  KSP            ksp;
  PC             pc;
  Mat            A,A_levels[MG_DEPTH],Mpscaled,Mpscaled_coarse;
  MatNullSpace   nullspace = NULL;
  IS             *is_saddle_field;

  PetscBool flg;
  char      options_file_yaml[PETSC_MAX_PATH_LEN];

  PetscInt  mx                           = 4;
  PetscInt  my;
  PetscInt  mz;                                 /* Ignored if NSD == 2 */
  PetscReal size_x                       = 1.0;
  PetscReal size_y                       = 1.0;
  PetscReal size_z                       = 1.0; /* ignored if NSD == 2 */
  PetscBool fs                           = PETSC_FALSE;
  PetscBool mg                           = PETSC_FALSE;
  PetscBool fs_coarse                    = PETSC_FALSE;
  PetscBool set_ksp_dm                   = PETSC_FALSE;
  PetscInt  nlevels                      = 1;
  PetscInt  refinefactor                 = 2;
  PetscBool view_fields                  = PETSC_FALSE;
  PetscBool diagnostics                  = PETSC_FALSE;
  PetscBool dump_operator                = PETSC_FALSE;
  PetscBool dump_solution                = PETSC_FALSE;
  PetscBool dump_preconditioner          = PETSC_FALSE;
  PetscBool dump_preconditioned_operator = PETSC_FALSE;
  PetscBool dump_smoother                = PETSC_FALSE;
  PetscBool dump_scaled_mass_matrix      = PETSC_FALSE;
  PetscBool twosolves                    = PETSC_FALSE;
  PetscBool pv_monitor                   = PETSC_FALSE;
  PetscBool check_solution               = PETSC_FALSE;
  PetscBool constant_pressure_nullspace  = PETSC_FALSE;


  PetscFunctionBeginUser;

  /* For convenience, this code accepts a .yml file with options */
  ierr = PetscOptionsGetString(NULL,NULL,"-options_file_yaml",options_file_yaml,PETSC_MAX_PATH_LEN-1, &flg);CHKERRQ(ierr);
  if (flg) {
#ifdef PETSC_HAVE_YAML
    const PetscBool require = PETSC_TRUE;
    ierr = PetscOptionsInsertFileYAML(PETSC_COMM_WORLD,options_file_yaml,require);CHKERRQ(ierr);
#else
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"You must configure PETSc with YAML (e.g. --download-yaml) to use -options_file_yaml");
#endif 
  }
  ierr = PetscOptionsGetInt (NULL,NULL,"-mx",                          &mx,                          NULL);CHKERRQ(ierr);
  my = mx;
  mz = mx;
  ierr = PetscOptionsGetInt (NULL,NULL,"-my",                          &my,                          NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (NULL,NULL,"-mz",                          &mz,                          NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-size_x",                      &size_x,                      NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-size_y",                      &size_y,                      NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-size_z",                      &size_z,                      NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-fs",                          &fs,                          NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-mg",                          &mg,                          NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-fs_coarse",                   &fs_coarse,                   NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-set_ksp_dm",                  &set_ksp_dm,                  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (NULL,NULL,"-refinefactor",                &refinefactor,                NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (NULL,NULL,"-nlevels",                     &nlevels,                     NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-view_fields",                 &view_fields,                 NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-diagnostics",                 &diagnostics,                 NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_operator",               &dump_operator,               NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_solution",               &dump_solution,               NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_preconditioner",         &dump_preconditioner,         NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_preconditioned_operator",&dump_preconditioned_operator,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_smoother",               &dump_smoother,               NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_scaled_mass_matrix",     &dump_scaled_mass_matrix,     NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-twosolves",                   &twosolves,                   NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-pv_monitor",                  &pv_monitor,                  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-check_solution",              &check_solution,              NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-constant_pressure_nullspace", &constant_pressure_nullspace, NULL);CHKERRQ(ierr);

  if (fs && mg)                       SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"both -fs and -mg supplied"                       );
  if (nlevels < 1)                    SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-nlevels < 1 supplied"                           );
  if (nlevels > 1 && fs)              SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-nlevels > 1 specified with -fs"                 );
  if (nlevels > 1 && !mg)             SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-nlevels > 1 specified without -mg"              );
  if (nlevels < 2 && mg)              SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-nlevels < 2 specified with -mg"                 );
  if (fs_coarse && !mg)               SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fs_coarse supplied without -mg"                 );
  if (nlevels > MG_DEPTH)             SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"MG levels must be less than %d",MG_DEPTH         );
  if (set_ksp_dm && (mg || fs) )      SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-set_ksp_dm not intended for use with -mg or -fs");
  if (dump_scaled_mass_matrix && !fs) SETERRQ (PETSC_COMM_WORLD,PETSC_ERR_SUP,"-dump_scaled_mass_matrix without -fs"            );

  /* Set up DMs and FE data for all levels */
  {    
    PetscInt mxcoarse=mx,mycoarse=my,mzcoarse=mz;
    if (nlevels > 1) {
      const PetscInt finetocoarseratio = (PetscInt)(PetscPowReal((PetscReal)refinefactor,nlevels-1) + 1e-6);
      if (finetocoarseratio > mx || finetocoarseratio > my || finetocoarseratio > mz) SETERRQ6(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,"Too much refinement %d ^ %d = %d requested for the given problem size (%d x %d x %d elements)",refinefactor,nlevels-1,finetocoarseratio,mx,my,mz);
      if ( mx % finetocoarseratio || my % finetocoarseratio || mz % finetocoarseratio) SETERRQ6(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,"Coarsening ratio of %d ^ %d = %d is incompatible with problem size (%d x %d x %d elements)",refinefactor,nlevels-1,finetocoarseratio,mx,my,mz);
      mxcoarse = mx/finetocoarseratio; 
      mycoarse = my/finetocoarseratio; 
      mzcoarse = mz/finetocoarseratio; 
    }
    for (k=0; k<nlevels; ++k) {
      const PetscInt factor = (PetscInt)(PetscPowReal((PetscReal)refinefactor,k) + 1.0e-6); 
      DM             dmv,dmp;

      ierr = FEMixedSpaceCreate(&fespace_levels[k]);CHKERRQ(ierr);
      ierr = DMCreate_SaddleQ2Q1(PETSC_COMM_WORLD,mxcoarse*factor,mycoarse*factor,mzcoarse*factor,&dm_saddle_levels[k],fespace_levels[k]);CHKERRQ(ierr);
      ierr = DMSetApplicationContext(dm_saddle_levels[k],(void*)fespace_levels[k]);CHKERRQ(ierr);
      dm_saddle_levels[k]->ops->createdomaindecomposition = DMCreateDomainDecomposition_DMDAFEQ2Q1;
      ierr = DMDASetUniformCoordinates_Saddle(dm_saddle_levels[k],0.0,size_x,0.0,size_y,0.0,size_z);CHKERRQ(ierr);
      ierr = FEMixedSpaceQuadratureCreate(fespace_levels[k],PETSC_FALSE,PETSC_TRUE);CHKERRQ(ierr);
      ierr = FEMixedSpaceBCISCreate(fespace_levels[k],dm_saddle_levels[k]);CHKERRQ(ierr);
      ierr = DMCompositeGetEntries(dm_saddle_levels[k],&dmv,&dmp);CHKERRQ(ierr);
      ierr = FEMixedSpaceDefineQPwiseProperties(fespace_levels[k],dmv);CHKERRQ(ierr); 
    }
  }
  fespace   = fespace_levels[nlevels-1];
  dm_saddle =  dm_saddle_levels[nlevels-1];

  {
    DM dmp,dmscalar[MG_DEPTH];

    for (k=0; k<nlevels; ++k) {
      ierr = DMCompositeGetEntries(dm_saddle_levels[k],NULL,&dmp);CHKERRQ(ierr);
      dmscalar[k] = dmp;
    }
    ierr = FEMixedSpaceDefineQPwiseProperties_Q1Projection(nlevels,fespace_levels,dmscalar);CHKERRQ(ierr);
  }

  /* Create System 
    
     We create the system, and additional rhs term due to non-zero
     Dirichlet boundary conditions, which we then add to the 
     usual RHS.
   */
  {
    Vec rhs_diri;
    ierr = DMCreateGlobalVector(dm_saddle,&X);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm_saddle,&F);CHKERRQ(ierr);
    ierr = VecDuplicate(F,&rhs_diri);CHKERRQ(ierr);
    for (k=0; k<nlevels; ++k) {
      Vec rhs_diri_k = k==nlevels-1 ? rhs_diri : NULL ; /* NULL if k<nlevels-10 */
      ierr = DMCreateMatrix(dm_saddle_levels[k],&A_levels[k]);CHKERRQ(ierr); 
      ierr = MatAssemble_Saddle_NULL(fespace_levels[k],dm_saddle_levels[k],A_levels[k]);CHKERRQ(ierr);
      ierr = MatAssemble_Saddle(fespace_levels[k],dm_saddle_levels[k],A_levels[k],rhs_diri_k);CHKERRQ(ierr);
    }
    A = A_levels[nlevels-1];
    ierr = PetscObjectSetName((PetscObject)A,"Asaddle");CHKERRQ(ierr);
    {
      Vec Fu,Fp;
      ierr = DMCompositeGetAccess(dm_saddle,F,&Fu,&Fp);CHKERRQ(ierr);
      ierr = VecAssemble_F1_qp(fespace,Fu);CHKERRQ(ierr);
      ierr = VecAssemble_F2_qp(fespace,Fp);CHKERRQ(ierr);
      ierr = ImposeDirichletValuesIS(Fu,fespace->u_is_global,fespace->u_bc_global);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dm_saddle,F,&Fu,&Fp);CHKERRQ(ierr);
    }
    ierr = VecAXPY(F,1.0,rhs_diri);CHKERRQ(ierr); /* This is adding zero for zero Diri cases, so an optimization could be to skip all of this logic for certain choices of boundary conditions, like SolCx */
    ierr = VecDestroy(&rhs_diri);CHKERRQ(ierr);
  }

  /* Define a constant-pressure null space, if requested (more properly, this would
     be a function of the model, as the BC's determine if you want this). This is NOT expected to work for anything other
     than a solve on the top level (that is, I don't know if it propagates to coarsened or sub-matrices) */
  if (constant_pressure_nullspace) {
      PetscInt     Np;
      Vec          Xnull,Xnullp;  

      ierr = DMCreateGlobalVector(dm_saddle,&Xnull);CHKERRQ(ierr);  
      ierr = VecZeroEntries(Xnull);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(dm_saddle,Xnull,NULL,&Xnullp);CHKERRQ(ierr);
      ierr = VecGetSize(Xnullp,&Np);CHKERRQ(ierr);
      ierr = VecSet(Xnullp,-1.0/(PetscSqrtScalar((PetscScalar)Np)));CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dm_saddle,X,NULL,&Xnullp);CHKERRQ(ierr);
      ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,1,&Xnull,&nullspace);CHKERRQ(ierr);
      ierr = MatSetNullSpace(A,nullspace);CHKERRQ(ierr); 
      ierr = VecDestroy(&Xnull);CHKERRQ(ierr);
  }

  /* Create Solver */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(ksp,"saddle_");CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  /* Fieldsplit-specific setup */
  if (fs) {
    ierr = DMCompositeGetEntries(dm_saddle,&dmv,&dmp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
    ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR);CHKERRQ(ierr);
    ierr = PCFieldSplitSetSchurFactType(pc,PC_FIELDSPLIT_SCHUR_FACT_UPPER);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dmp,&Mpscaled);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)Mpscaled,"Mpscaled");CHKERRQ(ierr);
    ierr = MatAssemble_Schur(fespace,dmp,Mpscaled);CHKERRQ(ierr);
    ierr = PCFieldSplitSetSchurPre(pc,PC_FIELDSPLIT_SCHUR_PRE_USER,Mpscaled);CHKERRQ(ierr);
    ierr = DMCompositeGetGlobalISs(fespace->dm_saddle,&is_saddle_field);CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(pc,"u",is_saddle_field[0]);CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(pc,"p",is_saddle_field[1]);CHKERRQ(ierr);
  }

  /* Setup for 1-level DM-aware PCs. This won't be used much outside of testing,
     but this is need, for instance to test using ASM as the PC */
      if(set_ksp_dm){
        ierr = KSPSetDM(ksp,dm_saddle);CHKERRQ(ierr);
        ierr = KSPSetDMActive(ksp,PETSC_FALSE);CHKERRQ(ierr);
      }

  /* Additional MG-specific setup */
  /* Note: By default, the smoothers will be essentially useless! They must be set from the command line */
  if (mg) {
    ierr = KSPSetDM(ksp,dm_saddle);CHKERRQ(ierr);
    ierr = KSPSetDMActive(ksp,PETSC_FALSE);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCMG);CHKERRQ(ierr);
    ierr = PCMGSetLevels(pc,nlevels,NULL);CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc,PC_MG_GALERKIN_NONE);CHKERRQ(ierr);
    for (k=1; k<nlevels; k++) {
      KSP ksp_smooth;
      Mat P;
      
      ierr = PCMGGetSmoother(pc,k,&ksp_smooth);CHKERRQ(ierr);
      ierr = KSPSetOperators(ksp_smooth,A_levels[k],A_levels[k]);CHKERRQ(ierr);
      ierr = KSPSetDM(ksp_smooth,dm_saddle_levels[k]);CHKERRQ(ierr);
      ierr = KSPSetDMActive(ksp_smooth,PETSC_FALSE);CHKERRQ(ierr);
      ierr = DMCreateInterpolation(dm_saddle_levels[k-1],dm_saddle_levels[k],&P,NULL);CHKERRQ(ierr);
      ierr = PCMGSetInterpolation(pc,k,P);CHKERRQ(ierr);
      ierr = MatDestroy(&P);
    }
    {
      KSP ksp_coarse;
      
      ierr = PCMGGetCoarseSolve(pc,&ksp_coarse);CHKERRQ(ierr);
      ierr = KSPSetOperators(ksp_coarse,A_levels[0],A_levels[0]);CHKERRQ(ierr);
      ierr = KSPSetDM(ksp_coarse,dm_saddle_levels[0]);CHKERRQ(ierr);
      ierr = KSPSetDMActive(ksp_coarse,PETSC_FALSE);CHKERRQ(ierr);
      /* If a command line flag is specified, set up a PCFIELDSPLIT-based coarse grid solver */
      /* GOTCHA: with PCMG many things will be set here, including insidiously that the convergence test will be skipped.
                You may need a flag like -saddle_mg_coarse_ksp_convergence_test default to get the behavior you want */
      if (fs_coarse) {
        DM dmv_coarse,dmp_coarse; 
        PC pc_coarse;
        IS *is_stokes_field_coarse;

        ierr = DMCompositeGetEntries(dm_saddle_levels[0],&dmv_coarse,&dmp_coarse);CHKERRQ(ierr);
        ierr = DMCreateMatrix(dmp_coarse,&Mpscaled_coarse);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)Mpscaled_coarse,"Mpscaled_coarse");CHKERRQ(ierr);
        ierr = MatAssemble_Schur(fespace_levels[0],dmp_coarse,Mpscaled_coarse);CHKERRQ(ierr);

        ierr = KSPGetPC(ksp_coarse,&pc_coarse);CHKERRQ(ierr);
        ierr = PCSetType(pc_coarse,PCFIELDSPLIT);CHKERRQ(ierr);
        ierr = PCFieldSplitSetType(pc_coarse,PC_COMPOSITE_SCHUR);CHKERRQ(ierr);
        ierr = PCFieldSplitSetSchurFactType(pc_coarse,PC_FIELDSPLIT_SCHUR_FACT_UPPER);CHKERRQ(ierr);
        ierr = PCFieldSplitSetSchurPre(pc_coarse,PC_FIELDSPLIT_SCHUR_PRE_USER,Mpscaled_coarse);CHKERRQ(ierr); 
        ierr = DMCompositeGetGlobalISs(fespace_levels[0]->dm_saddle,&is_stokes_field_coarse);CHKERRQ(ierr);
        ierr = PCFieldSplitSetIS(pc_coarse,"u",is_stokes_field_coarse[0]);CHKERRQ(ierr);
        ierr = PCFieldSplitSetIS(pc_coarse,"p",is_stokes_field_coarse[1]);CHKERRQ(ierr);

        ierr = KSPSetFromOptions(ksp_coarse);CHKERRQ(ierr);
        {
          KSP       *sub_ksp_coarse;
          PetscInt  nsplits_coarse;
          PetscBool same = PETSC_FALSE;

          ierr = KSPSetUp(ksp_coarse);CHKERRQ(ierr);
          ierr = KSPGetPC(ksp_coarse,&pc_coarse);CHKERRQ(ierr);
          ierr = PetscObjectTypeCompare((PetscObject)pc_coarse,PCFIELDSPLIT,&same);CHKERRQ(ierr);
          if (same) {
            ierr = PCFieldSplitGetSubKSP(pc_coarse,&nsplits_coarse,&sub_ksp_coarse);CHKERRQ(ierr);
            ierr = KSPSetDM(sub_ksp_coarse[0],dmv_coarse);CHKERRQ(ierr);
            ierr = KSPSetDMActive(sub_ksp_coarse[0],PETSC_FALSE);CHKERRQ(ierr);
            PetscFree(sub_ksp_coarse);
          }
        }
        ierr = ISDestroy(&is_stokes_field_coarse[0]);CHKERRQ(ierr);
        ierr = ISDestroy(&is_stokes_field_coarse[1]);CHKERRQ(ierr);
        ierr = PetscFree(is_stokes_field_coarse);CHKERRQ(ierr);
      }
    }
  }
 
  /* Set From Options */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* Additional Setup for fieldsplit */
  if (fs) {
    KSP       *sub_ksp;
    PetscInt  nsplits;
    PetscBool same = PETSC_FALSE;
    
    ierr = KSPSetUp(ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pc,PCFIELDSPLIT,&same);CHKERRQ(ierr);
    if (same) {
      ierr = PCFieldSplitGetSubKSP(pc,&nsplits,&sub_ksp);CHKERRQ(ierr);
      ierr = KSPSetDM(sub_ksp[0],dmv);CHKERRQ(ierr);
      ierr = KSPSetDMActive(sub_ksp[0],PETSC_FALSE);CHKERRQ(ierr);
      PetscFree(sub_ksp);
    }
  }

  /* Solve */
  ierr = KSPSolve(ksp,F,X);CHKERRQ(ierr);
  if (twosolves) {
    ierr = ExtraSolves(ksp,F,X);CHKERRQ(ierr);
  }

  /* Check Solution, if one exists */
  if (check_solution) {
    ierr = ComputeReferenceSolution(dm_saddle,&Xref);CHKERRQ(ierr);

    if (Xref) {
      Vec Xerror,Xerroru,Xu,Xrefu;
      PetscReal Xrefnorm,abs_err,rel_err;
      PetscReal Xrefnormu,abs_erru,rel_erru;

      if (constant_pressure_nullspace) {
        ierr = MatNullSpaceRemove(nullspace,Xref);CHKERRQ(ierr);
      }
      ierr = VecDuplicate(Xref,&Xerror);CHKERRQ(ierr);
      ierr = VecCopy(Xref,Xerror);CHKERRQ(ierr);
      ierr = VecAXPY(Xerror,-1.0,X);CHKERRQ(ierr);
      ierr = VecNorm(Xref,NORM_2,&Xrefnorm);CHKERRQ(ierr);
      ierr = VecNorm(Xerror,NORM_2,&abs_err);CHKERRQ(ierr);
      rel_err = abs_err/Xrefnorm;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Error in solution:\n");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  abs %g\n",abs_err);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  rel %g\n",rel_err);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");CHKERRQ(ierr);
      ierr = VecDestroy(&Xerror);CHKERRQ(ierr);

      ierr = DMCompositeGetAccess(dm_saddle,X,&Xu,NULL);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(dm_saddle,Xref,&Xrefu,NULL);CHKERRQ(ierr);
      ierr = VecDuplicate(Xrefu,&Xerroru);CHKERRQ(ierr);
      ierr = VecCopy(Xrefu,Xerroru);CHKERRQ(ierr);
      ierr = VecAXPY(Xerroru,-1.0,Xu);CHKERRQ(ierr);
      ierr = VecNorm(Xrefu,NORM_2,&Xrefnormu);CHKERRQ(ierr);
      ierr = VecNorm(Xerroru,NORM_2,&abs_erru);CHKERRQ(ierr);
      rel_erru = abs_erru/Xrefnormu;
      ierr = DMCompositeRestoreAccess(dm_saddle,X,&Xu,NULL);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(dm_saddle,Xref,&Xrefu,NULL);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Error in velocity solution:\n");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  abs %g\n",abs_erru);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  rel %g\n",rel_erru);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");CHKERRQ(ierr);
      ierr = VecDestroy(&Xerroru);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -check_solution supplied but no reference solution available\n");
    }
  }

  /* Post-solve output */
  if (diagnostics) {
    ierr = SaddleReportSolutionDiagnostics(dm_saddle,X);CHKERRQ(ierr);
  }

  if (view_fields) {
    ierr = ViewFields(dm_saddle,X,"");CHKERRQ(ierr);
    if (Xref) {
      ierr = ViewFields(dm_saddle,Xref,"ref_");CHKERRQ(ierr);
    }
  }

  if (dump_solution) {
    ierr = DumpSolution(X,"solution.petscbin");CHKERRQ(ierr);
    if (Xref) {
      ierr = DumpSolution(Xref,"ref_solution.petscbin");CHKERRQ(ierr);
    }
  }

  if (dump_operator) {
    for (k=0; k<nlevels; k++) {
      char name[PETSC_MAX_PATH_LEN];
      ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"operator_%D.petscbin",k);CHKERRQ(ierr);
      ierr = DumpOperator(A_levels[k],name);CHKERRQ(ierr);
    }
  }

  if (dump_preconditioner) {
    ierr = DumpPreconditioner(ksp,"preconditioner.petscbin");CHKERRQ(ierr);
  }

  if (dump_preconditioned_operator) {
    ierr = DumpPreconditionedOperator(ksp,"preconditioned_operator_out.petscbin");CHKERRQ(ierr);
  }

  if (dump_smoother) {
    // TODO move to a function in exSaddle_io.c
    PC        pc;
    PCType    pctype;
    PetscBool ismg;
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCGetType(pc,&pctype);CHKERRQ(ierr);
    ierr = PetscStrcmp(pctype,PCMG,&ismg);CHKERRQ(ierr);
    if (ismg) {
      for (k=1; k<nlevels; ++k) { /* don't dump the coarse grid */
        Mat  S;
        KSP  ksp_smoother;
        char name[PETSC_MAX_PATH_LEN];
        ierr = PCMGGetSmoother(pc,k,&ksp_smoother);CHKERRQ(ierr);
        ierr = PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"smoother_%D",k);CHKERRQ(ierr);
        ierr = KSPComputeOperator(ksp_smoother,NULL,&S);CHKERRQ(ierr);
        ierr = DumpOperator(S,name);CHKERRQ(ierr);
        ierr = MatDestroy(&S);CHKERRQ(ierr);
      }
    } else {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Smoother dump requires PC type PCMG");
    }
  }

  if(dump_scaled_mass_matrix) {
    ierr = DumpOperator(Mpscaled,"mpscaled.petscbin");CHKERRQ(ierr);
  }

  /* Cleanup */
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  if (Xref) {
    ierr = VecDestroy(&Xref);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  for (k=0; k<nlevels; ++k) {
    ierr = MatDestroy(&A_levels[k]);CHKERRQ(ierr);
  }
  if (nullspace) {
    ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  for (k=0; k<nlevels; ++k) {
    ierr = DMDestroy(&dm_saddle_levels[k]);CHKERRQ(ierr);
    ierr = FEMixedSpaceDestroy(&fespace_levels[k]);CHKERRQ(ierr);
  }
  if (fs) {
    ierr = MatDestroy(&Mpscaled);CHKERRQ(ierr);
    ierr = ISDestroy(&is_saddle_field[0]);CHKERRQ(ierr);
    ierr = ISDestroy(&is_saddle_field[1]);CHKERRQ(ierr);
    ierr = PetscFree(is_saddle_field);CHKERRQ(ierr);
  }
  if (fs_coarse){
      ierr = MatDestroy(&Mpscaled_coarse);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*****************************************************************************/
static PetscErrorCode ExtraSolves(KSP ksp, Vec F, Vec X)
{
  PetscErrorCode     ierr;
  const PetscInt     num_extra_solves = 1;
  PetscInt           i;
  PetscLogStage      stage;
  KSPConvergedReason reason;
  PetscBool          nonzero_initial_guess;

  PetscFunctionBeginUser;

  ierr = KSPMonitorCancel(ksp);CHKERRQ(ierr);
  ierr = KSPGetInitialGuessNonzero(ksp,&nonzero_initial_guess);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
  "  Commencing with %D additional solves. This will cancel a KSP monitor set on\n\
   saddle_, but no nested output. You should ensure that there is no output between\n\
   this output and the output which indicates the extra solves are completed. That\n\
   is, you should not use any ksp_view, ksp_converged_reason, or nested ksp_monitor\n\
   options if you want the results in this test to be meaningful.\n",
   num_extra_solves);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------------------\n");CHKERRQ(ierr);

  ierr = PetscLogStageRegister("Extra Solves",&stage);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stage);CHKERRQ(ierr);
  for (i=0;i<num_extra_solves;++i) {
    ierr = KSPSolve(ksp,F,X);CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
  if (reason<0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------------------------------------------\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n ERROR: EXTRA SOLVES(S) DIVERGED!\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  } else {
    PetscInt its;
    PetscReal rnorm;
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp,&rnorm);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n------------------------------------------------------------------------------\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  %D extra solve(s) succeeded with %D iterations and residual norm %1.6e \n",num_extra_solves,its,rnorm);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------------------\n");CHKERRQ(ierr);
  }
  ierr = KSPSetInitialGuessNonzero(ksp,nonzero_initial_guess);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
