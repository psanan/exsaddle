#ifndef __PCILDL_h__
#define __PCILDL_h__

#include <petscpc.h>   /*I "petscpc.h" I*/

#define PCILDL "ildl"

typedef enum {ILDL_ORDERING_METISN,ILDL_ORDERING_METISE,ILDL_ORDERING_RCM,ILDL_ORDERING_AMD} ILDLOrderingType;
PetscErrorCode PCCreate_ILDL(PC pc);

#endif
