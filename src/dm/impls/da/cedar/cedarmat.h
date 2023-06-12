#ifndef _CEDARMAT_H
#define _CEDARMAT_H

#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <cedar/capi.h>

typedef struct {
  MPI_Comm            hcomm;
  DM                  da;

  cedar_config* c_config;
  cedar_topo* c_topo;
  cedar_mat* c_matrix;

  PetscBool needsinitialization;
} Mat_CedarMatrix;

#endif
