#ifndef _CEDARMAT_H
#define _CEDARMAT_H

#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <cedar/capi.h>

typedef struct {
  MPI_Comm            hcomm;
  DM                  da;

  cedar_mat c_matrix;
  cedar_config c_config;

  PetscBool is_initialized;
} Mat_CedarMatrix;

#endif
