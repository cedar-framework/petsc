#ifndef _CEDARMAT_H
#define _CEDARMAT_H

#include <petscdmda.h> /*I "petscdmda.h" I*/
#define ENABLE_3D
#include <cedar/capi.h>
#undef ENABLE_3D

typedef struct {
  MPI_Comm            hcomm;
  DM                  da;

  PetscInt dim;
  PetscInt Nx, Ny, Nz;

  cedar_mat c_matrix;
  cedar_config c_config;
  cedar_topo c_topo;

  PetscBool is_initialized;
} Mat_CedarMatrix;

#endif
