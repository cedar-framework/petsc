#ifndef _CEDARMAT_H
#define _CEDARMAT_H

#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <petscvec.h>

#define ENABLE_3D
#include <cedar/capi.h>
#undef ENABLE_3D

typedef struct {
  MPI_Comm hcomm;
  DM da;

  PetscInt dim;
  PetscInt gnx, gny, gnz; /* global dimensions */
  PetscInt lx, ly, lz, lm, ln, lp; /* local coordinates: starting positions (x,y,z) and sizes (m,n,p) */
  PetscInt lgx, lgy, lgz, lgm, lgn, lgp; /* local ghost coordinates: starting positions (x,y,z) and sizes (m,n,p) */

  PetscInt row_start;
  const PetscInt* global_indices;

  cedar_mat c_matrix;
  cedar_config c_config;
  cedar_topo c_topo;

  union {
    cedar_stencil_2d c_stencil_2d;
    cedar_stencil_3d c_stencil_3d;
  };

  PetscBool is_initialized;
} Mat_CedarMatrix;

PetscErrorCode cedar_vec_copy_to(const PetscScalar* source, cedar_vec dest);
PetscErrorCode cedar_vec_copy_from(cedar_vec source, PetscScalar* dest);
PetscErrorCode petsc_vec_to_cedar_vec(Vec source, cedar_vec dest);
PetscErrorCode cedar_vec_to_petsc_vec(cedar_vec source, Vec dest);

#endif
