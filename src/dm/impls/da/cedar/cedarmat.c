#include <petscsys.h>
#include <petsc/private/matimpl.h>
#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <../src/dm/impls/da/cedar/cedarmat.h>

PetscErrorCode MatMult_CedarMatrix(Mat A, Vec x, Vec y)
{
  PetscFunctionBegin;
  /* Initialize values in the Cedar matrix object */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatAssemblyEnd_CedarMatrix(Mat mat, MatAssemblyType mode)
{
  PetscFunctionBegin;
  /* Any callbacks to be run after the matrix is assembled */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatZeroEntries_CedarMatrix(Mat mat)
{
  PetscFunctionBegin;
  /* Zero out entries in the Cedar matrix */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatDestroy_CedarMatrix(Mat mat)
{
  PetscFunctionBegin;
  /* Free the allocated Cedar data from MatCreate */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN PetscErrorCode MatCreate_CedarMatrix(Mat B)
{
  PetscFunctionBegin;
  /* Initialize Cedar, allocate the data for the matrix, etc. */
  PetscFunctionReturn(PETSC_SUCCESS);
}
