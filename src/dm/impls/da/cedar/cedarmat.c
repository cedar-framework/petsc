#include <petscsys.h>
#include <petsc/private/matimpl.h>
#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <../src/dm/impls/da/cedar/cedarmat.h>

static PetscBool cedar_initialized = PETSC_FALSE;
static cedar_config cedar_configuration = PETSC_NULLPTR;

const char* config_fpaths[] = {
  "cedar-config.json",
  "$HOME/cedar-config.json"
};

static PetscErrorCode cedar_initialize(cedar_config* config)
{
  size_t i;

  PetscFunctionBegin;
  if ((void*) cedar_configuration != PETSC_NULLPTR) {
    *config = cedar_configuration;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* Try all the paths for a valid Cedar configuration */
  for (i = 0; i < sizeof(config_fpaths) / sizeof(char*); ++i) {
    if (cedar_config_create(config_fpaths[i], &cedar_configuration) == CEDAR_SUCCESS) {
      cedar_log_init(cedar_configuration);
      cedar_initialized = PETSC_TRUE;
      *config = cedar_configuration;
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }

  /* Could not find any configuration */
  PetscFunctionReturn(PETSC_ERR_USER);
}

PetscErrorCode MatMult_CedarMatrix(Mat A, Vec x, Vec y)
{
  PetscFunctionBegin;
  /* y = Ax */
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

PetscErrorCode MatSetUp_CedarMatrix(Mat mat)
{
  PetscFunctionBegin;
  /* Initialize values in the Cedar matrix object */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN PetscErrorCode MatCreate_CedarMatrix(Mat B)
{
  Mat_CedarMatrix* cmat;

  PetscFunctionBegin;
  /* Initialize Cedar, allocate the data for the matrix, etc. */
  PetscCall(PetscNew(&cmat));
  B->data = (void*) cmat;
  B->rmap->bs = 1;
  B->assembled = PETSC_FALSE;

  B->insertmode = NOT_SET_VALUES;

  B->ops->assemblyend = MatAssemblyEnd_CedarMatrix;
  B->ops->mult = MatMult_CedarMatrix;
  B->ops->zeroentries = MatZeroEntries_CedarMatrix;
  B->ops->setup = MatSetUp_CedarMatrix;

  cmat->is_initialized = PETSC_FALSE;

  PetscCallMPI(MPI_Comm_dup(PetscObjectComm((PetscObject)B), &(cmat->hcomm)));
  PetscCall(PetscObjectChangeTypeName((PetscObject)B, MATCEDAR));
  PetscCall(cedar_initialize(&cmat->c_config));
  PetscFunctionReturn(PETSC_SUCCESS);
}
