#include <petscpkg_version.h>
#include <petsc/private/pcimpl.h> /*I "petscpc.h" I*/
#include <petsc/private/matimpl.h>
#include <petsc/private/vecimpl.h>
#include <../src/dm/impls/da/cedar/cedarmat.h>

typedef struct {
  cedar_config c_config;
  cedar_topo c_topo;
  cedar_mat c_matrix;
  cedar_solver c_solver;

  Mat_CedarMatrix* mat_handle;
  PetscInt dim;

  PetscBool initialized;
} PC_Cedar;

static PetscErrorCode PCSetUpFromMatrix_Cedar(PC_Cedar* pc, Mat_CedarMatrix* mat) {
  PetscFunctionBegin;
  if (pc->initialized) {
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  pc->c_config = mat->c_config;
  pc->c_topo = mat->c_topo;
  pc->c_matrix = mat->c_matrix;
  pc->dim = mat->dim;
  PetscCall(cedar_solver_create(pc->c_matrix, &pc->c_solver));

  pc->initialized = PETSC_TRUE;

  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCSetUp_Cedar(PC pc) {
  PC_Cedar* c_pc = (PC_Cedar*) pc->data;
  PetscFunctionBegin;

  if (c_pc->initialized) {
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  if (pc->pmat->data == PETSC_NULLPTR) {
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscCall(PCSetUpFromMatrix_Cedar(c_pc, (Mat_CedarMatrix*) pc->pmat->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCApply_Cedar(PC pc, Vec b, Vec x) {
  PC_Cedar* c_pc = (PC_Cedar*) pc->data;
  cedar_vec c_x, c_b;

  PetscFunctionBegin;
  PetscCall(PCSetUpFromMatrix_Cedar(c_pc, (Mat_CedarMatrix*) pc->pmat->data));

  if (c_pc->dim == 2) {
    PetscCall(cedar_vec_create2d(c_pc->c_topo, &c_x));
    PetscCall(cedar_vec_create2d(c_pc->c_topo, &c_b));
  } else {
    PetscCall(cedar_vec_create3d(c_pc->c_topo, &c_x));
    PetscCall(cedar_vec_create3d(c_pc->c_topo, &c_b));
  }
  PetscCall(petsc_vec_to_cedar_vec(b, c_b));
  PetscCall(petsc_vec_to_cedar_vec(x, c_x));

  PetscCall(cedar_solver_run(c_pc->c_solver, c_x, c_b));

  PetscCall(cedar_vec_to_petsc_vec(c_x, x));
  PetscCall(cedar_vec_free(&c_x));
  PetscCall(cedar_vec_free(&c_b));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN PetscErrorCode PCCreate_Cedar(PC pc) {
  PC_Cedar* c_pc;

  PetscFunctionBegin;
  PetscCall(PetscNew(&c_pc));

  pc->data = c_pc;
  pc->ops->setup = PCSetUp_Cedar;
  pc->ops->apply = PCApply_Cedar;

  c_pc->c_config = CEDAR_CONFIG_NULL;
  c_pc->c_topo = CEDAR_TOPO_NULL;
  c_pc->c_matrix = CEDAR_MAT_NULL;
  c_pc->c_solver = CEDAR_SOLVER_NULL;
  c_pc->initialized = PETSC_FALSE;
  c_pc->mat_handle = PETSC_NULLPTR;

  PetscFunctionReturn(PETSC_SUCCESS);
}
