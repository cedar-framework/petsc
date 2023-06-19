#include <petscsys.h>
#include <petsc/private/matimpl.h>
#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <../src/dm/impls/da/cedar/cedarmat.h>

static PetscBool cedar_initialized = PETSC_FALSE;
static cedar_config cedar_configuration = (cedar_config) PETSC_NULLPTR;

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
  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Could not find valid Cedar configuration.");
}

PetscErrorCode MatMult_Cedar(Mat A, Vec x, Vec y)
{
  PetscFunctionBegin;
  /* y = Ax */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatAssemblyEnd_Cedar(Mat mat, MatAssemblyType mode)
{
  PetscFunctionBegin;
  /* Any callbacks to be run after the matrix is assembled */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatZeroEntries_Cedar(Mat mat)
{
  PetscFunctionBegin;
  /* Zero out entries in the Cedar matrix */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatDestroy_Cedar(Mat mat)
{
  PetscFunctionBegin;
  /* Free the allocated Cedar data from MatCreate */
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatSetValuesLocal_Cedar2D(Mat mat, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode ins_mode) {
  PetscInt i, j, offset;
  cedar_coord_2d stencil_coords[9];
  cedar_real stencil_values[9];
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) mat->data;
  const PetscInt stride = ex->Nx;

  PetscFunctionBegin;
  PetscCheck(ins_mode == INSERT_VALUES, PETSC_COMM_SELF, PETSC_ERR_SUP, "Only setting values is supported (no updates).");
  PetscCheck(ncol <= 9, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Maximum stencil size for 2D is 9 points, got %" PetscInt_FMT " columns", ncol);

  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      offset = icol[j] - irow[i];
      stencil_coords[j].i = irow[i] / stride;
      stencil_coords[j].j = icol[j] % stride;
      stencil_values[j] = (cedar_real) y[j];

      if (offset == 0) {
        stencil_coords[j].dir = BMG2_C;
      } else if (offset == stride) {
        stencil_coords[j].dir = BMG2_N;
      } else if (offset == stride + 1) {
        stencil_coords[j].dir = BMG2_NE;
      } else if (offset == 1) {
        stencil_coords[j].dir = BMG2_E;
      } else if (offset == -stride + 1) {
        stencil_coords[j].dir = BMG2_SE;
      } else if (offset == -stride) {
        stencil_coords[j].dir = BMG2_S;
      } else if (offset == -stride - 1) {
        stencil_coords[j].dir = BMG2_SW;
      } else if (offset == -1) {
        stencil_coords[j].dir = BMG2_W;
      } else if (offset == stride - 1) {
        stencil_coords[j].dir = BMG2_NW;
      } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Local row %" PetscInt_FMT " local column %" PetscInt_FMT " have bad stencil offset %" PetscInt_FMT, irow[i], icol[j], offset);
      }
    }
    PetscCall(cedar_mat_set2d(ex->c_matrix, ncol, stencil_coords, stencil_values));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatSetValuesLocal_Cedar3D(Mat mat, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode ins_mode) {
  PetscInt i, j, offset;
  cedar_coord_2d stencil_coords[9];
  cedar_real stencil_values[9];
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) mat->data;
  const PetscInt stride1 = ex->Nx;
  const PetscInt stride2 = ex->Nx * ex->Ny;

  PetscFunctionBegin;
  PetscCheck(ins_mode == INSERT_VALUES, PETSC_COMM_SELF, PETSC_ERR_SUP, "Only setting values is supported (no updates).");
  PetscCheck(ncol <= 9, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Maximum stencil size for 2D is 9 points, got %" PetscInt_FMT " columns", ncol);

  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "todo: set values for 3d matrices.");
}

PetscErrorCode MatSetUp_Cedar(Mat mat)
{
  DM da;
  PetscInt dim, Nx, Ny, Nz, Px, Py, Pz, dof, psw;
  const PetscInt* Lx = PETSC_NULLPTR, *Ly = PETSC_NULLPTR, *Lz = PETSC_NULLPTR;
  DMBoundaryType Bx, By, Bz;
  DMDAStencilType stencil_type;
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) mat->data;

  PetscFunctionBegin;

  /* Grab information about the DM */
  PetscCall(MatGetDM(mat, (DM*) &da));
  PetscCall(DMDAGetInfo(da, &dim, &Nx, &Ny, &Nz, &Px, &Py, &Pz, &dof, &psw, &Bx, &By, &Bz, &stencil_type));
  PetscCall(DMDAGetOwnershipRanges(da, &Lx, &Ly, &Lz));

  ex->dim = dim;
  ex->Nx = Nx; ex->Ny = Ny; ex->Nz = Nz;

  /* Initialize Cedar */
  MPI_Comm comm = PetscObjectComm((PetscObject) da);
  PetscCheck(dim == 2 || dim == 3, PetscObjectComm((PetscObject) da), PETSC_ERR_SUP, "Cedar supports only 2D or 3D problems.");
  PetscCheck(dof == 1, PetscObjectComm((PetscObject) da), PETSC_ERR_SUP, "Cedar supports only scalar problems.");

  if (dim == 2) {
    cedar_stencil_2d cedar_stencil_type = (stencil_type == DMDA_STENCIL_STAR ? CEDAR_STENCIL_FIVE_PT : CEDAR_STENCIL_NINE_PT);
    PetscCall(cedar_topo_create2d(comm, Nx, Ny, Lx, Ly, Px, Py, &ex->c_topo));
    PetscCall(cedar_mat_create2d(ex->c_config, ex->c_topo, cedar_stencil_type, &ex->c_matrix));
    mat->ops->setvalueslocal = MatSetValuesLocal_Cedar2D;
  } else {
    cedar_stencil_3d cedar_stencil_type = (stencil_type == DMDA_STENCIL_STAR ? CEDAR_STENCIL_SEVEN_PT : CEDAR_STENCIL_XXVII_PT);
    PetscCall(cedar_topo_create3d(comm, Nx, Ny, Nz, Lx, Ly, Lz, Px, Py, Pz, &ex->c_topo));
    PetscCall(cedar_mat_create3d(ex->c_config, ex->c_topo, cedar_stencil_type, &ex->c_matrix));
    mat->ops->setvalueslocal = MatSetValuesLocal_Cedar3D;
  }

  ex->is_initialized = PETSC_TRUE;
  mat->assembled = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_EXTERN PetscErrorCode MatCreate_Cedar(Mat B)
{
  Mat_CedarMatrix* cmat;

  PetscFunctionBegin;
  /* Initialize Cedar, allocate the data for the matrix, etc. */
  PetscCall(PetscNew(&cmat));
  B->data = (void*) cmat;
  B->rmap->bs = 1;
  B->assembled = PETSC_FALSE;

  B->insertmode = NOT_SET_VALUES;

  B->ops->assemblyend = MatAssemblyEnd_Cedar;
  B->ops->mult = MatMult_Cedar;
  B->ops->zeroentries = MatZeroEntries_Cedar;
  B->ops->setup = MatSetUp_Cedar;

  cmat->is_initialized = PETSC_FALSE;

  PetscCall(cedar_initialize(&cmat->c_config));
  PetscCallMPI(MPI_Comm_dup(PetscObjectComm((PetscObject)B), &(cmat->hcomm)));
  PetscCall(PetscObjectChangeTypeName((PetscObject)B, MATCEDAR));
  PetscFunctionReturn(PETSC_SUCCESS);
}
