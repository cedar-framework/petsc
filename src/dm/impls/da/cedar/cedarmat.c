#include <petscsys.h>
#include <petsc/private/matimpl.h>
#include <petscdmda.h> /*I "petscdmda.h" I*/
#include <../src/dm/impls/da/cedar/cedarmat.h>

static PetscBool cedar_initialized = PETSC_FALSE;
static cedar_config cedar_configuration = CEDAR_CONFIG_NULL;

const char* config_fpaths[] = {
  "cedar-config.json",
  "$HOME/cedar-config.json"
};

static PetscErrorCode cedar_initialize(cedar_config* config)
{
  size_t i;

  PetscFunctionBegin;
  if (cedar_configuration != CEDAR_CONFIG_NULL) {
    *config = cedar_configuration;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* Try all the paths for a valid Cedar configuration */
  for (i = 0; i < sizeof(config_fpaths) / sizeof(char*); ++i) {
    if (cedar_config_create(config_fpaths[i], &cedar_configuration) == CEDAR_SUCCESS) {
      PetscCall(cedar_log_init(cedar_configuration));
      cedar_initialized = PETSC_TRUE;
      *config = cedar_configuration;
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }

  /* Could not find any configuration */
  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Could not find valid Cedar configuration.");
}

/* These are based on the Stella implementations: https://github.com/andrewreisner/stella/blob/master/src/stella_cedar.c */
PetscErrorCode cedar_vec_copy_to(const PetscScalar* source, cedar_vec dest) {
  cedar_real* dest_ptr;
  int dim;
  cedar_len i, j, k;
  cedar_len i_len, j_len, k_len;
  cedar_len src_ptr = 0;

  PetscFunctionBegin;
  PetscCall(cedar_vec_baseptr(dest, &dest_ptr));
  PetscCall(cedar_vec_getdim(dest, &dim));

  if (dim == 2) {
    PetscCall(cedar_vec_len2d(dest, &i_len, &j_len));

    for (j = 1; j < j_len - 1; ++j) {
      for (i = 1; i < i_len - 1; ++i) {
        dest_ptr[j * i_len + i] = source[src_ptr++];
      }
    }
  } else {
    PetscCall(cedar_vec_len3d(dest, &i_len, &j_len, &k_len));

    for (k = 1; k < k_len - 1; ++k) {
      for (j = 1; j < j_len - 1; ++j) {
        for (i = 1; i < i_len - 1; ++i) {
          dest_ptr[k * (i_len * j_len) + j * i_len + i] = source[src_ptr++];
        }
      }
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode petsc_vec_to_cedar_vec(Vec source, cedar_vec dest) {
  const PetscScalar* source_ptr;

  PetscFunctionBegin;
  PetscCall(VecGetArrayRead(source, &source_ptr));
  PetscCall(cedar_vec_copy_to(source_ptr, dest));
  PetscCall(VecRestoreArrayRead(source, &source_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode cedar_vec_copy_from(cedar_vec source, PetscScalar* dest) {
  cedar_real* source_ptr;
  int dim;
  cedar_len i, j, k;
  cedar_len i_len, j_len, k_len;
  cedar_len src_ptr = 0;

  PetscFunctionBegin;
  PetscCall(cedar_vec_baseptr(source, &source_ptr));
  PetscCall(cedar_vec_getdim(source, &dim));

  if (dim == 2) {
    PetscCall(cedar_vec_len2d(source, &i_len, &j_len));

    for (j = 1; j < j_len - 1; ++j) {
      for (i = 1; i < i_len - 1; ++i) {
        dest[src_ptr++] = source_ptr[j * i_len + i];
      }
    }
  } else {
    PetscCall(cedar_vec_len3d(source, &i_len, &j_len, &k_len));

    for (k = 1; k < k_len - 1; ++k) {
      for (j = 1; j < j_len - 1; ++j) {
        for (i = 1; i < i_len - 1; ++i) {
          dest[src_ptr++] = source_ptr[k * (i_len * j_len) + j * i_len + i];
        }
      }
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode cedar_vec_to_petsc_vec(cedar_vec source, Vec dest) {
  PetscScalar* dest_ptr;

  PetscFunctionBegin;
  PetscCall(VecGetArrayWrite(dest, &dest_ptr));
  PetscCall(cedar_vec_copy_from(source, dest_ptr));
  PetscCall(VecRestoreArrayWrite(dest, &dest_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_Cedar(Mat A, Vec x, Vec y)
{
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) A->data;
  cedar_vec c_x = CEDAR_VEC_NULL, c_y = CEDAR_VEC_NULL;
  cedar_topo c_topo = ex->c_topo;

  PetscFunctionBegin;

  /* y = Ax */
  if (ex->dim == 2) {
    PetscCall(cedar_vec_create2d(c_topo, &c_x));
    PetscCall(cedar_vec_create2d(c_topo, &c_y));
  } else {
    PetscCall(cedar_vec_create3d(c_topo, &c_x));
    PetscCall(cedar_vec_create3d(c_topo, &c_y));
  }

  /* Perform the mat-vec, copying data to and from Cedar/Petsc */
  PetscCall(petsc_vec_to_cedar_vec(x, c_x));
  PetscCall(cedar_matvec(ex->c_matrix, c_x, c_y));
  PetscCall(cedar_vec_to_petsc_vec(c_y, y));

  /* Free scratch space for Cedar */
  PetscCall(cedar_vec_free(&c_x));
  PetscCall(cedar_vec_free(&c_y));

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
  PetscInt row, col, i, j, colptr, offset;
  PetscBool valid;
  cedar_coord_2d stencil_coords[9];
  cedar_real stencil_values[9];
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) mat->data;
  const PetscInt Nx = ex->Nx;
  const PetscInt Ny = ex->Ny;
  const PetscInt stride = Nx;

  PetscFunctionBegin;
  PetscCheck(ins_mode == INSERT_VALUES, PETSC_COMM_SELF, PETSC_ERR_SUP, "Only setting values is supported (no updates).");
  PetscCheck(ncol <= 9, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Maximum stencil size for 2D is 9 points, got %" PetscInt_FMT " columns", ncol);

  for (row = 0; row < nrow; ++row) {
    colptr = 0;
    for (col = 0; col < ncol; ++col) {
      offset = icol[col] - irow[row];
      i = irow[row] % stride;
      j = irow[row] / stride;

      valid = PETSC_FALSE;

      /* Gracefully handle boundary conditions: we silently ignore parts of the stencil that will
        refer to out-of-bounds nodes */

      if (offset == 0) {
        stencil_coords[colptr].dir = BMG2_C;
        valid = PETSC_TRUE;
      } else if (offset == stride && j < Ny) {
        stencil_coords[colptr].dir = BMG2_N;
        valid = PETSC_TRUE;
      } else if (offset == stride + 1 && j < Ny && i < Nx) {
        stencil_coords[colptr].dir = BMG2_NE;
        valid = PETSC_TRUE;
      } else if (offset == 1 && i < Nx) {
        stencil_coords[colptr].dir = BMG2_E;
        valid = PETSC_TRUE;
      } else if (offset == -stride + 1 && j > 0 && i < Nx) {
        stencil_coords[colptr].dir = BMG2_SE;
        valid = PETSC_TRUE;
      } else if (offset == -stride && j > 0) {
        stencil_coords[colptr].dir = BMG2_S;
        valid = PETSC_TRUE;
      } else if (offset == -stride - 1 && j > 0 && i > 0) {
        stencil_coords[colptr].dir = BMG2_SW;
        valid = PETSC_TRUE;
      } else if (offset == -1 && i > 0) {
        stencil_coords[colptr].dir = BMG2_W;
        valid = PETSC_TRUE;
      } else if (offset == stride - 1 && j < Ny && i > 0) {
        stencil_coords[colptr].dir = BMG2_NW;
        valid = PETSC_TRUE;
      }

      if (valid) {
        stencil_coords[colptr].i = i;
        stencil_coords[colptr].j = j;
        stencil_values[colptr] = (cedar_real) y[col];
        colptr++;
      }
    }
    PetscCall(cedar_mat_set2d(ex->c_matrix, colptr, stencil_coords, stencil_values));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatZeroRowsLocal_Cedar2D(Mat mat, PetscInt nrow, const PetscInt irow[], PetscScalar d, Vec x, Vec b) {
  cedar_coord_2d stencil_coords[9];
  cedar_real stencil_values[9] = {0};
  PetscInt i, j;
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) mat->data;
  const PetscInt stride = ex->Nx;

  PetscFunctionBegin;
  PetscCheck(!x || !b, PetscObjectComm((PetscObject)mat), PETSC_ERR_SUP, "No support");

  stencil_coords[0].dir = BMG2_C;
  stencil_coords[1].dir = BMG2_N;
  stencil_coords[2].dir = BMG2_E;
  stencil_coords[3].dir = BMG2_S;
  stencil_coords[4].dir = BMG2_W;
  stencil_coords[5].dir = BMG2_NE;
  stencil_coords[6].dir = BMG2_SE;
  stencil_coords[7].dir = BMG2_SW;
  stencil_coords[8].dir = BMG2_NW;

  const int n_stencil = (ex->c_stencil_2d == CEDAR_STENCIL_FIVE_PT ? 5 : 9);

  for (i = 0; i < nrow; ++i) {
      for (j = 0; j < n_stencil; ++j) {
          stencil_coords[j].i = irow[i] / stride;
          stencil_coords[j].j = irow[i] % stride;
      }
      PetscCall(cedar_mat_set2d(ex->c_matrix, n_stencil, stencil_coords, stencil_values));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#pragma GCC diagnostic ignored "-Wunused-variable"
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
#pragma GCC diagnostic pop

PetscErrorCode MatSetUp_Cedar(Mat mat)
{
  DM da;
  PetscInt dim, Nx, Ny, Nz, Px, Py, Pz, dof, psw;
  const PetscInt* Lx = PETSC_NULLPTR, *Ly = PETSC_NULLPTR, *Lz = PETSC_NULLPTR;
  unsigned int* c_Lx = PETSC_NULLPTR, *c_Ly = PETSC_NULLPTR, *c_Lz = PETSC_NULLPTR;
  DMBoundaryType Bx, By, Bz;
  DMDAStencilType stencil_type;
  Mat_CedarMatrix* ex = (Mat_CedarMatrix*) mat->data;
  size_t i = 0;

  PetscFunctionBegin;

  /* Grab information about the DM */
  PetscCall(MatGetDM(mat, (DM*) &da));
  PetscCall(DMDAGetInfo(da, &dim, &Nx, &Ny, &Nz, &Px, &Py, &Pz, &dof, &psw, &Bx, &By, &Bz, &stencil_type));
  PetscCall(DMDAGetOwnershipRanges(da, &Lx, &Ly, &Lz));

  ex->dim = dim;
  ex->Nx = Nx; ex->Ny = Ny; ex->Nz = Nz;

  /* Convert Lx, Ly, Lz to correct Cedar types */
  PetscCall(PetscMalloc1(sizeof(unsigned int) * Px, &c_Lx));
  for (i = 0; i < Px; ++i) {
    c_Lx[i] = (unsigned int) Lx[i];
  }

  PetscCall(PetscMalloc1(sizeof(unsigned int) * Py, &c_Ly));
  for (i = 0; i < Py; ++i) {
    c_Ly[i] = (unsigned int) Ly[i];
  }

  if (dim == 3) {
    PetscCall(PetscMalloc1(sizeof(unsigned int) * Pz, &c_Lz));
    for (i = 0; i < Pz; ++i) {
      c_Lz[i] = (unsigned int) Lz[i];
    }
  }

  /* Initialize Cedar */
  MPI_Comm comm = PetscObjectComm((PetscObject) da);
  PetscCheck(dim == 2 || dim == 3, PetscObjectComm((PetscObject) da), PETSC_ERR_SUP, "Cedar supports only 2D or 3D problems.");
  PetscCheck(dof == 1, PetscObjectComm((PetscObject) da), PETSC_ERR_SUP, "Cedar supports only scalar problems.");

  if (dim == 2) {
    cedar_stencil_2d cedar_stencil_type = (stencil_type == DMDA_STENCIL_STAR ? CEDAR_STENCIL_FIVE_PT : CEDAR_STENCIL_NINE_PT);
    ex->c_stencil_2d = cedar_stencil_type;
    PetscCall(cedar_topo_create2d(comm, Nx, Ny, c_Lx, c_Ly, Px, Py, &ex->c_topo));
    PetscCall(cedar_mat_create2d(ex->c_config, ex->c_topo, cedar_stencil_type, &ex->c_matrix));
    mat->ops->setvalueslocal = MatSetValuesLocal_Cedar2D;
    mat->ops->zerorowslocal  = MatZeroRowsLocal_Cedar2D;
  } else {
    cedar_stencil_3d cedar_stencil_type = (stencil_type == DMDA_STENCIL_STAR ? CEDAR_STENCIL_SEVEN_PT : CEDAR_STENCIL_XXVII_PT);
    ex->c_stencil_3d = cedar_stencil_type;
    PetscCall(cedar_topo_create3d(comm, Nx, Ny, Nz, c_Lx, c_Ly, c_Lz, Px, Py, Pz, &ex->c_topo));
    PetscCall(cedar_mat_create3d(ex->c_config, ex->c_topo, cedar_stencil_type, &ex->c_matrix));
    mat->ops->setvalueslocal = MatSetValuesLocal_Cedar3D;
  }

  PetscCall(PetscFree(c_Lx));
  PetscCall(PetscFree(c_Ly));
  if (dim == 3) {
    PetscCall(PetscFree(c_Lz));
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
