
#include <petsc/private/matimpl.h> /*I "petscmat.h" I*/

typedef struct {
  Vec              diag;
  PetscBool        diag_valid;
  Vec              inv_diag;
  PetscBool        inv_diag_valid;
  PetscObjectState diag_state, inv_diag_state;
} Mat_Diagonal;

static PetscErrorCode MatDiagonalSetUpDiagonal(Mat A)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  if (!ctx->diag_valid) {
    PetscAssert(ctx->inv_diag_valid, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Neither diagonal nor inverse diagonal is in a valid state");
    PetscCall(VecCopy(ctx->inv_diag, ctx->diag));
    PetscCall(VecReciprocal(ctx->diag));
    ctx->diag_valid = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalSetUpInverseDiagonal(Mat A)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  if (!ctx->inv_diag_valid) {
    PetscAssert(ctx->diag_valid, PetscObjectComm((PetscObject)A), PETSC_ERR_PLIB, "Neither diagonal nor inverse diagonal is in a valid state");
    PetscCall(VecCopy(ctx->diag, ctx->inv_diag));
    PetscCall(VecReciprocal(ctx->inv_diag));
    ctx->inv_diag_valid = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatAXPY_Diagonal(Mat Y, PetscScalar a, Mat X, MatStructure str)
{
  Mat_Diagonal *yctx = (Mat_Diagonal *)Y->data;
  Mat_Diagonal *xctx = (Mat_Diagonal *)X->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(Y));
  PetscCall(MatDiagonalSetUpDiagonal(X));
  PetscCall(VecAXPY(yctx->diag, a, xctx->diag));
  yctx->inv_diag_valid = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_Diagonal(Mat A, Vec x, Vec y)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(A));
  PetscCall(VecPointwiseMult(y, ctx->diag, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMultAdd_Diagonal(Mat mat, Vec v1, Vec v2, Vec v3)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)mat->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(mat));
  if (v2 != v3) {
    PetscCall(VecPointwiseMult(v3, ctx->diag, v1));
    PetscCall(VecAXPY(v3, 1.0, v2));
  } else {
    Vec w;
    PetscCall(VecDuplicate(v3, &w));
    PetscCall(VecPointwiseMult(w, ctx->diag, v1));
    PetscCall(VecAXPY(v3, 1.0, w));
    PetscCall(VecDestroy(&w));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatNorm_Diagonal(Mat A, NormType type, PetscReal *nrm)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(A));
  type = (type == NORM_FROBENIUS) ? NORM_2 : NORM_INFINITY;
  PetscCall(VecNorm(ctx->diag, type, nrm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDuplicate_Diagonal(Mat A, MatDuplicateOption op, Mat *B)
{
  Mat_Diagonal *actx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  PetscCall(MatCreate(PetscObjectComm((PetscObject)A), B));
  PetscCall(MatSetSizes(*B, A->rmap->n, A->cmap->n, A->rmap->N, A->cmap->N));
  PetscCall(MatSetBlockSizesFromMats(*B, A, A));
  PetscCall(MatSetType(*B, MATDIAGONAL));
  PetscCall(PetscLayoutReference(A->rmap, &(*B)->rmap));
  PetscCall(PetscLayoutReference(A->cmap, &(*B)->cmap));
  if (op == MAT_COPY_VALUES) {
    Mat_Diagonal *bctx = (Mat_Diagonal *)(*B)->data;

    PetscCall(MatSetUp(A));
    PetscCall(MatSetUp(*B));
    PetscCall(VecCopy(actx->diag, bctx->diag));
    PetscCall(VecCopy(actx->inv_diag, bctx->inv_diag));
    bctx->diag_valid     = actx->diag_valid;
    bctx->inv_diag_valid = actx->inv_diag_valid;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
  MatDiagonalGetDiagonal - Get the diagonal of a `MATDIAGONAL`

  Input Parameter:
. A - the `MATDIAGONAL`

  Output Parameter:
. diag - the `Vec` that defines the diagonal

  Level: developer

  Note:
  The user must call
  `MatDiagonalRestoreDiagonal()` before using the matrix again.

  For a copy of the diagonal values, rather than a reference, use `MatGetDiagonal()`

  Any changes to the obtained vector immediately change the action of the `Mat`. The matrix can be changed more efficiently by accessing this vector and changing its values, instead of filling a work vector and using `MatDiagonalSet()`

.seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalRestoreDiagonal()`, `MatDiagonalGetInverseDiagonal()`, `MatGetDiagonal()`
@*/
PetscErrorCode MatDiagonalGetDiagonal(Mat A, Vec *diag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(A, MAT_CLASSID, 1);
  PetscValidPointer(diag, 2);
  *diag = NULL;
  PetscUseMethod((PetscObject)A, "MatDiagonalGetDiagonal_C", (Mat, Vec *), (A, diag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalGetDiagonal_Diagonal(Mat A, Vec *diag)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  PetscCall(MatSetUp(A));
  PetscCall(MatDiagonalSetUpDiagonal(A));
  *diag = ctx->diag;
  PetscCall(PetscObjectStateGet((PetscObject)*diag, &ctx->diag_state));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
  MatDiagonalRestoreDiagonal - Restore the diagonal of a `MATDIAGONAL`

  Input Parameters:
+ A - the `MATDIAGONAL`
- diag - the `Vec` obtained from `MatDiagonalGetDiagonal()`

  Level: developer

  Note:
  Use `MatDiagonalSet()` to change the values by copy, rather than reference.

.seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalGetDiagonal()`
@*/
PetscErrorCode MatDiagonalRestoreDiagonal(Mat A, Vec *diag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(A, MAT_CLASSID, 1);
  PetscValidPointer(diag, 2);
  PetscUseMethod((PetscObject)A, "MatDiagonalRestoreDiagonal_C", (Mat, Vec *), (A, diag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalRestoreDiagonal_Diagonal(Mat A, Vec *diag)
{
  Mat_Diagonal    *ctx = (Mat_Diagonal *)A->data;
  PetscObjectState diag_state;

  PetscFunctionBegin;
  PetscCheck(ctx->diag == *diag, PetscObjectComm((PetscObject)A), PETSC_ERR_ARG_WRONG, "Restored a different diagonal vector");
  ctx->diag_valid = PETSC_TRUE;
  PetscCall(PetscObjectStateGet((PetscObject)*diag, &diag_state));
  if (ctx->diag_state != diag_state) {
    PetscCall(PetscObjectStateIncrease((PetscObject)A));
    ctx->inv_diag_valid = PETSC_FALSE;
  }
  *diag = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
  MatDiagonalGetInverseDiagonal - Get the inverse diagonal of a `MATDIAGONAL`

  Input Parameter:
. A - the `MATDIAGONAL`

  Output Parameter:
. inv_diag - the `Vec` that defines the inverse diagonal

  Level: developer

  Note:
   The user must call
  `MatDiagonalRestoreInverseDiagonal()` before using the matrix again.

  If a matrix is created only to call `MatSolve()` (which happens for `MATLMVMDIAGBROYDEN`),
  using `MatDiagonalGetInverseDiagonal()` and `MatDiagonalRestoreInverseDiagonal()` avoids copies
  and avoids any call to `VecReciprocal()`.

.seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalRestoreInverseDiagonal()`, `MatDiagonalGetDiagonal()`, `MATLMVMBROYDEN`, `MatSolve()`
@*/
PetscErrorCode MatDiagonalGetInverseDiagonal(Mat A, Vec *inv_diag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(A, MAT_CLASSID, 1);
  PetscValidPointer(inv_diag, 2);
  *inv_diag = NULL;
  PetscUseMethod((PetscObject)A, "MatDiagonalGetInverseDiagonal_C", (Mat, Vec *), (A, inv_diag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalGetInverseDiagonal_Diagonal(Mat A, Vec *inv_diag)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  PetscCall(MatSetUp(A));
  PetscCall(MatDiagonalSetUpInverseDiagonal(A));
  *inv_diag = ctx->inv_diag;
  PetscCall(PetscObjectStateGet((PetscObject)*inv_diag, &ctx->inv_diag_state));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
  MatDiagonalRestoreInverseDiagonal - Restore the inverse diagonal of a `MATDIAGONAL`

  Input Parameters:
+ A - the `MATDIAGONAL`
- inv_diag - the `Vec` obtained from `MatDiagonalGetInverseDiagonal()`

  Level: developer

.seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalGetInverseDiagonal()`
@*/
PetscErrorCode MatDiagonalRestoreInverseDiagonal(Mat A, Vec *inv_diag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(A, MAT_CLASSID, 1);
  PetscValidPointer(inv_diag, 2);
  PetscUseMethod((PetscObject)A, "MatDiagonalRestoreInverseDiagonal_C", (Mat, Vec *), (A, inv_diag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalRestoreInverseDiagonal_Diagonal(Mat A, Vec *inv_diag)
{
  Mat_Diagonal    *ctx = (Mat_Diagonal *)A->data;
  PetscObjectState inv_diag_state;

  PetscFunctionBegin;
  PetscCheck(ctx->inv_diag == *inv_diag, PetscObjectComm((PetscObject)A), PETSC_ERR_ARG_WRONG, "Restored a different diagonal vector");
  ctx->inv_diag_valid = PETSC_TRUE;
  PetscCall(PetscObjectStateGet((PetscObject)*inv_diag, &inv_diag_state));
  if (ctx->inv_diag_state != inv_diag_state) {
    PetscCall(PetscObjectStateIncrease((PetscObject)A));
    ctx->diag_valid = PETSC_FALSE;
  }
  *inv_diag = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_Diagonal(Mat mat)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)mat->data;

  PetscFunctionBegin;
  PetscCall(VecDestroy(&ctx->diag));
  PetscCall(VecDestroy(&ctx->inv_diag));
  PetscCall(PetscObjectComposeFunction((PetscObject)mat, "MatDiagonalGetDiagonal_C", NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)mat, "MatDiagonalRestoreDiagonal_C", NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)mat, "MatDiagonalGetInverseDiagonal_C", NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)mat, "MatDiagonalRestoreInverseDiagonal_C", NULL));
  PetscCall(PetscFree(mat->data));
  mat->structural_symmetry_eternal = PETSC_FALSE;
  mat->symmetry_eternal            = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatView_Diagonal(Mat J, PetscViewer viewer)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)J->data;
  PetscBool     iascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &iascii));
  if (iascii) {
    PetscViewerFormat format;

    PetscCall(PetscViewerGetFormat(viewer, &format));
    if (format == PETSC_VIEWER_ASCII_FACTOR_INFO || format == PETSC_VIEWER_ASCII_INFO) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(MatDiagonalSetUpDiagonal(J));
    PetscCall(VecView(ctx->diag, viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatGetDiagonal_Diagonal(Mat J, Vec x)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)J->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(J));
  PetscCall(VecCopy(ctx->diag, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalSet_Diagonal(Mat J, Vec D, InsertMode is)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)J->data;

  PetscFunctionBegin;
  switch (is) {
  case ADD_VALUES:
  case ADD_ALL_VALUES:
  case ADD_BC_VALUES:
    PetscCall(MatDiagonalSetUpDiagonal(J));
    PetscCall(VecAXPY(ctx->diag, 1.0, D));
    ctx->inv_diag_valid = PETSC_FALSE;
    break;
  case INSERT_VALUES:
  case INSERT_BC_VALUES:
  case INSERT_ALL_VALUES:
    PetscCall(MatSetUp(J));
    PetscCall(VecCopy(D, ctx->diag));
    ctx->diag_valid     = PETSC_TRUE;
    ctx->inv_diag_valid = PETSC_FALSE;
    break;
  case MAX_VALUES:
    PetscCall(MatDiagonalSetUpDiagonal(J));
    PetscCall(VecPointwiseMax(ctx->diag, D, ctx->diag));
    ctx->inv_diag_valid = PETSC_FALSE;
    break;
  case MIN_VALUES:
    PetscCall(MatDiagonalSetUpDiagonal(J));
    PetscCall(VecPointwiseMin(ctx->diag, D, ctx->diag));
    ctx->inv_diag_valid = PETSC_FALSE;
    break;
  case NOT_SET_VALUES:
    break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatShift_Diagonal(Mat Y, PetscScalar a)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)Y->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(Y));
  PetscCall(VecShift(ctx->diag, a));
  ctx->inv_diag_valid = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatScale_Diagonal(Mat Y, PetscScalar a)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)Y->data;

  PetscFunctionBegin;
  PetscCall(MatSetUp(Y));
  PetscCall(MatDiagonalSetUpDiagonal(Y));
  PetscCall(VecScale(ctx->diag, a));
  ctx->inv_diag_valid = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDiagonalScale_Diagonal(Mat Y, Vec l, Vec r)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)Y->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpDiagonal(Y));
  if (l) {
    PetscCall(VecPointwiseMult(ctx->diag, ctx->diag, l));
    ctx->inv_diag_valid = PETSC_FALSE;
  }
  if (r) {
    PetscCall(VecPointwiseMult(ctx->diag, ctx->diag, r));
    ctx->inv_diag_valid = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatSetUp_Diagonal(Mat A)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)A->data;

  PetscFunctionBegin;
  if (!ctx->diag) {
    PetscCall(PetscLayoutSetUp(A->cmap));
    PetscCall(PetscLayoutSetUp(A->rmap));
    PetscCall(MatCreateVecs(A, &ctx->diag, NULL));
    PetscCall(VecDuplicate(ctx->diag, &ctx->inv_diag));
    PetscCall(VecZeroEntries(ctx->diag));
    ctx->diag_valid = PETSC_TRUE;
  }
  A->assembled = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatZeroEntries_Diagonal(Mat Y)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)Y->data;

  PetscFunctionBegin;
  PetscCall(MatSetUp(Y));
  PetscCall(VecZeroEntries(ctx->diag));
  ctx->diag_valid     = PETSC_TRUE;
  ctx->inv_diag_valid = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatSolve_Diagonal(Mat matin, Vec b, Vec x)
{
  Mat_Diagonal *ctx = (Mat_Diagonal *)matin->data;

  PetscFunctionBegin;
  PetscCall(MatDiagonalSetUpInverseDiagonal(matin));
  PetscCall(VecPointwiseMult(x, b, ctx->inv_diag));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetInfo_Diagonal(Mat A, MatInfoType flag, MatInfo *info)
{
  PetscFunctionBegin;
  info->block_size        = 1.0;
  info->nz_allocated      = A->cmap->N;
  info->nz_used           = A->cmap->N;
  info->nz_unneeded       = 0.0;
  info->assemblies        = A->num_ass;
  info->mallocs           = 0.0;
  info->memory            = 0;
  info->fill_ratio_given  = 0;
  info->fill_ratio_needed = 0;
  info->factor_mallocs    = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MatCreateDiagonal - Creates a matrix defined by a given vector along its diagonal.

   Collective

   Input Parameter:
.  diag - vector for the diagonal

   Output Parameter:
.  J - the diagonal matrix

   Level: advanced

   Notes:
    Only supports square matrices with the same number of local rows and columns.

    The input vector `diag` will be referenced internally: any changes to `diag`
    will affect the matrix `J`.

.seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MATCONSTANTDIAGONAL`, `MatScale()`, `MatShift()`, `MatMult()`, `MatGetDiagonal()`, `MatSolve()`
          `MatDiagonalRestoreInverseDiagonal()`, `MatDiagonalGetDiagonal()`, `MatDiagonalRestoreDiagonal()`, `MatDiagonalGetInverseDiagonal()`
@*/
PetscErrorCode MatCreateDiagonal(Vec diag, Mat *J)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(diag, VEC_CLASSID, 1);
  PetscCall(MatCreate(PetscObjectComm((PetscObject)diag), J));
  PetscInt m, M;
  PetscCall(VecGetLocalSize(diag, &m));
  PetscCall(VecGetSize(diag, &M));
  PetscCall(MatSetSizes(*J, m, m, M, M));
  PetscCall(MatSetType(*J, MATDIAGONAL));

  PetscLayout map;
  PetscCall(VecGetLayout(diag, &map));
  PetscCall(MatSetLayouts(*J, map, map));
  Mat_Diagonal *ctx = (Mat_Diagonal *)(*J)->data;
  PetscCall(PetscObjectReference((PetscObject)diag));
  PetscCall(VecDestroy(&ctx->diag));
  PetscCall(VecDestroy(&ctx->inv_diag));
  ctx->diag           = diag;
  ctx->diag_valid     = PETSC_TRUE;
  ctx->inv_diag_valid = PETSC_FALSE;
  VecType type;
  PetscCall(VecDuplicate(ctx->diag, &ctx->inv_diag));
  PetscCall(VecGetType(diag, &type));
  PetscCall(PetscFree((*J)->defaultvectype));
  PetscCall(PetscStrallocpy(type, &(*J)->defaultvectype));
  PetscCall(MatSetUp(*J));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   MATDIAGONAL - MATDIAGONAL = "diagonal" - A diagonal matrix type with the diagonal implemented as a `Vec`.  Useful for
   cases where `VecPointwiseMult()` or `VecPointwiseDivide()` should be thought of as the actions of a linear operator.

  Level: advanced

.seealso: [](ch_matrices), `Mat`, `MatCreateDiagonal()`, `MatDiagonalRestoreInverseDiagonal()`, `MatDiagonalGetDiagonal()`, `MatDiagonalRestoreDiagonal()`, `MatDiagonalGetInverseDiagonal()`
M*/
PETSC_INTERN PetscErrorCode MatCreate_Diagonal(Mat A)
{
  Mat_Diagonal *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  A->data = (void *)ctx;

  A->structurally_symmetric      = PETSC_BOOL3_TRUE;
  A->structural_symmetry_eternal = PETSC_TRUE;
  A->symmetry_eternal            = PETSC_TRUE;
  A->symmetric                   = PETSC_BOOL3_TRUE;
  if (!PetscDefined(USE_COMPLEX)) A->hermitian = PETSC_BOOL3_TRUE;

  A->ops->mult             = MatMult_Diagonal;
  A->ops->multadd          = MatMultAdd_Diagonal;
  A->ops->multtranspose    = MatMult_Diagonal;
  A->ops->multtransposeadd = MatMultAdd_Diagonal;
  A->ops->norm             = MatNorm_Diagonal;
  A->ops->duplicate        = MatDuplicate_Diagonal;
  A->ops->solve            = MatSolve_Diagonal;
  A->ops->solvetranspose   = MatSolve_Diagonal;
  A->ops->shift            = MatShift_Diagonal;
  A->ops->scale            = MatScale_Diagonal;
  A->ops->diagonalscale    = MatDiagonalScale_Diagonal;
  A->ops->getdiagonal      = MatGetDiagonal_Diagonal;
  A->ops->diagonalset      = MatDiagonalSet_Diagonal;
  A->ops->view             = MatView_Diagonal;
  A->ops->zeroentries      = MatZeroEntries_Diagonal;
  A->ops->destroy          = MatDestroy_Diagonal;
  A->ops->getinfo          = MatGetInfo_Diagonal;
  A->ops->axpy             = MatAXPY_Diagonal;
  A->ops->setup            = MatSetUp_Diagonal;

  PetscCall(PetscObjectComposeFunction((PetscObject)A, "MatDiagonalGetDiagonal_C", MatDiagonalGetDiagonal_Diagonal));
  PetscCall(PetscObjectComposeFunction((PetscObject)A, "MatDiagonalRestoreDiagonal_C", MatDiagonalRestoreDiagonal_Diagonal));
  PetscCall(PetscObjectComposeFunction((PetscObject)A, "MatDiagonalGetInverseDiagonal_C", MatDiagonalGetInverseDiagonal_Diagonal));
  PetscCall(PetscObjectComposeFunction((PetscObject)A, "MatDiagonalRestoreInverseDiagonal_C", MatDiagonalRestoreInverseDiagonal_Diagonal));
  PetscCall(PetscObjectChangeTypeName((PetscObject)A, MATDIAGONAL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
