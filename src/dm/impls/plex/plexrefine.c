#include <petsc-private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include <petscsf.h>

#undef __FUNCT__
#define __FUNCT__ "GetDepthStart_Private"
PETSC_STATIC_INLINE PetscErrorCode GetDepthStart_Private(PetscInt depth, PetscInt depthSize[], PetscInt *cStart, PetscInt *fStart, PetscInt *eStart, PetscInt *vStart)
{
  PetscFunctionBegin;
  if (cStart) *cStart = 0;
  if (vStart) *vStart = depthSize[depth];
  if (fStart) *fStart = depthSize[depth] + depthSize[0];
  if (eStart) *eStart = depthSize[depth] + depthSize[0] + depthSize[depth-1];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetDepthEnd_Private"
PETSC_STATIC_INLINE PetscErrorCode GetDepthEnd_Private(PetscInt depth, PetscInt depthSize[], PetscInt *cEnd, PetscInt *fEnd, PetscInt *eEnd, PetscInt *vEnd)
{
  PetscFunctionBegin;
  if (cEnd) *cEnd = depthSize[depth];
  if (vEnd) *vEnd = depthSize[depth] + depthSize[0];
  if (fEnd) *fEnd = depthSize[depth] + depthSize[0] + depthSize[depth-1];
  if (eEnd) *eEnd = depthSize[depth] + depthSize[0] + depthSize[depth-1] + depthSize[1];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellRefinerGetSizes"
PetscErrorCode CellRefinerGetSizes(CellRefiner refiner, DM dm, PetscInt depthSize[])
{
  PetscInt       cStart, cEnd, cMax, vStart, vEnd, vMax, fStart, fEnd, fMax, eStart, eEnd, eMax;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cMax, &fMax, &eMax, &vMax);CHKERRQ(ierr);
  switch (refiner) {
  case 1:
    /* Simplicial 2D */
    depthSize[0] = vEnd - vStart + fEnd - fStart;         /* Add a vertex on every face */
    depthSize[1] = 2*(fEnd - fStart) + 3*(cEnd - cStart); /* Every face is split into 2 faces and 3 faces are added for each cell */
    depthSize[2] = 4*(cEnd - cStart);                     /* Every cell split into 4 cells */
    break;
  case 3:
    /* Hybrid 2D */
    if (cMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No cell maximum specified in hybrid mesh");
    cMax = PetscMin(cEnd, cMax);
    if (fMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No face maximum specified in hybrid mesh");
    fMax         = PetscMin(fEnd, fMax);
    depthSize[0] = vEnd - vStart + fMax - fStart;                                         /* Add a vertex on every face, but not hybrid faces */
    depthSize[1] = 2*(fMax - fStart) + 3*(cMax - cStart) + (fEnd - fMax) + (cEnd - cMax); /* Every interior face is split into 2 faces, 3 faces are added for each interior cell, and one in each hybrid cell */
    depthSize[2] = 4*(cMax - cStart) + 2*(cEnd - cMax);                                   /* Interior cells split into 4 cells, Hybrid cells split into 2 cells */
    break;
  case 2:
    /* Hex 2D */
    depthSize[0] = vEnd - vStart + cEnd - cStart + fEnd - fStart; /* Add a vertex on every face and cell */
    depthSize[1] = 2*(fEnd - fStart) + 4*(cEnd - cStart);         /* Every face is split into 2 faces and 4 faces are added for each cell */
    depthSize[2] = 4*(cEnd - cStart);                             /* Every cell split into 4 cells */
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellRefinerSetConeSizes"
PetscErrorCode CellRefinerSetConeSizes(CellRefiner refiner, DM dm, PetscInt depthSize[], DM rdm)
{
  PetscInt       depth, cStart, cStartNew, cEnd, cMax, c, vStart, vStartNew, vEnd, vMax, v, fStart, fStartNew, fEnd, fMax, f, eStart, eStartNew, eEnd, eMax, r;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cMax, &fMax, &eMax, &vMax);CHKERRQ(ierr);
  ierr = GetDepthStart_Private(depth, depthSize, &cStartNew, &fStartNew, &eStartNew, &vStartNew);CHKERRQ(ierr);
  switch (refiner) {
  case 1:
    /* Simplicial 2D */
    /* All cells have 3 faces */
    for (c = cStart; c < cEnd; ++c) {
      for (r = 0; r < 4; ++r) {
        const PetscInt newp = (c - cStart)*4 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 3);CHKERRQ(ierr);
      }
    }
    /* Split faces have 2 vertices and the same cells as the parent */
    for (f = fStart; f < fEnd; ++f) {
      for (r = 0; r < 2; ++r) {
        const PetscInt newp = fStartNew + (f - fStart)*2 + r;
        PetscInt       size;

        ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
        ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
        ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
      }
    }
    /* Interior faces have 2 vertices and 2 cells */
    for (c = cStart; c < cEnd; ++c) {
      for (r = 0; r < 3; ++r) {
        const PetscInt newp = fStartNew + (fEnd - fStart)*2 + (c - cStart)*3 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
        ierr = DMPlexSetSupportSize(rdm, newp, 2);CHKERRQ(ierr);
      }
    }
    /* Old vertices have identical supports */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt newp = vStartNew + (v - vStart);
      PetscInt       size;

      ierr = DMPlexGetSupportSize(dm, v, &size);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
    }
    /* Face vertices have 2 + cells*2 supports */
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt newp = vStartNew + (vEnd - vStart) + (f - fStart);
      PetscInt       size;

      ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, 2 + size*2);CHKERRQ(ierr);
    }
    break;
  case 2:
    /* Hex 2D */
    /* All cells have 4 faces */
    for (c = cStart; c < cEnd; ++c) {
      for (r = 0; r < 4; ++r) {
        const PetscInt newp = (c - cStart)*4 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 4);CHKERRQ(ierr);
      }
    }
    /* Split faces have 2 vertices and the same cells as the parent */
    for (f = fStart; f < fEnd; ++f) {
      for (r = 0; r < 2; ++r) {
        const PetscInt newp = fStartNew + (f - fStart)*2 + r;
        PetscInt       size;

        ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
        ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
        ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
      }
    }
    /* Interior faces have 2 vertices and 2 cells */
    for (c = cStart; c < cEnd; ++c) {
      for (r = 0; r < 4; ++r) {
        const PetscInt newp = fStartNew + (fEnd - fStart)*2 + (c - cStart)*4 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
        ierr = DMPlexSetSupportSize(rdm, newp, 2);CHKERRQ(ierr);
      }
    }
    /* Old vertices have identical supports */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt newp = vStartNew + (v - vStart);
      PetscInt       size;

      ierr = DMPlexGetSupportSize(dm, v, &size);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
    }
    /* Face vertices have 2 + cells supports */
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt newp = vStartNew + (vEnd - vStart) + (f - fStart);
      PetscInt       size;

      ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, 2 + size);CHKERRQ(ierr);
    }
    /* Cell vertices have 4 supports */
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt newp = vStartNew + (vEnd - vStart) + (fEnd - fStart) + (c - cStart);

      ierr = DMPlexSetSupportSize(rdm, newp, 4);CHKERRQ(ierr);
    }
    break;
  case 3:
    /* Hybrid 2D */
    if (cMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No cell maximum specified in hybrid mesh");
    cMax = PetscMin(cEnd, cMax);
    if (fMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No face maximum specified in hybrid mesh");
    fMax = PetscMin(fEnd, fMax);
    ierr = DMPlexSetHybridBounds(rdm, cStartNew + (cMax - cStart)*4, fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3, PETSC_DETERMINE, PETSC_DETERMINE);CHKERRQ(ierr);
    /* Interior cells have 3 faces */
    for (c = cStart; c < cMax; ++c) {
      for (r = 0; r < 4; ++r) {
        const PetscInt newp = cStartNew + (c - cStart)*4 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 3);CHKERRQ(ierr);
      }
    }
    /* Hybrid cells have 4 faces */
    for (c = cMax; c < cEnd; ++c) {
      for (r = 0; r < 2; ++r) {
        const PetscInt newp = cStartNew + (cMax - cStart)*4 + (c - cMax)*2 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 4);CHKERRQ(ierr);
      }
    }
    /* Interior split faces have 2 vertices and the same cells as the parent */
    for (f = fStart; f < fMax; ++f) {
      for (r = 0; r < 2; ++r) {
        const PetscInt newp = fStartNew + (f - fStart)*2 + r;
        PetscInt       size;

        ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
        ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
        ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
      }
    }
    /* Interior cell faces have 2 vertices and 2 cells */
    for (c = cStart; c < cMax; ++c) {
      for (r = 0; r < 3; ++r) {
        const PetscInt newp = fStartNew + (fMax - fStart)*2 + (c - cStart)*3 + r;

        ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
        ierr = DMPlexSetSupportSize(rdm, newp, 2);CHKERRQ(ierr);
      }
    }
    /* Hybrid faces have 2 vertices and the same cells */
    for (f = fMax; f < fEnd; ++f) {
      const PetscInt newp = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (f - fMax);
      PetscInt       size;

      ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
      ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
    }
    /* Hybrid cell faces have 2 vertices and 2 cells */
    for (c = cMax; c < cEnd; ++c) {
      const PetscInt newp = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (fEnd - fMax) + (c - cMax);

      ierr = DMPlexSetConeSize(rdm, newp, 2);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, 2);CHKERRQ(ierr);
    }
    /* Old vertices have identical supports */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt newp = vStartNew + (v - vStart);
      PetscInt       size;

      ierr = DMPlexGetSupportSize(dm, v, &size);CHKERRQ(ierr);
      ierr = DMPlexSetSupportSize(rdm, newp, size);CHKERRQ(ierr);
    }
    /* Face vertices have 2 + (2 interior, 1 hybrid) supports */
    for (f = fStart; f < fMax; ++f) {
      const PetscInt newp = vStartNew + (vEnd - vStart) + (f - fStart);
      const PetscInt *support;
      PetscInt       size, newSize = 2, s;

      ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
      for (s = 0; s < size; ++s) {
        if (support[s] >= cMax) newSize += 1;
        else newSize += 2;
      }
      ierr = DMPlexSetSupportSize(rdm, newp, newSize);CHKERRQ(ierr);
    }
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellRefinerSetCones"
PetscErrorCode CellRefinerSetCones(CellRefiner refiner, DM dm, PetscInt depthSize[], DM rdm)
{
  PetscInt       depth, cStart, cEnd, cMax, cStartNew, cEndNew, c, vStart, vEnd, vMax, vStartNew, vEndNew, v, fStart, fEnd, fMax, fStartNew, fEndNew, f, eStart, eEnd, eMax, eStartNew, eEndNew, r, p;
  PetscInt       maxSupportSize, *supportRef;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cMax, &fMax, &eMax, &vMax);CHKERRQ(ierr);
  ierr = GetDepthStart_Private(depth, depthSize, &cStartNew, &fStartNew, &eStartNew, &vStartNew);CHKERRQ(ierr);
  ierr = GetDepthEnd_Private(depth, depthSize, &cEndNew, &fEndNew, &eEndNew, &vEndNew);CHKERRQ(ierr);
  switch (refiner) {
  case 1:
    /* Simplicial 2D */
    /*
     2
     |\
     | \
     |  \
     |   \
     | C  \
     |     \
     |      \
     2---1---1
     |\  D  / \
     | 2   0   \
     |A \ /  B  \
     0---0-------1
     */
    /* All cells have 3 faces */
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt  newp = cStartNew + (c - cStart)*4;
      const PetscInt *cone, *ornt;
      PetscInt        coneNew[3], orntNew[3];

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      ierr = DMPlexGetConeOrientation(dm, c, &ornt);CHKERRQ(ierr);
      /* A triangle */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 1 : 0);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*3 + 2;
      orntNew[1] = -2;
      coneNew[2] = fStartNew + (cone[2] - fStart)*2 + (ornt[2] < 0 ? 0 : 1);
      orntNew[2] = ornt[2];
      ierr       = DMPlexSetCone(rdm, newp+0, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+0, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+0 < cStartNew) || (newp+0 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+0, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* B triangle */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 0 : 1);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 1 : 0);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*3 + 0;
      orntNew[2] = -2;
      ierr       = DMPlexSetCone(rdm, newp+1, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+1, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+1 < cStartNew) || (newp+1 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+1, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* C triangle */
      coneNew[0] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*3 + 1;
      orntNew[0] = -2;
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 0 : 1);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (cone[2] - fStart)*2 + (ornt[2] < 0 ? 1 : 0);
      orntNew[2] = ornt[2];
      ierr       = DMPlexSetCone(rdm, newp+2, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+2, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+2 < cStartNew) || (newp+2 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+2, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* D triangle */
      coneNew[0] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*3 + 0;
      orntNew[0] = 0;
      coneNew[1] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*3 + 1;
      orntNew[1] = 0;
      coneNew[2] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*3 + 2;
      orntNew[2] = 0;
      ierr       = DMPlexSetCone(rdm, newp+3, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+3, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+3 < cStartNew) || (newp+3 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+3, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Split faces have 2 vertices and the same cells as the parent */
    ierr = DMPlexGetMaxSizes(dm, NULL, &maxSupportSize);CHKERRQ(ierr);
    ierr = PetscMalloc((2 + maxSupportSize*2) * sizeof(PetscInt), &supportRef);CHKERRQ(ierr);
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt newv = vStartNew + (vEnd - vStart) + (f - fStart);

      for (r = 0; r < 2; ++r) {
        const PetscInt  newp = fStartNew + (f - fStart)*2 + r;
        const PetscInt *cone, *support;
        PetscInt        coneNew[2], coneSize, c, supportSize, s;

        ierr             = DMPlexGetCone(dm, f, &cone);CHKERRQ(ierr);
        coneNew[0]       = vStartNew + (cone[0] - vStart);
        coneNew[1]       = vStartNew + (cone[1] - vStart);
        coneNew[(r+1)%2] = newv;
        ierr             = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
        }
#endif
        ierr = DMPlexGetSupportSize(dm, f, &supportSize);CHKERRQ(ierr);
        ierr = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
        for (s = 0; s < supportSize; ++s) {
          ierr = DMPlexGetConeSize(dm, support[s], &coneSize);CHKERRQ(ierr);
          ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
          for (c = 0; c < coneSize; ++c) {
            if (cone[c] == f) break;
          }
          supportRef[s] = cStartNew + (support[s] - cStart)*4 + (c+r)%3;
        }
        ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < supportSize; ++p) {
          if ((supportRef[p] < cStartNew) || (supportRef[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportRef[p], cStartNew, cEndNew);
        }
#endif
      }
    }
    /* Interior faces have 2 vertices and 2 cells */
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt *cone;

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      for (r = 0; r < 3; ++r) {
        const PetscInt newp = fStartNew + (fEnd - fStart)*2 + (c - cStart)*3 + r;
        PetscInt       coneNew[2];
        PetscInt       supportNew[2];

        coneNew[0] = vStartNew + (vEnd - vStart) + (cone[r]       - fStart);
        coneNew[1] = vStartNew + (vEnd - vStart) + (cone[(r+1)%3] - fStart);
        ierr       = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
        }
#endif
        supportNew[0] = (c - cStart)*4 + (r+1)%3;
        supportNew[1] = (c - cStart)*4 + 3;
        ierr          = DMPlexSetSupport(rdm, newp, supportNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((supportNew[p] < cStartNew) || (supportNew[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportNew[p], cStartNew, cEndNew);
        }
#endif
      }
    }
    /* Old vertices have identical supports */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt  newp = vStartNew + (v - vStart);
      const PetscInt *support, *cone;
      PetscInt        size, s;

      ierr = DMPlexGetSupportSize(dm, v, &size);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, v, &support);CHKERRQ(ierr);
      for (s = 0; s < size; ++s) {
        PetscInt r = 0;

        ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
        if (cone[1] == v) r = 1;
        supportRef[s] = fStartNew + (support[s] - fStart)*2 + r;
      }
      ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
      if ((newp < vStartNew) || (newp >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", newp, vStartNew, vEndNew);
      for (p = 0; p < size; ++p) {
        if ((supportRef[p] < fStartNew) || (supportRef[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", supportRef[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Face vertices have 2 + cells*2 supports */
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt  newp = vStartNew + (vEnd - vStart) + (f - fStart);
      const PetscInt *cone, *support;
      PetscInt        size, s;

      ierr          = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr          = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
      supportRef[0] = fStartNew + (f - fStart)*2 + 0;
      supportRef[1] = fStartNew + (f - fStart)*2 + 1;
      for (s = 0; s < size; ++s) {
        PetscInt r = 0;

        ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
        if      (cone[1] == f) r = 1;
        else if (cone[2] == f) r = 2;
        supportRef[2+s*2+0] = fStartNew + (fEnd - fStart)*2 + (support[s] - cStart)*3 + (r+2)%3;
        supportRef[2+s*2+1] = fStartNew + (fEnd - fStart)*2 + (support[s] - cStart)*3 + r;
      }
      ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
      if ((newp < vStartNew) || (newp >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", newp, vStartNew, vEndNew);
      for (p = 0; p < 2+size*2; ++p) {
        if ((supportRef[p] < fStartNew) || (supportRef[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", supportRef[p], fStartNew, fEndNew);
      }
#endif
    }
    ierr = PetscFree(supportRef);CHKERRQ(ierr);
    break;
  case 2:
    /* Hex 2D */
    /*
     3---------2---------2
     |         |         |
     |    D    2    C    |
     |         |         |
     3----3----0----1----1
     |         |         |
     |    A    0    B    |
     |         |         |
     0---------0---------1
     */
    /* All cells have 4 faces */
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt  newp = (c - cStart)*4;
      const PetscInt *cone, *ornt;
      PetscInt        coneNew[4], orntNew[4];

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      ierr = DMPlexGetConeOrientation(dm, c, &ornt);CHKERRQ(ierr);
      /* A quad */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 1 : 0);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 0;
      orntNew[1] = 0;
      coneNew[2] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 3;
      orntNew[2] = -2;
      coneNew[3] = fStartNew + (cone[3] - fStart)*2 + (ornt[3] < 0 ? 0 : 1);
      orntNew[3] = ornt[3];
      ierr       = DMPlexSetCone(rdm, newp+0, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+0, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+0 < cStartNew) || (newp+0 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+0, cStartNew, cEndNew);
      for (p = 0; p < 4; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* B quad */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 0 : 1);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 1 : 0);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 1;
      orntNew[2] = 0;
      coneNew[3] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 0;
      orntNew[3] = -2;
      ierr       = DMPlexSetCone(rdm, newp+1, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+1, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+1 < cStartNew) || (newp+1 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+1, cStartNew, cEndNew);
      for (p = 0; p < 4; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* C quad */
      coneNew[0] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 1;
      orntNew[0] = -2;
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 0 : 1);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (cone[2] - fStart)*2 + (ornt[2] < 0 ? 1 : 0);
      orntNew[2] = ornt[2];
      coneNew[3] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 2;
      orntNew[3] = 0;
      ierr       = DMPlexSetCone(rdm, newp+2, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+2, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+2 < cStartNew) || (newp+2 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+2, cStartNew, cEndNew);
      for (p = 0; p < 4; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* D quad */
      coneNew[0] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 3;
      orntNew[0] = 0;
      coneNew[1] = fStartNew + (fEnd    - fStart)*2 + (c - cStart)*4 + 2;
      orntNew[1] = -2;
      coneNew[2] = fStartNew + (cone[2] - fStart)*2 + (ornt[2] < 0 ? 0 : 1);
      orntNew[2] = ornt[2];
      coneNew[3] = fStartNew + (cone[3] - fStart)*2 + (ornt[3] < 0 ? 1 : 0);
      orntNew[3] = ornt[3];
      ierr       = DMPlexSetCone(rdm, newp+3, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+3, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+3 < cStartNew) || (newp+3 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+3, cStartNew, cEndNew);
      for (p = 0; p < 4; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Split faces have 2 vertices and the same cells as the parent */
    ierr = DMPlexGetMaxSizes(dm, NULL, &maxSupportSize);CHKERRQ(ierr);
    ierr = PetscMalloc((2 + maxSupportSize*2) * sizeof(PetscInt), &supportRef);CHKERRQ(ierr);
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt newv = vStartNew + (vEnd - vStart) + (f - fStart);

      for (r = 0; r < 2; ++r) {
        const PetscInt  newp = fStartNew + (f - fStart)*2 + r;
        const PetscInt *cone, *support;
        PetscInt        coneNew[2], coneSize, c, supportSize, s;

        ierr             = DMPlexGetCone(dm, f, &cone);CHKERRQ(ierr);
        coneNew[0]       = vStartNew + (cone[0] - vStart);
        coneNew[1]       = vStartNew + (cone[1] - vStart);
        coneNew[(r+1)%2] = newv;
        ierr             = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
        }
#endif
        ierr = DMPlexGetSupportSize(dm, f, &supportSize);CHKERRQ(ierr);
        ierr = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
        for (s = 0; s < supportSize; ++s) {
          ierr = DMPlexGetConeSize(dm, support[s], &coneSize);CHKERRQ(ierr);
          ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
          for (c = 0; c < coneSize; ++c) {
            if (cone[c] == f) break;
          }
          supportRef[s] = cStartNew + (support[s] - cStart)*4 + (c+r)%4;
        }
        ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < supportSize; ++p) {
          if ((supportRef[p] < cStartNew) || (supportRef[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportRef[p], cStartNew, cEndNew);
        }
#endif
      }
    }
    /* Interior faces have 2 vertices and 2 cells */
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt *cone;
      PetscInt        coneNew[2], supportNew[2];

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      for (r = 0; r < 4; ++r) {
        const PetscInt newp = fStartNew + (fEnd - fStart)*2 + (c - cStart)*4 + r;

        coneNew[0] = vStartNew + (vEnd - vStart) + (cone[r] - fStart);
        coneNew[1] = vStartNew + (vEnd - vStart) + (fEnd    - fStart) + (c - cStart);
        ierr       = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
        }
#endif
        supportNew[0] = (c - cStart)*4 + r;
        supportNew[1] = (c - cStart)*4 + (r+1)%4;
        ierr          = DMPlexSetSupport(rdm, newp, supportNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((supportNew[p] < cStartNew) || (supportNew[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportNew[p], cStartNew, cEndNew);
        }
#endif
      }
    }
    /* Old vertices have identical supports */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt  newp = vStartNew + (v - vStart);
      const PetscInt *support, *cone;
      PetscInt        size, s;

      ierr = DMPlexGetSupportSize(dm, v, &size);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, v, &support);CHKERRQ(ierr);
      for (s = 0; s < size; ++s) {
        PetscInt r = 0;

        ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
        if (cone[1] == v) r = 1;
        supportRef[s] = fStartNew + (support[s] - fStart)*2 + r;
      }
      ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
      if ((newp < vStartNew) || (newp >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", newp, vStartNew, vEndNew);
      for (p = 0; p < size; ++p) {
        if ((supportRef[p] < fStartNew) || (supportRef[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", supportRef[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Face vertices have 2 + cells supports */
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt  newp = vStartNew + (vEnd - vStart) + (f - fStart);
      const PetscInt *cone, *support;
      PetscInt        size, s;

      ierr          = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr          = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
      supportRef[0] = fStartNew + (f - fStart)*2 + 0;
      supportRef[1] = fStartNew + (f - fStart)*2 + 1;
      for (s = 0; s < size; ++s) {
        PetscInt r = 0;

        ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
        if      (cone[1] == f) r = 1;
        else if (cone[2] == f) r = 2;
        else if (cone[3] == f) r = 3;
        supportRef[2+s] = fStartNew + (fEnd - fStart)*2 + (support[s] - cStart)*4 + r;
      }
      ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
      if ((newp < vStartNew) || (newp >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", newp, vStartNew, vEndNew);
      for (p = 0; p < 2+size; ++p) {
        if ((supportRef[p] < fStartNew) || (supportRef[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", supportRef[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Cell vertices have 4 supports */
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt newp = vStartNew + (vEnd - vStart) + (fEnd - fStart) + (c - cStart);
      PetscInt       supportNew[4];

      for (r = 0; r < 4; ++r) {
        supportNew[r] = fStartNew + (fEnd - fStart)*2 + (c - cStart)*4 + r;
      }
      ierr = DMPlexSetSupport(rdm, newp, supportNew);CHKERRQ(ierr);
    }
    break;
  case 3:
    if (cMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No cell maximum specified in hybrid mesh");
    cMax = PetscMin(cEnd, cMax);
    if (fMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No face maximum specified in hybrid mesh");
    fMax = PetscMin(fEnd, fMax);
    /* Interior cells have 3 faces */
    for (c = cStart; c < cMax; ++c) {
      const PetscInt  newp = cStartNew + (c - cStart)*4;
      const PetscInt *cone, *ornt;
      PetscInt        coneNew[3], orntNew[3];

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      ierr = DMPlexGetConeOrientation(dm, c, &ornt);CHKERRQ(ierr);
      /* A triangle */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 1 : 0);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (fMax    - fStart)*2 + (c - cStart)*3 + 2;
      orntNew[1] = -2;
      coneNew[2] = fStartNew + (cone[2] - fStart)*2 + (ornt[2] < 0 ? 0 : 1);
      orntNew[2] = ornt[2];
      ierr       = DMPlexSetCone(rdm, newp+0, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+0, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+0 < cStartNew) || (newp+0 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+0, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* B triangle */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 0 : 1);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 1 : 0);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (fMax    - fStart)*2 + (c - cStart)*3 + 0;
      orntNew[2] = -2;
      ierr       = DMPlexSetCone(rdm, newp+1, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+1, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+1 < cStartNew) || (newp+1 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+1, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* C triangle */
      coneNew[0] = fStartNew + (fMax    - fStart)*2 + (c - cStart)*3 + 1;
      orntNew[0] = -2;
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 0 : 1);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (cone[2] - fStart)*2 + (ornt[2] < 0 ? 1 : 0);
      orntNew[2] = ornt[2];
      ierr       = DMPlexSetCone(rdm, newp+2, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+2, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+2 < cStartNew) || (newp+2 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+2, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* D triangle */
      coneNew[0] = fStartNew + (fMax    - fStart)*2 + (c - cStart)*3 + 0;
      orntNew[0] = 0;
      coneNew[1] = fStartNew + (fMax    - fStart)*2 + (c - cStart)*3 + 1;
      orntNew[1] = 0;
      coneNew[2] = fStartNew + (fMax    - fStart)*2 + (c - cStart)*3 + 2;
      orntNew[2] = 0;
      ierr       = DMPlexSetCone(rdm, newp+3, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+3, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+3 < cStartNew) || (newp+3 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+3, cStartNew, cEndNew);
      for (p = 0; p < 3; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
    }
    /*
     2----3----3
     |         |
     |    B    |
     |         |
     0----4--- 1
     |         |
     |    A    |
     |         |
     0----2----1
     */
    /* Hybrid cells have 4 faces */
    for (c = cMax; c < cEnd; ++c) {
      const PetscInt  newp = cStartNew + (cMax - cStart)*4 + (c - cMax)*2;
      const PetscInt *cone, *ornt;
      PetscInt        coneNew[4], orntNew[4];

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      ierr = DMPlexGetConeOrientation(dm, c, &ornt);CHKERRQ(ierr);
      /* A quad */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 1 : 0);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 1 : 0);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (fMax    - fStart)*2 + (cMax - cStart)*3 + (cone[2] - fMax);
      orntNew[2] = 0;
      coneNew[3] = fStartNew + (fMax    - fStart)*2 + (cMax - cStart)*3 + (fEnd    - fMax) + (c - cMax);
      orntNew[3] = 0;
      ierr       = DMPlexSetCone(rdm, newp+0, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+0, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+0 < cStartNew) || (newp+0 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+0, cStartNew, cEndNew);
      for (p = 0; p < 4; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
      /* B quad */
      coneNew[0] = fStartNew + (cone[0] - fStart)*2 + (ornt[0] < 0 ? 0 : 1);
      orntNew[0] = ornt[0];
      coneNew[1] = fStartNew + (cone[1] - fStart)*2 + (ornt[1] < 0 ? 0 : 1);
      orntNew[1] = ornt[1];
      coneNew[2] = fStartNew + (fMax    - fStart)*2 + (cMax - cStart)*3 + (fEnd    - fMax) + (c - cMax);
      orntNew[2] = 0;
      coneNew[3] = fStartNew + (fMax    - fStart)*2 + (cMax - cStart)*3 + (cone[3] - fMax);
      orntNew[3] = 0;
      ierr       = DMPlexSetCone(rdm, newp+1, coneNew);CHKERRQ(ierr);
      ierr       = DMPlexSetConeOrientation(rdm, newp+1, orntNew);CHKERRQ(ierr);
#if 1
      if ((newp+1 < cStartNew) || (newp+1 >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", newp+1, cStartNew, cEndNew);
      for (p = 0; p < 4; ++p) {
        if ((coneNew[p] < fStartNew) || (coneNew[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", coneNew[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Interior split faces have 2 vertices and the same cells as the parent */
    ierr = DMPlexGetMaxSizes(dm, NULL, &maxSupportSize);CHKERRQ(ierr);
    ierr = PetscMalloc((2 + maxSupportSize*2) * sizeof(PetscInt), &supportRef);CHKERRQ(ierr);
    for (f = fStart; f < fMax; ++f) {
      const PetscInt newv = vStartNew + (vEnd - vStart) + (f - fStart);

      for (r = 0; r < 2; ++r) {
        const PetscInt  newp = fStartNew + (f - fStart)*2 + r;
        const PetscInt *cone, *support;
        PetscInt        coneNew[2], coneSize, c, supportSize, s;

        ierr             = DMPlexGetCone(dm, f, &cone);CHKERRQ(ierr);
        coneNew[0]       = vStartNew + (cone[0] - vStart);
        coneNew[1]       = vStartNew + (cone[1] - vStart);
        coneNew[(r+1)%2] = newv;
        ierr             = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
        }
#endif
        ierr = DMPlexGetSupportSize(dm, f, &supportSize);CHKERRQ(ierr);
        ierr = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
        for (s = 0; s < supportSize; ++s) {
          if (support[s] >= cMax) {
            supportRef[s] = cStartNew + (cMax - cStart)*4 + (support[s] - cMax)*2 + r;
          } else {
            ierr = DMPlexGetConeSize(dm, support[s], &coneSize);CHKERRQ(ierr);
            ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
            for (c = 0; c < coneSize; ++c) {
              if (cone[c] == f) break;
            }
            supportRef[s] = cStartNew + (support[s] - cStart)*4 + (c+r)%3;
          }
        }
        ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < supportSize; ++p) {
          if ((supportRef[p] < cStartNew) || (supportRef[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportRef[p], cStartNew, cEndNew);
        }
#endif
      }
    }
    /* Interior cell faces have 2 vertices and 2 cells */
    for (c = cStart; c < cMax; ++c) {
      const PetscInt *cone;

      ierr = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      for (r = 0; r < 3; ++r) {
        const PetscInt newp = fStartNew + (fMax - fStart)*2 + (c - cStart)*3 + r;
        PetscInt       coneNew[2];
        PetscInt       supportNew[2];

        coneNew[0] = vStartNew + (vEnd - vStart) + (cone[r]       - fStart);
        coneNew[1] = vStartNew + (vEnd - vStart) + (cone[(r+1)%3] - fStart);
        ierr       = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
        }
#endif
        supportNew[0] = (c - cStart)*4 + (r+1)%3;
        supportNew[1] = (c - cStart)*4 + 3;
        ierr          = DMPlexSetSupport(rdm, newp, supportNew);CHKERRQ(ierr);
#if 1
        if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
        for (p = 0; p < 2; ++p) {
          if ((supportNew[p] < cStartNew) || (supportNew[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportNew[p], cStartNew, cEndNew);
        }
#endif
      }
    }
    /* Interior hybrid faces have 2 vertices and the same cells */
    for (f = fMax; f < fEnd; ++f) {
      const PetscInt  newp = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (f - fMax);
      const PetscInt *cone;
      const PetscInt *support;
      PetscInt        coneNew[2];
      PetscInt        supportNew[2];
      PetscInt        size, s, r;

      ierr       = DMPlexGetCone(dm, f, &cone);CHKERRQ(ierr);
      coneNew[0] = vStartNew + (cone[0] - vStart);
      coneNew[1] = vStartNew + (cone[1] - vStart);
      ierr       = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
      if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
      for (p = 0; p < 2; ++p) {
        if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
      }
#endif
      ierr = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
      for (s = 0; s < size; ++s) {
        ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
        for (r = 0; r < 2; ++r) {
          if (cone[r+2] == f) break;
        }
        supportNew[s] = (cMax - cStart)*4 + (support[s] - cMax)*2 + r;
      }
      ierr = DMPlexSetSupport(rdm, newp, supportNew);CHKERRQ(ierr);
#if 1
      if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
      for (p = 0; p < size; ++p) {
        if ((supportNew[p] < cStartNew) || (supportNew[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportNew[p], cStartNew, cEndNew);
      }
#endif
    }
    /* Cell hybrid faces have 2 vertices and 2 cells */
    for (c = cMax; c < cEnd; ++c) {
      const PetscInt  newp = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (fEnd - fMax) + (c - cMax);
      const PetscInt *cone;
      PetscInt        coneNew[2];
      PetscInt        supportNew[2];

      ierr       = DMPlexGetCone(dm, c, &cone);CHKERRQ(ierr);
      coneNew[0] = vStartNew + (vEnd - vStart) + (cone[0] - fStart);
      coneNew[1] = vStartNew + (vEnd - vStart) + (cone[1] - fStart);
      ierr       = DMPlexSetCone(rdm, newp, coneNew);CHKERRQ(ierr);
#if 1
      if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
      for (p = 0; p < 2; ++p) {
        if ((coneNew[p] < vStartNew) || (coneNew[p] >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", coneNew[p], vStartNew, vEndNew);
      }
#endif
      supportNew[0] = (cMax - cStart)*4 + (c - cMax)*2 + 0;
      supportNew[1] = (cMax - cStart)*4 + (c - cMax)*2 + 1;
      ierr          = DMPlexSetSupport(rdm, newp, supportNew);CHKERRQ(ierr);
#if 1
      if ((newp < fStartNew) || (newp >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", newp, fStartNew, fEndNew);
      for (p = 0; p < 2; ++p) {
        if ((supportNew[p] < cStartNew) || (supportNew[p] >= cEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a cell [%d, %d)", supportNew[p], cStartNew, cEndNew);
      }
#endif
    }
    /* Old vertices have identical supports */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt  newp = vStartNew + (v - vStart);
      const PetscInt *support, *cone;
      PetscInt        size, s;

      ierr = DMPlexGetSupportSize(dm, v, &size);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, v, &support);CHKERRQ(ierr);
      for (s = 0; s < size; ++s) {
        if (support[s] >= fMax) {
          supportRef[s] = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (support[s] - fMax);
        } else {
          PetscInt r = 0;

          ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
          if (cone[1] == v) r = 1;
          supportRef[s] = fStartNew + (support[s] - fStart)*2 + r;
        }
      }
      ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
      if ((newp < vStartNew) || (newp >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", newp, vStartNew, vEndNew);
      for (p = 0; p < size; ++p) {
        if ((supportRef[p] < fStartNew) || (supportRef[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", supportRef[p], fStartNew, fEndNew);
      }
#endif
    }
    /* Face vertices have 2 + (2 interior, 1 hybrid) supports */
    for (f = fStart; f < fMax; ++f) {
      const PetscInt  newp = vStartNew + (vEnd - vStart) + (f - fStart);
      const PetscInt *cone, *support;
      PetscInt        size, newSize = 2, s;

      ierr          = DMPlexGetSupportSize(dm, f, &size);CHKERRQ(ierr);
      ierr          = DMPlexGetSupport(dm, f, &support);CHKERRQ(ierr);
      supportRef[0] = fStartNew + (f - fStart)*2 + 0;
      supportRef[1] = fStartNew + (f - fStart)*2 + 1;
      for (s = 0; s < size; ++s) {
        PetscInt r = 0;

        ierr = DMPlexGetCone(dm, support[s], &cone);CHKERRQ(ierr);
        if (support[s] >= cMax) {
          supportRef[newSize+0] = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (fEnd - fMax) + (support[s] - cMax);

          newSize += 1;
        } else {
          if      (cone[1] == f) r = 1;
          else if (cone[2] == f) r = 2;
          supportRef[newSize+0] = fStartNew + (fMax - fStart)*2 + (support[s] - cStart)*3 + (r+2)%3;
          supportRef[newSize+1] = fStartNew + (fMax - fStart)*2 + (support[s] - cStart)*3 + r;

          newSize += 2;
        }
      }
      ierr = DMPlexSetSupport(rdm, newp, supportRef);CHKERRQ(ierr);
#if 1
      if ((newp < vStartNew) || (newp >= vEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a vertex [%d, %d)", newp, vStartNew, vEndNew);
      for (p = 0; p < newSize; ++p) {
        if ((supportRef[p] < fStartNew) || (supportRef[p] >= fEndNew)) SETERRQ3(PetscObjectComm((PetscObject)dm), PETSC_ERR_PLIB, "Point %d is not a face [%d, %d)", supportRef[p], fStartNew, fEndNew);
      }
#endif
    }
    ierr = PetscFree(supportRef);CHKERRQ(ierr);
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellRefinerSetCoordinates"
PetscErrorCode CellRefinerSetCoordinates(CellRefiner refiner, DM dm, PetscInt depthSize[], DM rdm)
{
  PetscSection   coordSection, coordSectionNew;
  Vec            coordinates, coordinatesNew;
  PetscScalar   *coords, *coordsNew;
  PetscInt       dim, depth, coordSizeNew, cStart, cEnd, c, vStart, vStartNew, vEnd, v, fStart, fEnd, fMax, f;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, NULL, &fMax, NULL, NULL);CHKERRQ(ierr);
  ierr = GetDepthStart_Private(depth, depthSize, NULL, NULL, NULL, &vStartNew);CHKERRQ(ierr);
  ierr = DMPlexGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
  ierr = PetscSectionCreate(PetscObjectComm((PetscObject)dm), &coordSectionNew);CHKERRQ(ierr);
  ierr = PetscSectionSetNumFields(coordSectionNew, 1);CHKERRQ(ierr);
  ierr = PetscSectionSetFieldComponents(coordSectionNew, 0, dim);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(coordSectionNew, vStartNew, vStartNew+depthSize[0]);CHKERRQ(ierr);
  if (fMax < 0) fMax = fEnd;
  switch (refiner) {
  case 1:
  case 2:
  case 3:
    /* Simplicial and Hex 2D */
    /* All vertices have the dim coordinates */
    for (v = vStartNew; v < vStartNew+depthSize[0]; ++v) {
      ierr = PetscSectionSetDof(coordSectionNew, v, dim);CHKERRQ(ierr);
      ierr = PetscSectionSetFieldDof(coordSectionNew, v, 0, dim);CHKERRQ(ierr);
    }
    ierr = PetscSectionSetUp(coordSectionNew);CHKERRQ(ierr);
    ierr = DMPlexSetCoordinateSection(rdm, coordSectionNew);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
    ierr = PetscSectionGetStorageSize(coordSectionNew, &coordSizeNew);CHKERRQ(ierr);
    ierr = VecCreate(PetscObjectComm((PetscObject)dm), &coordinatesNew);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) coordinatesNew, "coordinates");CHKERRQ(ierr);
    ierr = VecSetSizes(coordinatesNew, coordSizeNew, PETSC_DETERMINE);CHKERRQ(ierr);
    ierr = VecSetFromOptions(coordinatesNew);CHKERRQ(ierr);
    ierr = VecGetArray(coordinates, &coords);CHKERRQ(ierr);
    ierr = VecGetArray(coordinatesNew, &coordsNew);CHKERRQ(ierr);
    /* Old vertices have the same coordinates */
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt newv = vStartNew + (v - vStart);
      PetscInt       off, offnew, d;

      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      ierr = PetscSectionGetOffset(coordSectionNew, newv, &offnew);CHKERRQ(ierr);
      for (d = 0; d < dim; ++d) {
        coordsNew[offnew+d] = coords[off+d];
      }
    }
    /* Face vertices have the average of endpoint coordinates */
    for (f = fStart; f < fMax; ++f) {
      const PetscInt  newv = vStartNew + (vEnd - vStart) + (f - fStart);
      const PetscInt *cone;
      PetscInt        coneSize, offA, offB, offnew, d;

      ierr = DMPlexGetConeSize(dm, f, &coneSize);CHKERRQ(ierr);
      if (coneSize != 2) SETERRQ2(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "Face %d cone should have two vertices, not %d", f, coneSize);
      ierr = DMPlexGetCone(dm, f, &cone);CHKERRQ(ierr);
      ierr = PetscSectionGetOffset(coordSection, cone[0], &offA);CHKERRQ(ierr);
      ierr = PetscSectionGetOffset(coordSection, cone[1], &offB);CHKERRQ(ierr);
      ierr = PetscSectionGetOffset(coordSectionNew, newv, &offnew);CHKERRQ(ierr);
      for (d = 0; d < dim; ++d) {
        coordsNew[offnew+d] = 0.5*(coords[offA+d] + coords[offB+d]);
      }
    }
    /* Just Hex 2D */
    if (refiner == 2) {
      /* Cell vertices have the average of corner coordinates */
      for (c = cStart; c < cEnd; ++c) {
        const PetscInt newv = vStartNew + (vEnd - vStart) + (fEnd - fStart) + (c - cStart);
        PetscInt      *cone = NULL;
        PetscInt       closureSize, coneSize = 0, offA, offB, offC, offD, offnew, p, d;

        ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &cone);CHKERRQ(ierr);
        for (p = 0; p < closureSize*2; p += 2) {
          const PetscInt point = cone[p];
          if ((point >= vStart) && (point < vEnd)) cone[coneSize++] = point;
        }
        if (coneSize != 4) SETERRQ2(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "Quad %d cone should have four vertices, not %d", c, coneSize);
        ierr = PetscSectionGetOffset(coordSection, cone[0], &offA);CHKERRQ(ierr);
        ierr = PetscSectionGetOffset(coordSection, cone[1], &offB);CHKERRQ(ierr);
        ierr = PetscSectionGetOffset(coordSection, cone[2], &offC);CHKERRQ(ierr);
        ierr = PetscSectionGetOffset(coordSection, cone[3], &offD);CHKERRQ(ierr);
        ierr = PetscSectionGetOffset(coordSectionNew, newv, &offnew);CHKERRQ(ierr);
        for (d = 0; d < dim; ++d) {
          coordsNew[offnew+d] = 0.25*(coords[offA+d] + coords[offB+d] + coords[offC+d] + coords[offD+d]);
        }
        ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &cone);CHKERRQ(ierr);
      }
    }
    ierr = VecRestoreArray(coordinates, &coords);CHKERRQ(ierr);
    ierr = VecRestoreArray(coordinatesNew, &coordsNew);CHKERRQ(ierr);
    ierr = DMSetCoordinatesLocal(rdm, coordinatesNew);CHKERRQ(ierr);
    ierr = VecDestroy(&coordinatesNew);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&coordSectionNew);CHKERRQ(ierr);
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexCreateProcessSF"
PetscErrorCode DMPlexCreateProcessSF(DM dm, PetscSF sfPoint, IS *processRanks, PetscSF *sfProcess)
{
  PetscInt           numRoots, numLeaves, l;
  const PetscInt    *localPoints;
  const PetscSFNode *remotePoints;
  PetscInt          *localPointsNew;
  PetscSFNode       *remotePointsNew;
  PetscInt          *ranks, *ranksNew;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = PetscSFGetGraph(sfPoint, &numRoots, &numLeaves, &localPoints, &remotePoints);CHKERRQ(ierr);
  ierr = PetscMalloc(numLeaves * sizeof(PetscInt), &ranks);CHKERRQ(ierr);
  for (l = 0; l < numLeaves; ++l) {
    ranks[l] = remotePoints[l].rank;
  }
  ierr = PetscSortRemoveDupsInt(&numLeaves, ranks);CHKERRQ(ierr);
  ierr = PetscMalloc(numLeaves * sizeof(PetscInt),    &ranksNew);CHKERRQ(ierr);
  ierr = PetscMalloc(numLeaves * sizeof(PetscInt),    &localPointsNew);CHKERRQ(ierr);
  ierr = PetscMalloc(numLeaves * sizeof(PetscSFNode), &remotePointsNew);CHKERRQ(ierr);
  for (l = 0; l < numLeaves; ++l) {
    ranksNew[l]              = ranks[l];
    localPointsNew[l]        = l;
    remotePointsNew[l].index = 0;
    remotePointsNew[l].rank  = ranksNew[l];
  }
  ierr = PetscFree(ranks);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm), numLeaves, ranksNew, PETSC_OWN_POINTER, processRanks);CHKERRQ(ierr);
  ierr = PetscSFCreate(PetscObjectComm((PetscObject)dm), sfProcess);CHKERRQ(ierr);
  ierr = PetscSFSetFromOptions(*sfProcess);CHKERRQ(ierr);
  ierr = PetscSFSetGraph(*sfProcess, 1, numLeaves, localPointsNew, PETSC_OWN_POINTER, remotePointsNew, PETSC_OWN_POINTER);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellRefinerCreateSF"
PetscErrorCode CellRefinerCreateSF(CellRefiner refiner, DM dm, PetscInt depthSize[], DM rdm)
{
  PetscSF            sf, sfNew, sfProcess;
  IS                 processRanks;
  MPI_Datatype       depthType;
  PetscInt           numRoots, numLeaves, numLeavesNew = 0, l, m;
  const PetscInt    *localPoints, *neighbors;
  const PetscSFNode *remotePoints;
  PetscInt          *localPointsNew;
  PetscSFNode       *remotePointsNew;
  PetscInt          *depthSizeOld, *rdepthSize, *rdepthSizeOld, *rdepthMaxOld, *rvStart, *rvStartNew, *reStart, *reStartNew, *rfStart, *rfStartNew, *rcStart, *rcStartNew;
  PetscInt           depth, numNeighbors, pStartNew, pEndNew, cStart, cStartNew, cEnd, cMax, vStart, vStartNew, vEnd, vMax, fStart, fStartNew, fEnd, fMax, eStart, eStartNew, eEnd, eMax, r, n;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetChart(rdm, &pStartNew, &pEndNew);CHKERRQ(ierr);
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cMax, &fMax, &eMax, &vMax);CHKERRQ(ierr);
  ierr = GetDepthStart_Private(depth, depthSize, &cStartNew, &fStartNew, &eStartNew, &vStartNew);CHKERRQ(ierr);
  switch (refiner) {
  case 3:
    if (cMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No cell maximum specified in hybrid mesh");
    cMax = PetscMin(cEnd, cMax);
    if (fMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No face maximum specified in hybrid mesh");
    fMax = PetscMin(fEnd, fMax);
  }
  ierr = DMGetPointSF(dm, &sf);CHKERRQ(ierr);
  ierr = DMGetPointSF(rdm, &sfNew);CHKERRQ(ierr);
  /* Caculate size of new SF */
  ierr = PetscSFGetGraph(sf, &numRoots, &numLeaves, &localPoints, &remotePoints);CHKERRQ(ierr);
  if (numRoots < 0) PetscFunctionReturn(0);
  for (l = 0; l < numLeaves; ++l) {
    const PetscInt p = localPoints[l];

    switch (refiner) {
    case 1:
      /* Simplicial 2D */
      if ((p >= vStart) && (p < vEnd)) {
        /* Old vertices stay the same */
        ++numLeavesNew;
      } else if ((p >= fStart) && (p < fEnd)) {
        /* Old faces add new faces and vertex */
        numLeavesNew += 1 + 2;
      } else if ((p >= cStart) && (p < cEnd)) {
        /* Old cells add new cells and interior faces */
        numLeavesNew += 4 + 3;
      }
      break;
    case 2:
      /* Hex 2D */
      if ((p >= vStart) && (p < vEnd)) {
        /* Old vertices stay the same */
        ++numLeavesNew;
      } else if ((p >= fStart) && (p < fEnd)) {
        /* Old faces add new faces and vertex */
        numLeavesNew += 1 + 2;
      } else if ((p >= cStart) && (p < cEnd)) {
        /* Old cells add new cells and interior faces */
        numLeavesNew += 4 + 4;
      }
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
    }
  }
  /* Communicate depthSizes for each remote rank */
  ierr = DMPlexCreateProcessSF(dm, sf, &processRanks, &sfProcess);CHKERRQ(ierr);
  ierr = ISGetLocalSize(processRanks, &numNeighbors);CHKERRQ(ierr);
  ierr = PetscMalloc5((depth+1)*numNeighbors,PetscInt,&rdepthSize,numNeighbors,PetscInt,&rvStartNew,numNeighbors,PetscInt,&reStartNew,numNeighbors,PetscInt,&rfStartNew,numNeighbors,PetscInt,&rcStartNew);CHKERRQ(ierr);
  ierr = PetscMalloc7(depth+1,PetscInt,&depthSizeOld,(depth+1)*numNeighbors,PetscInt,&rdepthSizeOld,(depth+1)*numNeighbors,PetscInt,&rdepthMaxOld,numNeighbors,PetscInt,&rvStart,numNeighbors,PetscInt,&reStart,numNeighbors,PetscInt,&rfStart,numNeighbors,PetscInt,&rcStart);CHKERRQ(ierr);
  ierr = MPI_Type_contiguous(depth+1, MPIU_INT, &depthType);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&depthType);CHKERRQ(ierr);
  ierr = PetscSFBcastBegin(sfProcess, depthType, depthSize, rdepthSize);CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfProcess, depthType, depthSize, rdepthSize);CHKERRQ(ierr);
  for (n = 0; n < numNeighbors; ++n) {
    ierr = GetDepthStart_Private(depth, &rdepthSize[n*(depth+1)], &rcStartNew[n], &rfStartNew[n], &reStartNew[n], &rvStartNew[n]);CHKERRQ(ierr);
  }
  depthSizeOld[depth]   = cMax;
  depthSizeOld[0]       = vMax;
  depthSizeOld[depth-1] = fMax;
  depthSizeOld[1]       = eMax;

  ierr = PetscSFBcastBegin(sfProcess, depthType, depthSizeOld, rdepthMaxOld);CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfProcess, depthType, depthSizeOld, rdepthMaxOld);CHKERRQ(ierr);

  depthSizeOld[depth]   = cEnd - cStart;
  depthSizeOld[0]       = vEnd - vStart;
  depthSizeOld[depth-1] = fEnd - fStart;
  depthSizeOld[1]       = eEnd - eStart;

  ierr = PetscSFBcastBegin(sfProcess, depthType, depthSizeOld, rdepthSizeOld);CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfProcess, depthType, depthSizeOld, rdepthSizeOld);CHKERRQ(ierr);
  for (n = 0; n < numNeighbors; ++n) {
    ierr = GetDepthStart_Private(depth, &rdepthSizeOld[n*(depth+1)], &rcStart[n], &rfStart[n], &reStart[n], &rvStart[n]);CHKERRQ(ierr);
  }
  ierr = MPI_Type_free(&depthType);CHKERRQ(ierr);
  ierr = PetscSFDestroy(&sfProcess);CHKERRQ(ierr);
  /* Calculate new point SF */
  ierr = PetscMalloc(numLeavesNew * sizeof(PetscInt),    &localPointsNew);CHKERRQ(ierr);
  ierr = PetscMalloc(numLeavesNew * sizeof(PetscSFNode), &remotePointsNew);CHKERRQ(ierr);
  ierr = ISGetIndices(processRanks, &neighbors);CHKERRQ(ierr);
  for (l = 0, m = 0; l < numLeaves; ++l) {
    PetscInt    p     = localPoints[l];
    PetscInt    rp    = remotePoints[l].index, n;
    PetscMPIInt rrank = remotePoints[l].rank;

    ierr = PetscFindInt(rrank, numNeighbors, neighbors, &n);CHKERRQ(ierr);
    if (n < 0) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Could not locate remote rank %d", rrank);
    switch (refiner) {
    case 1:
      /* Simplicial 2D */
      if ((p >= vStart) && (p < vEnd)) {
        /* Old vertices stay the same */
        localPointsNew[m]        = vStartNew     + (p  - vStart);
        remotePointsNew[m].index = rvStartNew[n] + (rp - rvStart[n]);
        remotePointsNew[m].rank  = rrank;
        ++m;
      } else if ((p >= fStart) && (p < fEnd)) {
        /* Old faces add new faces and vertex */
        localPointsNew[m]        = vStartNew     + (vEnd - vStart)              + (p  - fStart);
        remotePointsNew[m].index = rvStartNew[n] + rdepthSizeOld[n*(depth+1)+0] + (rp - rfStart[n]);
        remotePointsNew[m].rank  = rrank;
        ++m;
        for (r = 0; r < 2; ++r, ++m) {
          localPointsNew[m]        = fStartNew     + (p  - fStart)*2     + r;
          remotePointsNew[m].index = rfStartNew[n] + (rp - rfStart[n])*2 + r;
          remotePointsNew[m].rank  = rrank;
        }
      } else if ((p >= cStart) && (p < cEnd)) {
        /* Old cells add new cells and interior faces */
        for (r = 0; r < 4; ++r, ++m) {
          localPointsNew[m]        = cStartNew     + (p  - cStart)*4     + r;
          remotePointsNew[m].index = rcStartNew[n] + (rp - rcStart[n])*4 + r;
          remotePointsNew[m].rank  = rrank;
        }
        for (r = 0; r < 3; ++r, ++m) {
          localPointsNew[m]        = fStartNew     + (fEnd - fStart)*2                    + (p  - cStart)*3     + r;
          remotePointsNew[m].index = rfStartNew[n] + rdepthSizeOld[n*(depth+1)+depth-1]*2 + (rp - rcStart[n])*3 + r;
          remotePointsNew[m].rank  = rrank;
        }
      }
      break;
    case 2:
      /* Hex 2D */
      if ((p >= vStart) && (p < vEnd)) {
        /* Old vertices stay the same */
        localPointsNew[m]        = vStartNew     + (p  - vStart);
        remotePointsNew[m].index = rvStartNew[n] + (rp - rvStart[n]);
        remotePointsNew[m].rank  = rrank;
        ++m;
      } else if ((p >= fStart) && (p < fEnd)) {
        /* Old faces add new faces and vertex */
        localPointsNew[m]        = vStartNew     + (vEnd - vStart)              + (p  - fStart);
        remotePointsNew[m].index = rvStartNew[n] + rdepthSizeOld[n*(depth+1)+0] + (rp - rfStart[n]);
        remotePointsNew[m].rank  = rrank;
        ++m;
        for (r = 0; r < 2; ++r, ++m) {
          localPointsNew[m]        = fStartNew     + (p  - fStart)*2     + r;
          remotePointsNew[m].index = rfStartNew[n] + (rp - rfStart[n])*2 + r;
          remotePointsNew[m].rank  = rrank;
        }
      } else if ((p >= cStart) && (p < cEnd)) {
        /* Old cells add new cells and interior faces */
        for (r = 0; r < 4; ++r, ++m) {
          localPointsNew[m]        = cStartNew     + (p  - cStart)*4     + r;
          remotePointsNew[m].index = rcStartNew[n] + (rp - rcStart[n])*4 + r;
          remotePointsNew[m].rank  = rrank;
        }
        for (r = 0; r < 4; ++r, ++m) {
          localPointsNew[m]        = fStartNew     + (fEnd - fStart)*2                    + (p  - cStart)*4     + r;
          remotePointsNew[m].index = rfStartNew[n] + rdepthSizeOld[n*(depth+1)+depth-1]*2 + (rp - rcStart[n])*4 + r;
          remotePointsNew[m].rank  = rrank;
        }
      }
      break;
    case 3:
      /* Hybrid simplicial 2D */
      if ((p >= vStart) && (p < vEnd)) {
        /* Old vertices stay the same */
        localPointsNew[m]        = vStartNew     + (p  - vStart);
        remotePointsNew[m].index = rvStartNew[n] + (rp - rvStart[n]);
        remotePointsNew[m].rank  = rrank;
        ++m;
      } else if ((p >= fStart) && (p < fMax)) {
        /* Old interior faces add new faces and vertex */
        localPointsNew[m]        = vStartNew     + (vEnd - vStart)              + (p  - fStart);
        remotePointsNew[m].index = rvStartNew[n] + rdepthSizeOld[n*(depth+1)+0] + (rp - rfStart[n]);
        remotePointsNew[m].rank  = rrank;
        ++m;
        for (r = 0; r < 2; ++r, ++m) {
          localPointsNew[m]        = fStartNew     + (p  - fStart)*2     + r;
          remotePointsNew[m].index = rfStartNew[n] + (rp - rfStart[n])*2 + r;
          remotePointsNew[m].rank  = rrank;
        }
      } else if ((p >= fMax) && (p < fEnd)) {
        /* Old hybrid faces stay the same */
        localPointsNew[m]        = fStartNew     + (fMax                              - fStart)*2     + (p  - fMax);
        remotePointsNew[m].index = rfStartNew[n] + (rdepthMaxOld[n*(depth+1)+depth-1] - rfStart[n])*2 + (rp - rdepthMaxOld[n*(depth+1)+depth-1]);
        remotePointsNew[m].rank  = rrank;
        ++m;
      } else if ((p >= cStart) && (p < cMax)) {
        /* Old interior cells add new cells and interior faces */
        for (r = 0; r < 4; ++r, ++m) {
          localPointsNew[m]        = cStartNew     + (p  - cStart)*4     + r;
          remotePointsNew[m].index = rcStartNew[n] + (rp - rcStart[n])*4 + r;
          remotePointsNew[m].rank  = rrank;
        }
        for (r = 0; r < 3; ++r, ++m) {
          localPointsNew[m]        = fStartNew     + (fMax                              - fStart)*2     + (p  - cStart)*3     + r;
          remotePointsNew[m].index = rfStartNew[n] + (rdepthMaxOld[n*(depth+1)+depth-1] - rfStart[n])*2 + (rp - rcStart[n])*3 + r;
          remotePointsNew[m].rank  = rrank;
        }
      } else if ((p >= cStart) && (p < cMax)) {
        /* Old hybrid cells add new cells and hybrid face */
        for (r = 0; r < 2; ++r, ++m) {
          localPointsNew[m]        = cStartNew     + (p  - cStart)*4     + r;
          remotePointsNew[m].index = rcStartNew[n] + (rp - rcStart[n])*4 + r;
          remotePointsNew[m].rank  = rrank;
        }
        localPointsNew[m]        = fStartNew     + (fMax                              - fStart)*2     + (cMax                            - cStart)*3     + (p  - cMax);
        remotePointsNew[m].index = rfStartNew[n] + (rdepthMaxOld[n*(depth+1)+depth-1] - rfStart[n])*2 + (rdepthMaxOld[n*(depth+1)+depth] - rcStart[n])*3 + (rp - rdepthMaxOld[n*(depth+1)+depth]);
        remotePointsNew[m].rank  = rrank;
        ++m;
      }
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
    }
  }
  ierr = ISRestoreIndices(processRanks, &neighbors);CHKERRQ(ierr);
  ierr = ISDestroy(&processRanks);CHKERRQ(ierr);
  ierr = PetscSFSetGraph(sfNew, pEndNew-pStartNew, numLeavesNew, localPointsNew, PETSC_OWN_POINTER, remotePointsNew, PETSC_OWN_POINTER);CHKERRQ(ierr);
  ierr = PetscFree5(rdepthSize,rvStartNew,reStartNew,rfStartNew,rcStartNew);CHKERRQ(ierr);
  ierr = PetscFree6(depthSizeOld,rdepthSizeOld,rvStart,reStart,rfStart,rcStart);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CellRefinerCreateLabels"
PetscErrorCode CellRefinerCreateLabels(CellRefiner refiner, DM dm, PetscInt depthSize[], DM rdm)
{
  PetscInt       numLabels, l;
  PetscInt       newp, cStart, cStartNew, cEnd, cMax, vStart, vStartNew, vEnd, vMax, fStart, fStartNew, fEnd, fMax, eStart, eEnd, eMax, r;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);

  cStartNew = 0;
  vStartNew = depthSize[2];
  fStartNew = depthSize[2] + depthSize[0];

  ierr = DMPlexGetNumLabels(dm, &numLabels);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cMax, &fMax, &eMax, &vMax);CHKERRQ(ierr);
  switch (refiner) {
  case 3:
    if (cMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No cell maximum specified in hybrid mesh");
    cMax = PetscMin(cEnd, cMax);
    if (fMax < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No face maximum specified in hybrid mesh");
    fMax = PetscMin(fEnd, fMax);
  }
  for (l = 0; l < numLabels; ++l) {
    DMLabel         label, labelNew;
    const char     *lname;
    PetscBool       isDepth;
    IS              valueIS;
    const PetscInt *values;
    PetscInt        numValues, val;

    ierr = DMPlexGetLabelName(dm, l, &lname);CHKERRQ(ierr);
    ierr = PetscStrcmp(lname, "depth", &isDepth);CHKERRQ(ierr);
    if (isDepth) continue;
    ierr = DMPlexCreateLabel(rdm, lname);CHKERRQ(ierr);
    ierr = DMPlexGetLabel(dm, lname, &label);CHKERRQ(ierr);
    ierr = DMPlexGetLabel(rdm, lname, &labelNew);CHKERRQ(ierr);
    ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
    ierr = ISGetLocalSize(valueIS, &numValues);CHKERRQ(ierr);
    ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
    for (val = 0; val < numValues; ++val) {
      IS              pointIS;
      const PetscInt *points;
      PetscInt        numPoints, n;

      ierr = DMLabelGetStratumIS(label, values[val], &pointIS);CHKERRQ(ierr);
      ierr = ISGetLocalSize(pointIS, &numPoints);CHKERRQ(ierr);
      ierr = ISGetIndices(pointIS, &points);CHKERRQ(ierr);
      for (n = 0; n < numPoints; ++n) {
        const PetscInt p = points[n];
        switch (refiner) {
        case 1:
          /* Simplicial 2D */
          if ((p >= vStart) && (p < vEnd)) {
            /* Old vertices stay the same */
            newp = vStartNew + (p - vStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
          } else if ((p >= fStart) && (p < fEnd)) {
            /* Old faces add new faces and vertex */
            newp = vStartNew + (vEnd - vStart) + (p - fStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            for (r = 0; r < 2; ++r) {
              newp = fStartNew + (p - fStart)*2 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
          } else if ((p >= cStart) && (p < cEnd)) {
            /* Old cells add new cells and interior faces */
            for (r = 0; r < 4; ++r) {
              newp = cStartNew + (p - cStart)*4 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
            for (r = 0; r < 3; ++r) {
              newp = fStartNew + (fEnd - fStart)*2 + (p - cStart)*3 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
          }
          break;
        case 2:
          /* Hex 2D */
          if ((p >= vStart) && (p < vEnd)) {
            /* Old vertices stay the same */
            newp = vStartNew + (p - vStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
          } else if ((p >= fStart) && (p < fEnd)) {
            /* Old faces add new faces and vertex */
            newp = vStartNew + (vEnd - vStart) + (p - fStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            for (r = 0; r < 2; ++r) {
              newp = fStartNew + (p - fStart)*2 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
          } else if ((p >= cStart) && (p < cEnd)) {
            /* Old cells add new cells and interior faces and vertex */
            for (r = 0; r < 4; ++r) {
              newp = cStartNew + (p - cStart)*4 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
            for (r = 0; r < 4; ++r) {
              newp = fStartNew + (fEnd - fStart)*2 + (p - cStart)*4 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
            newp = vStartNew + (vEnd - vStart) + (fEnd - fStart) + (p - cStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
          }
          break;
        case 3:
          /* Hybrid simplicial 2D */
          if ((p >= vStart) && (p < vEnd)) {
            /* Old vertices stay the same */
            newp = vStartNew + (p - vStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
          } else if ((p >= fStart) && (p < fMax)) {
            /* Old interior faces add new faces and vertex */
            newp = vStartNew + (vEnd - vStart) + (p - fStart);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            for (r = 0; r < 2; ++r) {
              newp = fStartNew + (p - fStart)*2 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
          } else if ((p >= fMax) && (p < fEnd)) {
            /* Old hybrid faces stay the same */
            newp = fStartNew + (fMax - fStart)*2 + (p - fMax);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
          } else if ((p >= cStart) && (p < cMax)) {
            /* Old interior cells add new cells and interior faces */
            for (r = 0; r < 4; ++r) {
              newp = cStartNew + (p - cStart)*4 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
            for (r = 0; r < 3; ++r) {
              newp = fStartNew + (fEnd - fStart)*2 + (p - cStart)*3 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
          } else if ((p >= cMax) && (p < cEnd)) {
            /* Old hybrid cells add new cells and hybrid face */
            for (r = 0; r < 2; ++r) {
              newp = cStartNew + (cMax - cStart)*4 + (p - cMax)*2 + r;
              ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
            }
            newp = fStartNew + (fMax - fStart)*2 + (cMax - cStart)*3 + (p - cMax);
            ierr = DMLabelSetValue(labelNew, newp, values[val]);CHKERRQ(ierr);
          }
          break;
        default:
          SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown cell refiner %d", refiner);
        }
      }
      ierr = ISRestoreIndices(pointIS, &points);CHKERRQ(ierr);
      ierr = ISDestroy(&pointIS);CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(valueIS, &values);CHKERRQ(ierr);
    ierr = ISDestroy(&valueIS);CHKERRQ(ierr);
    if (0) {
      ierr = PetscViewerASCIISynchronizedAllow(PETSC_VIEWER_STDOUT_WORLD, PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMLabelView(labelNew, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = PetscViewerFlush(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexRefineUniform_Internal"
/* This will only work for interpolated meshes */
PetscErrorCode DMPlexRefineUniform_Internal(DM dm, CellRefiner cellRefiner, DM *dmRefined)
{
  DM             rdm;
  PetscInt      *depthSize;
  PetscInt       dim, depth = 0, d, pStart = 0, pEnd = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreate(PetscObjectComm((PetscObject)dm), &rdm);CHKERRQ(ierr);
  ierr = DMSetType(rdm, DMPLEX);CHKERRQ(ierr);
  ierr = DMPlexGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexSetDimension(rdm, dim);CHKERRQ(ierr);
  /* Calculate number of new points of each depth */
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = PetscMalloc((depth+1) * sizeof(PetscInt), &depthSize);CHKERRQ(ierr);
  ierr = PetscMemzero(depthSize, (depth+1) * sizeof(PetscInt));CHKERRQ(ierr);
  ierr = CellRefinerGetSizes(cellRefiner, dm, depthSize);CHKERRQ(ierr);
  /* Step 1: Set chart */
  for (d = 0; d <= depth; ++d) pEnd += depthSize[d];
  ierr = DMPlexSetChart(rdm, pStart, pEnd);CHKERRQ(ierr);
  /* Step 2: Set cone/support sizes */
  ierr = CellRefinerSetConeSizes(cellRefiner, dm, depthSize, rdm);CHKERRQ(ierr);
  /* Step 3: Setup refined DM */
  ierr = DMSetUp(rdm);CHKERRQ(ierr);
  /* Step 4: Set cones and supports */
  ierr = CellRefinerSetCones(cellRefiner, dm, depthSize, rdm);CHKERRQ(ierr);
  /* Step 5: Stratify */
  ierr = DMPlexStratify(rdm);CHKERRQ(ierr);
  /* Step 6: Set coordinates for vertices */
  ierr = CellRefinerSetCoordinates(cellRefiner, dm, depthSize, rdm);CHKERRQ(ierr);
  /* Step 7: Create pointSF */
  ierr = CellRefinerCreateSF(cellRefiner, dm, depthSize, rdm);CHKERRQ(ierr);
  /* Step 8: Create labels */
  ierr = CellRefinerCreateLabels(cellRefiner, dm, depthSize, rdm);CHKERRQ(ierr);
  ierr = PetscFree(depthSize);CHKERRQ(ierr);

  *dmRefined = rdm;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexSetRefinementUniform"
PetscErrorCode DMPlexSetRefinementUniform(DM dm, PetscBool refinementUniform)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  mesh->refinementUniform = refinementUniform;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexGetRefinementUniform"
PetscErrorCode DMPlexGetRefinementUniform(DM dm, PetscBool *refinementUniform)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(refinementUniform,  2);
  *refinementUniform = mesh->refinementUniform;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexSetRefinementLimit"
PetscErrorCode DMPlexSetRefinementLimit(DM dm, PetscReal refinementLimit)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  mesh->refinementLimit = refinementLimit;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexGetRefinementLimit"
PetscErrorCode DMPlexGetRefinementLimit(DM dm, PetscReal *refinementLimit)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(refinementLimit,  2);
  /* if (mesh->refinementLimit < 0) = getMaxVolume()/2.0; */
  *refinementLimit = mesh->refinementLimit;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexGetCellRefiner_Internal"
PetscErrorCode DMPlexGetCellRefiner_Internal(DM dm, CellRefiner *cellRefiner)
{
  PetscInt       dim, cStart, coneSize, cMax;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, NULL);CHKERRQ(ierr);
  ierr = DMPlexGetConeSize(dm, cStart, &coneSize);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cMax, NULL, NULL, NULL);CHKERRQ(ierr);
  switch (dim) {
  case 2:
    switch (coneSize) {
    case 3:
      if (cMax >= 0) *cellRefiner = 3; /* Hybrid */
      else *cellRefiner = 1; /* Triangular */
      break;
    case 4:
      if (cMax >= 0) *cellRefiner = 4; /* Hybrid */
      else *cellRefiner = 2; /* Quadrilateral */
      break;
    default:
      SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown coneSize %d in dimension %d for cell refiner", coneSize, dim);
    }
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unknown dimension %d for cell refiner", dim);
  }
  PetscFunctionReturn(0);
}
