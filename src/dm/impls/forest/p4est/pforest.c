#include <petscds.h>
#include <petsc/private/dmimpl.h>
#include <petsc/private/dmforestimpl.h>
#include <petsc/private/dmpleximpl.h>
#include <petsc/private/dmlabelimpl.h>
#include <petsc/private/viewerimpl.h>
#include <../src/sys/classes/viewer/impls/vtk/vtkvimpl.h>
#include "petsc_p4est_package.h"

#if defined(PETSC_HAVE_P4EST)

/* we need two levels of macros to stringify the results of macro expansion */
#define _pforest_string(a) _pforest_string_internal(a)
#define _pforest_string_internal(a) #a

#if !defined(P4_TO_P8)
#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_geometry.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#include <p4est_plex.h>
#include <p4est_bits.h>
#include <p4est_algorithms.h>
#else
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_geometry.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#include <p8est_plex.h>
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#endif

typedef enum {PATTERN_HASH,PATTERN_FRACTAL,PATTERN_CORNER,PATTERN_CENTER,PATTERN_COUNT} DMRefinePattern;
static const char *DMRefinePatternName[PATTERN_COUNT] = {"hash","fractal","corner","center"};

typedef struct _DMRefinePatternCtx
{
  PetscInt  corner;
  PetscBool fractal[P4EST_CHILDREN];
  PetscReal hashLikelihood;
  PetscInt  maxLevel;
  p4est_refine_t refine_fn;
}
DMRefinePatternCtx;

static int DMRefinePattern_Corner(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  p4est_quadrant_t   root, rootcorner;
  DMRefinePatternCtx *ctx;

  ctx = (DMRefinePatternCtx *) p4est->user_pointer;
  if (quadrant->level >= ctx->maxLevel) {
    return 0;
  }

  memset(&root,0,sizeof(p4est_quadrant_t));
  p4est_quadrant_corner_descendant(&root,&rootcorner,ctx->corner,quadrant->level);
  if (p4est_quadrant_is_equal(quadrant,&rootcorner)) {
    return 1;
  }
  return 0;
}

static int DMRefinePattern_Center(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  int                cid;
  p4est_quadrant_t   ancestor, ancestorcorner;
  DMRefinePatternCtx *ctx;

  ctx = (DMRefinePatternCtx *) p4est->user_pointer;
  if (quadrant->level >= ctx->maxLevel) {
    return 0;
  }
  if (quadrant->level <= 1) {
    return 1;
  }

  p4est_quadrant_ancestor(quadrant,1,&ancestor);
  cid = p4est_quadrant_child_id(&ancestor);
  p4est_quadrant_corner_descendant(&ancestor,&ancestorcorner,P4EST_CHILDREN - 1 - cid,quadrant->level);
  if (p4est_quadrant_is_equal(quadrant,&ancestorcorner)) {
    return 1;
  }
  return 0;
}

static int DMRefinePattern_Fractal(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  int                cid;
  DMRefinePatternCtx *ctx;

  ctx = (DMRefinePatternCtx *) p4est->user_pointer;
  if (quadrant->level >= ctx->maxLevel) {
    return 0;
  }
  if (!quadrant->level) {
    return 1;
  }
  cid = p4est_quadrant_child_id(quadrant);
  if (ctx->fractal[cid ^ ((int) (quadrant->level % P4EST_CHILDREN))]) {
    return 1;
  }
  return 0;
}

/* simplified from MurmurHash3 by Austin Appleby */
#define DMPROT32(x, y) ((x << y) | (x >> (32 - y)))
static uint32_t DMPforestHash(const uint32_t *blocks, uint32_t nblocks) {
  uint32_t c1 = 0xcc9e2d51;
  uint32_t c2 = 0x1b873593;
  uint32_t r1 = 15;
  uint32_t r2 = 13;
  uint32_t m = 5;
  uint32_t n = 0xe6546b64;
  uint32_t hash = 0;
  int      len = nblocks * 4;
  int      i;

  for (i = 0; i < nblocks; i++) {
    uint32_t k;

    k = blocks[i];
    k *= c1;
    k = DMPROT32(k, r1);
    k *= c2;

    hash ^= k;
    hash = DMPROT32(hash, r2) * m + n;
  }

  hash ^= len;
  hash ^= (hash >> 16);
  hash *= 0x85ebca6b;
  hash ^= (hash >> 13);
  hash *= 0xc2b2ae35;
  hash ^= (hash >> 16);

  return hash;
}

#if defined(UINT32_MAX)
#define DMP4EST_HASH_MAX UINT32_MAX
#else
#define DMP4EST_HASH_MAX ((uint32_t) 0xffffffff)
#endif

static int DMRefinePattern_Hash(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  uint32_t           data[5];
  uint32_t           result;
  DMRefinePatternCtx *ctx;

  ctx = (DMRefinePatternCtx *) p4est->user_pointer;
  if (quadrant->level >= ctx->maxLevel) {
    return 0;
  }
  data[0] = ((uint32_t) quadrant->level) << 24;
  data[1] = (uint32_t) which_tree;
  data[2] = (uint32_t) quadrant->x;
  data[3] = (uint32_t) quadrant->y;
#if defined(P4_TO_P8)
  data[4] = (uint32_t) quadrant->z;
#endif

  result = DMPforestHash(data,2+P4EST_DIM);
  if (((double) result / (double) DMP4EST_HASH_MAX) < ctx->hashLikelihood) {
    return 1;
  }
  return 0;
}

#define DMConvert_pforest_plex _infix_pforest(DMConvert,_plex)
static PetscErrorCode DMConvert_pforest_plex(DM,DMType,DM*);

#define DMFTopology_pforest _append_pforest(DMFTopology)
typedef struct {
  PetscInt             refct;
  p4est_connectivity_t *conn;
  p4est_geometry_t     *geom;
  PetscInt             *tree_face_to_uniq; /* p4est does not explicitly enumerate facets, but we must to keep track of labels */
} DMFTopology_pforest;

#define DM_Forest_pforest _append_pforest(DM_Forest)
typedef struct {
  DMFTopology_pforest *topo;
  p4est_t             *forest;
  p4est_ghost_t       *ghost;
  p4est_lnodes_t      *lnodes;
  PetscBool            partition_for_coarsening;
  PetscBool            coarsen_hierarchy;
  PetscBool            labelsFinalized;
  PetscInt             cLocalStart;
  PetscInt             cLocalEnd;
  DM                   plex;
  char                *ghostName;
} DM_Forest_pforest;

#define DMFTopologyDestroy_pforest _append_pforest(DMFTopologyDestroy)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMFTopologyDestroy_pforest)
static PetscErrorCode DMFTopologyDestroy_pforest(DMFTopology_pforest **topo)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!(*topo)) {
    PetscFunctionReturn(0);
  }
  if (--((*topo)->refct) > 0) {
    *topo = NULL;
    PetscFunctionReturn(0);
  }
  PetscStackCallP4est(p4est_geometry_destroy,((*topo)->geom));
  PetscStackCallP4est(p4est_connectivity_destroy,((*topo)->conn));
  ierr = PetscFree((*topo)->tree_face_to_uniq);CHKERRQ(ierr);
  ierr = PetscFree(*topo);CHKERRQ(ierr);
  *topo = NULL;
  PetscFunctionReturn(0);
}

static PetscErrorCode PforestConnectivityEnumerateFacets(p4est_connectivity_t*,PetscInt **);

#define DMFTopologyCreateBrick_pforest _append_pforest(DMFTopologyCreateBrick)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMFTopologyCreateBrick_pforest)
static PetscErrorCode DMFTopologyCreateBrick_pforest(DM dm,PetscInt N[], PetscInt P[], PetscReal B[],DMFTopology_pforest **topo, PetscBool useMorton)
{
  double         *vertices;
  PetscInt       i, numVerts;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!useMorton) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Lexicographic ordering not implemented yet");
  ierr = PetscNewLog(dm,topo);CHKERRQ(ierr);

  (*topo)->refct = 1;
#if !defined(P4_TO_P8)
  PetscStackCallP4estReturn((*topo)->conn,p4est_connectivity_new_brick,((int) N[0], (int) N[1], (P[0] == DM_BOUNDARY_NONE) ? 0 : 1, (P[1] == DM_BOUNDARY_NONE) ? 0 : 1));
#else
  PetscStackCallP4estReturn((*topo)->conn,p8est_connectivity_new_brick,((int) N[0], (int) N[1], (int) N[2], (P[0] == DM_BOUNDARY_NONE) ? 0 : 1, (P[1] == DM_BOUNDARY_NONE) ? 0 : 1, (P[2] == DM_BOUNDARY_NONE) ? 0 : 1));
#endif
  numVerts = (*topo)->conn->num_vertices;
  vertices = (*topo)->conn->vertices;
  for (i = 0; i < 3 * numVerts; i++) {
    PetscInt j = i % 3;

    vertices[i] = B[2 * j] + (vertices[i]/N[j]) * (B[2 * j + 1] - B[2 * j]);
  }
  PetscStackCallP4estReturn((*topo)->geom,p4est_geometry_new_connectivity,((*topo)->conn));
  ierr = PforestConnectivityEnumerateFacets((*topo)->conn,&(*topo)->tree_face_to_uniq);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMFTopologyCreate_pforest _append_pforest(DMFTopologyCreate)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMFTopologyCreate_pforest)
static PetscErrorCode DMFTopologyCreate_pforest(DM dm, DMForestTopology topologyName, DMFTopology_pforest **topo)
{
  DM_Forest  *forest = (DM_Forest *) dm->data;
  const char *name   = (const char *) topologyName;
  const char *prefix;
  PetscBool  isBrick, isShell, isSphere;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidCharPointer(name,2);
  PetscValidPointer(topo,3);
  ierr = PetscStrcmp(name,"brick",&isBrick);CHKERRQ(ierr);
  ierr = PetscStrcmp(name,"shell",&isShell);CHKERRQ(ierr);
  ierr = PetscStrcmp(name,"sphere",&isSphere);CHKERRQ(ierr);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)dm,&prefix);CHKERRQ(ierr);
  if (isBrick) {
    PetscBool  flgN, flgP, flgM, flgB, useMorton = PETSC_TRUE;
    PetscInt   N[3] = {2,2,2}, P[3] = {0,0,0}, nretN = P4EST_DIM, nretP = P4EST_DIM, nretB = 2 * P4EST_DIM, i;
    PetscReal  B[6] = {0.0,1.0,0.0,1.0,0.0,1.0};

    if (forest->setfromoptionscalled) {
      ierr = PetscOptionsGetIntArray(((PetscObject)dm)->options,prefix,"-dm_p4est_brick_size",N,&nretN,&flgN);CHKERRQ(ierr);
      ierr = PetscOptionsGetIntArray(((PetscObject)dm)->options,prefix,"-dm_p4est_brick_periodicity",P,&nretP,&flgP);CHKERRQ(ierr);
      ierr = PetscOptionsGetRealArray(((PetscObject)dm)->options,prefix,"-dm_p4est_brick_bounds",B,&nretB,&flgB);CHKERRQ(ierr);
      ierr = PetscOptionsGetBool(((PetscObject)dm)->options,prefix,"-dm_p4est_brick_use_morton_curve",&useMorton,&flgM);CHKERRQ(ierr);
      if (flgN && nretN != P4EST_DIM) SETERRQ2(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_SIZ,"Need to give %d sizes in -dm_p4est_brick_size, gave %d",P4EST_DIM,nretN);
      if (flgP && nretP != P4EST_DIM) SETERRQ2(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_SIZ,"Need to give %d periodicities in -dm_p4est_brick_periodicity, gave %d",P4EST_DIM,nretP);
      if (flgB && nretB != 2 * P4EST_DIM) SETERRQ2(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_SIZ,"Need to give %d bounds in -dm_p4est_brick_bounds, gave %d",P4EST_DIM,nretP);
    }
    for (i = 0; i < P4EST_DIM; i++) {
      P[i] = (P[i] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE);
      if (!flgB) {
        B[2 * i + 1] = N[i];
      }
    }
    ierr = DMFTopologyCreateBrick_pforest(dm,N,P,B,topo,useMorton);CHKERRQ(ierr);
  }
  else {
    ierr = PetscNewLog(dm,topo);CHKERRQ(ierr);

    (*topo)->refct = 1;
    PetscStackCallP4estReturn((*topo)->conn,p4est_connectivity_new_byname,(name));
#if !defined(P4_TO_P8)
    PetscStackCallP4estReturn((*topo)->geom,p4est_geometry_new_connectivity,((*topo)->conn));
#else
    if (isShell) {
      PetscReal R2 = 1., R1 = .55;

      if (forest->setfromoptionscalled) {
        ierr = PetscOptionsGetReal(((PetscObject)dm)->options,prefix,"-dm_p4est_shell_outer_radius",&R2,NULL);CHKERRQ(ierr);
        ierr = PetscOptionsGetReal(((PetscObject)dm)->options,prefix,"-dm_p4est_shell_inner_radius",&R1,NULL);CHKERRQ(ierr);
      }
      PetscStackCallP4estReturn((*topo)->geom,p8est_geometry_new_shell,((*topo)->conn,R2,R1));
    }
    else if (isSphere) {
      PetscReal R2 = 1., R1 = 0.191728, R0 = 0.039856;

      if (forest->setfromoptionscalled) {
        ierr = PetscOptionsGetReal(((PetscObject)dm)->options,prefix,"-dm_p4est_sphere_outer_radius",&R2,NULL);CHKERRQ(ierr);
        ierr = PetscOptionsGetReal(((PetscObject)dm)->options,prefix,"-dm_p4est_sphere_inner_radius",&R1,NULL);CHKERRQ(ierr);
        ierr = PetscOptionsGetReal(((PetscObject)dm)->options,prefix,"-dm_p4est_sphere_core_radius",&R0,NULL);CHKERRQ(ierr);
      }
      PetscStackCallP4estReturn((*topo)->geom,p8est_geometry_new_sphere,((*topo)->conn,R2,R1,R0));
    }
    else {
      PetscStackCallP4estReturn((*topo)->geom,p4est_geometry_new_connectivity,((*topo)->conn));
    }
#endif
    ierr = PforestConnectivityEnumerateFacets((*topo)->conn,&(*topo)->tree_face_to_uniq);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#define DMConvert_plex_pforest _append_pforest(DMConvert_plex)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMConvert_plex_pforest)
static PetscErrorCode DMConvert_plex_pforest(DM dm, DMType newtype, DM *pforest)
{
  MPI_Comm       comm;
  PetscBool      isPlex;
  PetscInt       dim;
  void           *ctx;
  PetscDS        ds;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  comm = PetscObjectComm((PetscObject)dm);
  ierr = PetscObjectTypeCompare((PetscObject)dm,DMPLEX,&isPlex);CHKERRQ(ierr);
  if (!isPlex) SETERRQ2(comm,PETSC_ERR_ARG_WRONG,"Expected DM type %s, got %s\n",DMPLEX,((PetscObject)dm)->type_name);
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  if (dim != P4EST_DIM) SETERRQ2(comm,PETSC_ERR_ARG_WRONG,"Expected DM dimension %d, got %d\n",P4EST_DIM,dim);
  ierr = DMCreate(comm,pforest);CHKERRQ(ierr);
  ierr = DMSetType(*pforest,DMPFOREST);CHKERRQ(ierr);
  ierr = DMForestSetBaseDM(*pforest,dm);CHKERRQ(ierr);
  ierr = DMGetApplicationContext(dm,&ctx);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(*pforest,ctx);CHKERRQ(ierr);
  ierr = DMGetDS(dm,&ds);CHKERRQ(ierr);
  ierr = DMSetDS(*pforest,ds);CHKERRQ(ierr);
  if (dm->maxCell) {
    const PetscReal *maxCell, *L;
    const DMBoundaryType *bd;

    ierr = DMGetPeriodicity(dm,&maxCell,&L,&bd);CHKERRQ(ierr);
    ierr = DMSetPeriodicity(*pforest,maxCell,L,bd);CHKERRQ(ierr);
  }
  ierr = DMBoundaryDestroy(&(*pforest)->boundary);CHKERRQ(ierr);
  ierr = DMBoundaryDuplicate(dm->boundary, &(*pforest)->boundary);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMForestDestroy_pforest _append_pforest(DMForestDestroy)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMForestDestroy_pforest)
static PetscErrorCode DMForestDestroy_pforest(DM dm)
{
  DM_Forest         *forest  = (DM_Forest *) dm->data;
  DM_Forest_pforest *pforest = (DM_Forest_pforest *) forest->data;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  if (pforest->lnodes) PetscStackCallP4est(p4est_lnodes_destroy,(pforest->lnodes));
  pforest->lnodes = NULL;
  if (pforest->ghost)  PetscStackCallP4est(p4est_ghost_destroy,(pforest->ghost));
  pforest->ghost = NULL;
  if (pforest->forest) PetscStackCallP4est(p4est_destroy,(pforest->forest));
  pforest->forest = NULL;
  ierr = DMFTopologyDestroy_pforest(&pforest->topo);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)dm,_pforest_string(DMConvert_plex_pforest) "_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)dm,_pforest_string(DMConvert_pforest_plex) "_C",NULL);CHKERRQ(ierr);
  ierr = PetscFree(pforest->ghostName);CHKERRQ(ierr);
  ierr = DMDestroy(&pforest->plex);CHKERRQ(ierr);
  ierr = PetscFree(forest->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMForestTemplate_pforest _append_pforest(DMForestTemplate)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMForestTemplate_pforest)
static PetscErrorCode DMForestTemplate_pforest(DM dm, DM tdm)
{
  DM_Forest_pforest *pforest  = (DM_Forest_pforest *) ((DM_Forest *) dm->data)->data;
  DM_Forest_pforest *tpforest = (DM_Forest_pforest *) ((DM_Forest *) tdm->data)->data;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  if (pforest->topo) pforest->topo->refct++;
  ierr = DMFTopologyDestroy_pforest(&(tpforest->topo));CHKERRQ(ierr);
  tpforest->topo = pforest->topo;
  PetscFunctionReturn(0);
}

#define DMPlexCreateConnectivity_pforest _append_pforest(DMPlexCreateConnectivity)
static PetscErrorCode DMPlexCreateConnectivity_pforest(DM,p4est_connectivity_t**,PetscInt**);

static int pforest_coarsen_uniform (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrants[])
{
  PetscInt minLevel = *((PetscInt *) p4est->user_pointer);

  if ((PetscInt) quadrants[0]->level > minLevel) {
    return 1;
  }
  return 0;
}

static int pforest_refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  PetscInt maxLevel = *((PetscInt *) p4est->user_pointer);

  if ((PetscInt) quadrant->level < maxLevel) {
    return 1;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestComputeLocalCellTransferSF_loop"
static PetscErrorCode DMPforestComputeLocalCellTransferSF_loop(p4est_t *p4estFrom, PetscInt FromOffset, p4est_t *p4estTo, PetscInt ToOffset, p4est_topidx_t flt, p4est_topidx_t llt, PetscInt *toFineLeavesCount, PetscInt *toLeaves, PetscSFNode *fromRoots, PetscInt *fromFineLeavesCount, PetscInt *fromLeaves, PetscSFNode *toRoots)
{
  PetscMPIInt rank = p4estFrom->mpirank;
  p4est_topidx_t t;
  PetscInt toFineLeaves = 0, fromFineLeaves = 0;

  PetscFunctionBegin;
  for (t = flt; t <= llt; t++) { /* count roots and leaves */
    p4est_tree_t *treeFrom = &(((p4est_tree_t *) p4estFrom->trees->array)[t]);
    p4est_tree_t *treeTo   = &(((p4est_tree_t *) p4estTo->trees->array)[t]);
    p4est_quadrant_t *firstFrom = &treeFrom->first_desc;
    p4est_quadrant_t *firstTo   = &treeTo->first_desc;
    PetscInt numFrom = (PetscInt) &treeFrom->quadrants.elem_count;
    PetscInt numTo   = (PetscInt) &treeTo->quadrants.elem_count;
    p4est_quadrant_t *quadsFrom = (p4est_quadrant_t *) treeFrom->quadrants.array;
    p4est_quadrant_t *quadsTo   = (p4est_quadrant_t *) treeTo->quadrants.array;
    PetscInt currentFrom, currentTo;
    PetscInt treeOffsetFrom = (PetscInt) treeFrom->quadrants_offset;
    PetscInt treeOffsetTo   = (PetscInt) treeTo->quadrants_offset;
    int comp;

    PetscStackCallP4estReturn(comp,p4est_quadrant_is_equal,(firstFrom,firstTo));
    if (!comp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"non-matching partitions");

    for (currentFrom = 0, currentTo = 0; currentFrom < numFrom && currentTo < numTo;) {
      p4est_quadrant_t *quadFrom = &quadsFrom[currentFrom];
      p4est_quadrant_t *quadTo = &quadsTo[currentTo];

      if (quadFrom->level == quadTo->level) {
        if (toLeaves) {
          toLeaves[toFineLeaves] = currentTo + treeOffsetTo + ToOffset;
          fromRoots[toFineLeaves].rank = rank;
          fromRoots[toFineLeaves].index = currentFrom + treeOffsetFrom + FromOffset;
        }
        toFineLeaves++;
        currentFrom++;
        currentTo++;
      }
      else {
        int fromIsAncestor;

        PetscStackCallP4estReturn(fromIsAncestor,p4est_quadrant_is_ancestor,(quadFrom,quadTo));
        if (fromIsAncestor) {
          p4est_quadrant_t lastDesc;

          if (toLeaves) {
            toLeaves[toFineLeaves] = currentTo + treeOffsetTo + ToOffset;
            fromRoots[toFineLeaves].rank = rank;
            fromRoots[toFineLeaves].index = currentFrom + treeOffsetFrom + FromOffset;
          }
          toFineLeaves++;
          currentTo++;
          PetscStackCallP4est(p4est_quadrant_last_descendant,(quadFrom,&lastDesc,quadTo->level));
          PetscStackCallP4estReturn(comp,p4est_quadrant_is_equal,(quadTo,&lastDesc));
          if (comp) {
            currentFrom++;
          }
        }
        else {
          p4est_quadrant_t lastDesc;

          if (fromLeaves) {
            fromLeaves[fromFineLeaves] = currentFrom + treeOffsetFrom + FromOffset;
            toRoots[fromFineLeaves].rank = rank;
            toRoots[fromFineLeaves].index = currentTo + treeOffsetTo + ToOffset;
          }
          fromFineLeaves++;
          currentFrom++;
          PetscStackCallP4est(p4est_quadrant_last_descendant,(quadTo,&lastDesc,quadFrom->level));
          PetscStackCallP4estReturn(comp,p4est_quadrant_is_equal,(quadFrom,&lastDesc));
          if (comp) {
            currentTo++;
          }
        }
      }
    }
  }
  *toFineLeavesCount = toFineLeaves;
  *fromFineLeavesCount = fromFineLeaves;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestComputeLocalCellTransferSF"
/* Puts identity in coarseToFine */
/* assumes a matching partition */
static PetscErrorCode DMPforestComputeLocalCellTransferSF(MPI_Comm comm, p4est_t *p4estFrom, PetscInt FromOffset, p4est_t *p4estTo, PetscInt ToOffset, PetscSF *fromCoarseToFine, PetscSF *toCoarseFromFine)
{
  p4est_connectivity_t * conn;
  p4est_topidx_t flt, llt;
  PetscSF fromCoarse, toCoarse;
  PetscInt numRootsFrom, numRootsTo, numLeavesFrom, numLeavesTo;
  PetscInt *fromLeaves = NULL, *toLeaves = NULL;
  PetscSFNode *fromRoots = NULL, *toRoots = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  conn = p4estFrom->connectivity;
  flt  = p4estFrom->first_local_tree;
  llt  = p4estFrom->last_local_tree;
  ierr = PetscSFCreate(comm,&fromCoarse);CHKERRQ(ierr);
  if (toCoarseFromFine) {
    ierr = PetscSFCreate(comm,&toCoarse);CHKERRQ(ierr);
  }
  numRootsFrom = p4estFrom->local_num_quadrants + FromOffset;
  numRootsTo = p4estTo->local_num_quadrants + ToOffset;
  ierr = DMPforestComputeLocalCellTransferSF_loop(p4estFrom,FromOffset,p4estTo,ToOffset,flt,llt,&numLeavesTo,NULL,NULL,&numLeavesFrom,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeavesTo,&toLeaves);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeavesTo,&fromRoots);CHKERRQ(ierr);
  if (toCoarseFromFine) {
    ierr = PetscMalloc1(numLeavesFrom,&fromLeaves);CHKERRQ(ierr);
    ierr = PetscMalloc1(numLeavesFrom,&fromRoots);CHKERRQ(ierr);
  }
  ierr = DMPforestComputeLocalCellTransferSF_loop(p4estFrom,FromOffset,p4estTo,ToOffset,flt,llt,&numLeavesTo,toLeaves,fromRoots,&numLeavesFrom,fromLeaves,toRoots);CHKERRQ(ierr);
  if (!ToOffset && (numLeavesTo == numRootsTo)) { /* compress */
    ierr = PetscFree(toLeaves);CHKERRQ(ierr);
    ierr = PetscSFSetGraph(fromCoarse,numRootsFrom,numLeavesTo,NULL,PETSC_OWN_POINTER,fromRoots,PETSC_OWN_POINTER);CHKERRQ(ierr);
  }
  else { /* generic */
    ierr = PetscSFSetGraph(fromCoarse,numRootsFrom,numLeavesTo,toLeaves,PETSC_OWN_POINTER,fromRoots,PETSC_OWN_POINTER);CHKERRQ(ierr);
  }
  *fromCoarseToFine = fromCoarse;
  if (toCoarseFromFine) {
    ierr = PetscSFSetGraph(toCoarse,numRootsTo,numLeavesFrom,fromLeaves,PETSC_OWN_POINTER,toRoots,PETSC_OWN_POINTER);CHKERRQ(ierr);
    *toCoarseFromFine = toCoarse;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestComputeOverlappingRanks"
/* range of processes whose B sections overlap this ranks A section */
static PetscErrorCode DMPforestComputeOverlappingRanks(PetscMPIInt size, PetscMPIInt rank, p4est_t *p4estA, p4est_t *p4estB, PetscInt *startB, PetscInt *endB)
{
  p4est_quadrant_t * myCoarseStart = &(p4estA->global_first_position[rank]);
  p4est_quadrant_t * myCoarseEnd   = &(p4estA->global_first_position[rank+1]);
  p4est_quadrant_t * globalFirstB  = p4estB->global_first_position;

  PetscFunctionBegin;
  *startB = -1;
  *endB = -1;
  if (p4estA->local_num_quadrants) {
    PetscInt lo, hi, guess;
    /* binary search to find interval containing myCoarseStart */
    lo    = 0;
    hi    = size;
    guess = rank;
    while (1) {
      int startCompMy, myCompEnd;

      PetscStackCallP4estReturn(startCompMy,p4est_quadrant_compare_piggy,(&globalFirstB[guess],myCoarseStart));
      PetscStackCallP4estReturn(myCompEnd,p4est_quadrant_compare_piggy,(myCoarseStart,&globalFirstB[guess+1]));
      if (startCompMy <= 0 && myCompEnd < 0) {
        *startB = guess;
        break;
      }
      else if (startCompMy > 0) { /* guess is to high */
        hi = guess;
      }
      else { /* guess is to low */
        lo = guess + 1;
      }
      guess = lo + (hi - lo) / 2;
    }
    /* reset bounds, but not guess */
    lo = 0;
    hi = size;
    while (1) {
      int startCompMy, myCompEnd;

      PetscStackCallP4estReturn(startCompMy,p4est_quadrant_compare_piggy,(&globalFirstB[guess],myCoarseEnd));
      PetscStackCallP4estReturn(myCompEnd,p4est_quadrant_compare_piggy,(myCoarseEnd,&globalFirstB[guess+1]));
      if (startCompMy < 0 && myCompEnd <= 0) { /* notice that the comparison operators are different from above */
        *endB = guess + 1;
        break;
      }
      else if (startCompMy >= 0) { /* guess is to high */
        hi = guess;
      }
      else { /* guess is to low */
        lo = guess + 1;
      }
      guess = lo + (hi - lo) / 2;
    }
  }
  PetscFunctionReturn(0);
}

#define DMSetUp_pforest _append_pforest(DMSetUp)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMSetUp_pforest)
static PetscErrorCode DMSetUp_pforest(DM dm)
{
  DM_Forest         *forest  = (DM_Forest *) dm->data;
  DM_Forest_pforest *pforest = (DM_Forest_pforest *) forest->data;
  DM                base, adaptFrom;
  DMForestTopology  topoName;
  PetscSF           preCoarseToFine = NULL, coarseToPreFine = NULL;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  /* sanity check */
  ierr = DMForestGetAdaptivityForest(dm,&adaptFrom);CHKERRQ(ierr);
  ierr = DMForestGetBaseDM(dm,&base);CHKERRQ(ierr);
  ierr = DMForestGetTopology(dm,&topoName);CHKERRQ(ierr);
  if (!adaptFrom && !base && !topoName) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_WRONGSTATE,"A forest needs a topology, a base DM, or a DM to adapt from");

  /* === Step 1: DMFTopology === */
  if (adaptFrom) { /* reference already created topology */
    PetscBool         ispforest;
    DM_Forest         *aforest  = (DM_Forest *) adaptFrom->data;
    DM_Forest_pforest *apforest = (DM_Forest_pforest *) aforest->data;

    ierr = PetscObjectTypeCompare((PetscObject)adaptFrom,DMPFOREST,&ispforest);CHKERRQ(ierr);
    if (!ispforest) SETERRQ2(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_NOTSAMETYPE,"Trying to adapt from %s, which is not %s\n",((PetscObject)adaptFrom)->type_name,DMPFOREST);
    if (!apforest->topo) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_WRONGSTATE,"The pre-adaptation forest must have a topology");
    ierr = DMSetUp(adaptFrom);CHKERRQ(ierr);
    ierr = DMForestGetBaseDM(dm,&base);CHKERRQ(ierr);
    ierr = DMForestGetTopology(dm,&topoName);CHKERRQ(ierr);
  }
  else if (base) { /* construct a connectivity from base */
    PetscBool isPlex, isDA;

    ierr = PetscObjectGetName((PetscObject)base,&topoName);CHKERRQ(ierr);
    ierr = DMForestSetTopology(dm,topoName);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)base,DMPLEX,&isPlex);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)base,DMDA,&isDA);CHKERRQ(ierr);
    if (isPlex) {
      MPI_Comm             comm = PetscObjectComm((PetscObject)dm);
      PetscInt             depth;
      PetscMPIInt          size;
      p4est_connectivity_t *conn = NULL;
      DMFTopology_pforest  *topo;
      PetscInt             *tree_face_to_uniq;
      PetscErrorCode       ierr;

      ierr = DMPlexGetDepth(base,&depth);CHKERRQ(ierr);
      if (depth == 1) {
        DM connDM;

        ierr = DMPlexInterpolate(base,&connDM);CHKERRQ(ierr);
        base = connDM;
        ierr = DMForestSetBaseDM(dm,base);CHKERRQ(ierr);
        ierr = DMDestroy(&connDM);CHKERRQ(ierr);
      }
      else if (depth != P4EST_DIM) {
        SETERRQ2(comm,PETSC_ERR_ARG_WRONG,"Base plex is neither interpolated nor uninterpolated? depth %d, expected 2 or %d\n",depth,P4EST_DIM + 1);
      }
      ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
      if (size > 1) {
        DM dmRedundant;

        ierr = DMPlexGetRedundantDM(base,&dmRedundant);CHKERRQ(ierr);
        if (!dmRedundant) SETERRQ(comm,PETSC_ERR_PLIB,"Could not create redundant DM\n");
        base = dmRedundant;
        ierr = DMForestSetBaseDM(dm,base);CHKERRQ(ierr);
        ierr = DMDestroy(&dmRedundant);CHKERRQ(ierr);
      }
      ierr = DMPlexCreateConnectivity_pforest(base,&conn,&tree_face_to_uniq);CHKERRQ(ierr);
      ierr = PetscNewLog(dm,&topo);CHKERRQ(ierr);
      topo->refct             = 1;
      topo->conn              = conn;
      PetscStackCallP4estReturn(topo->geom,p4est_geometry_new_connectivity,(conn));
      topo->tree_face_to_uniq = tree_face_to_uniq;
      pforest->topo           = topo;
    }
    else if (isDA) {
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_PLIB,"Not implemented yet");
#if 0
      PetscInt N[3], P[3];

      /* get the sizes, periodicities */
      /* ... */
                                                                  /* don't use Morton order */
      ierr = DMFTopologyCreateBrick_pforest(dm,N,P,&pforest->topo,PETSC_FALSE);CHKERRQ(ierr);
#endif
    }
    {
      PetscInt numLabels, l;

      ierr = DMGetNumLabels(base,&numLabels);CHKERRQ(ierr);
      for (l = 0; l < numLabels; l++) {
        PetscBool  isDepth, isGhost, isVTK;
        DMLabel    label, labelNew;
        PetscInt   defVal;
        const char *name;

        ierr = DMGetLabelName(base, l, &name);CHKERRQ(ierr);
        ierr = DMGetLabelByNum(base, l, &label);CHKERRQ(ierr);
        ierr = PetscStrcmp(name,"depth",&isDepth);CHKERRQ(ierr);
        if (isDepth) continue;
        ierr = PetscStrcmp(name,"ghost",&isGhost);CHKERRQ(ierr);
        if (isGhost) continue;
        ierr = PetscStrcmp(name,"vtk",&isVTK);CHKERRQ(ierr);
        if (isVTK) continue;
        ierr = DMCreateLabel(dm,name);CHKERRQ(ierr);
        ierr = DMGetLabel(dm,name,&labelNew);CHKERRQ(ierr);
        ierr = DMLabelGetDefaultValue(label,&defVal);CHKERRQ(ierr);
        ierr = DMLabelSetDefaultValue(labelNew,defVal);CHKERRQ(ierr);
      }
    }
  }
  else { /* construct from topology name */
    DMFTopology_pforest *topo;

    ierr = DMFTopologyCreate_pforest(dm,topoName,&topo);CHKERRQ(ierr);
    pforest->topo = topo;
    /* TODO: construct base? */
  }

  /* === Step 2: get the leaves of the forest === */
  if (adaptFrom) { /* start with the old forest */
    const char        *adaptName;
    DMLabel           adaptLabel;
    PetscInt          defaultValue;
    PetscInt          numValues;
    DM_Forest         *aforest  = (DM_Forest *) adaptFrom->data;
    DM_Forest_pforest *apforest = (DM_Forest_pforest *) aforest->data;
    PetscBool         computeAdaptSF;

    ierr = DMForestGetComputeAdaptivitySF(dm,&computeAdaptSF);CHKERRQ(ierr);
    PetscStackCallP4estReturn(pforest->forest,p4est_copy,(apforest->forest, 0)); /* 0 indicates no data copying */
    ierr = DMForestGetAdaptivityLabel(dm,&adaptName);CHKERRQ(ierr);
    if (adaptName) {
      ierr = DMGetLabel(adaptFrom,adaptName,&adaptLabel);CHKERRQ(ierr);
      if (!adaptLabel) SETERRQ(PetscObjectComm((PetscObject)adaptFrom),PETSC_ERR_USER,"No adaptivity label found in pre-adaptation forest");
      /* apply the refinement/coarsening by flags, plus minimum/maximum refinement */
      ierr = DMLabelGetNumValues(adaptLabel,&numValues);CHKERRQ(ierr);
      ierr = DMLabelGetDefaultValue(adaptLabel,&defaultValue);CHKERRQ(ierr);
      if (!numValues && defaultValue == DM_FOREST_COARSEN) { /* uniform coarsen */
        PetscInt minLevel;

        ierr = DMForestGetMinimumRefinement(dm,&minLevel);CHKERRQ(ierr);
        pforest->forest->user_pointer = (void *) &minLevel;
        PetscStackCallP4est(p4est_coarsen,(pforest->forest,0,pforest_coarsen_uniform,NULL));
        pforest->forest->user_pointer = (void *) dm;
        PetscStackCallP4est(p4est_balance,(pforest->forest,P4EST_CONNECT_FULL,NULL));
        /* we will have to change the offset after we compute the overlap */
        if (computeAdaptSF) {
          ierr = DMPforestComputeLocalCellTransferSF(PetscObjectComm((PetscObject)dm),pforest->forest,0,apforest->forest,apforest->cLocalStart,&coarseToPreFine,NULL);CHKERRQ(ierr);
        }
      }
      else if (!numValues && defaultValue == DM_FOREST_REFINE) { /* uniform refine */
        PetscInt maxLevel;

        ierr = DMForestGetMaximumRefinement(dm,&maxLevel);CHKERRQ(ierr);
        pforest->forest->user_pointer = (void *) &maxLevel;
        PetscStackCallP4est(p4est_refine,(pforest->forest,0,pforest_refine_uniform,NULL));
        pforest->forest->user_pointer = (void *) dm;
        PetscStackCallP4est(p4est_balance,(pforest->forest,P4EST_CONNECT_FULL,NULL));
        /* we will have to change the offset after we compute the overlap */
        if (computeAdaptSF) {
          ierr = DMPforestComputeLocalCellTransferSF(PetscObjectComm((PetscObject)dm),apforest->forest,apforest->cLocalStart,pforest->forest,0,&preCoarseToFine,NULL);CHKERRQ(ierr);
        }
      }
      else if (numValues) {
        SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"non-uniform adaptivity not implemented");
        PetscStackCallP4est(p4est_balance,(pforest->forest,P4EST_CONNECT_FULL,NULL));
        if (computeAdaptSF) {
          ierr = DMPforestComputeLocalCellTransferSF(PetscObjectComm((PetscObject)dm),apforest->forest,apforest->cLocalStart,pforest->forest,0,&preCoarseToFine,&coarseToPreFine);CHKERRQ(ierr);
        }
      }
      {
        PetscInt numLabels, l;

        ierr = DMGetNumLabels(adaptFrom,&numLabels);CHKERRQ(ierr);
        for (l = 0; l < numLabels; l++) {
          PetscBool  isDepth, isGhost, isVTK;
          DMLabel    label, labelNew;
          PetscInt   defVal;
          const char *name;

          ierr = DMGetLabelName(adaptFrom, l, &name);CHKERRQ(ierr);
          ierr = DMGetLabelByNum(adaptFrom, l, &label);CHKERRQ(ierr);
          ierr = PetscStrcmp(name,"depth",&isDepth);CHKERRQ(ierr);
          if (isDepth) continue;
          ierr = PetscStrcmp(name,"ghost",&isGhost);CHKERRQ(ierr);
          if (isGhost) continue;
          ierr = PetscStrcmp(name,"vtk",&isVTK);CHKERRQ(ierr);
          if (isVTK) continue;
          ierr = DMCreateLabel(dm,name);CHKERRQ(ierr);
          ierr = DMGetLabel(dm,name,&labelNew);CHKERRQ(ierr);
          ierr = DMLabelGetDefaultValue(label,&defVal);CHKERRQ(ierr);
          ierr = DMLabelSetDefaultValue(labelNew,defVal);CHKERRQ(ierr);
        }
      }
    }
  }
  else { /* initial */
    PetscInt initLevel, minLevel;

    ierr = DMForestGetInitialRefinement(dm,&initLevel);CHKERRQ(ierr);
    ierr = DMForestGetMinimumRefinement(dm,&minLevel);CHKERRQ(ierr);
    PetscStackCallP4estReturn(pforest->forest,p4est_new_ext,(PetscObjectComm((PetscObject)dm),pforest->topo->conn,
                                                             0,           /* minimum number of quadrants per processor */
                                                             initLevel,   /* level of refinement */
                                                             1,           /* uniform refinement */
                                                             0,           /* we don't allocate any per quadrant data */
                                                             NULL,        /* there is no special quadrant initialization */
                                                             (void *)dm)); /* this dm is the user context */
    if (forest->setfromoptionscalled) {
      PetscBool       flgPattern, flgFractal;
      PetscInt        corner = 0;
      PetscInt        corners[P4EST_CHILDREN], ncorner = P4EST_CHILDREN;
      PetscReal       likelihood = 1./ P4EST_DIM;
      PetscInt        pattern;
      const char      *prefix;

      ierr = PetscObjectGetOptionsPrefix((PetscObject)dm,&prefix);CHKERRQ(ierr);
      ierr = PetscOptionsGetEList(((PetscObject)dm)->options,prefix,"-dm_p4est_refine_pattern",DMRefinePatternName,PATTERN_COUNT,&pattern,&flgPattern);CHKERRQ(ierr);
      ierr = PetscOptionsGetInt(((PetscObject)dm)->options,prefix,"-dm_p4est_refine_corner",&corner,NULL);CHKERRQ(ierr);
      ierr = PetscOptionsGetIntArray(((PetscObject)dm)->options,prefix,"-dm_p4est_refine_fractal_corners",corners,&ncorner,&flgFractal);CHKERRQ(ierr);
      ierr = PetscOptionsGetReal(((PetscObject)dm)->options,prefix,"-dm_p4est_refine_hash_likelihood",&likelihood,NULL);CHKERRQ(ierr);

      if (flgPattern) {
        DMRefinePatternCtx *ctx;
        PetscInt           maxLevel;

        ierr = DMForestGetMaximumRefinement(dm,&maxLevel);CHKERRQ(ierr);
        ierr = PetscNewLog(dm,&ctx);CHKERRQ(ierr);
        ctx->maxLevel = PetscMin(maxLevel,P4EST_QMAXLEVEL);
        switch (pattern) {
        case PATTERN_HASH:
          ctx->refine_fn      = DMRefinePattern_Hash;
          ctx->hashLikelihood = likelihood;
          break;
        case PATTERN_CORNER:
          ctx->corner    = corner;
          ctx->refine_fn = DMRefinePattern_Corner;
          break;
        case PATTERN_CENTER:
          ctx->refine_fn = DMRefinePattern_Center;
          break;
        case PATTERN_FRACTAL:
          if (flgFractal) {
            PetscInt i;

            for (i = 0; i < ncorner; i++) {
              ctx->fractal[corners[i]] = PETSC_TRUE;
            }
          }
          else {
#if !defined (P4_TO_P8)
            ctx->fractal[0] = ctx->fractal[1] = ctx->fractal[2] = PETSC_TRUE;
#else
            ctx->fractal[0] = ctx->fractal[3] = ctx->fractal[5] = ctx->fractal[6] = PETSC_TRUE;
#endif
          }
          ctx->refine_fn = DMRefinePattern_Fractal;
          break;
        default:
          SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Not a valid refinement pattern\n");
        }

        pforest->forest->user_pointer = (void *) ctx;
        PetscStackCallP4est(p4est_refine,(pforest->forest,1,ctx->refine_fn,NULL));
        PetscStackCallP4est(p4est_balance,(pforest->forest,P4EST_CONNECT_FULL,NULL));
        ierr = PetscFree(ctx);CHKERRQ(ierr);
        pforest->forest->user_pointer = (void *) dm;
      }
    }
    if (initLevel > minLevel) {
      pforest->coarsen_hierarchy = PETSC_TRUE;
    }
  }
  if (pforest->coarsen_hierarchy) {
    PetscInt initLevel, minLevel;

    ierr = DMForestGetInitialRefinement(dm,&initLevel);CHKERRQ(ierr);
    ierr = DMForestGetMinimumRefinement(dm,&minLevel);CHKERRQ(ierr);
    if (initLevel > minLevel) {
      DM_Forest_pforest *coarse_pforest;
      DMLabel           coarsen;
      DM                coarseDM;

      ierr = DMForestTemplate(dm,MPI_COMM_NULL,&coarseDM);CHKERRQ(ierr);
      ierr = DMForestSetAdaptivityPurpose(coarseDM,DM_FOREST_COARSEN);CHKERRQ(ierr);
      ierr = DMCreateLabel(dm,"coarsen");CHKERRQ(ierr);
      ierr = DMGetLabel(dm,"coarsen",&coarsen);CHKERRQ(ierr);
      ierr = DMLabelSetDefaultValue(coarsen,DM_FOREST_COARSEN);CHKERRQ(ierr);
      ierr = DMForestSetAdaptivityLabel(coarseDM,"coarsen");CHKERRQ(ierr);
      ierr = DMSetCoarseDM(dm,coarseDM);CHKERRQ(ierr);
      if (forest->setfromoptionscalled) {
        ierr = DMSetFromOptions(coarseDM);CHKERRQ(ierr);
      }
      ierr = DMForestSetInitialRefinement(coarseDM,initLevel - 1);CHKERRQ(ierr);
      coarse_pforest = (DM_Forest_pforest *) ((DM_Forest *) coarseDM->data)->data;
      coarse_pforest->coarsen_hierarchy = PETSC_TRUE;
      ierr = DMDestroy(&coarseDM);CHKERRQ(ierr);
    }
  }
  { /* repartitioning and overlap */
    PetscMPIInt size, rank;

    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
    if ((size > 1) && (pforest->partition_for_coarsening || forest->cellWeights || forest->weightCapacity != 1. || forest->weightsFactor != 1.)) {
      PetscBool copyForest = PETSC_FALSE;
      p4est_t *forest_copy = NULL;
      if (preCoarseToFine || coarseToPreFine) {
        copyForest = PETSC_TRUE;
      }
      if (copyForest) {
        PetscStackCallP4estReturn(forest_copy,p4est_copy,(pforest->forest,0));
      }

      if (!forest->cellWeights && forest->weightCapacity == 1. && forest->weightsFactor == 1.) {
        PetscStackCallP4est(p4est_partition,(pforest->forest,(int)pforest->partition_for_coarsening,NULL));
      }
      else {
        /* TODO: handle non-uniform partition cases */
        SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_PLIB,"Not implemented yet");
      }
      if (forest_copy) {
        if (preCoarseToFine || coarseToPreFine) {
          PetscSF repartSF; /* repartSF has roots in the old partition */
          PetscInt pStart = -1, pEnd = -1, p;
          PetscInt numRoots, numLeaves;
          PetscSFNode *repartRoots;
          p4est_gloidx_t postStart  = pforest->forest->global_first_quadrant[rank];
          p4est_gloidx_t postEnd    = pforest->forest->global_first_quadrant[rank+1];
          p4est_gloidx_t partOffset = postStart;

          numRoots  = (PetscInt) (forest_copy->global_first_quadrant[rank + 1] - forest_copy->global_first_quadrant[rank]);
          numLeaves = (PetscInt) (postEnd - postStart);
          ierr = DMPforestComputeOverlappingRanks(size,rank,pforest->forest,forest_copy,&pStart,&pEnd);CHKERRQ(ierr);
          ierr = PetscMalloc1((PetscInt) pforest->forest->local_num_quadrants,&repartRoots);CHKERRQ(ierr);
          for (p = pStart; p < pEnd; p++) {
            p4est_gloidx_t preStart = forest_copy->global_first_quadrant[p];
            p4est_gloidx_t preEnd   = forest_copy->global_first_quadrant[p+1];
            PetscInt q;

            if (preEnd == preStart) continue;
            if (preStart > postStart) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Bad partition overlap computation");
            preEnd = preEnd > postEnd ? postEnd : preEnd;
            for (q = partOffset; q < preEnd; q++) {
              repartRoots[q - postStart].rank = p;
              repartRoots[q - postStart].index = partOffset - preStart;
            }
            partOffset = preEnd;
          }
          ierr = PetscSFCreate(PetscObjectComm((PetscObject)dm),&repartSF);CHKERRQ(ierr);
          ierr = PetscSFSetGraph(repartSF,numRoots,numLeaves,NULL,PETSC_OWN_POINTER,repartRoots,PETSC_OWN_POINTER);CHKERRQ(ierr);
          ierr = PetscSFSetUp(repartSF);CHKERRQ(ierr);
          if (preCoarseToFine) {
            PetscSF repartSFembed, preCoarseToFineNew;
            PetscInt nleaves;
            const PetscInt *leaves;

            ierr = PetscSFSetUp(preCoarseToFine);CHKERRQ(ierr);
            ierr = PetscSFGetGraph(preCoarseToFine,NULL,&nleaves,&leaves,NULL);CHKERRQ(ierr);
            if (leaves) {
              ierr = PetscSFCreateEmbeddedSF(repartSF,nleaves,leaves,&repartSFembed);CHKERRQ(ierr);
            }
            else {
              repartSFembed = repartSF;
              ierr = PetscObjectReference((PetscObject)repartSFembed);CHKERRQ(ierr);
            }
            ierr = PetscSFCompose(preCoarseToFine,repartSFembed,&preCoarseToFineNew);CHKERRQ(ierr);
            ierr = PetscSFDestroy(&preCoarseToFine);CHKERRQ(ierr);
            ierr = PetscSFDestroy(&repartSFembed);CHKERRQ(ierr);
            preCoarseToFine = preCoarseToFineNew;
          }
          if (coarseToPreFine) {
            PetscSF repartSFinv, coarseToPreFineNew;

            ierr = PetscSFCreateInverseSF(repartSF,&repartSFinv);CHKERRQ(ierr);
            ierr = PetscSFCompose(repartSFinv,coarseToPreFine,&coarseToPreFineNew);CHKERRQ(ierr);
            ierr = PetscSFDestroy(&coarseToPreFine);CHKERRQ(ierr);
            ierr = PetscSFDestroy(&repartSFinv);CHKERRQ(ierr);
            coarseToPreFine = coarseToPreFineNew;
          }
          ierr = PetscSFDestroy(&repartSF);CHKERRQ(ierr);
        }
        PetscStackCallP4est(p4est_destroy,(forest_copy));
      }
    }
    if (size > 1) {
      PetscInt overlap;

      ierr = DMForestGetPartitionOverlap(dm,&overlap);CHKERRQ(ierr);

      if (overlap > 0) {
        PetscInt i, cLocalStart, cLocalEnd;
        PetscInt cEnd;
        PetscSF  preCellSF, cellSF;

        PetscStackCallP4estReturn(pforest->ghost,p4est_ghost_new,(pforest->forest,P4EST_CONNECT_FULL));
        PetscStackCallP4estReturn(pforest->lnodes,p4est_lnodes_new,(pforest->forest,pforest->ghost,-P4EST_DIM));
        PetscStackCallP4est(p4est_ghost_support_lnodes,(pforest->forest,pforest->lnodes,pforest->ghost));
        for (i = 1; i < overlap; i++) {
          PetscStackCallP4est(p4est_ghost_expand_by_lnodes,(pforest->forest,pforest->lnodes,pforest->ghost));
        }

        cLocalStart = pforest->cLocalStart = pforest->ghost->proc_offsets[rank];
        cLocalEnd = pforest->cLocalEnd = cLocalStart + pforest->forest->local_num_quadrants;
        cEnd = pforest->forest->local_num_quadrants + pforest->ghost->proc_offsets[size];

        /* shift sfs by cLocalStart, expand by cell SFs */
        ierr = DMForestGetCellSF(adaptFrom,&preCellSF);CHKERRQ(ierr);
        ierr = DMForestGetCellSF(dm,&cellSF);CHKERRQ(ierr);
        if (preCoarseToFine) {
          PetscSF preCoarseToFineNew;
          PetscInt nleaves, nroots, *leavesNew, i, nleavesNew;
          const PetscInt *leaves;
          const PetscSFNode *remotes;
          PetscSFNode *remotesAll;

          ierr = PetscSFSetUp(preCoarseToFine);CHKERRQ(ierr);
          ierr = PetscSFGetGraph(preCoarseToFine,&nroots,&nleaves,&leaves,&remotes);CHKERRQ(ierr);
          ierr = PetscMalloc1(cEnd,&remotesAll);CHKERRQ(ierr);
          for (i = 0; i < cEnd; i++) {
            remotesAll[i].rank  = -1;
            remotesAll[i].index = -1;
          }
          for (i = 0; i < nleaves; i++) {
            remotesAll[(leaves ? leaves[i] : i) + cLocalStart] = remotes[i];
          }
          ierr = PetscSFSetUp(cellSF);CHKERRQ(ierr);
          ierr = PetscSFBcastBegin(cellSF,MPIU_2INT,remotesAll,remotesAll);CHKERRQ(ierr);
          ierr = PetscSFBcastEnd(cellSF,MPIU_2INT,remotesAll,remotesAll);CHKERRQ(ierr);
          nleavesNew = 0;
          for (i = 0; i < nleaves; i++) {
            if (remotesAll[i].rank >= 0) {
              nleavesNew++;
            }
          }
          ierr = PetscMalloc1(nleavesNew,&leavesNew);CHKERRQ(ierr);
          nleavesNew = 0;
          for (i = 0; i < nleaves; i++) {
            if (remotesAll[i].rank >= 0) {
              leavesNew[nleavesNew] = i;
              if (i > nleavesNew) {
                remotesAll[nleavesNew] = remotesAll[i];
              }
              nleavesNew++;
            }
          }
          ierr = PetscSFCreate(PetscObjectComm((PetscObject)dm),&preCoarseToFineNew);CHKERRQ(ierr);
          if (nleavesNew < cEnd) {
            ierr = PetscSFSetGraph(preCoarseToFineNew,nroots,nleavesNew,leavesNew,PETSC_OWN_POINTER,remotesAll,PETSC_COPY_VALUES);CHKERRQ(ierr);
          }
          else { /* all cells are leaves */
            ierr = PetscFree(leavesNew);CHKERRQ(ierr);
            ierr = PetscSFSetGraph(preCoarseToFineNew,nroots,nleavesNew,NULL,PETSC_OWN_POINTER,remotesAll,PETSC_COPY_VALUES);CHKERRQ(ierr);
          }
          ierr = PetscFree(remotesAll);CHKERRQ(ierr);
          ierr = PetscSFDestroy(&preCoarseToFine);CHKERRQ(ierr);
          preCoarseToFine = preCoarseToFineNew;
          preCoarseToFine = preCoarseToFineNew;
        }
        if (coarseToPreFine) {
          PetscSF coarseToPreFineNew;
          PetscInt nleaves, nroots, i, nleavesCellSF, nleavesExpanded, *leavesNew;
          const PetscInt *leaves;
          const PetscSFNode *remotes;
          PetscSFNode *remotesNew, *remotesNewRoot, *remotesExpanded;

          ierr = PetscSFSetUp(coarseToPreFine);CHKERRQ(ierr);
          ierr = PetscSFGetGraph(coarseToPreFine,&nroots,&nleaves,&leaves,&remotes);CHKERRQ(ierr);
          ierr = PetscSFGetGraph(preCellSF,NULL,&nleavesCellSF,NULL,NULL);CHKERRQ(ierr);
          ierr = PetscMalloc1(nroots,&remotesNewRoot);CHKERRQ(ierr);
          ierr = PetscMalloc1(nleaves,&remotesNew);CHKERRQ(ierr);
          for (i = 0; i < nroots; i++) {
            remotesNewRoot[i].rank = rank;
            remotesNewRoot[i].index = i + cLocalStart;
          }
          ierr = PetscSFBcastBegin(coarseToPreFine,MPIU_2INT,remotesNewRoot,remotesNew);CHKERRQ(ierr);
          ierr = PetscSFBcastEnd(coarseToPreFine,MPIU_2INT,remotesNewRoot,remotesNew);CHKERRQ(ierr);
          ierr = PetscFree(remotesNewRoot);CHKERRQ(ierr);
          ierr = PetscMalloc1(nleavesCellSF,&remotesExpanded);CHKERRQ(ierr);
          for (i = 0; i < nleavesCellSF; i++) {
            remotesExpanded[i].rank = -1;
            remotesExpanded[i].index = -1;
          }
          for (i = 0; i < nleaves; i++) {
            remotesExpanded[leaves ? leaves[i] : i] = remotesNew[i];
          }
          ierr = PetscFree(remotesNew);CHKERRQ(ierr);
          ierr = PetscSFBcastBegin(preCellSF,MPIU_2INT,remotesExpanded,remotesExpanded);CHKERRQ(ierr);
          ierr = PetscSFBcastEnd(preCellSF,MPIU_2INT,remotesExpanded,remotesExpanded);CHKERRQ(ierr);

          nleavesExpanded = 0;
          for (i = 0; i < nleavesCellSF; i++) {
            if (remotesExpanded[i].rank >= 0) {
              nleavesExpanded++;
            }
          }
          ierr = PetscMalloc1(nleavesExpanded,&leavesNew);CHKERRQ(ierr);
          nleavesExpanded = 0;
          for (i = 0; i < nleavesCellSF; i++) {
            if (remotesExpanded[i].rank >= 0) {
              leavesNew[nleavesExpanded] = i;
              if (i > nleavesExpanded) {
                remotesExpanded[nleavesExpanded] = remotes[i];
              }
              nleavesExpanded++;
            }
          }
          ierr = PetscSFCreate(PetscObjectComm((PetscObject)dm),&coarseToPreFineNew);CHKERRQ(ierr);
          if (nleavesExpanded < nleavesCellSF) {
            ierr = PetscSFSetGraph(coarseToPreFineNew,cEnd,nleavesExpanded,leavesNew,PETSC_OWN_POINTER,remotesExpanded,PETSC_COPY_VALUES);CHKERRQ(ierr);
          }
          else {
            ierr = PetscFree(leavesNew);CHKERRQ(ierr);
            ierr = PetscSFSetGraph(coarseToPreFineNew,cEnd,nleavesExpanded,NULL,PETSC_OWN_POINTER,remotesExpanded,PETSC_COPY_VALUES);CHKERRQ(ierr);
          }
          ierr = PetscFree(remotesExpanded);CHKERRQ(ierr);
          ierr = PetscSFDestroy(&coarseToPreFine);CHKERRQ(ierr);
          coarseToPreFine = coarseToPreFineNew;
        }
      }
    }
  }
  forest->preCoarseToFine = preCoarseToFine;
  forest->coarseToPreFine = coarseToPreFine;
  PetscFunctionReturn(0);
}

#define DMView_ASCII_pforest _append_pforest(DMView_ASCII)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMView_ASCII_pforest)
static PetscErrorCode DMView_ASCII_pforest(PetscObject odm, PetscViewer viewer)
{
  DM                dm       = (DM) odm;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidHeaderSpecific(viewer, PETSC_VIEWER_CLASSID, 2);
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  switch (viewer->format) {
  case PETSC_VIEWER_DEFAULT:
  case PETSC_VIEWER_ASCII_INFO:
  case PETSC_VIEWER_ASCII_INFO_DETAIL:
  {
    PetscInt    dim;
    const char *name;

    ierr = PetscObjectGetName((PetscObject) dm, &name);CHKERRQ(ierr);
    ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
    if (name) {ierr = PetscViewerASCIIPrintf(viewer, "Forest %s in %D dimensions:\n", name, dim);CHKERRQ(ierr);}
    else      {ierr = PetscViewerASCIIPrintf(viewer, "Forest in %D dimensions:\n", dim);CHKERRQ(ierr);}
  }
    break;
  default: SETERRQ1(PetscObjectComm((PetscObject) dm), PETSC_ERR_SUP, "No support for format '%s'", PetscViewerFormats[viewer->format]);
  }
  PetscFunctionReturn(0);
}

#define DMView_VTK_pforest _append_pforest(DMView_VTK)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMView_VTK_pforest)
static PetscErrorCode DMView_VTK_pforest(PetscObject odm, PetscViewer viewer)
{
  DM                dm       = (DM) odm;
  DM_Forest         *forest  = (DM_Forest *) dm->data;
  DM_Forest_pforest *pforest = (DM_Forest_pforest *) forest->data;
  PetscBool         isvtk;
  PetscReal         vtkScale = 1. - PETSC_MACHINE_EPSILON;
  PetscViewer_VTK   *vtk = (PetscViewer_VTK*)viewer->data;
  const char        *name;
  char              *filenameStrip = NULL;
  PetscBool         hasExt;
  size_t            len;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidHeaderSpecific(viewer, PETSC_VIEWER_CLASSID, 2);
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERVTK, &isvtk);CHKERRQ(ierr);
  if (!isvtk) SETERRQ1(PetscObjectComm((PetscObject)viewer), PETSC_ERR_ARG_INCOMP, "Cannot use viewer type %s", ((PetscObject)viewer)->type_name);
  switch (viewer->format) {
  case PETSC_VIEWER_VTK_VTU:
    if (!pforest->forest) SETERRQ (PetscObjectComm(odm),PETSC_ERR_ARG_WRONG,"DM has not been setup with a valid forest");
    name = vtk->filename;
    ierr = PetscStrlen(name,&len);CHKERRQ(ierr);
    ierr = PetscStrcasecmp(name+len-4,".vtu",&hasExt);CHKERRQ(ierr);
    if (hasExt) {
      ierr = PetscStrallocpy(name,&filenameStrip);CHKERRQ(ierr);
      filenameStrip[len-4]='\0';
      name = filenameStrip;
    }
    PetscStackCallP4est(p4est_vtk_write_all,(pforest->forest,pforest->topo->geom,(double)vtkScale,
                                             1, /* write tree */
                                             1, /* write level */
                                             1, /* write rank */
                                             0, /* do not wrap rank */
                                             0, /* no scalar fields */
                                             0, /* no vector fields */
                                             name));
    ierr = PetscFree(filenameStrip);CHKERRQ(ierr);
    break;
  default: SETERRQ1(PetscObjectComm((PetscObject)dm), PETSC_ERR_SUP, "No support for format '%s'", PetscViewerFormats[viewer->format]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPforestGetPlex(DM,DM*);

#define DMView_HDF5_pforest _append_pforest(DMView_HDF5)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMView_HDF5_pforest)
static PetscErrorCode DMView_HDF5_pforest(DM dm, PetscViewer viewer)
{
  DM             plex;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  ierr = DMPforestGetPlex(dm, &plex);CHKERRQ(ierr);
  ierr = DMView(plex, viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMView_pforest _append_pforest(DMView)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMView_pforest)
static PetscErrorCode DMView_pforest(DM dm, PetscViewer viewer)
{
  PetscBool      isascii, isvtk, ishdf5;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidHeaderSpecific(viewer, PETSC_VIEWER_CLASSID, 2);
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERASCII, &isascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERVTK,   &isvtk);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERHDF5,  &ishdf5);CHKERRQ(ierr);
  if (isascii) {
    ierr = DMView_ASCII_pforest((PetscObject) dm,viewer);CHKERRQ(ierr);
  } else if (isvtk) {
    ierr = DMView_VTK_pforest((PetscObject) dm,viewer);CHKERRQ(ierr);
  } else if (ishdf5) {
    ierr = DMView_HDF5_pforest(dm, viewer);CHKERRQ(ierr);
  } else {
    SETERRQ(PetscObjectComm((PetscObject) dm),PETSC_ERR_SUP,"Viewer not supported (not VTK or HDF5)");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PforestConnectivityEnumerateFacets"
static PetscErrorCode PforestConnectivityEnumerateFacets(p4est_connectivity_t *conn, PetscInt **tree_face_to_uniq)
{
  PetscInt       *ttf, f, t, g, count;
  PetscInt       numFacets;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  numFacets = conn->num_trees * P4EST_FACES;
  ierr = PetscMalloc1(numFacets,&ttf);CHKERRQ(ierr);
  for (f = 0; f < numFacets; f++) ttf[f] = -1;
  for (g = 0, count = 0, t = 0; t < conn->num_trees; t++) {
    for (f = 0; f < P4EST_FACES; f++, g++) {
      if (ttf[g] == -1) {
        PetscInt ng;

        ttf[g] = count++;
        ng = conn->tree_to_tree[g] * P4EST_FACES + (conn->tree_to_face[g] % P4EST_FACES);
        ttf[ng] = ttf[g];
      }
    }
  }
  *tree_face_to_uniq = ttf;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMPlexCreateConnectivity_pforest)
static PetscErrorCode DMPlexCreateConnectivity_pforest(DM dm, p4est_connectivity_t **connOut, PetscInt **tree_face_to_uniq)
{
  p4est_topidx_t       numTrees, numVerts, numCorns, numCtt;
  PetscSection         ctt;
#if defined(P4_TO_P8)
  p4est_topidx_t       numEdges, numEtt;
  PetscSection         ett;
  PetscInt             eStart, eEnd, e, ettSize;
  PetscInt             vertOff = 1 + P4EST_FACES + P8EST_EDGES;
  PetscInt             edgeOff = 1 + P4EST_FACES;
#else
  PetscInt             vertOff = 1 + P4EST_FACES;
#endif
  p4est_connectivity_t *conn;
  PetscInt             cStart, cEnd, cEndInterior, c, vStart, vEnd, vEndInterior, v, fStart, fEnd, fEndInterior, f, eEndInterior;
  PetscInt             *star = NULL, *closure = NULL, closureSize, starSize, cttSize;
  PetscInt             *ttf;
  PetscErrorCode       ierr;

  PetscFunctionBegin;
  /* 1: count objects, allocate */
  ierr = DMPlexGetHeightStratum(dm,0,&cStart,&cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(dm, &cEndInterior, &fEndInterior, &eEndInterior, &vEndInterior);CHKERRQ(ierr);
  cEnd = cEndInterior < 0 ? cEnd : cEndInterior;
  if (cEndInterior < 0) {
    cEndInterior = cEnd;
  }
  ierr = P4estTopidxCast(cEnd-cStart,&numTrees);CHKERRQ(ierr);
  numVerts = P4EST_CHILDREN * numTrees;
  ierr = DMPlexGetDepthStratum(dm,0,&vStart,&vEnd);CHKERRQ(ierr);
  vEnd = vEndInterior < 0 ? vEnd : vEndInterior;
  if (vEndInterior < 0) {
    vEndInterior = vEnd;
  }
  ierr = P4estTopidxCast(vEnd-vStart,&numCorns);CHKERRQ(ierr);
  ierr = PetscSectionCreate(PETSC_COMM_SELF,&ctt);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(ctt,vStart,vEnd);CHKERRQ(ierr);
  for (v = vStart; v < vEnd; v++) {
    PetscInt s;

    ierr = DMPlexGetTransitiveClosure(dm,v,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
    for (s = 0; s < starSize; s++) {
      PetscInt p = star[2*s];

      if (p >= cStart && p < cEnd) {
        /* we want to count every time cell p references v, so we see how many times it comes up in the closure.  This
         * only protects against periodicity problems */
        ierr = DMPlexGetTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
        if (closureSize != P4EST_INSUL) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Cell with wrong closure size");
        for (c = 0; c < P4EST_CHILDREN; c++) {
          PetscInt cellVert = closure[2 * (c + vertOff)];

          if (cellVert < vStart || cellVert >= vEnd) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Non-standard closure: vertices\n");
          if (cellVert == v) {
            ierr = PetscSectionAddDof(ctt,v,1);CHKERRQ(ierr);
          }
        }
        ierr = DMPlexRestoreTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
      }
    }
    ierr = DMPlexRestoreTransitiveClosure(dm,v,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(ctt);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(ctt,&cttSize);CHKERRQ(ierr);
  ierr = P4estTopidxCast(cttSize,&numCtt);CHKERRQ(ierr);
#if defined(P4_TO_P8)
  ierr = DMPlexGetDepthStratum(dm,1,&eStart,&eEnd);CHKERRQ(ierr);
  eEnd = eEndInterior < 0 ? eEnd : eEndInterior;
  if (eEndInterior < 0) {
    eEndInterior = eEnd;
  }
  ierr = P4estTopidxCast(eEnd-eStart,&numEdges);CHKERRQ(ierr);
  ierr = PetscSectionCreate(PETSC_COMM_SELF,&ett);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(ett,eStart,eEnd);CHKERRQ(ierr);
  for (e = eStart; e < eEnd; e++) {
    PetscInt s;

    ierr = DMPlexGetTransitiveClosure(dm,e,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
    for (s = 0; s < starSize; s++) {
      PetscInt p = star[2*s];

      if (p >= cStart && p < cEnd) {
        /* we want to count every time cell p references e, so we see how many times it comes up in the closure.  This
         * only protects against periodicity problems */
        ierr = DMPlexGetTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
        if (closureSize != P4EST_INSUL) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Cell with wrong closure size");
        for (c = 0; c < P8EST_EDGES; c++) {
          PetscInt cellEdge = closure[2 * (c + edgeOff)];

          if (cellEdge < eStart || cellEdge >= eEnd) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Non-standard closure: edges\n");
          if (cellEdge == e) {
            ierr = PetscSectionAddDof(ett,e,1);CHKERRQ(ierr);
          }
        }
        ierr = DMPlexRestoreTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
      }
    }
    ierr = DMPlexRestoreTransitiveClosure(dm,e,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(ett);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(ett,&ettSize);CHKERRQ(ierr);
  ierr = P4estTopidxCast(ettSize,&numEtt);CHKERRQ(ierr);

  /* This routine allocates space for the arrays, which we fill below */
  PetscStackCallP4estReturn(conn,p8est_connectivity_new,(numVerts,numTrees,numEdges,numEtt,numCorns,numCtt));
#else
  PetscStackCallP4estReturn(conn,p4est_connectivity_new,(numVerts,numTrees,numCorns,numCtt));
#endif

  /* 2: visit every face, determine neighboring cells(trees) */
  ierr = DMPlexGetHeightStratum(dm,1,&fStart,&fEnd);CHKERRQ(ierr);
  fEnd = fEndInterior < 0 ? fEnd : fEndInterior;
  if (fEndInterior < 0) {
    fEndInterior = fEnd;
  }
  ierr = PetscMalloc1((cEnd-cStart) * P4EST_FACES,&ttf);CHKERRQ(ierr);
  for (f = fStart; f < fEnd; f++) {
    PetscInt numSupp, s;
    PetscInt myFace[2] = {-1, -1};
    PetscInt myOrnt[2] = {PETSC_MIN_INT, PETSC_MIN_INT};
    const PetscInt *supp;

    ierr = DMPlexGetSupportSize(dm, f, &numSupp);CHKERRQ(ierr);
    if (numSupp != 1 && numSupp != 2) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"point %d has facet with %d sides: must be 1 or 2 (boundary or conformal)\n",f,numSupp);
    ierr = DMPlexGetSupport(dm, f, &supp);CHKERRQ(ierr);

    for (s = 0; s < numSupp; s++) {
      PetscInt p = supp[s];

      if (p >= cEnd) {
        numSupp--;
        if (s) {
          supp = &supp[1 - s];
        }
        break;
      }
    }
    for (s = 0; s < numSupp; s++) {
      PetscInt p = supp[s], i;
      PetscInt numCone;
      const PetscInt *cone;
      const PetscInt *ornt;
      PetscInt orient = PETSC_MIN_INT;

      ierr = DMPlexGetConeSize(dm, p, &numCone);CHKERRQ(ierr);
      if (numCone != P4EST_FACES) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"cell %d has %d facets, expect %d\n",p,numCone,P4EST_FACES);
      ierr = DMPlexGetCone(dm, p, &cone);CHKERRQ(ierr);
      ierr = DMPlexGetConeOrientation(dm, p, &ornt);CHKERRQ(ierr);
      for (i = 0; i < P4EST_FACES; i++) {
        if (cone[i] == f) {
          orient = ornt[i];
          break;
        }
      }
      if (i >= P4EST_FACES) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"cell %d faced %d mismatch\n",p,f);
      ttf[P4EST_FACES * (p - cStart) + PetscFaceToP4estFace[i]] = f - fStart;
      if (numSupp == 1) {
        /* boundary faces indicated by self reference */
        conn->tree_to_tree[P4EST_FACES * (p - cStart) + PetscFaceToP4estFace[i]] = p - cStart;
        conn->tree_to_face[P4EST_FACES * (p - cStart) + PetscFaceToP4estFace[i]] = (int8_t) PetscFaceToP4estFace[i];
      }
      else {
        const PetscInt N = P4EST_CHILDREN / 2;

        conn->tree_to_tree[P4EST_FACES * (p - cStart) + PetscFaceToP4estFace[i]] = supp[1 - s] - cStart;
        myFace[s] = PetscFaceToP4estFace[i];
        /* get the orientation of cell p in p4est-type closure to facet f, by composing the p4est-closure to
         * petsc-closure permutation and the petsc-closure to facet orientation */
        myOrnt[s] = DihedralCompose(N,orient,P4estFaceToPetscOrnt[myFace[s]]);
      }
    }
    if (numSupp == 2) {
      for (s = 0; s < numSupp; s++) {
        PetscInt p = supp[s];
        PetscInt orntAtoB;
        PetscInt p4estOrient;
        const PetscInt N = P4EST_CHILDREN / 2;

        /* composing the forward permutation with the other cell's inverse permutation gives the self-to-neighbor
         * permutation of this cell-facet's cone */
        orntAtoB = DihedralCompose(N,DihedralInvert(N,myOrnt[1-s]),myOrnt[s]);

        /* convert cone-description permutation (i.e., edges around facet) to cap-description permutation (i.e.,
         * vertices around facet) */
#if !defined(P4_TO_P8)
        p4estOrient = orntAtoB < 0 ? -(orntAtoB + 1) : orntAtoB;
#else
        {
          PetscInt firstVert = orntAtoB < 0 ? ((-orntAtoB) % N): orntAtoB;
          PetscInt p4estFirstVert = firstVert < 2 ? firstVert : (firstVert ^ 1);

                                                                                           /* swap bits */
          p4estOrient = ((myFace[s] <= myFace[1 - s]) || (orntAtoB < 0)) ? p4estFirstVert : ((p4estFirstVert >> 1) | ((p4estFirstVert & 1) << 1));
        }
#endif
        /* encode neighbor face and orientation in tree_to_face per p4est_connectivity standard (see
         * p4est_connectivity.h, p8est_connectivity.h) */
        conn->tree_to_face[P4EST_FACES * (p - cStart) + myFace[s]] = (int8_t) myFace[1 - s] + p4estOrient * P4EST_FACES;
      }
    }
  }

#if defined(P4_TO_P8)
  /* 3: visit every edge */
  conn->ett_offset[0] = 0;
  for (e = eStart; e < eEnd; e++) {
    PetscInt off, s;

    ierr = PetscSectionGetOffset(ett,e,&off);CHKERRQ(ierr);
    conn->ett_offset[e - eStart] = (p4est_topidx_t) off;
    ierr = DMPlexGetTransitiveClosure(dm,e,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
    for (s = 0; s < starSize; s++) {
      PetscInt p = star[2 * s];

      if (p >= cStart && p < cEnd) {
        ierr = DMPlexGetTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
        if (closureSize != P4EST_INSUL) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Non-standard closure\n");
        for (c = 0; c < P8EST_EDGES; c++) {
          PetscInt cellEdge = closure[2 * (c + edgeOff)];
          PetscInt cellOrnt = closure[2 * (c + edgeOff) + 1];

          if (cellEdge == e) {
            PetscInt p4estEdge = PetscEdgeToP4estEdge[c];
            PetscInt totalOrient;

            /* compose p4est-closure to petsc-closure permutation and petsc-closure to edge orientation */
            totalOrient = DihedralCompose(2,cellOrnt,P4estEdgeToPetscOrnt[p4estEdge]);
            /* p4est orientations are positive: -2 => 1, -1 => 0 */
            totalOrient = (totalOrient < 0) ? -(totalOrient + 1) : totalOrient;
            conn->edge_to_tree[off] = (p4est_locidx_t) (p - cStart);
            /* encode cell-edge and orientation in edge_to_edge per p8est_connectivity standart (see
             * p8est_connectivity.h) */
            conn->edge_to_edge[off++] = (int8_t) p4estEdge + P8EST_EDGES * totalOrient;
            conn->tree_to_edge[P8EST_EDGES * (p - cStart) + p4estEdge] = e - eStart;
          }
        }
        ierr = DMPlexRestoreTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
      }
    }
    ierr = DMPlexRestoreTransitiveClosure(dm,e,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
  }
  ierr = PetscSectionDestroy(&ett);CHKERRQ(ierr);
#endif

  /* 4: visit every vertex */
  conn->ctt_offset[0] = 0;
  for (v = vStart; v < vEnd; v++) {
    PetscInt off, s;

    ierr = PetscSectionGetOffset(ctt,v,&off);CHKERRQ(ierr);
    conn->ctt_offset[v - vStart] = (p4est_topidx_t) off;
    ierr = DMPlexGetTransitiveClosure(dm,v,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
    for (s = 0; s < starSize; s++) {
      PetscInt p = star[2 * s];

      if (p >= cStart && p < cEnd) {
        ierr = DMPlexGetTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
        if (closureSize != P4EST_INSUL) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Non-standard closure\n");
        for (c = 0; c < P4EST_CHILDREN; c++) {
          PetscInt cellVert = closure[2 * (c + vertOff)];

          if (cellVert == v) {
            PetscInt p4estVert = PetscVertToP4estVert[c];

            conn->corner_to_tree[off] = (p4est_locidx_t) (p - cStart);
            conn->corner_to_corner[off++] = (int8_t) p4estVert;
            conn->tree_to_corner[P4EST_CHILDREN * (p - cStart) + p4estVert] = v - vStart;
          }
        }
        ierr = DMPlexRestoreTransitiveClosure(dm,p,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
      }
    }
    ierr = DMPlexRestoreTransitiveClosure(dm,v,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
  }
  ierr = PetscSectionDestroy(&ctt);CHKERRQ(ierr);

  /* 5: Compute the coordinates */
  {
    PetscInt     coordDim;
    Vec          coordVec;
    PetscSection coordSec;

    ierr = DMGetCoordinateDim(dm,&coordDim);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm,&coordVec);CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(dm,&coordSec);CHKERRQ(ierr);
    for (c = cStart; c < cEnd; c++) {
      PetscInt    dof;
      PetscScalar *cellCoords = NULL;

      ierr = DMPlexVecGetClosure(dm, coordSec, coordVec, c, &dof, &cellCoords);CHKERRQ(ierr);
      if (!dm->maxCell && dof != P4EST_CHILDREN * coordDim) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Need coordinates at the corners\n");
      for (v = 0; v < P4EST_CHILDREN; v++) {
        PetscInt i, lim = PetscMin(3, coordDim);
        PetscInt p4estVert = PetscVertToP4estVert[v];

        conn->tree_to_vertex[P4EST_CHILDREN * (c - cStart) + v] = P4EST_CHILDREN * (c - cStart) + v;
        /* p4est vertices are always embedded in R^3 */
        for (i = 0; i < 3; i++) {
          conn->vertices[3 * (P4EST_CHILDREN * (c - cStart) + p4estVert) + i] = 0.;
        }
        for (i = 0; i < lim; i++) {
          conn->vertices[3 * (P4EST_CHILDREN * (c - cStart) + p4estVert) + i] = PetscRealPart(cellCoords[v * coordDim + i]);
        }
      }
      ierr = DMPlexVecRestoreClosure(dm, coordSec, coordVec, c, &dof, &cellCoords);CHKERRQ(ierr);
    }
  }

#if defined(P4EST_ENABLE_DEBUG)
  if (!p4est_connectivity_is_valid(conn)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Plex to p4est conversion failed\n");
#endif

  *connOut = conn;

  *tree_face_to_uniq = ttf;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "locidx_to_PetscInt"
static PetscErrorCode locidx_to_PetscInt (sc_array_t * array)
{
  sc_array_t         *newarray;
  size_t              zz, count = array->elem_count;

  PetscFunctionBegin;
  if (array->elem_size != sizeof (p4est_locidx_t)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Wrong locidx size");

  if (sizeof (p4est_locidx_t) == sizeof (PetscInt)) {
    PetscFunctionReturn(0);
  }

  newarray = sc_array_new_size (sizeof (PetscInt), array->elem_count);
  for (zz = 0; zz < count; zz++) {
    p4est_locidx_t      il = *((p4est_locidx_t *) sc_array_index (array, zz));
    PetscInt           *ip = (PetscInt *) sc_array_index (newarray, zz);

    *ip = (PetscInt) il;
  }

  sc_array_reset (array);
  sc_array_init_size (array, sizeof (PetscInt), count);
  sc_array_copy (array, newarray);
  sc_array_destroy (newarray);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "coords_double_to_PetscScalar"
static PetscErrorCode coords_double_to_PetscScalar (sc_array_t * array, PetscInt dim)
{
  sc_array_t         *newarray;
  size_t              zz, count = array->elem_count;

  PetscFunctionBegin;
  if (array->elem_size != 3 * sizeof (double)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Wrong coordinate size");

  if (sizeof (double) == sizeof (PetscScalar) && dim == 3) {
    PetscFunctionReturn(0);
  }

  newarray = sc_array_new_size (dim * sizeof (PetscScalar), array->elem_count);
  for (zz = 0; zz < count; zz++) {
    int                 i;
    double             *id = (double *) sc_array_index (array, zz);
    PetscScalar        *ip = (PetscScalar *) sc_array_index (newarray, zz);

    for (i = 0; i < dim; i++) {
      ip[i] = 0.;
    }
    for (i = 0; i < PetscMin(dim,3); i++) {
      ip[i] = (PetscScalar) id[i];
    }
  }

  sc_array_reset (array);
  sc_array_init_size (array, dim * sizeof (PetscScalar), count);
  sc_array_copy (array, newarray);
  sc_array_destroy (newarray);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "locidx_pair_to_PetscSFNode"
static PetscErrorCode locidx_pair_to_PetscSFNode (sc_array_t * array)
{
  sc_array_t         *newarray;
  size_t              zz, count = array->elem_count;

  PetscFunctionBegin;
  if (array->elem_size != 2 * sizeof (p4est_locidx_t)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Wrong locidx size");

  newarray = sc_array_new_size (sizeof (PetscSFNode), array->elem_count);
  for (zz = 0; zz < count; zz++) {
    p4est_locidx_t     *il = (p4est_locidx_t *) sc_array_index (array, zz);
    PetscSFNode        *ip = (PetscSFNode *) sc_array_index (newarray, zz);

    ip->rank = (PetscInt) il[0];
    ip->index = (PetscInt) il[1];
  }

  sc_array_reset (array);
  sc_array_init_size (array, sizeof (PetscSFNode), count);
  sc_array_copy (array, newarray);
  sc_array_destroy (newarray);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "P4estToPlex_Local"
static PetscErrorCode P4estToPlex_Local(p4est_t *p4est, DM * plex)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  {
    sc_array_t     *points_per_dim    = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *cone_sizes        = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *cones             = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *cone_orientations = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *coords            = sc_array_new(3 * sizeof (double));
    sc_array_t     *children          = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *parents           = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *childids          = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *leaves            = sc_array_new(sizeof (p4est_locidx_t));
    sc_array_t     *remotes           = sc_array_new(2 * sizeof (p4est_locidx_t));
    p4est_locidx_t first_local_quad;

    PetscStackCallP4est(p4est_get_plex_data,(p4est,P4EST_CONNECT_FULL,0,&first_local_quad,points_per_dim,cone_sizes,cones,cone_orientations,coords,children,parents,childids,leaves,remotes));

    ierr = locidx_to_PetscInt(points_per_dim);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(cone_sizes);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(cones);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(cone_orientations);CHKERRQ(ierr);
    ierr = coords_double_to_PetscScalar(coords, P4EST_DIM);CHKERRQ(ierr);

    ierr = DMPlexCreate(PETSC_COMM_SELF,plex);CHKERRQ(ierr);
    ierr = DMSetDimension(*plex,P4EST_DIM);CHKERRQ(ierr);
    ierr = DMPlexCreateFromDAG(*plex,P4EST_DIM,(PetscInt *)points_per_dim->array,(PetscInt *)cone_sizes->array,(PetscInt *)cones->array,(PetscInt *)cone_orientations->array,(PetscScalar *)coords->array);CHKERRQ(ierr);
    sc_array_destroy (points_per_dim);
    sc_array_destroy (cone_sizes);
    sc_array_destroy (cones);
    sc_array_destroy (cone_orientations);
    sc_array_destroy (coords);
    sc_array_destroy (children);
    sc_array_destroy (parents);
    sc_array_destroy (childids);
    sc_array_destroy (leaves);
    sc_array_destroy (remotes);
  }
  PetscFunctionReturn(0);
}

#define DMReferenceTreeGetChildSymmetry_pforest _append_pforest(DMReferenceTreeGetChildSymmetry)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMReferenceTreeGetChildSymmetry_pforest)
static PetscErrorCode DMReferenceTreeGetChildSymmetry_pforest(DM dm, PetscInt parent, PetscInt parentOrientA, PetscInt childOrientA, PetscInt childA, PetscInt parentOrientB, PetscInt *childOrientB,
                                                              PetscInt *childB)
{
  PetscInt       coneSize, dStart, dEnd, vStart, vEnd, dim, ABswap, oAvert, oBvert, ABswapVert;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (parentOrientA == parentOrientB) {
    if (childOrientB) *childOrientB = childOrientA;
    if (childB) *childB = childA;
    PetscFunctionReturn(0);
  }
  ierr = DMPlexGetDepthStratum(dm,0,&vStart,&vEnd);CHKERRQ(ierr);
  if (childA >= vStart && childA < vEnd) { /* vertices (always in the middle) are invarient under rotation */
    if (childOrientB) *childOrientB = 0;
    if (childB) *childB = childA;
    PetscFunctionReturn(0);
  }
  for (dim = 0; dim < 3; dim++) {
    ierr = DMPlexGetDepthStratum(dm,dim,&dStart,&dEnd);CHKERRQ(ierr);
    if (parent >= dStart && parent <= dEnd) {
      break;
    }
  }
  if (dim > 2) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot perform child symmetry for %d-cells",dim);
  if (!dim) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"A vertex has no children");
  if (childA < dStart || childA >= dEnd) { /* a 1-cell in a 2-cell */
    /* this is a lower-dimensional child: bootstrap */
    PetscInt size, i, sA = -1, sB, sOrientB, sConeSize;
    const PetscInt *supp, *coneA, *coneB, *oA, *oB;

    ierr = DMPlexGetSupportSize(dm,childA,&size);CHKERRQ(ierr);
    ierr = DMPlexGetSupport(dm,childA,&supp);CHKERRQ(ierr);

    /* find a point sA in supp(childA) that has the same parent */
    for (i = 0; i < size; i++) {
      PetscInt sParent;

      sA   = supp[i];
      if (sA == parent) continue;
      ierr = DMPlexGetTreeParent(dm,sA,&sParent,NULL);CHKERRQ(ierr);
      if (sParent == parent) {
        break;
      }
    }
    if (i == size) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"could not find support in children");
    /* find out which point sB is in an equivalent position to sA under
     * parentOrientB */
    ierr = DMReferenceTreeGetChildSymmetry_pforest(dm,parent,parentOrientA,0,sA,parentOrientB,&sOrientB,&sB);CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm,sA,&sConeSize);CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm,sA,&coneA);CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm,sB,&coneB);CHKERRQ(ierr);
    ierr = DMPlexGetConeOrientation(dm,sA,&oA);CHKERRQ(ierr);
    ierr = DMPlexGetConeOrientation(dm,sB,&oB);CHKERRQ(ierr);
    /* step through the cone of sA in natural order */
    for (i = 0; i < sConeSize; i++) {
      if (coneA[i] == childA) {
        /* if childA is at position i in coneA,
         * then we want the point that is at sOrientB*i in coneB */
        PetscInt j = (sOrientB >= 0) ? ((sOrientB + i) % sConeSize) : ((sConeSize -(sOrientB+1) - i) % sConeSize);
        if (childB) *childB = coneB[j];
        if (childOrientB) {
          PetscInt oBtrue;

          ierr          = DMPlexGetConeSize(dm,childA,&coneSize);CHKERRQ(ierr);
          /* compose sOrientB and oB[j] */
          if (coneSize != 0 && coneSize != 2) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Expected a vertex or an edge");
          /* we may have to flip an edge */
          oBtrue        = coneSize ? ((sOrientB >= 0) ? oB[j] : -(oB[j] + 2)) : 0;
          ABswap        = DihedralSwap(coneSize,oA[i],oBtrue);CHKERRQ(ierr);
          *childOrientB = DihedralCompose(coneSize,childOrientA,ABswap);
        }
        break;
      }
    }
    if (i == sConeSize) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"support cone mismatch");
    PetscFunctionReturn(0);
  }
  /* get the cone size and symmetry swap */
  ierr   = DMPlexGetConeSize(dm,parent,&coneSize);CHKERRQ(ierr);
  ABswap = DihedralSwap(coneSize, parentOrientA, parentOrientB);
  if (dim == 2) {
    /* orientations refer to cones: we want them to refer to vertices:
     * if it's a rotation, they are the same, but if the order is reversed, a
     * permutation that puts side i first does *not* put vertex i first */
    oAvert     = (parentOrientA >= 0) ? parentOrientA : -((-parentOrientA % coneSize) + 1);
    oBvert     = (parentOrientB >= 0) ? parentOrientB : -((-parentOrientB % coneSize) + 1);
    ABswapVert = DihedralSwap(coneSize, oAvert, oBvert);
  }
  else {
    oAvert     = parentOrientA;
    oBvert     = parentOrientB;
    ABswapVert = ABswap;
  }
  if (childB) {
    /* assume that each child corresponds to a vertex, in the same order */
    PetscInt p, posA = -1, numChildren, i;
    const PetscInt *children;

    /* count which position the child is in */
    ierr = DMPlexGetTreeChildren(dm,parent,&numChildren,&children);CHKERRQ(ierr);
    for (i = 0; i < numChildren; i++) {
      p = children[i];
      if (p == childA) {
        if (dim == 1) {
          posA = i;
        }
        else { /* 2D Morton to rotation */
          posA = (i & 2) ? (i ^ 1) : i;
        }
        break;
      }
    }
    if (posA >= coneSize) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Could not find childA in children of parent");
    }
    else {
      /* figure out position B by applying ABswapVert */
      PetscInt posB, childIdB;

      posB = (ABswapVert >= 0) ? ((ABswapVert + posA) % coneSize) : ((coneSize -(ABswapVert + 1) - posA) % coneSize);
      if (dim == 1) {
        childIdB = posB;
      }
      else { /* 2D rotation to Morton */
        childIdB = (posB & 2) ? (posB ^ 1) : posB;
      }
      if (childB) *childB = children[childIdB];
    }
  }
  if (childOrientB) *childOrientB = DihedralCompose(coneSize,childOrientA,ABswap);
  PetscFunctionReturn(0);
}

#define DMCreateReferenceTree_pforest _append_pforest(DMCreateReferenceTree)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateReferenceTree_pforest)
static PetscErrorCode DMCreateReferenceTree_pforest(MPI_Comm comm, DM *dm)
{
  p4est_connectivity_t *refcube;
  p4est_t        *root, *refined;
  DM             dmRoot, dmRefined;
  DM_Plex        *mesh;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscStackCallP4estReturn(refcube,p4est_connectivity_new_byname,("unit"));
  { /* [-1,1]^d geometry */
    PetscInt i, j;

    for(i = 0; i < P4EST_CHILDREN; i++) {
      for (j = 0; j < 3; j++) {
        refcube->vertices[3 * i + j] *= 2.;
        refcube->vertices[3 * i + j] -= 1.;
      }
    }
  }
  PetscStackCallP4estReturn(root,p4est_new,(PETSC_COMM_SELF,refcube,0,NULL,NULL));
  PetscStackCallP4estReturn(refined,p4est_new_ext,(PETSC_COMM_SELF,refcube,0,1,1,0,NULL,NULL));
  ierr = P4estToPlex_Local(root,&dmRoot);CHKERRQ(ierr);
  ierr = P4estToPlex_Local(refined,&dmRefined);CHKERRQ(ierr);
  {
#if !defined(P4_TO_P8)
    PetscInt nPoints = 25;
    PetscInt perm[25] = {0, 1, 2, 3,
                          4, 12, 8, 14,
                              6, 9, 15,
                          5, 13,    10,
                              7,    11,
                         16, 22, 20, 24,
                             17,     21,
                                 18, 23,
                                     19};
    PetscInt ident[25] = {0, 0, 0, 0,
                          1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 0, 0,
                          5, 6, 7, 8, 1, 2, 3, 4, 0};
#else
    PetscInt nPoints = 125;
    PetscInt perm[125] = {0, 1, 2, 3, 4, 5, 6, 7,
                           8, 32, 16, 36, 24, 40,
                              12, 17, 37, 25, 41,
                           9, 33,     20, 26, 42,
                              13,     21, 27, 43,
                          10, 34, 18, 38,     28,
                              14, 19, 39,     29,
                          11, 35,     22,     30,
                              15,     23,     31,
                          44, 84, 76, 92, 52, 86, 68, 94, 60, 78, 70, 96,
                          45, 85, 77, 93,     54,     72,     62,     74,
                              46,     80, 53, 87, 69, 95,         64, 82,
                              47,     81,     55,     73,             66,
                                  48, 88,         56, 90, 61, 79, 71, 97,
                                  49, 89,             58,     63,     75,
                                      50,         57, 91,         65, 83,
                                      51,             59,             67,
                           98, 106, 110, 122, 114, 120, 118, 124,
                                99,      111,      115,      119,
                                    100, 107,           116, 121,
                                         101,                117,
                                              102, 108, 112, 123,
                                                   103,      113,
                                                        104, 109,
                                                             105};
    PetscInt ident[125] = {0, 0, 0, 0, 0, 0, 0, 0,
                           1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18,
                           1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6,
                           0, 0, 0, 0, 0, 0,
                           19, 20, 21, 22, 23, 24, 25, 26,
                           7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                           1, 2, 3, 4, 5, 6,
                           0};

#endif
    IS permIS;
    DM dmPerm;

    ierr = ISCreateGeneral(PETSC_COMM_SELF,nPoints,perm,PETSC_USE_POINTER,&permIS);CHKERRQ(ierr);
    ierr = DMPlexPermute(dmRefined,permIS,&dmPerm);CHKERRQ(ierr);
    if (dmPerm) {
      ierr = DMDestroy(&dmRefined);CHKERRQ(ierr);
      dmRefined = dmPerm;
    }
    ierr = ISDestroy(&permIS);CHKERRQ(ierr);
    {
      PetscInt p;
      ierr = DMCreateLabel(dmRoot,"identity");CHKERRQ(ierr);
      ierr = DMCreateLabel(dmRefined,"identity");CHKERRQ(ierr);
      for (p = 0; p < P4EST_INSUL; p++) {
        DMSetLabelValue(dmRoot,"identity",p,p);CHKERRQ(ierr);
      }
      for (p = 0; p < nPoints; p++) {
        DMSetLabelValue(dmRefined,"identity",p,ident[p]);CHKERRQ(ierr);
      }
    }
  }
  ierr = DMPlexCreateReferenceTree_Union(dmRoot,dmRefined,"identity",dm);CHKERRQ(ierr);
  mesh = (DM_Plex *) (*dm)->data;
  mesh->getchildsymmetry = DMReferenceTreeGetChildSymmetry_pforest;
  ierr = DMDestroy(&dmRefined);CHKERRQ(ierr);
  ierr = DMDestroy(&dmRoot);CHKERRQ(ierr);
  PetscStackCallP4est(p4est_destroy,(refined));
  PetscStackCallP4est(p4est_destroy,(root));
  PetscStackCallP4est(p4est_connectivity_destroy,(refcube));
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMShareDiscretization"
static PetscErrorCode DMShareDiscretization(DM dmA, DM dmB)
{
  PetscDS        ds, dsB;
  PetscBool      newDS;
  void           *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dmA,&ctx);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dmB,ctx);CHKERRQ(ierr);
  ierr = DMGetDS(dmA,&ds);CHKERRQ(ierr);
  ierr = DMGetDS(dmB,&dsB);CHKERRQ(ierr);
  newDS = (PetscBool) (ds != dsB);
  ierr = DMSetDS(dmB,ds);CHKERRQ(ierr);
  if (newDS) {
    ierr = DMClearGlobalVectors(dmB);CHKERRQ(ierr);
    ierr = DMClearLocalVectors(dmB);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)dmA->defaultSection);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&(dmB->defaultSection));CHKERRQ(ierr);
    dmB->defaultSection = dmA->defaultSection;
    ierr = PetscObjectReference((PetscObject)dmA->defaultConstraintSection);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&(dmB->defaultConstraintSection));CHKERRQ(ierr);
    dmB->defaultConstraintSection = dmA->defaultConstraintSection;
    ierr = PetscObjectReference((PetscObject)dmA->defaultConstraintMat);CHKERRQ(ierr);
    ierr = MatDestroy(&(dmB->defaultConstraintMat));CHKERRQ(ierr);
    dmB->defaultConstraintMat = dmA->defaultConstraintMat;
    ierr = PetscObjectReference((PetscObject)dmA->defaultGlobalSection);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&(dmB->defaultGlobalSection));CHKERRQ(ierr);
    dmB->defaultGlobalSection = dmA->defaultGlobalSection;
    ierr = PetscObjectReference((PetscObject)dmA->defaultSF);CHKERRQ(ierr);
    ierr = PetscSFDestroy(&dmB->defaultSF);CHKERRQ(ierr);
    dmB->defaultSF = dmA->defaultSF;
    if (dmA->map) {ierr = PetscLayoutReference(dmA->map,&dmB->map);CHKERRQ(ierr);}
  }
  dmA->boundary->refct++;
  ierr = DMBoundaryDestroy(&(dmB->boundary));CHKERRQ(ierr);
  dmB->boundary = dmA->boundary;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestGetCellCoveringSF"
/* Get an SF that broadcasts a coarse-cell covering of the local fine cells */
static PetscErrorCode DMPforestGetCellCoveringSF(MPI_Comm comm,p4est_t *p4estC, p4est_t *p4estF, PetscInt cStart, PetscInt cEnd, PetscInt cLocalStart, PetscSF *coveringSF)
{
  PetscInt       startF, endF, startC, endC, p, nLeaves;
  PetscSFNode    *leaves;
  PetscSF        sf;
  PetscInt       *recv, *send;
  PetscMPIInt    tag;
  MPI_Request    *recvReqs, *sendReqs;
  PetscSection   section;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPforestComputeOverlappingRanks(p4estC->mpisize,p4estC->mpirank,p4estF,p4estC,&startC,&endC);CHKERRQ(ierr);
  ierr = PetscMalloc2(2*(endC-startC),&recv,endC-startC,&recvReqs);CHKERRQ(ierr);
  ierr = PetscCommGetNewTag(comm,&tag);CHKERRQ(ierr);
  for (p = startC; p < endC; p++) {
    recvReqs[p-startC] = MPI_REQUEST_NULL; /* just in case we don't initiate a receive */
    if (p4estC->global_first_quadrant[p] == p4estC->global_first_quadrant[p+1]) { /* empty coarse partition */
      recv[2*(p-startC)] = 0;
      recv[2*(p-startC)+1] = 0;
      continue;
    }

    ierr = MPI_Irecv(&recv[2*(p-startC)],2,MPIU_INT,p,tag,comm,&recvReqs[p-startC]);CHKERRQ(ierr);
  }
  ierr = DMPforestComputeOverlappingRanks(p4estC->mpisize,p4estC->mpirank,p4estC,p4estF,&startF,&endF);CHKERRQ(ierr);
  ierr = PetscMalloc2(2*(endF-startF),&send,endF-startF,&sendReqs);CHKERRQ(ierr);
  ierr = PetscSFCreate(p4estC->mpicomm,&sf);CHKERRQ(ierr);
  /* count the quadrants rank will send to each of [startF,endF) */
  for (p = startF; p < endF; p++) {
    p4est_quadrant_t *myFineStart = &p4estF->global_first_position[p];
    p4est_quadrant_t *myFineEnd   = &p4estF->global_first_position[p+1];
    PetscInt tStart = (PetscInt) myFineStart->p.which_tree;
    PetscInt tEnd   = (PetscInt) myFineEnd->p.which_tree;
    PetscInt firstCell = -1, lastCell = -1;
    p4est_tree_t *treeStart = &(((p4est_tree_t *) p4estC->trees->array)[tStart]);
    p4est_tree_t *treeEnd   = (size_t) tEnd < p4estC->trees->elem_count ? &(((p4est_tree_t *) p4estC->trees->array)[tEnd]) : NULL;
    ssize_t overlapIndex;

    sendReqs[p-startF] = MPI_REQUEST_NULL; /* just in case we don't initiate a send */
    if (p4estF->global_first_quadrant[p] == p4estF->global_first_quadrant[p+1]) {
      continue;
    }

    /* locate myFineStart in (or before) a cell */
    if (treeStart->quadrants.elem_count) {
      PetscStackCallP4estReturn(overlapIndex,sc_array_bsearch,(&(treeStart->quadrants),myFineStart,p4est_quadrant_disjoint));
      if (overlapIndex < 0) {
        firstCell = 0;
      }
      else {
        firstCell = treeStart->quadrants_offset + overlapIndex;
      }
    }
    else {
      firstCell = 0;
    }
    if (treeEnd && treeEnd->quadrants.elem_count) {
      PetscStackCallP4estReturn(overlapIndex,sc_array_bsearch,(&(treeEnd->quadrants),myFineEnd,p4est_quadrant_disjoint));
      if (overlapIndex < 0) { /* all of this local section is overlapped */
        lastCell = p4estC->local_num_quadrants;
      }
      else {
        p4est_quadrant_t *container = &(((p4est_quadrant_t *) treeEnd->quadrants.array)[overlapIndex]);
        p4est_quadrant_t first_desc;
        int equal;

        PetscStackCallP4est(p4est_quadrant_first_descendant,(container,&first_desc,P4EST_QMAXLEVEL));
        PetscStackCallP4estReturn(equal,p4est_quadrant_is_equal,(myFineEnd,&first_desc));
        if (equal) {
          lastCell = treeEnd->quadrants_offset + overlapIndex;
        }
        else {
          lastCell = treeEnd->quadrants_offset + overlapIndex + 1;
        }
      }
    }
    else {
      lastCell = p4estC->local_num_quadrants;
    }
    send[2*(p-startF)] = firstCell + cLocalStart;
    send[2*(p-startF)+1] = lastCell - firstCell;
    ierr = MPI_Isend(&send[2*(p-startF)],2,MPIU_INT,p,tag,comm,&sendReqs[p-startF]);CHKERRQ(ierr);
  }
  ierr = MPI_Waitall((PetscMPIInt)(endC-startC),recvReqs,NULL);CHKERRQ(ierr);
  ierr = PetscSectionCreate(PETSC_COMM_SELF,&section);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(section,startC,endC);CHKERRQ(ierr);
  for (p = startC; p < endC; p++) {
    PetscInt numCells = recv[2*(p-startC)+1];
    ierr = PetscSectionSetDof(section,p,numCells);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(section);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(section,&nLeaves);CHKERRQ(ierr);
  ierr = PetscMalloc1(nLeaves,&leaves);CHKERRQ(ierr);
  for (p = startC; p < endC; p++) {
    PetscInt firstCell = recv[2*(p-startC)];
    PetscInt numCells = recv[2*(p-startC)+1];
    PetscInt off, i;

    ierr = PetscSectionGetOffset(section,p,&off);CHKERRQ(ierr);
    for (i = 0; i < numCells; i++) {
      leaves[off+i].rank = p;
      leaves[off+i].index = firstCell + i;
    }
  }
  ierr = PetscSFCreate(comm,&sf);CHKERRQ(ierr);
  ierr = PetscSFSetGraph(sf,cEnd-cStart,nLeaves,NULL,PETSC_OWN_POINTER,leaves,PETSC_OWN_POINTER);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&section);CHKERRQ(ierr);
  ierr = MPI_Waitall((PetscMPIInt)(endF-startF),sendReqs,NULL);CHKERRQ(ierr);
  ierr = PetscFree2(send,sendReqs);CHKERRQ(ierr);
  ierr = PetscFree2(recv,recvReqs);CHKERRQ(ierr);
  *coveringSF = sf;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestGetCellSFNodes"
/* closure points for locally-owned cells */
static PetscErrorCode DMPforestGetCellSFNodes(DM dm, PetscInt numClosureIndices, const PetscInt closureIndices[], PetscInt *numClosurePoints, PetscSFNode **closurePoints,PetscBool redirect)
{
  PetscInt          cStart, cEnd, cEndInterior;
  PetscInt          count, c;
  PetscMPIInt       rank;
  PetscInt          closureSize = -1;
  PetscInt          *closure = NULL;
  PetscSF           pointSF;
  PetscInt          nleaves, nroots;
  const PetscInt    *ilocal;
  const PetscSFNode *iremote;
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(plex,0,&cStart,&cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHybridBounds(plex,&cEndInterior,NULL,NULL,NULL);CHKERRQ(ierr);
  cEnd = cEndInterior < 0 ? cEnd : cEndInterior;
  ierr = DMGetPointSF(dm,&pointSF);CHKERRQ(ierr);
  ierr = PetscSFGetGraph(pointSF,&nroots,&nleaves,&ilocal,&iremote);CHKERRQ(ierr);
  nleaves = PetscMax(0,nleaves);
  nroots = PetscMax(0,nroots);
  *numClosurePoints = numClosureIndices * (cEnd - cStart);
  ierr = PetscMalloc1(*numClosurePoints,closurePoints);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
  for (c = cStart, count = 0; c < cEnd; c++) {
    PetscInt i;
    ierr = DMPlexGetTransitiveClosure(plex,c,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);

    for (i = 0; i < numClosureIndices; i++, count++) {
      PetscInt j = closureIndices[i];
      PetscInt p = closure[2 * j];
      PetscInt loc = -1;

      ierr = PetscFindInt(p,nleaves,ilocal,&loc);CHKERRQ(ierr);
      if (redirect && loc >= 0) {
        (*closurePoints)[count].rank  = iremote[loc].rank;
        (*closurePoints)[count].index = iremote[loc].index;
      }
      else {
        (*closurePoints)[count].rank  = rank;
        (*closurePoints)[count].index = p;
      }
    }
    ierr = DMPlexRestoreTransitiveClosure(plex,c,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static void DMPforestMaxSFNode(void *a, void *b, PetscMPIInt *len, MPI_Datatype *type)
{
  PetscMPIInt i;

  for (i = 0; i < *len; i++) {
    PetscSFNode *A = (PetscSFNode *)a;
    PetscSFNode *B = (PetscSFNode *)b;

    if (B->rank < 0) {
      *B = *A;
    }
  }
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestGetTransferSF_Internal"
/* children are sf leaves of parents */
static PetscErrorCode DMPforestGetTransferSF_Internal(DM coarse, DM fine, const PetscInt dofPerDim[], PetscSF *sf, PetscBool transferIdent, PetscInt *childIds[])
{
  MPI_Comm          comm;
  PetscMPIInt       rank, size;
  DM_Forest_pforest *pforestC, *pforestF;
  p4est_t           *p4estC, *p4estF;
  PetscInt          numClosureIndices, *closureIndices;
  PetscInt          numClosurePointsC, numClosurePointsF;
  PetscSFNode       *closurePointsC, *closurePointsF;
  p4est_quadrant_t  *coverQuads = NULL;
  p4est_quadrant_t  **treeQuads;
  PetscInt          *treeQuadCounts;
  MPI_Datatype      nodeType;
  MPI_Datatype      nodeClosureType;
  MPI_Op            sfNodeReduce;
  p4est_topidx_t    fltF, lltF, t;
  DM                plexC;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  pforestC = (DM_Forest_pforest *) ((DM_Forest *) coarse->data)->data;
  pforestF = (DM_Forest_pforest *) ((DM_Forest *) fine->data)->data;
  p4estC   = pforestC->forest;
  p4estF   = pforestF->forest;
  if (pforestC->topo != pforestF->topo) SETERRQ(PetscObjectComm((PetscObject)coarse),PETSC_ERR_ARG_INCOMP,"DM's must have the same base DM");
  comm = PetscObjectComm((PetscObject)coarse);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

  /* count the number of closure points that have dofs and create a list */
  numClosureIndices = 0;
  if (dofPerDim[P4EST_DIM]     > 0) numClosureIndices += 1;
  if (dofPerDim[P4EST_DIM - 1] > 0) numClosureIndices += P4EST_FACES;
#if defined(P4_TO_P8)
  if (dofPerDim[P4EST_DIM - 2] > 0) numClosureIndices += P8EST_EDGES;
#endif
  if (dofPerDim[0]             > 0) numClosureIndices += P4EST_CHILDREN;
  ierr = PetscMalloc1(numClosureIndices,&closureIndices);CHKERRQ(ierr);
  {
    PetscInt count = 0, offset = 0;
    if (dofPerDim[P4EST_DIM] > 0) {
      closureIndices[count++] = offset;
    }
    offset++;
    if (dofPerDim[P4EST_DIM - 1] > 0) {
      PetscInt i;

      for (i = 0; i < P4EST_FACES; i++) {
        closureIndices[count + i] = offset + i;
      }
      count += P4EST_FACES;
    }
    offset += P4EST_FACES;
#if defined(P4_TO_P8)
    if (dofPerDim[P4EST_DIM - 2] > 0) {
      PetscInt i;

      for (i = 0; i < P8EST_EDGES; i++) {
        closureIndices[count + i] = offset + i;
      }
      count += P8EST_EDGES;
    }
    offset += P8EST_EDGES;
#endif
    if (dofPerDim[0] > 0) {
      PetscInt i;

      for (i = 0; i < P4EST_CHILDREN; i++) {
        closureIndices[count + i] = offset + i;
      }
    }
  }
  /* create the datatype */
  ierr = MPI_Type_contiguous(2,MPIU_INT,&nodeType);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&nodeType);CHKERRQ(ierr);
  ierr = MPI_Op_create(DMPforestMaxSFNode,PETSC_FALSE,&sfNodeReduce);CHKERRQ(ierr);
  ierr = MPI_Type_contiguous(numClosureIndices*2,MPIU_INT,&nodeClosureType);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&nodeClosureType);CHKERRQ(ierr);
  /* everything has to go through cells: for each cell, create a list of the sfnodes in its closure */
  /* get lists of closure point SF nodes for every cell */
  ierr = DMPforestGetCellSFNodes(coarse,numClosureIndices,closureIndices,&numClosurePointsC,&closurePointsC,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DMPforestGetCellSFNodes(fine  ,numClosureIndices,closureIndices,&numClosurePointsF,&closurePointsF,PETSC_FALSE);CHKERRQ(ierr);
  /* create pointers for tree lists */
  fltF = p4estF->first_local_tree;
  lltF = p4estF->last_local_tree;
  ierr = PetscCalloc2(lltF + 1  - fltF, &treeQuads, lltF + 1 - fltF, &treeQuadCounts);CHKERRQ(ierr);
  ierr = DMPforestGetPlex(coarse,&plexC);CHKERRQ(ierr);
  /* if the partitions don't match, ship the coarse to cover the fine */
  if (size > 1) {
    PetscInt p;
    for (p = 0; p < size; p++) {
      int equal;

      PetscStackCallP4estReturn(equal,p4est_quadrant_is_equal_piggy,(&p4estC->global_first_position[p],&p4estF->global_first_position[p]));
      if (!equal) {
        break;
      }
    }
    if (p < size) { /* non-matching distribution: send the coarse to cover the fine */
      PetscInt cStartC, cEndC, cEndCInterior;
      PetscSF coveringSF;
      PetscInt nleaves;
      PetscInt count;
      PetscSFNode *newClosurePointsC;
      p4est_quadrant_t *coverQuadsSend;
      p4est_topidx_t fltC = p4estC->first_local_tree;
      p4est_topidx_t lltC = p4estC->last_local_tree;
      p4est_topidx_t t;
      PetscMPIInt blockSizes[5] = {P4EST_DIM,2,1,1,1};
      MPI_Aint    blockOffsets[5] = {offsetof(p4est_quadrant_t,x),
                                    offsetof(p4est_quadrant_t,level),
                                    offsetof(p4est_quadrant_t,pad16),
                                    offsetof(p4est_quadrant_t,p),
                                    sizeof (p4est_quadrant_t)};
      MPI_Datatype blockTypes[5] = {MPI_INT32_T,MPI_INT8_T,MPI_INT16_T,MPI_INT32_T,MPI_UB};
      MPI_Datatype quadType;
      ierr = DMPlexGetHeightStratum(plexC,0,&cStartC,&cEndC);CHKERRQ(ierr);
      ierr = DMPlexGetHybridBounds(plexC,&cEndCInterior,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = DMPforestGetCellCoveringSF(comm,p4estC,p4estF,cStartC,cEndC,pforestC->cLocalStart,&coveringSF);CHKERRQ(ierr);
      ierr = PetscSFGetGraph(coveringSF,NULL,&nleaves,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscMalloc1(numClosureIndices*nleaves,&newClosurePointsC);CHKERRQ(ierr);
      ierr = PetscMalloc1(nleaves,&coverQuads);CHKERRQ(ierr);
      ierr = PetscMalloc1(cEndC-cStartC,&coverQuadsSend);CHKERRQ(ierr);
      count = 0;
      for (t = fltC; t <= lltC; t++) { /* unfortunately, we need to pack a send array, since quads are not stored packed in p4est */
        p4est_tree_t *tree = &(((p4est_tree_t *) p4estC->trees->array)[t]);
        PetscInt q;

        ierr = PetscMemcpy(&coverQuadsSend[count],tree->quadrants.array,tree->quadrants.elem_count * sizeof(p4est_quadrant_t));CHKERRQ(ierr);
        for (q = 0; q < tree->quadrants.elem_count; q++) {
          coverQuadsSend[count+q].p.which_tree = t;
        }
        count += tree->quadrants.elem_count;
      }
      ierr = MPI_Type_create_struct(5,blockSizes,blockOffsets,blockTypes,&quadType);CHKERRQ(ierr);
      ierr = MPI_Type_commit(&quadType);CHKERRQ(ierr);
      ierr = PetscSFBcastBegin(coveringSF,nodeClosureType,closurePointsC,newClosurePointsC);CHKERRQ(ierr);
      ierr = PetscSFBcastBegin(coveringSF,quadType,coverQuadsSend,coverQuads);CHKERRQ(ierr);
      ierr = PetscSFBcastEnd(coveringSF,nodeClosureType,closurePointsC,newClosurePointsC);CHKERRQ(ierr);
      ierr = PetscSFBcastEnd(coveringSF,quadType,coverQuadsSend,coverQuads);CHKERRQ(ierr);
      ierr = MPI_Type_free(&quadType);CHKERRQ(ierr);
      ierr = PetscFree(coverQuadsSend);CHKERRQ(ierr);
      ierr = PetscFree(closurePointsC);CHKERRQ(ierr);
      ierr = PetscSFDestroy(&coveringSF);CHKERRQ(ierr);
      closurePointsC = newClosurePointsC;

      /* assign tree quads based on locations in coverQuads */
      {
        PetscInt q;
        for (q = 0; q < nleaves; q++) {
          p4est_locidx_t t = coverQuads[q].p.which_tree;
          if (!treeQuadCounts[t-fltF]++) {
            treeQuads[t-fltF] = &coverQuads[q];
          }
        }
      }
    }
  }
  if (!coverQuads) { /* matching partitions: assign tree quads based on locations in p4est native arrays */
    for (t = fltF; t <= lltF; t++) {
      p4est_tree_t *tree = &(((p4est_tree_t *) p4estC->trees->array)[t]);

      treeQuadCounts[t - fltF] = tree->quadrants.elem_count;
      treeQuads[t - fltF] = (p4est_quadrant_t *) tree->quadrants.array;
    }
  }

  {
    PetscInt pStartF, pEndF, p;
    PetscInt cLocalStartF;
    PetscSF  pointSF;
    PetscSFNode *roots;
    PetscInt *rootType;
    DM       plexF, refTree = NULL;
    DMLabel  canonical;
    PetscInt *childClosures[P4EST_CHILDREN] = {NULL};
    PetscInt *rootClosure = NULL;
    PetscInt *cids = NULL;
    PetscInt coarseOffset;
    PetscInt numCoarseQuads;

    ierr = DMPforestGetPlex(fine,&plexF);CHKERRQ(ierr);
    ierr = DMPlexGetChart(plexF,&pStartF,&pEndF);CHKERRQ(ierr);
    ierr = PetscMalloc1(pEndF-pStartF,&roots);CHKERRQ(ierr);
    ierr = PetscMalloc1(pEndF-pStartF,&rootType);CHKERRQ(ierr);
    ierr = DMGetPointSF(fine,&pointSF);CHKERRQ(ierr);
    for (p = pStartF; p < pEndF; p++) {
      roots[p-pStartF].rank  = -1;
      roots[p-pStartF].index = -1;
      rootType[p-pStartF] = -1;
    }
    if (childIds) {
      PetscInt child;

      ierr = PetscMalloc1(pEndF-pStartF,&cids);CHKERRQ(ierr);
      for (p = pStartF; p < pEndF; p++) {
        cids[p - pStartF] = -1;
      }
      ierr = DMPlexGetReferenceTree(plexF,&refTree);CHKERRQ(ierr);
      ierr = DMPlexGetTransitiveClosure(refTree,0,PETSC_TRUE,NULL,&rootClosure);CHKERRQ(ierr);
      for (child = 0; child < P4EST_CHILDREN; child++) { /* get the closures of the child cells in the reference tree */
        ierr = DMPlexGetTransitiveClosure(refTree,child+1,PETSC_TRUE,NULL,&childClosures[child]);CHKERRQ(ierr);
      }
      ierr = DMGetLabel(refTree,"canonical",&canonical);CHKERRQ(ierr);
    }
    cLocalStartF = pforestF->cLocalStart;
    for (t = fltF, coarseOffset = 0, numCoarseQuads = 0; t <= lltF; t++, coarseOffset += numCoarseQuads) {
      p4est_tree_t *tree = &(((p4est_tree_t *) p4estF->trees->array)[t]);
      PetscInt numFineQuads = tree->quadrants.elem_count;
      p4est_quadrant_t *coarseQuads = treeQuads[t - fltF];
      p4est_quadrant_t *fineQuads = (p4est_quadrant_t *) tree->quadrants.array;
      PetscInt i, coarseCount = 0;
      PetscInt offset = cLocalStartF + tree->quadrants_offset;
      sc_array_t coarseQuadsArray;

      numCoarseQuads = treeQuadCounts[t - fltF];
      PetscStackCallP4est(sc_array_init_data,(&coarseQuadsArray,coarseQuads,sizeof(p4est_quadrant_t),(size_t) numCoarseQuads));
      for (i = 0; i < numFineQuads; i++) {
        PetscInt c = i + offset;
        p4est_quadrant_t *quad = &fineQuads[i];
        p4est_quadrant_t *quadCoarse;
        ssize_t overlap = -1;

        while (overlap < 0 && coarseCount < numCoarseQuads) {
          quadCoarse = &coarseQuads[coarseCount];
          PetscStackCallP4estReturn(overlap,p4est_quadrant_disjoint,(quadCoarse,quad));
          if (overlap < 0) {
            coarseCount++;
          }
        }
        if (overlap != 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"did not find overlapping coarse quad");
        if (quadCoarse->level > quad->level || (quadCoarse->level == quad->level && !transferIdent)) { /* the "coarse" mesh is finer than the fine mesh at the point: continue */
          continue;
        }
        if (quadCoarse->level == quad->level) { /* same quad present in coarse and fine mesh */
          PetscInt j;
          for (j = 0; j < numClosureIndices; j++) {
            PetscInt p = closurePointsF[numClosureIndices * c + j].index;

            if (j > rootType[p-pStartF]) {
              roots[p-pStartF] = closurePointsC[numClosureIndices * (coarseCount + coarseOffset) + j];
              rootType[p-pStartF] = j;
            }
          }
        }
        else {
          if (childIds) {
            PetscInt cl;
            PetscInt *pointClosure = NULL;
            int cid;

            if (quadCoarse->level < quad->level - 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Recursive child ids not implemented");
            PetscStackCallP4estReturn(cid,p4est_quadrant_child_id,(quad));
            ierr = DMPlexGetTransitiveClosure(plexF,c,PETSC_TRUE,NULL,&pointClosure);CHKERRQ(ierr);
            for (cl = 0; cl < P4EST_INSUL; cl++) {
              PetscInt point = childClosures[cid][2 * cl];
              PetscInt ornt  = childClosures[cid][2 * cl + 1];
              PetscInt newcid = -1;
              if (!cl) {
                newcid = cid + 1;
              } else {
                PetscInt rcl, parent, parentOrnt;

                ierr = DMPlexGetTreeParent(refTree,point,&parent,NULL);CHKERRQ(ierr);
                if (parent == point) {
                  newcid = -1;
                }
                else if (!parent) { /* in the root */
                  newcid = point;
                }
                else {
                  for (rcl = 1; rcl < P4EST_INSUL; rcl++) {
                    if (rootClosure[2 * rcl] == parent) {
                      parentOrnt = rootClosure[2 * rcl + 1];
                      break;
                    }
                  }
                  if (rcl >= P4EST_INSUL) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Couldn't find parent in root closure");
                  ierr = DMPlexReferenceTreeGetChildSymmetry(refTree,parent,parentOrnt,ornt,point,pointClosure[2 * rcl + 1],NULL,&newcid);CHKERRQ(ierr);
                }
              }
              if (newcid >= 0) {
                if (canonical) {
                  ierr = DMLabelGetValue(canonical,newcid,&newcid);CHKERRQ(ierr);
                }
                cids[pointClosure[2 * cl] - pStartF] = newcid;
              }
            }
            ierr = DMPlexRestoreTransitiveClosure(plexF,c,PETSC_TRUE,NULL,&pointClosure);CHKERRQ(ierr);
          }
          p4est_qcoord_t coarseBound[2][P4EST_DIM] = {{quadCoarse->x,quadCoarse->y,
#if defined(P4_TO_P8)
                                                   quadCoarse->z
#endif
          },{0}};
          p4est_qcoord_t fineBound[2][P4EST_DIM] = {{quad->x,quad->y,
#if defined(P4_TO_P8)
                                                 quad->z
#endif
          },{0}};
          PetscInt j;
          for (j = 0; j < P4EST_DIM; j++) { /* get the coordinates of cell boundaries in each direction */
            coarseBound[1][j] = coarseBound[0][j] + P4EST_QUADRANT_LEN(quadCoarse->level);
            fineBound[1][j]   = fineBound[0][j]   + P4EST_QUADRANT_LEN(quad->level);
          }
          for (j = 0; j < numClosureIndices; j++) {
            PetscInt k = closureIndices[j];
            PetscInt l, m;

            if (k == 0) { /* volume: ancestor is volume */
              l = 0;
            }
            else if (k < 1 + P4EST_FACES) { /* facet */
              PetscInt face = PetscFaceToP4estFace[k - 1];
              PetscInt direction = face / 2;
              PetscInt coarseFace = -1;

              if (coarseBound[face % 2][direction] == fineBound[face % 2][direction]) {
                coarseFace = face;
                l = 1 + P4estFaceToPetscFace[coarseFace];
              }
              else {
                l = 0;
              }
            }
#if defined(P4_TO_P8)
            else if (k < 1 + P4EST_FACES + P8EST_EDGES) {
              PetscInt edge = PetscEdgeToP4estEdge[k - (1 + P4EST_FACES)];
              PetscInt direction = edge / 4;
              PetscInt mod = edge % 4;
              PetscInt coarseEdge = -1, coarseFace = -1;
              PetscInt minDir = PetscMin((direction + 1) % 3,(direction + 2) % 3);
              PetscInt maxDir = PetscMax((direction + 1) % 3,(direction + 2) % 3);
              PetscBool dirTest[2];

              dirTest[0] = (PetscBool) (coarseBound[mod % 2][minDir] == fineBound[mod % 2][minDir]);
              dirTest[1] = (PetscBool) (coarseBound[mod / 2][maxDir] == fineBound[mod / 2][maxDir]);

              if (dirTest[0] && dirTest[1]) { /* fine edge falls on coarse edge */
                coarseEdge = edge;
                l = 1 + P4EST_FACES + P4estEdgeToPetscEdge[coarseEdge];
              }
              else if (dirTest[0]) { /* fine edge falls on a coarse face in the minDir direction */
                coarseFace = 2 * minDir + (mod % 2);
                l = 1 + P4estFaceToPetscFace[coarseFace];
              }
              else if (dirTest[1]) { /* fine edge falls on a coarse face in the maxDir direction */
                coarseFace = 2 * maxDir + (mod / 2);
                l = 1 + P4estFaceToPetscFace[coarseFace];
              }
              else {
                l = 0;
              }
            }
#endif
            else {
              PetscInt vertex = PetscVertToP4estVert[P4EST_CHILDREN - (P4EST_INSUL - k)];
              PetscBool dirTest[P4EST_DIM];
              PetscInt m;
              PetscInt numMatch = 0;
              PetscInt coarseVertex = -1, coarseFace = -1;
#if defined(P4_TO_P8)
              PetscInt coarseEdge = -1;
#endif

              for (m = 0; m < P4EST_DIM; m++) {
                dirTest[m] = (PetscBool) (coarseBound[(vertex >> m) & 1][m] == fineBound[(vertex >> m) & 1][m]);
                if (dirTest[m]) {
                  numMatch++;
                }
              }
              if (numMatch == P4EST_DIM) { /* vertex on vertex */
                coarseVertex = vertex;
                l = P4EST_INSUL - (P4EST_CHILDREN - P4estVertToPetscVert[coarseVertex]);
              }
              else if (numMatch == 1) { /* vertex on face */
                for (m = 0; m < P4EST_DIM; m++) {
                  if (dirTest[m]) {
                    coarseFace = 2 * m + ((vertex >> m) & 1);
                    break;
                  }
                }
                l = 1 + P4estFaceToPetscFace[coarseFace];
              }
#if defined(P4_TO_P8)
              else if (numMatch == 2) { /* vertex on edge */
                for (m = 0; m < P4EST_DIM; m++) {
                  if (!dirTest[m]) {
                    PetscInt otherDir1 = (m + 1) % 3;
                    PetscInt otherDir2 = (m + 2) % 3;
                    PetscInt minDir = PetscMin(otherDir1,otherDir2);
                    PetscInt maxDir = PetscMax(otherDir1,otherDir2);

                    coarseEdge = m * 4 + maxDir * 2 * ((vertex >> maxDir) & 1) + minDir * ((vertex >> minDir) & 1);
                    break;
                  }
                }
                l = 1 + P4EST_FACES + P4estEdgeToPetscEdge[coarseEdge];
              }
#endif
              else { /* volume */
                l = 0;
              }
            }
            for (m = 0; m < numClosureIndices; m++) {
              if (closureIndices[m] == l) {
                PetscInt p = closurePointsF[numClosureIndices * c + j].index;

                if (m > rootType[p-pStartF]) {
                  roots[p-pStartF] = closurePointsC[numClosureIndices * (coarseCount + coarseOffset) + m];
                  rootType[p-pStartF] = m;
                }
                break;
              }
            }
          }
        }
      }
    }
    ierr = PetscFree(rootType);CHKERRQ(ierr);

    /* now every cell has labeled the points in its closure, so we first make sure everyone agrees by reducing to roots, and the broadcast the agreements */
    if (size > 1) {
      ierr = PetscSFReduceBegin(pointSF,nodeType,roots,roots,sfNodeReduce);CHKERRQ(ierr);
      ierr = PetscSFReduceEnd(pointSF,nodeType,roots,roots,sfNodeReduce);CHKERRQ(ierr);
      ierr = PetscSFBcastBegin(pointSF,nodeType,roots,roots);CHKERRQ(ierr);
      ierr = PetscSFBcastEnd(pointSF,nodeType,roots,roots);CHKERRQ(ierr);
    }

    {
      PetscInt pStartC, pEndC;
      PetscInt numRoots;
      PetscInt numLeaves;
      PetscInt *leaves;
      PetscInt d;
      PetscSFNode *iremote;
      PetscSF  pointTransferSF;
      PetscSection leafSection, rootSection;
      PetscInt endInterior[4];
      /* count leaves */

      ierr = DMPlexGetChart(plexC,&pStartC,&pEndC);CHKERRQ(ierr);
      numRoots = pEndC - pStartC;

      numLeaves = 0;
      for (p = pStartF; p < pEndF; p++) {
        if (roots[p-pStartF].index >= 0) {
          numLeaves++;
        }
      }
      ierr = PetscMalloc1(numLeaves,&leaves);CHKERRQ(ierr);
      ierr = PetscMalloc1(numLeaves,&iremote);CHKERRQ(ierr);
      numLeaves = 0;
      for (p = pStartF; p < pEndF; p++) {
        if (roots[p-pStartF].index >= 0) {
          leaves[numLeaves] = p-pStartF;
          iremote[numLeaves] = roots[p-pStartF];
          numLeaves++;
        }
      }
      ierr = PetscFree(roots);CHKERRQ(ierr);
      ierr = PetscSFCreate(comm,&pointTransferSF);CHKERRQ(ierr);
      ierr = PetscSFSetGraph(pointTransferSF,numRoots,numLeaves,leaves,PETSC_OWN_POINTER,iremote,PETSC_OWN_POINTER);CHKERRQ(ierr);
      ierr = PetscSectionCreate(PETSC_COMM_SELF,&rootSection);CHKERRQ(ierr);
      ierr = PetscSectionCreate(PETSC_COMM_SELF,&leafSection);CHKERRQ(ierr);
      ierr = PetscSectionSetChart(rootSection,pStartC,pEndC);CHKERRQ(ierr);
      ierr = PetscSectionSetChart(leafSection,pStartF,pEndF);CHKERRQ(ierr);

      ierr = DMPlexGetHybridBounds(plexC,&endInterior[P4EST_DIM],&endInterior[P4EST_DIM - 1],&endInterior[1],&endInterior[0]);CHKERRQ(ierr);
      for (d = 0; d <= P4EST_DIM; d++) {
        PetscInt startC, endC, e;

        ierr = DMPlexGetDepthStratum(plexC,d,&startC,&endC);CHKERRQ(ierr);
        endC = endInterior[d] < 0 ? endC : endInterior[d];
        for (e = startC; e < endC; e++) {
          ierr = PetscSectionSetDof(rootSection,e,dofPerDim[d]);CHKERRQ(ierr);
        }
      }

      ierr = DMPlexGetHybridBounds(plexF,&endInterior[P4EST_DIM],&endInterior[P4EST_DIM - 1],&endInterior[1],&endInterior[0]);CHKERRQ(ierr);
      for (d = 0; d <= P4EST_DIM; d++) {
        PetscInt startF, endF, e;

        ierr = DMPlexGetDepthStratum(plexF,d,&startF,&endF);CHKERRQ(ierr);
        endF = endInterior[d] < 0 ? endF : endInterior[d];
        for (e = startF; e < endF; e++) {
          ierr = PetscSectionSetDof(leafSection,e,dofPerDim[d]);CHKERRQ(ierr);
        }
      }

      ierr = PetscSectionSetUp(rootSection);CHKERRQ(ierr);
      ierr = PetscSectionSetUp(leafSection);CHKERRQ(ierr);
      {
        PetscInt nroots, nleaves;
        PetscInt *mine, i, p;
        PetscInt *offsets, *offsetsRoot;
        PetscSFNode *remote;

        ierr = PetscMalloc1(pEndF-pStartF,&offsets);CHKERRQ(ierr);
        ierr = PetscMalloc1(pEndC-pStartC,&offsetsRoot);CHKERRQ(ierr);
        for (p = pStartC; p < pEndC; p++) {
          ierr = PetscSectionGetOffset(rootSection,p,&offsetsRoot[p-pStartC]);CHKERRQ(ierr);
        }
        ierr = PetscSFBcastBegin(pointTransferSF,MPIU_INT,offsetsRoot,offsets);CHKERRQ(ierr);
        ierr = PetscSFBcastEnd(pointTransferSF,MPIU_INT,offsetsRoot,offsets);CHKERRQ(ierr);
        ierr = PetscSectionGetStorageSize(rootSection,&nroots);CHKERRQ(ierr);
        nleaves = 0;
        for (i = 0; i < numLeaves; i++) {
          PetscInt leaf = leaves[i];
          PetscInt dof;

          ierr = PetscSectionGetDof(leafSection,leaf,&dof);CHKERRQ(ierr);
          nleaves += dof;
        }
        ierr = PetscMalloc1(nleaves,&mine);CHKERRQ(ierr);
        ierr = PetscMalloc1(nleaves,&remote);CHKERRQ(ierr);
        nleaves = 0;
        for (i = 0; i < numLeaves; i++) {
          PetscInt leaf = leaves[i];
          PetscInt dof;
          PetscInt off, j;

          ierr = PetscSectionGetDof(leafSection,leaf,&dof);CHKERRQ(ierr);
          ierr = PetscSectionGetOffset(leafSection,leaf,&off);CHKERRQ(ierr);
          for (j = 0; j < dof; j++) {
            remote[nleaves].rank = iremote[i].rank;
            remote[nleaves].index = offsets[leaf] + j;
            mine[nleaves++] = off + j;
          }
        }
        ierr = PetscFree(offsetsRoot);CHKERRQ(ierr);
        ierr = PetscFree(offsets);CHKERRQ(ierr);
        ierr = PetscSFCreate(comm,sf);CHKERRQ(ierr);
        ierr = PetscSFSetGraph(*sf,nroots,nleaves,mine,PETSC_OWN_POINTER,remote,PETSC_OWN_POINTER);CHKERRQ(ierr);
      }
      ierr = PetscSectionDestroy(&leafSection);CHKERRQ(ierr);
      ierr = PetscSectionDestroy(&rootSection);CHKERRQ(ierr);
      ierr = PetscSFDestroy(&pointTransferSF);CHKERRQ(ierr);
    }
    if (childIds) {
      PetscSF  sf;
      PetscInt child;

      ierr = DMPlexGetReferenceTree(plexF,&refTree);CHKERRQ(ierr);
      ierr = DMGetPointSF(plexF,&sf);CHKERRQ(ierr);
      ierr = PetscSFReduceBegin(sf,MPIU_INT,cids,cids,MPIU_MAX);CHKERRQ(ierr);
      ierr = PetscSFReduceEnd(sf,MPIU_INT,cids,cids,MPIU_MAX);CHKERRQ(ierr);
      *childIds = cids;
      for (child = 0; child < P4EST_CHILDREN; child++) {
        ierr = DMPlexRestoreTransitiveClosure(refTree,child+1,PETSC_TRUE,NULL,&childClosures[child]);CHKERRQ(ierr);
      }
      ierr = DMPlexRestoreTransitiveClosure(refTree,0,PETSC_TRUE,NULL,&rootClosure);CHKERRQ(ierr);
    }
  }
  ierr = PetscFree2(treeQuads,treeQuadCounts);CHKERRQ(ierr);
  ierr = PetscFree(coverQuads);CHKERRQ(ierr);
  ierr = PetscFree(closurePointsC);CHKERRQ(ierr);
  ierr = PetscFree(closurePointsF);CHKERRQ(ierr);
  ierr = PetscFree(closureIndices);CHKERRQ(ierr);
  ierr = MPI_Type_free(&nodeClosureType);CHKERRQ(ierr);
  ierr = MPI_Op_free(&sfNodeReduce);CHKERRQ(ierr);
  ierr = MPI_Type_free(&nodeType);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestGetTransferSF"
static PetscErrorCode DMPforestGetTransferSF(DM dmA, DM dmB, const PetscInt dofPerDim[], PetscSF *sfAtoB, PetscSF *sfBtoA)
{
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  if (sfAtoB) {
    ierr = DMPforestGetTransferSF_Internal(dmA,dmB,dofPerDim,sfAtoB,PETSC_TRUE,NULL);CHKERRQ(ierr);
  }
  if (sfBtoA) {
    ierr = DMPforestGetTransferSF_Internal(dmB,dmA,dofPerDim,sfBtoA,(PetscBool) (sfAtoB == NULL),NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestLabelsInitialize"
static PetscErrorCode DMPforestLabelsInitialize(DM dm, DM plex)
{
  DM_Forest         *forest  = (DM_Forest *) dm->data;
  DM_Forest_pforest *pforest = (DM_Forest_pforest *) forest->data;
  PetscInt          cLocalStart, cLocalEnd, cStart, cEnd, fStart, fEnd, eStart, eEnd, vStart, vEnd;
  PetscInt          cStartBase, cEndBase, fStartBase, fEndBase, vStartBase, vEndBase, eStartBase, eEndBase;
  PetscInt          cEndBaseInterior, fEndBaseInterior, vEndBaseInterior, eEndBaseInterior;
  PetscInt          cEndInterior, fEndInterior, vEndInterior, eEndInterior;
  PetscInt          pStart, pEnd, pStartBase, pEndBase, p;
  DM                base;
  PetscInt          *star = NULL, starSize;
  DMLabelLink       next = dm->labels->next;
  PetscInt          guess = 0;
  p4est_topidx_t    num_trees = pforest->topo->conn->num_trees;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  pforest->labelsFinalized = PETSC_TRUE;
  cLocalStart = pforest->cLocalStart;
  cLocalEnd   = pforest->cLocalEnd;
  ierr = DMForestGetBaseDM(dm,&base);CHKERRQ(ierr);
  if (!base) {
    if (pforest->ghostName) { /* insert a label to make the boundaries, with stratum values that are denote which face of the element touches the boundary */
      p4est_connectivity_t *conn  = pforest->topo->conn;
      p4est_t              *p4est = pforest->forest;
      p4est_tree_t         *trees = (p4est_tree_t *) p4est->trees->array;
      p4est_topidx_t t, flt = p4est->first_local_tree;
      p4est_topidx_t llt = pforest->forest->last_local_tree;
      DMLabel ghostLabel;
      PetscInt c;

      ierr = DMCreateLabel(plex,pforest->ghostName);CHKERRQ(ierr);
      ierr = DMGetLabel(plex,pforest->ghostName,&ghostLabel);CHKERRQ(ierr);
      for (c = cLocalStart, t = flt; t <= llt; t++) {
        p4est_tree_t     *tree     = &trees[t];
        p4est_quadrant_t *quads    = (p4est_quadrant_t *) tree->quadrants.array;
        PetscInt          numQuads = (PetscInt) tree->quadrants.elem_count;
        PetscInt          q;

        for (q = 0; q < numQuads; q++, c++) {
          p4est_quadrant_t *quad = &quads[q];
          PetscInt f;

          for (f = 0; f < P4EST_FACES; f++) {
            p4est_quadrant_t neigh;
            int              isOutside;

            PetscStackCallP4est(p4est_quadrant_face_neighbor,(quad,f,&neigh));
            PetscStackCallP4estReturn(isOutside,p4est_quadrant_is_outside_face,(&neigh));
            if (isOutside) {
              p4est_topidx_t nt;
              PetscInt nf;

              nt = conn->tree_to_tree[t * P4EST_FACES + f];
              nf = (PetscInt) conn->tree_to_face[t * P4EST_FACES + f];
              nf = nf % P4EST_FACES;
              if (nt == t && nf == f) {
                PetscInt plexF = P4estFaceToPetscFace[f];
                const PetscInt *cone;

                ierr = DMPlexGetCone(plex,c,&cone);CHKERRQ(ierr);
                ierr = DMLabelSetValue(ghostLabel,cone[plexF],plexF+1);CHKERRQ(ierr);
              }
            }
          }
        }
      }
    }
    PetscFunctionReturn(0);
  }
  ierr = DMPlexGetHybridBounds(base,&cEndBaseInterior,&fEndBaseInterior,&eEndBaseInterior,&vEndBaseInterior);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(base,0,&cStartBase,&cEndBase);CHKERRQ(ierr);
  cEndBase = cEndBaseInterior < 0 ? cEndBase : cEndBaseInterior;
  ierr = DMPlexGetHeightStratum(base,1,&fStartBase,&fEndBase);CHKERRQ(ierr);
  fEndBase = fEndBaseInterior < 0 ? fEndBase : fEndBaseInterior;
  ierr = DMPlexGetDepthStratum(base,1,&eStartBase,&eEndBase);CHKERRQ(ierr);
  eEndBase = eEndBaseInterior < 0 ? eEndBase : eEndBaseInterior;
  ierr = DMPlexGetDepthStratum(base,0,&vStartBase,&vEndBase);CHKERRQ(ierr);
  vEndBase = vEndBaseInterior < 0 ? vEndBase : vEndBaseInterior;

  ierr = DMPlexGetHybridBounds(plex,&cEndInterior,&fEndInterior,&eEndInterior,&vEndInterior);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(plex,0,&cStart,&cEnd);CHKERRQ(ierr);
  cEnd = cEndInterior < 0 ? cEnd : cEndInterior;
  ierr = DMPlexGetHeightStratum(plex,1,&fStart,&fEnd);CHKERRQ(ierr);
  fEnd = fEndInterior < 0 ? fEnd : fEndInterior;
  ierr = DMPlexGetDepthStratum(plex,1,&eStart,&eEnd);CHKERRQ(ierr);
  eEnd = eEndInterior < 0 ? eEnd : eEndInterior;
  ierr = DMPlexGetDepthStratum(plex,0,&vStart,&vEnd);CHKERRQ(ierr);
  vEnd = vEndInterior < 0 ? vEnd : vEndInterior;

  ierr = DMPlexGetChart(plex,&pStart,&pEnd);CHKERRQ(ierr);
  ierr = DMPlexGetChart(base,&pStartBase,&pEndBase);CHKERRQ(ierr);
  /* go through the mesh: use star to find a quadrant that borders a point.  Use the closure to determine the
   * orientation of the quadrant relative to that point.  Use that to relate the point to the numbering in the base
   * mesh, and extract a label value (since the base mesh is redundantly distributed, can be found locally). */
  while (next) {
    DMLabel baseLabel;
    DMLabel label = next->label;
    PetscBool isDepth, isGhost, isVTK;

    ierr = PetscStrcmp(label->name,"depth",&isDepth);CHKERRQ(ierr);
    if (isDepth) {
      next = next->next;
      continue;
    }
    ierr = PetscStrcmp(label->name,"ghost",&isGhost);CHKERRQ(ierr);
    if (isGhost) {
      next = next->next;
      continue;
    }
    ierr = PetscStrcmp(label->name,"vtk",&isVTK);CHKERRQ(ierr);
    if (isVTK) {
      next = next->next;
      continue;
    }
    ierr = DMGetLabel(base,label->name,&baseLabel);CHKERRQ(ierr);
    if (!baseLabel) {
      next = next->next;
      continue;
    }
    ierr = DMLabelCreateIndex(baseLabel,pStartBase,pEndBase);CHKERRQ(ierr);
    for (p = pStart; p < pEnd; p++) {
      PetscInt s, c = -1, l;
      PetscInt *closure = NULL, closureSize;
      p4est_quadrant_t * ghosts = (p4est_quadrant_t *) pforest->ghost->ghosts.array;
      p4est_tree_t *trees = (p4est_tree_t *) pforest->forest->trees->array;
      p4est_quadrant_t * q;
      PetscInt t, val;

      ierr = DMPlexGetTransitiveClosure(plex,p,PETSC_FALSE,&starSize,&star);CHKERRQ(ierr);
      for (s = 0; s < starSize; s++) {
        PetscInt point = star[2*s];

        if (cStart <= point && point < cEnd) {
          ierr = DMPlexGetTransitiveClosure(plex,point,PETSC_TRUE,&closureSize,&closure);CHKERRQ(ierr);
          for (l = 0; l < closureSize; l++) {
            PetscInt qParent = closure[2 * l], q;
            do {
              q = qParent;
              if (q == p) {
                c = point;
                break;
              }
              ierr = DMPlexGetTreeParent(plex,q,&qParent,NULL);CHKERRQ(ierr);
            } while (qParent != q);
            if (c != -1) {
              break;
            }
          }
          ierr = DMPlexRestoreTransitiveClosure(plex,point,PETSC_TRUE,NULL,&closure);CHKERRQ(ierr);
          if (l < closureSize) {
            break;
          }
        }
      }
      if (s == starSize) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Failed to find cell with point %d in its closure",p);
      ierr = DMPlexRestoreTransitiveClosure(plex,p,PETSC_FALSE,NULL,&star);CHKERRQ(ierr);

      if (c < cLocalStart) {
        /* get from the beginning of the ghost layer */
        q = &(ghosts[c]);
        t = (PetscInt) q->p.which_tree;
      }
      else if (c < cLocalEnd) {
        PetscInt lo = 0, hi = num_trees;
        /* get from local quadrants: have to find the right tree */

        c -= cLocalStart;

        do {
          p4est_tree_t *tree;

          if (guess < lo || guess >= num_trees || lo >= hi) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"failed binary search");
          tree = &trees[guess];
          if (c < tree->quadrants_offset) {
            hi = guess;
          }
          else if (c < tree->quadrants_offset + (PetscInt) tree->quadrants.elem_count) {
            q = &((p4est_quadrant_t *)tree->quadrants.array)[c - (PetscInt) tree->quadrants_offset];
            t = guess;
            break;
          }
          else {
            lo = guess + 1;
          }
          guess = lo + (hi - lo) / 2;
        } while (1);
      }
      else {
        /* get from the end of the ghost layer */
        c -= (cLocalEnd - cLocalEnd);

        q = &(ghosts[c]);
        t = (PetscInt) q->p.which_tree;
      }

      if (l == 0) { /* cell */
        ierr = DMLabelGetValue(baseLabel,t+cStartBase,&val);CHKERRQ(ierr);
        ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
      }
      else if (l >= 1 && l < 1 + P4EST_FACES) { /* facet */
        p4est_quadrant_t nq;
        int              isInside;

        l = PetscFaceToP4estFace[l - 1];
        PetscStackCallP4est(p4est_quadrant_face_neighbor,(q,l,&nq));
        PetscStackCallP4estReturn(isInside,p4est_quadrant_is_inside_root,(&nq));
        if (isInside) {
          /* this facet is in the interior of a tree, so it inherits the label of the tree */
          ierr = DMLabelGetValue(baseLabel,t+cStartBase,&val);CHKERRQ(ierr);
          ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
        }
        else {
          PetscInt f = pforest->topo->tree_face_to_uniq[P4EST_FACES * t + l];

          ierr = DMLabelGetValue(baseLabel,f+fStartBase,&val);CHKERRQ(ierr);
          ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
        }
      }
#if defined(P4_TO_P8)
      else if (l >= 1 + P4EST_FACES && l < 1 + P4EST_FACES + P8EST_EDGES) { /* edge */
        p4est_quadrant_t nq;
        int              isInside;

        l = PetscEdgeToP4estEdge[l - (1 + P4EST_FACES)];
        PetscStackCallP4est(p8est_quadrant_edge_neighbor,(q,l,&nq));
        PetscStackCallP4estReturn(isInside,p4est_quadrant_is_inside_root,(&nq));
        if (isInside) {
          /* this edge is in the interior of a tree, so it inherits the label of the tree */
          ierr = DMLabelGetValue(baseLabel,t+cStartBase,&val);CHKERRQ(ierr);
          ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
        }
        else {
          int isOutsideFace;

          PetscStackCallP4estReturn(isOutsideFace,p4est_quadrant_is_outside_face,(&nq));
          if (isOutsideFace) {
            PetscInt f;

            if (nq.x < 0) {
              f = 0;
            }
            else if (nq.x >= P4EST_ROOT_LEN) {
              f = 1;
            }
            else if (nq.y < 0) {
              f = 2;
            }
            else if (nq.y >= P4EST_ROOT_LEN) {
              f = 3;
            }
            else if (nq.z < 0) {
              f = 4;
            }
            else {
              f = 5;
            }
            f = pforest->topo->tree_face_to_uniq[P4EST_FACES * t + f];
            ierr = DMLabelGetValue(baseLabel,f+fStartBase,&val);CHKERRQ(ierr);
            ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
          }
          else { /* the quadrant edge corresponds to the tree edge */
            PetscInt e = pforest->topo->conn->tree_to_edge[P8EST_EDGES * t + l];

            ierr = DMLabelGetValue(baseLabel,e+eStartBase,&val);CHKERRQ(ierr);
            ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
          }
        }
      }
#endif
      else { /* vertex */
        p4est_quadrant_t nq;
        int              isInside;

#if defined(P4_TO_P8)
        l = PetscVertToP4estVert[l - (1 + P4EST_FACES + P8EST_EDGES)];
#else
        l = PetscVertToP4estVert[l - (1 + P4EST_FACES)];
#endif
        PetscStackCallP4est(p4est_quadrant_corner_neighbor,(q,l,&nq));
        PetscStackCallP4estReturn(isInside,p4est_quadrant_is_inside_root,(&nq));
        if (isInside) {
          ierr = DMLabelGetValue(baseLabel,t+cStartBase,&val);CHKERRQ(ierr);
          ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
        }
        else {
          int isOutside;

          PetscStackCallP4estReturn(isOutside,p4est_quadrant_is_outside_face,(&nq));
          if (isOutside) {
            PetscInt f = -1;

            if (nq.x < 0) {
              f = 0;
            }
            else if (nq.x >= P4EST_ROOT_LEN) {
              f = 1;
            }
            else if (nq.y < 0) {
              f = 2;
            }
            else if (nq.y >= P4EST_ROOT_LEN) {
              f = 3;
            }
#if defined(P4_TO_P8)
            else if (nq.z < 0) {
              f = 4;
            }
            else {
              f = 5;
            }
#endif
            f = pforest->topo->tree_face_to_uniq[P4EST_FACES * t + f];
            ierr = DMLabelGetValue(baseLabel,f+fStartBase,&val);CHKERRQ(ierr);
            ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
            continue;
          }
#if defined(P4_TO_P8)
          PetscStackCallP4estReturn(isOutside,p8est_quadrant_is_outside_edge,(&nq));
          if (isOutside) {
            /* outside edge */
            PetscInt e = -1;

            if (nq.x >= 0 && nq.x < P4EST_ROOT_LEN) {
              if (nq.z < 0) {
                if (nq.y < 0) {
                  e = 0;
                }
                else {
                  e = 1;
                }
              }
              else {
                if (nq.y < 0) {
                  e = 2;
                }
                else {
                  e = 3;
                }
              }
            }
            else if (nq.y >= 0 && nq.y < P4EST_ROOT_LEN) {
              if (nq.z < 0) {
                if (nq.x < 0) {
                  e = 4;
                }
                else {
                  e = 5;
                }
              }
              else {
                if (nq.x < 0) {
                  e = 6;
                }
                else {
                  e = 7;
                }
              }
            }
            else {
              if (nq.y < 0) {
                if (nq.x < 0) {
                  e = 8;
                }
                else {
                  e = 9;
                }
              }
              else {
                if (nq.x < 0) {
                  e = 10;
                }
                else {
                  e = 11;
                }
              }
            }

            e = pforest->topo->conn->tree_to_edge[P8EST_EDGES * t + e];
            ierr = DMLabelGetValue(baseLabel,e+eStartBase,&val);CHKERRQ(ierr);
            ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
            continue;
          }
#endif
          else {
            /* outside vertex: same corner as quadrant corner */
            PetscInt v = pforest->topo->conn->tree_to_corner[P4EST_CHILDREN * t + l];

            ierr = DMLabelGetValue(baseLabel,v+vStartBase,&val);CHKERRQ(ierr);
            ierr = DMLabelSetValue(label,p,val);CHKERRQ(ierr);
          }
        }
      }
    }
    next = next->next;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestLabelsFinalize"
static PetscErrorCode DMPforestLabelsFinalize(DM dm, DM plex)
{
  DM_Forest_pforest *pforest = (DM_Forest_pforest *) ((DM_Forest *) dm->data)->data;
  DM                adapt;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  if (pforest->labelsFinalized) PetscFunctionReturn(0);
  pforest->labelsFinalized = PETSC_TRUE;
  ierr = DMForestGetAdaptivityForest(dm,&adapt);CHKERRQ(ierr);
  if (!adapt) {
    /* Initialize labels from the base dm */
    ierr = DMPforestLabelsInitialize(dm,plex);CHKERRQ(ierr);
  }
  else {
    PetscInt    dofPerDim[4]={1, 1, 1, 1};
    PetscSF     transferForward, transferBackward, pointSF;
    PetscInt    pStart, pEnd, pStartA, pEndA;
    PetscInt    *values, *adaptValues;
    DMLabelLink next = adapt->labels->next;
    DM          adaptPlex;

    ierr = DMPforestGetPlex(adapt,&adaptPlex);CHKERRQ(ierr);
    ierr = DMPforestGetTransferSF(adapt,dm,dofPerDim,&transferForward,&transferBackward);CHKERRQ(ierr);
    ierr = DMPlexGetChart(plex,&pStart,&pEnd);CHKERRQ(ierr);
    ierr = DMPlexGetChart(adaptPlex,&pStartA,&pEndA);CHKERRQ(ierr);
    ierr = PetscMalloc2(pEnd-pStart,&values,pEndA-pStartA,&adaptValues);CHKERRQ(ierr);
    ierr = DMGetPointSF(plex,&pointSF);CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
    {
      PetscInt p;
      for (p = pStartA; p < pEndA; p++) {
        adaptValues[p-pStartA] = -1;
      }
      for (p = pStart; p < pEnd; p++) {
        values[p-pStart] = -2;
      }
      if (transferForward) {
        ierr = PetscSFBcastBegin(transferForward,MPIU_INT,adaptValues,values);CHKERRQ(ierr);
        ierr = PetscSFBcastEnd(transferForward,MPIU_INT,adaptValues,values);CHKERRQ(ierr);
      }
      if (transferBackward) {
        ierr = PetscSFReduceBegin(transferBackward,MPIU_INT,adaptValues,values,MPIU_MAX);CHKERRQ(ierr);
        ierr = PetscSFReduceEnd(transferBackward,MPIU_INT,adaptValues,values,MPIU_MAX);CHKERRQ(ierr);
      }
      ierr = PetscSFBcastBegin(pointSF,MPIU_INT,values,values);CHKERRQ(ierr);
      ierr = PetscSFBcastEnd(pointSF,MPIU_INT,values,values);CHKERRQ(ierr);
      for (p = pStart; p < pEnd; p++) {
        if (values[p-pStart] == -2) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_PLIB,"uncovered point %D",p);
      }
    }
#endif
    while (next) {
      DMLabel adaptLabel = next->label;
      const char *name = adaptLabel->name;
      PetscBool isDepth, isGhost, isVTK;
      DMLabel label;
      PetscInt  p;

      ierr = PetscStrcmp(name,"depth",&isDepth);CHKERRQ(ierr);
      if (isDepth) {
        next = next->next;
        continue;
      }
      ierr = PetscStrcmp(name,"ghost",&isGhost);CHKERRQ(ierr);
      if (isGhost) {
        next = next->next;
        continue;
      }
      ierr = PetscStrcmp(name,"vtk",&isVTK);CHKERRQ(ierr);
      if (isVTK) {
        next = next->next;
        continue;
      }
      /* label was created earlier */
      ierr = DMGetLabel(dm,name,&label);CHKERRQ(ierr);

      for (p = pStartA; p < pEndA; p++) {
        ierr = DMLabelGetValue(adaptLabel,p,&adaptValues[p]);CHKERRQ(ierr);
      }
      for (p = pStart; p < pEnd; p++) {
        values[p] = PETSC_MIN_INT;
      }

      if (transferForward) {
        ierr = PetscSFBcastBegin(transferForward,MPIU_INT,adaptValues,values);CHKERRQ(ierr);
      }
      if (transferBackward) {
        ierr = PetscSFReduceBegin(transferBackward,MPIU_INT,adaptValues,values,MPIU_MAX);CHKERRQ(ierr);
      }
      if (transferForward) {
        ierr = PetscSFBcastEnd(transferForward,MPIU_INT,adaptValues,values);CHKERRQ(ierr);
      }
      if (transferBackward) {
        ierr = PetscSFReduceEnd(transferBackward,MPIU_INT,adaptValues,values,MPIU_MAX);CHKERRQ(ierr);
      }
      ierr = PetscSFBcastBegin(pointSF,MPIU_INT,values,values);CHKERRQ(ierr);
      ierr = PetscSFBcastEnd(pointSF,MPIU_INT,values,values);CHKERRQ(ierr);

      for (p = pStart; p < pEnd; p++) {
        ierr = DMLabelSetValue(label,p,values[p]);CHKERRQ(ierr);
      }
      next = next->next;
    }
    ierr = PetscFree2(values,adaptValues);CHKERRQ(ierr);
    ierr = PetscSFDestroy(&transferForward);CHKERRQ(ierr);
    ierr = PetscSFDestroy(&transferBackward);CHKERRQ(ierr);
    pforest->labelsFinalized = PETSC_TRUE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMConvert_pforest_plex)
static PetscErrorCode DMConvert_pforest_plex(DM dm, DMType newtype, DM *plex)
{
  DM_Forest      *forest;
  DM_Forest_pforest *pforest;
  DM             refTree, newPlex, base;
  PetscInt       adjDim, adjCodim, coordDim;
  MPI_Comm       comm;
  PetscBool      isPforest;
  PetscInt       dim;
  PetscInt       overlap;
  p4est_connect_type_t ctype;
  p4est_locidx_t first_local_quad = -1;
  sc_array_t     *points_per_dim, *cone_sizes, *cones, *cone_orientations, *coords, *children, *parents, *childids, *leaves, *remotes;
  PetscSection   parentSection;
  PetscSF        pointSF;
  size_t         zz, count;
  PetscInt       pStart, pEnd;
  DMLabel        ghostLabelBase = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  comm = PetscObjectComm((PetscObject)dm);
  ierr = PetscObjectTypeCompare((PetscObject)dm,DMPFOREST,&isPforest);CHKERRQ(ierr);
  if (!isPforest) SETERRQ2(comm,PETSC_ERR_ARG_WRONG,"Expected DM type %s, got %s\n",DMPFOREST,((PetscObject)dm)->type_name);
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  if (dim != P4EST_DIM) SETERRQ2(comm,PETSC_ERR_ARG_WRONG,"Expected DM dimension %d, got %d\n",P4EST_DIM,dim);
  forest = (DM_Forest *) dm->data;
  pforest = (DM_Forest_pforest *) forest->data;
  ierr = DMForestGetBaseDM(dm,&base);CHKERRQ(ierr);
  if (base) {
    ierr = DMGetLabel(base,"ghost",&ghostLabelBase);CHKERRQ(ierr);
  }
  if (!pforest->plex) {
    ierr = DMCreate(comm,&newPlex);CHKERRQ(ierr);
    ierr = DMSetType(newPlex,DMPLEX);CHKERRQ(ierr);
    ierr = PetscFree(newPlex->labels);CHKERRQ(ierr); /* share labels */
    dm->labels->refct++;
    newPlex->labels = dm->labels;
    ierr = DMForestGetAdjacencyDimension(dm,&adjDim);CHKERRQ(ierr);
    ierr = DMForestGetAdjacencyCodimension(dm,&adjCodim);CHKERRQ(ierr);
    ierr = DMGetCoordinateDim(dm,&coordDim);CHKERRQ(ierr);
    if (adjDim == 0) {
      ctype = P4EST_CONNECT_FULL;
    }
    else if (adjCodim == 1) {
      ctype = P4EST_CONNECT_FACE;
    }
#if defined(P4_TO_P8)
    else if (adjDim == 1) {
      ctype = P8EST_CONNECT_EDGE;
    }
#endif
    else {
      SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_WRONG,"Invalid adjacency dimension %d",adjDim);
    }
    ierr = DMForestGetPartitionOverlap(dm,&overlap);CHKERRQ(ierr);

    points_per_dim    = sc_array_new(sizeof (p4est_locidx_t));
    cone_sizes        = sc_array_new(sizeof (p4est_locidx_t));
    cones             = sc_array_new(sizeof (p4est_locidx_t));
    cone_orientations = sc_array_new(sizeof (p4est_locidx_t));
    coords            = sc_array_new(3 * sizeof (double));
    children          = sc_array_new(sizeof (p4est_locidx_t));
    parents           = sc_array_new(sizeof (p4est_locidx_t));
    childids          = sc_array_new(sizeof (p4est_locidx_t));
    leaves            = sc_array_new(sizeof (p4est_locidx_t));
    remotes           = sc_array_new(2 * sizeof (p4est_locidx_t));

    PetscStackCallP4est(p4est_get_plex_data_ext,(pforest->forest,&pforest->ghost,&pforest->lnodes,ctype,(int)overlap,&first_local_quad,points_per_dim,cone_sizes,cones,cone_orientations,coords,children,parents,childids,leaves,remotes,1));

    pforest->cLocalStart = (PetscInt) first_local_quad;
    pforest->cLocalEnd   = pforest->cLocalStart + (PetscInt) pforest->forest->local_num_quadrants;
    ierr = locidx_to_PetscInt(points_per_dim);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(cone_sizes);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(cones);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(cone_orientations);CHKERRQ(ierr);
    ierr = coords_double_to_PetscScalar(coords, coordDim);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(children);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(parents);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(childids);CHKERRQ(ierr);
    ierr = locidx_to_PetscInt(leaves);CHKERRQ(ierr);
    ierr = locidx_pair_to_PetscSFNode(remotes);CHKERRQ(ierr);

    ierr = DMSetDimension(newPlex,P4EST_DIM);CHKERRQ(ierr);
    ierr = DMSetCoordinateDim(newPlex,coordDim);CHKERRQ(ierr);
    ierr = DMPlexCreateFromDAG(newPlex,P4EST_DIM,(PetscInt *)points_per_dim->array,(PetscInt *)cone_sizes->array,(PetscInt *)cones->array,(PetscInt *)cone_orientations->array,(PetscScalar *)coords->array);CHKERRQ(ierr);
    ierr = PetscSFCreate(comm,&pointSF);CHKERRQ(ierr);
    ierr = DMCreateReferenceTree_pforest(comm,&refTree);CHKERRQ(ierr);
    ierr = DMPlexSetReferenceTree(newPlex,refTree);CHKERRQ(ierr);
    ierr = DMDestroy(&refTree);CHKERRQ(ierr);
    ierr = PetscSectionCreate(comm,&parentSection);CHKERRQ(ierr);
    ierr = DMPlexGetChart(newPlex,&pStart,&pEnd);CHKERRQ(ierr);
    ierr = PetscSectionSetChart(parentSection,pStart,pEnd);CHKERRQ(ierr);
    count = children->elem_count;
    for(zz = 0;zz < count;zz++) {
      PetscInt            child = *((PetscInt *) sc_array_index(children,zz));

      ierr = PetscSectionSetDof(parentSection,child,1);CHKERRQ(ierr);
    }
    ierr = PetscSectionSetUp(parentSection);CHKERRQ(ierr);
    ierr = DMPlexSetTree(newPlex,parentSection,(PetscInt *)parents->array,(PetscInt *)childids->array);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&parentSection);CHKERRQ(ierr);
    ierr = PetscSFSetGraph(pointSF,pEnd - pStart,(PetscInt)leaves->elem_count,(PetscInt *)leaves->array,PETSC_COPY_VALUES,(PetscSFNode *)remotes->array,PETSC_COPY_VALUES);CHKERRQ(ierr);
    ierr = DMSetPointSF(newPlex,pointSF);CHKERRQ(ierr);
    ierr = DMSetPointSF(dm,pointSF);CHKERRQ(ierr);
    ierr = PetscSFDestroy(&pointSF);CHKERRQ(ierr);
    if (dm->maxCell) {
      const PetscReal *maxCell, *L;
      const DMBoundaryType *bd;

      ierr = DMGetPeriodicity(dm,&maxCell,&L,&bd);CHKERRQ(ierr);
      ierr = DMSetPeriodicity(newPlex,maxCell,L,bd);CHKERRQ(ierr);
      ierr = DMLocalizeCoordinates(newPlex);CHKERRQ(ierr);
    }
    sc_array_destroy (points_per_dim);
    sc_array_destroy (cone_sizes);
    sc_array_destroy (cones);
    sc_array_destroy (cone_orientations);
    sc_array_destroy (coords);
    sc_array_destroy (children);
    sc_array_destroy (parents);
    sc_array_destroy (childids);
    sc_array_destroy (leaves);
    sc_array_destroy (remotes);

    pforest->plex = newPlex;

    /* copy labels */
    ierr = DMPforestLabelsFinalize(dm,newPlex);CHKERRQ(ierr);

    if (ghostLabelBase || pforest->ghostName) { /* we have to do this after copying labels because the labels drive the construction of ghost cells */
      PetscInt numAdded;
      DM       newPlexGhosted;
      void     *ctx;

      ierr = DMPlexConstructGhostCells(newPlex,pforest->ghostName,&numAdded,&newPlexGhosted);CHKERRQ(ierr);
      ierr = DMGetApplicationContext(newPlex,&ctx);CHKERRQ(ierr);
      ierr = DMSetApplicationContext(newPlexGhosted,ctx);CHKERRQ(ierr);
      /* we want the sf for the ghost dm to be the one for the p4est dm as well */
      ierr = DMGetPointSF(newPlexGhosted,&pointSF);CHKERRQ(ierr);
      ierr = DMSetPointSF(dm,pointSF);CHKERRQ(ierr);
      ierr = DMDestroy(&newPlex);CHKERRQ(ierr);
      newPlex = newPlexGhosted;

      /* share the labels back */
      ierr = DMDestroyLabelLinkList(dm);CHKERRQ(ierr);
      newPlex->labels->refct++;
      dm->labels = newPlex->labels;

      pforest->plex = newPlex;
    }

    if (forest->setfromoptionscalled) {
      ierr = PetscObjectOptionsBegin((PetscObject)newPlex);CHKERRQ(ierr);
      ierr = DMSetFromOptions_NonRefinement_Plex(PetscOptionsObject,newPlex);CHKERRQ(ierr);
      ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject) newPlex);CHKERRQ(ierr);
      ierr = PetscOptionsEnd();CHKERRQ(ierr);
    }
    {
      Vec coords;

      ierr = DMGetCoordinatesLocal(newPlex,&coords);CHKERRQ(ierr);
      ierr = DMSetCoordinatesLocal(dm,coords);CHKERRQ(ierr);
    }
  }
  newPlex = pforest->plex;
  if (plex) {
    DM      coordDM;

    ierr = DMClone(newPlex,plex);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(newPlex,&coordDM);CHKERRQ(ierr);
    ierr = DMSetCoordinateDM(*plex,coordDM);CHKERRQ(ierr);

    ierr = DMShareDiscretization(dm,*plex);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMSetFromOptions_pforest)
static PetscErrorCode DMSetFromOptions_pforest(PetscOptionItems *PetscOptionsObject,DM dm)
{
  DM_Forest_pforest *pforest = (DM_Forest_pforest *) ((DM_Forest *) dm->data)->data;
  char              stringBuffer[256];
  PetscBool         flg;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = DMSetFromOptions_Forest(PetscOptionsObject,dm);CHKERRQ(ierr);
  ierr = PetscOptionsHead(PetscOptionsObject,"DM" P4EST_STRING " options");CHKERRQ(ierr);
  ierr = PetscOptionsBool("-dm_p4est_partition_for_coarsening","partition forest to allow for coarsening","DMP4estSetPartitionForCoarsening",pforest->partition_for_coarsening,&(pforest->partition_for_coarsening),NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-dm_p4est_ghost_label_name","the name of the ghost label when converting from a DMPlex",NULL,NULL,stringBuffer,256,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  if (flg) {
    ierr = PetscFree(pforest->ghostName);CHKERRQ(ierr);
    ierr = PetscStrallocpy(stringBuffer,&pforest->ghostName);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#if !defined(P4_TO_P8)
#define DMPforestGetPartitionForCoarsening DMP4estGetPartitionForCoarsening
#define DMPforestSetPartitionForCoarsening DMP4estSetPartitionForCoarsening
#else
#define DMPforestGetPartitionForCoarsening DMP8estGetPartitionForCoarsening
#define DMPforestSetPartitionForCoarsening DMP8estSetPartitionForCoarsening
#endif

#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMPforestGetPartitionForCoarsening)
PETSC_EXTERN PetscErrorCode DMPforestGetPartitionForCoarsening(DM dm, PetscBool *flg)
{
  DM_Forest_pforest *pforest;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  pforest = (DM_Forest_pforest *) ((DM_Forest *) dm->data)->data;
  *flg = pforest->partition_for_coarsening;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMPforestSetPartitionForCoarsening)
PETSC_EXTERN PetscErrorCode DMPforestSetPartitionForCoarsening(DM dm, PetscBool flg)
{
  DM_Forest_pforest *pforest;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  pforest = (DM_Forest_pforest *) ((DM_Forest *) dm->data)->data;
  pforest->partition_for_coarsening = flg;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPforestGetPlex"
static PetscErrorCode DMPforestGetPlex(DM dm,DM *plex)
{
  DM_Forest_pforest *pforest;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr    = DMSetUp(dm);CHKERRQ(ierr);
  pforest = (DM_Forest_pforest *) ((DM_Forest *) dm->data)->data;
  if (!pforest->plex) {
    ierr = DMConvert_pforest_plex(dm,DMPLEX,NULL);CHKERRQ(ierr);
  }
  ierr  = DMShareDiscretization(dm,pforest->plex);CHKERRQ(ierr);
  *plex = pforest->plex;
  PetscFunctionReturn(0);
}

#define DMCoarsen_pforest _append_pforest(DMCoarsen)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCoarsen_pforest)
static PetscErrorCode DMCoarsen_pforest(DM dm, MPI_Comm comm, DM *dmc)
{
  DM             coarseDM;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  {
    PetscMPIInt mpiComparison;
    MPI_Comm    dmcomm = PetscObjectComm((PetscObject)dm);

    ierr = MPI_Comm_compare(comm, dmcomm, &mpiComparison);CHKERRQ(ierr);
    if (mpiComparison != MPI_IDENT && mpiComparison != MPI_CONGRUENT) {
      SETERRQ(dmcomm,PETSC_ERR_SUP,"No support for different communicators yet");
    }
  }
  ierr = DMGetCoarseDM(dm,&coarseDM);CHKERRQ(ierr);
  if (!coarseDM) {
    DM coarseDM;

    ierr = DMCreate(PetscObjectComm((PetscObject)dm),&coarseDM);CHKERRQ(ierr);
  }
  if (coarseDM) {
    void *ctx;
    PetscDS ds;

    ierr = DMGetApplicationContext(dm,&ctx);CHKERRQ(ierr);
    ierr = DMSetApplicationContext(coarseDM,ctx);CHKERRQ(ierr);
    ierr = DMGetDS(dm,&ds);CHKERRQ(ierr);
    ierr = DMSetDS(coarseDM,ds);CHKERRQ(ierr);
    ierr = DMCopyBoundary(dm,coarseDM);CHKERRQ(ierr);
  }
  *dmc = coarseDM;
  PetscFunctionReturn(0);
}

#define DMCreateInterpolation_pforest _append_pforest(DMCreateInterpolation)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateInterpolation_pforest)
static PetscErrorCode DMCreateInterpolation_pforest (DM dmCoarse, DM dmFine, Mat *interpolation, Vec *scaling)
{
  PetscSection   gsc, gsf;
  PetscInt       m, n;
  DM             cdm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetDefaultGlobalSection(dmFine, &gsf);CHKERRQ(ierr);
  ierr = PetscSectionGetConstrainedStorageSize(gsf, &m);CHKERRQ(ierr);
  ierr = DMGetDefaultGlobalSection(dmCoarse, &gsc);CHKERRQ(ierr);
  ierr = PetscSectionGetConstrainedStorageSize(gsc, &n);CHKERRQ(ierr);

  ierr = MatCreate(PetscObjectComm((PetscObject) dmFine), interpolation);CHKERRQ(ierr);
  ierr = MatSetSizes(*interpolation, m, n, PETSC_DETERMINE, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetType(*interpolation, dmCoarse->mattype);CHKERRQ(ierr);

  ierr = DMGetCoarseDM(dmFine, &cdm);CHKERRQ(ierr);
  if (cdm != dmCoarse) {
    SETERRQ(PetscObjectComm((PetscObject)dmFine),PETSC_ERR_SUP,"Only interpolation to coarse DM for now");
  }
  {
    DM           plexF, plexC;
    PetscSF      sf;
    PetscInt     *cids;
    PetscInt     dofPerDim[4] = {1,1,1,1};

    ierr = DMPforestGetPlex(dmCoarse,&plexC);CHKERRQ(ierr);
    ierr = DMPforestGetPlex(dmFine,&plexF);CHKERRQ(ierr);
    ierr = DMPforestGetTransferSF_Internal(dmCoarse, dmFine, dofPerDim, &sf, PETSC_TRUE, &cids);CHKERRQ(ierr);
    ierr = PetscSFSetUp(sf);CHKERRQ(ierr);
    ierr = DMPlexComputeInterpolatorTree(plexC, plexF, sf, cids, *interpolation);CHKERRQ(ierr);
    ierr = PetscSFDestroy(&sf);CHKERRQ(ierr);
    ierr = PetscFree(cids);
  }
  ierr = MatViewFromOptions(*interpolation, NULL, "-interp_mat_view");CHKERRQ(ierr);
  /* Use naive scaling */
  ierr = DMCreateInterpolationScale(dmCoarse, dmFine, *interpolation, scaling);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMCreateCoordinateDM_pforest _append_pforest(DMCreateCoordinateDM)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateCoordinateDM_pforest)
static PetscErrorCode DMCreateCoordinateDM_pforest(DM dm,DM *cdm)
{
  DM                 plex;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(plex,cdm);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)*cdm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define VecView_pforest _append_pforest(VecView)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(VecView_pforest)
static PetscErrorCode VecView_pforest(Vec vec,PetscViewer viewer)
{
  DM             dm, plex;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = VecSetDM(vec,plex);CHKERRQ(ierr);
  ierr = VecView_Plex(vec,viewer);CHKERRQ(ierr);
  ierr = VecSetDM(vec,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define VecView_pforest_Native _infix_pforest(VecView,_Native)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(VecView_pforest_Native)
static PetscErrorCode VecView_pforest_Native(Vec vec,PetscViewer viewer)
{
  DM             dm, plex;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = VecSetDM(vec,plex);CHKERRQ(ierr);
  ierr = VecView_Plex_Native(vec,viewer);CHKERRQ(ierr);
  ierr = VecSetDM(vec,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define VecLoad_pforest _append_pforest(VecLoad)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(VecLoad_pforest)
static PetscErrorCode VecLoad_pforest(Vec vec,PetscViewer viewer)
{
  DM             dm, plex;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = VecSetDM(vec,plex);CHKERRQ(ierr);
  ierr = VecLoad_Plex(vec,viewer);CHKERRQ(ierr);
  ierr = VecSetDM(vec,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define VecLoad_pforest_Native _infix_pforest(VecLoad,_Native)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(VecLoad_pforest_Native)
static PetscErrorCode VecLoad_pforest_Native(Vec vec,PetscViewer viewer)
{
  DM             dm, plex;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetDM(vec,&dm);CHKERRQ(ierr);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = VecSetDM(vec,plex);CHKERRQ(ierr);
  ierr = VecLoad_Plex_Native(vec,viewer);CHKERRQ(ierr);
  ierr = VecSetDM(vec,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMCreateGlobalVector_pforest _append_pforest(DMCreateGlobalVector)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateGlobalVector_pforest)
static PetscErrorCode DMCreateGlobalVector_pforest(DM dm,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateGlobalVector_Section_Private(dm,vec);CHKERRQ(ierr);
  /* ierr = VecSetOperation(*vec, VECOP_DUPLICATE, (void(*)(void)) VecDuplicate_MPI_DM);CHKERRQ(ierr); */
  ierr = VecSetOperation(*vec, VECOP_VIEW, (void (*)(void)) VecView_pforest);CHKERRQ(ierr);
  ierr = VecSetOperation(*vec, VECOP_VIEWNATIVE, (void (*)(void)) VecView_pforest_Native);CHKERRQ(ierr);
  ierr = VecSetOperation(*vec, VECOP_LOAD, (void (*)(void)) VecLoad_pforest);CHKERRQ(ierr);
  ierr = VecSetOperation(*vec, VECOP_LOADNATIVE, (void (*)(void)) VecLoad_pforest_Native);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if 0
#define DMCreateLocalVector_pforest _append_pforest(DMCreateLocalVector)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateLocalVector_pforest)
static PetscErrorCode DMCreateLocalVector_pforest(DM dm,Vec *vec)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(plex,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

#define DMCreateMatrix_pforest _append_pforest(DMCreateMatrix)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateMatrix_pforest)
static PetscErrorCode DMCreateMatrix_pforest(DM dm,Mat *mat)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMCreateMatrix(plex,mat);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMProjectFunctionLocal_pforest _append_pforest(DMProjectFunctionLocal)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMProjectFunctionLocal_pforest)
static PetscErrorCode DMProjectFunctionLocal_pforest(DM dm, PetscReal time, PetscErrorCode (**funcs)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void **ctxs, InsertMode mode, Vec localX)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMProjectFunctionLocal(plex,time,funcs,ctxs,mode,localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMProjectFunctionLabelLocal_pforest _append_pforest(DMProjectFunctionLabelLocal)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMProjectFunctionLabelLocal_pforest)
static PetscErrorCode DMProjectFunctionLabelLocal_pforest(DM dm, PetscReal time, DMLabel label, PetscInt numIds, const PetscInt ids[], PetscErrorCode (**funcs)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void **ctxs, InsertMode mode, Vec localX)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMProjectFunctionLabelLocal(plex,time,label,numIds,ids,funcs,ctxs,mode,localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMProjectFieldLocal_pforest _append_pforest(DMProjectFieldLocal)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMProjectFieldLocal_pforest)
PetscErrorCode DMProjectFieldLocal_pforest(DM dm, Vec localU,
                                           void (**funcs)(PetscInt, PetscInt, PetscInt,
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                                                          PetscReal, const PetscReal[], PetscScalar[]),
                                           InsertMode mode, Vec localX)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMProjectFieldLocal(plex,localU,funcs,mode,localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMComputeL2Diff_pforest _append_pforest(DMComputeL2Diff)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMComputeL2Diff_pforest)
PetscErrorCode DMComputeL2Diff_pforest(DM dm, PetscReal time, PetscErrorCode (**funcs)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void **ctxs, Vec X, PetscReal *diff)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMComputeL2Diff(plex,time,funcs,ctxs,X,diff);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMComputeL2FieldDiff_pforest _append_pforest(DMComputeL2FieldDiff)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMComputeL2FieldDiff_pforest)
PetscErrorCode DMComputeL2FieldDiff_pforest(DM dm, PetscReal time, PetscErrorCode (**funcs)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void **ctxs, Vec X, PetscReal diff[])
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMComputeL2FieldDiff(plex,time,funcs,ctxs,X,diff);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMCreateDefaultSection_pforest _append_pforest(DMCreateDefaultSection)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateDefaultSection_pforest)
static PetscErrorCode DMCreateDefaultSection_pforest(DM dm)
{
  DM                plex;
  PetscSection      section;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMGetDefaultSection(plex,&section);CHKERRQ(ierr);
  ierr = DMSetDefaultSection(dm,section);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMCreateDefaultConstraints_pforest _append_pforest(DMCreateDefaultConstraints)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreateDefaultConstraints_pforest)
static PetscErrorCode DMCreateDefaultConstraints_pforest(DM dm)
{
  DM                plex;
  Mat               mat;
  PetscSection      section;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMGetDefaultConstraints(plex,&section,&mat);CHKERRQ(ierr);
  ierr = DMSetDefaultConstraints(dm,section,mat);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMGetDimPoints_pforest _append_pforest(DMGetDimPoints)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMGetDimPoints_pforest)
static PetscErrorCode DMGetDimPoints_pforest(DM dm, PetscInt dim, PetscInt *cStart, PetscInt *cEnd)
{
  DM                plex;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = DMPforestGetPlex(dm,&plex);CHKERRQ(ierr);
  ierr = DMGetDimPoints(plex,dim,cStart,cEnd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Need to forward declare */
#define DMInitialize_pforest _append_pforest(DMInitialize)
static PetscErrorCode DMInitialize_pforest(DM dm);

#define DMClone_pforest _append_pforest(DMClone)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMClone_pforest)
static PetscErrorCode DMClone_pforest(DM dm, DM *newdm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMClone_Forest(dm,newdm);CHKERRQ(ierr);
  ierr = DMInitialize_pforest(*newdm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMForestCreateCellChart_pforest _append_pforest(DMForestCreateCellChart)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMForestCreateCellChart_pforest)
static PetscErrorCode DMForestCreateCellChart_pforest(DM dm, PetscInt *cStart, PetscInt *cEnd)
{
  DM_Forest         *forest;
  DM_Forest_pforest *pforest;
  PetscInt          overlap;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  forest  = (DM_Forest *) dm->data;
  pforest = (DM_Forest_pforest *) forest->data;
  *cStart = 0;
  ierr = DMForestGetPartitionOverlap(dm,&overlap);CHKERRQ(ierr);
  if (overlap && pforest->ghost) {
    *cEnd = pforest->forest->local_num_quadrants + pforest->ghost->proc_offsets[pforest->forest->mpisize];
  }
  else {
    *cEnd = pforest->forest->local_num_quadrants;
  }
  PetscFunctionReturn(0);
}

#define DMForestCreateCellSF_pforest _append_pforest(DMForestCreateCellSF)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMForestCreateCellSF_pforest)
static PetscErrorCode DMForestCreateCellSF_pforest(DM dm, PetscSF *cellSF)
{
  DM_Forest         *forest;
  DM_Forest_pforest *pforest;
  PetscMPIInt       rank;
  PetscInt          overlap;
  PetscInt          cStart, cEnd, cLocalStart, cLocalEnd;
  PetscInt          nRoots, nLeaves;
  PetscSFNode       *remote;
  PetscSF           sf;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = DMForestGetCellChart(dm,&cStart,&cEnd);CHKERRQ(ierr);
  forest  = (DM_Forest *)         dm->data;
  pforest = (DM_Forest_pforest *) forest->data;
  nRoots  = nLeaves = cEnd - cStart;
  cLocalStart = pforest->cLocalStart;
  cLocalEnd   = pforest->cLocalEnd;
  ierr = DMForestGetPartitionOverlap(dm,&overlap);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
  ierr = PetscMalloc1(cEnd-cStart,&remote);CHKERRQ(ierr);
  if (overlap && pforest->ghost) {
    PetscSFNode *mirror, *ghost;
    p4est_quadrant_t *mirror_array;
    PetscInt nMirror, nGhost, nGhostPre, nGhostPost, nSelf, q;
    void **mirrorPtrs;

    nMirror    = (PetscInt) pforest->ghost->mirrors.elem_count;
    nSelf      = cLocalEnd - cLocalStart;
    nGhost     = (cEnd - cStart) - nSelf;
    nGhostPre  = (PetscInt) pforest->ghost->proc_offsets[rank];
    nGhostPost = nGhost - nGhostPre;
    ierr = PetscMalloc3(nMirror,&mirror,nMirror,&mirrorPtrs,nGhost,&ghost);CHKERRQ(ierr);
    mirror_array = (p4est_quadrant_t *) pforest->ghost->mirrors.array;
    for (q = 0; q < nMirror; q++) {
      p4est_quadrant_t *mir = &(mirror_array[q]);

      mirror[q].rank = rank;
      mirror[q].index = (PetscInt) mir->p.piggy3.local_num;
      mirrorPtrs[q] = (void *) &(mirror[q]);
    }
    PetscStackCallP4est(p4est_ghost_exchange_custom,(pforest->forest,pforest->ghost,sizeof(PetscSFNode),mirrorPtrs,ghost));
    ierr = PetscFree3(mirror,mirrorPtrs,ghost);CHKERRQ(ierr);
    ierr = PetscMemcpy(remote,ghost,nGhostPre * sizeof(PetscSFNode));CHKERRQ(ierr);
    for (q = cLocalStart; q < cLocalEnd; q++) {
      remote[q].rank  = rank;
      remote[q].index = q;
    }
    ierr = PetscMemcpy(&remote[cLocalEnd],&ghost[nGhostPre],nGhostPost * sizeof(PetscSFNode));CHKERRQ(ierr);
  }
  else {
    PetscInt q;

    for (q = 0; q < cEnd; q++) {
      remote[q].rank = rank;
      remote[q].index = q;
    }
  }
  ierr = PetscSFCreate(PetscObjectComm((PetscObject)dm),&sf);CHKERRQ(ierr);
  ierr = PetscSFSetGraph(sf,nRoots,nLeaves,NULL,PETSC_OWN_POINTER,remote,PETSC_OWN_POINTER);CHKERRQ(ierr);
  *cellSF = sf;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMInitialize_pforest)
static PetscErrorCode DMInitialize_pforest(DM dm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dm->ops->setup                     = DMSetUp_pforest;
  dm->ops->view                      = DMView_pforest;
  dm->ops->clone                     = DMClone_pforest;
  dm->ops->coarsen                   = DMCoarsen_pforest;
  dm->ops->createinterpolation       = DMCreateInterpolation_pforest;
  dm->ops->setfromoptions            = DMSetFromOptions_pforest;
  dm->ops->createcoordinatedm        = DMCreateCoordinateDM_pforest;
  dm->ops->createglobalvector        = DMCreateGlobalVector_pforest;
  dm->ops->createlocalvector         = DMCreateLocalVector_Section_Private;
  dm->ops->creatematrix              = DMCreateMatrix_pforest;
  dm->ops->projectfunctionlocal      = DMProjectFunctionLocal_pforest;
  dm->ops->projectfunctionlabellocal = DMProjectFunctionLabelLocal_pforest;
  dm->ops->projectfieldlocal         = DMProjectFieldLocal_pforest;
  dm->ops->createdefaultsection      = DMCreateDefaultSection_pforest;
  dm->ops->createdefaultconstraints  = DMCreateDefaultConstraints_pforest;
  dm->ops->computel2diff             = DMComputeL2Diff_pforest;
  dm->ops->computel2fielddiff        = DMComputeL2FieldDiff_pforest;
  dm->ops->getdimpoints              = DMGetDimPoints_pforest;
  ierr = PetscObjectComposeFunction((PetscObject)dm,_pforest_string(DMConvert_plex_pforest) "_C",DMConvert_plex_pforest);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)dm,_pforest_string(DMConvert_pforest_plex) "_C",DMConvert_pforest_plex);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define DMCreate_pforest _append_pforest(DMCreate)
#undef __FUNCT__
#define __FUNCT__ _pforest_string(DMCreate_pforest)
PETSC_EXTERN PetscErrorCode DMCreate_pforest(DM dm)
{
  DM_Forest         *forest;
  DM_Forest_pforest *pforest;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscP4estInitialize();CHKERRQ(ierr);
  ierr = DMCreate_Forest(dm);CHKERRQ(ierr);
  ierr = DMInitialize_pforest(dm);CHKERRQ(ierr);
  ierr = DMSetDimension(dm,P4EST_DIM);CHKERRQ(ierr);

  /* set forest defaults */
  ierr = DMForestSetTopology(dm,"unit");CHKERRQ(ierr);
  ierr = DMForestSetMinimumRefinement(dm,0);CHKERRQ(ierr);
  ierr = DMForestSetInitialRefinement(dm,0);CHKERRQ(ierr);
  ierr = DMForestSetMaximumRefinement(dm,P4EST_QMAXLEVEL);CHKERRQ(ierr);
  ierr = DMForestSetGradeFactor(dm,2);CHKERRQ(ierr);
  ierr = DMForestSetAdjacencyDimension(dm,0);CHKERRQ(ierr);
  ierr = DMForestSetPartitionOverlap(dm,0);CHKERRQ(ierr);

  /* create p4est data */
  ierr = PetscNewLog(dm,&pforest);CHKERRQ(ierr);

  forest                            = (DM_Forest *) dm->data;
  forest->data                      = pforest;
  forest->destroy                   = DMForestDestroy_pforest;
  forest->ftemplate                 = DMForestTemplate_pforest;
  forest->createcellchart           = DMForestCreateCellChart_pforest;
  forest->createcellsf              = DMForestCreateCellSF_pforest;
  pforest->topo                     = NULL;
  pforest->forest                   = NULL;
  pforest->ghost                    = NULL;
  pforest->lnodes                   = NULL;
  pforest->partition_for_coarsening = PETSC_TRUE;
  pforest->coarsen_hierarchy        = PETSC_FALSE;
  pforest->cLocalStart              = -1;
  pforest->cLocalEnd                = -1;
  pforest->labelsFinalized          = PETSC_FALSE;
  pforest->ghostName                = NULL;
  PetscFunctionReturn(0);
}

#endif /* defined(PETSC_HAVE_P4EST) */
