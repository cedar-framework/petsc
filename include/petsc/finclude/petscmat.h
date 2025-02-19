!
!
!  Include file for Fortran use of the Mat package in PETSc
!  Portions of this code are under:
!  Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
!
#if !defined (PETSCMATDEF_H)
#define PETSCMATDEF_H

#include "petsc/finclude/petscvec.h"

#define Mat type(tMat)
#define MatNullSpace type(tMatNullSpace)
#define MatFDColoring type(tMatFDColoring)

#define MatColoring PetscFortranAddr
#define MatPartitioning PetscFortranAddr
#define MatCoarsen PetscFortranAddr
#define MatAIJIndices PetscFortranAddr
#define MatType character*(80)
#define MatSolverType character*(80)
#define MatOption PetscEnum
#define MatCreateSubMatrixOption PetscEnum
#define MPChacoGlobalType PetscEnum
#define MPChacoLocalType PetscEnum
#define MPChacoEigenType PetscEnum
#define MPPTScotchStragegyType PetscEnum
#define MatAssemblyType PetscEnum
#define MatFactorType PetscEnum
#define MatFactorError PetscEnum
#define MatFactorShiftType PetscEnum
#define MatProductType PetscEnum
#define MatProductAlgorithm character*(80)
#define MatFactorSchurStatus PetscEnum
#define MatOrderingType character*(80)
#define MatSORType PetscEnum
#define MatInfoType PetscEnum
#define MatReuse PetscEnum
#define MatOperation PetscEnum
#define MatColoringType character*(80)
#define MatInfo PetscLogDouble
#define MatFactorInfo PetscReal
#define MatDuplicateOption PetscEnum
#define MatStructure PetscEnum
#define MatPartitioningType character*(80)
#define MatCoarsenType character*(80)
#define MatCompositeType PetscEnum
#define MatCompositeMergeType PetscEnum
#define MatStencil PetscInt
#define MatStencil_k 1
#define MatStencil_j 2
#define MatStencil_i 3
#define MatStencil_c 4

#define MATPARTITIONING_CURRENT 'current'
#define MATPARTITIONING_PARMETIS 'parmetis'

#define MATCOARSEN_MIS 'mis'

#define MATCOLORINGJP      'jp'
#define MATCOLORINGPOWER   'power'
#define MATCOLORINGNATURAL 'natural'
#define MATCOLORINGSL      'sl'
#define MATCOLORINGLF      'lf'
#define MATCOLORINGID      'id'
#define MATCOLORINGGREEDY  'greedy'

#define MATORDERINGNATURAL   'natural'
#define MATORDERINGNATURAL_OR_ND 'natural_or_nd'
#define MATORDERINGND        'nd'
#define MATORDERING1WD       '1wd'
#define MATORDERINGRCM       'rcm'
#define MATORDERINGQMD       'qmd'
#define MATORDERINGROWLENGTH 'rowlength'
#define MATORDERINGWBM       'wbm'
#define MATORDERINGSPECTRAL  'spectral'
#define MATORDERINGAMD       'amd'
#define MATORDERINGEXTERNAL  'external'
!
!  Matrix types
!
#define MATSAME            'same'
#define MATMAIJ            'maij'
#define MATSEQMAIJ         'seqmaij'
#define MATMPIMAIJ         'mpimaij'
#define MATKAIJ            'kaij'
#define MATSEQKAIJ         'seqkaij'
#define MATMPIKAIJ         'mpikaij'
#define MATIS              'is'
#define MATAIJ             'aij'
#define MATSEQAIJ          'seqaij'
#define MATMPIAIJ          'mpiaij'
#define MATAIJCRL          'aijcrl'
#define MATSEQAIJCRL       'seqaijcrl'
#define MATMPIAIJCRL       'mpiaijcrl'
#define MATAIJCUSPARSE     'aijcusparse'
#define MATSEQAIJCUSPARSE  'seqaijcusparse'
#define MATMPIAIJCUSPARSE  'mpiaijcusparse'
#define MATAIJHIPSPARSE    'aijhipsparse'
#define MATSEQAIJHIPSPARSE 'seqaijhipsparse'
#define MATMPIAIJHIPSPARSE 'mpiaijhipsparse'
#define MATAIJKOKKOS       'aijkokkos'
#define MATSEQAIJKOKKOS    'seqaijkokkos'
#define MATMPIAIJKOKKOS    'mpiaijkokkos'
#define MATAIJVIENNACL     'aijviennacl'
#define MATSEQAIJVIENNACL  'seqaijviennacl'
#define MATMPIAIJVIENNACL  'mpiaijviennacl'
#define MATAIJPERM         'aijperm'
#define MATSEQAIJPERM      'seqaijperm'
#define MATMPIAIJPERM      'mpiaijperm'
#define MATAIJSELL         'aijsell'
#define MATSEQAIJSELL      'seqaijsell'
#define MATMPIAIJSELL      'mpiaijsell'
#define MATAIJMKL          'aijmkl'
#define MATSEQAIJMKL       'seqaijmkl'
#define MATMPIAIJMKL       'mpiaijmkl'
#define MATBAIJMKL         'baijmkl'
#define MATSEQBAIJMKL      'seqbaijmkl'
#define MATMPIBAIJMKL      'mpibaijmkl'
#define MATSHELL           'shell'
#define MATCENTERING       'centering'
#define MATDENSE           'dense'
#define MATDENSECUDA       'densecuda'
#define MATDENSEHIP        'densehip'
#define MATSEQDENSE        'seqdense'
#define MATSEQDENSECUDA    'seqdensecuda'
#define MATSEQDENSEHIP     'seqdensehip'
#define MATMPIDENSE        'mpidense'
#define MATMPIDENSECUDA    'mpidensecuda'
#define MATMPIDENSEHIP     'mpidensehip'
#define MATELEMENTAL       'elemental'
#define MATSCALAPACK       'scalapack'
#define MATBAIJ            'baij'
#define MATSEQBAIJ         'seqbaij'
#define MATMPIBAIJ         'mpibaij'
#define MATMPIADJ          'mpiadj'
#define MATSBAIJ           'sbaij'
#define MATSEQSBAIJ        'seqsbaij'
#define MATMPISBAIJ        'mpisbaij'
#define MATMFFD            'mffd'
#define MATNORMAL          'normal'
#define MATNORMALHERMITIAN 'normalh'
#define MATLRC             'lrc'
#define MATSCATTER         'scatter'
#define MATBLOCKMAT        'blockmat'
#define MATCOMPOSITE       'composite'
#define MATFFT             'fft'
#define MATFFTW            'fftw'
#define MATSEQCUFFT        'seqcufft'
#define MATTRANSPOSEVIRTUAL       'transpose'
#define MATHERMITIANTRANSPOSEVIRTUAL 'hermitiantranspose'
#define MATSCHURCOMPLEMENT 'schurcomplement'
#define MATPYTHON          'python'
#define MATHYPRE           'hypre'
#define MATHYPRESTRUCT     'hyprestruct'
#define MATHYPRESSTRUCT    'hypresstruct'
#define MATSUBMATRIX       'submatrix'
#define MATLOCALREF        'localref'
#define MATNEST            'nest'
#define MATPREALLOCATOR    'preallocator'
#define MATSELL            'sell'
#define MATSEQSELL         'seqsell'
#define MATMPISELL         'mpisell'
#define MATDUMMY           'dummy'
#define MATLMVM            'lmvm'
#define MATLMVMDFP         'lmvmdfp'
#define MATLMVMBFGS        'lmvmbfgs'
#define MATLMVMSR1         'lmvmsr1'
#define MATLMVMBROYDEN     'lmvmbroyden'
#define MATLMVMBADBROYDEN  'lmvmbadbroyden'
#define MATLMVMSYMBROYDEN  'lmvmsymbroyden'
#define MATLMVMSYMBADBROYDEN 'lmvmsymbadbroyden'
#define MATLMVMDIAGBROYDEN 'lmvmdiagbroyden'
#define MATCONSTANTDIAGONAL 'constantdiagonal'
#define MATHTOOL           'htool'
#define MATH2OPUS          'h2opus'

!
! MatMFFDType values
!
#define MATMFFD_DS 'ds'
#define MATMFFD_WP 'wp'

!
! MatSolverTypes
!
#define MATSOLVERSUPERLU         'superlu'
#define MATSOLVERSUPERLU_DIST    'superlu_dist'
#define MATSOLVERSTRUMPACK       'strumpack'
#define MATSOLVERUMFPACK         'umfpack'
#define MATSOLVERCHOLMOD         'cholmod'
#define MATSOLVERKLU             'klu'
#define MATSOLVERELEMENTAL       'elemental'
#define MATSOLVERSCALAPACK       'scalapack'
#define MATSOLVERESSL            'essl'
#define MATSOLVERLUSOL           'lusol'
#define MATSOLVERMUMPS           'mumps'
#define MATSOLVERMKL_PARDISO     'mkl_pardiso'
#define MATSOLVERMKL_CPARDISO    'mkl_cpardiso'
#define MATSOLVERPASTIX          'pastix'
#define MATSOLVERMATLAB          'matlab'
#define MATSOLVERPETSC           'petsc'
#define MATSOLVERBAS             'bas'
#define MATSOLVERCUSPARSE        'cusparse'
#define MATSOLVERCUDA            'cuda'
#define MATSOLVERHIPSPARSE       'hipsparse'
#define MATSOLVERHIP             'hip'
#define MATSOLVERKOKKOS          'kokkos'
#define MATSOLVERSPQR            'spqr'

!
! GPU Storage Formats for CUSPARSE
!
#define MatCUSPARSEStorageFormat PetscEnum
#define MatCUSPARSEFormatOperation PetscEnum

!
! GPU Storage Formats for HIPSPARSE
!
#define MatHIPSPARSEStorageFormat PetscEnum
#define MatHIPSPARSEFormatOperation PetscEnum

!
! sparsity reducing ordering for STRUMPACK
!
#define MatSTRUMPACKReordering PetscEnum

#endif
