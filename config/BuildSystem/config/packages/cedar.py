import config.package
import os
import traceback

class Configure(config.package.CMakePackage):
  def __init__(self, framework):
    config.package.CMakePackage.__init__(self, framework)
    self.versionname = 'CEDAR_GIT_VERSION'
    self.gitcommit = '8642979198120b3546fe43d6eb74e3ee9bc18b41' # master commit, as of 12 June 2023
    self.download = ['git://https://github.com/cedar-framework/cedar.git']
    self.includes = ['cedar/capi.h']
    self.liblist = [['libcedar.a']]
    self.functions = ['cedar_config_create',
                      'cedar_log_init',
                      'cedar_topo_create2d',
                      'cedar_mat_create2d',
                      'cedar_topo_create3d',
                      'cedar_mat_create3d']
    self.hastests = 1
    self.hastestsdatafiles = 1
    self.precisions = ['double']
    self.buildLanguages = ['Cxx']
    self.minCxxVersion = 'c++14'
    self.builtafterpetsc = 0
    self.minCmakeVersion = (3,14,0)

  def setupHelp(self, help):
    import nargs
    config.package.CMakePackage.setupHelp(self, help)

  def setupDependencies(self, framework):
    config.package.CMakePackage.setupDependencies(self, framework)
    self.blasLapack     = framework.require('config.packages.BlasLapack',self)
    self.mpi            = framework.require('config.packages.MPI',self)
    self.kokkos         = framework.require('config.packages.kokkos',self)
    self.odeps          = [self.kokkos]
    self.deps           = [self.mpi,self.blasLapack]

  def formCMakeConfigureArgs(self):
    args = config.package.CMakePackage.formCMakeConfigureArgs(self)
    return args

  def configureLibrary(self):
    config.package.Package.configureLibrary(self)
    self.addDefine('HAVE_CEDAR', 1)
