#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
#          --------------------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, April 2014
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

# Makefile.  Generated from Makefile.in by configure.

#######################################

FC = mpif90
FCFLAGS = -g
FC_DEFINE = -D
MPIFC = mpif90
MPILIBS = 

FLAGS_CHECK =   -DFORCE_VECTORIZATION

FCFLAGS_f90 = -mod ./obj -I./obj -I. -I${SETUP}

FC_MODEXT = mod
FC_MODDIR = ./obj

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK)

MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK)

CC = gcc
CFLAGS = -g -O2 
CPPFLAGS = -I${SETUP}  -DFORCE_VECTORIZATION

FCLINK = $(MPIFCCOMPILE_CHECK)


#######################################
####
#### GPU
#### with configure: ./configure --with-cuda=cuda5 CUDA_FLAGS=.. CUDA_LIB=.. CUDA_INC=.. MPI_INC=.. ..
#### with configure: ./configure --with-opencl OCL_GPU_FLAGS=.. OCL_LIB=.. OCL_INC=.. MPI_INC=.. ..
####
#######################################

#CUDA = yes
CUDA = no

#CUDA5 = yes
CUDA5 = no

#OCL = yes
OCL = no

MPI_INCLUDES = 

CUDA_FLAGS = 
CUDA_INC = 
CUDA_LINK =  
CUDA_DEBUG = --cudart=shared

#NVCC = nvcc
NVCC = gcc

OCL_CPU_FLAGS = 
OCL_GPU_FLAGS = 

OCL_INC = 
OCL_LINK =  

ifeq ($(OCL), yes)
  ifeq ($(CUDA), yes)
    GPU_CUDA_AND_OCL = yes
  endif
endif

ifeq ($(OCL), no)
  ifeq ($(CUDA), no)
    NO_GPU = yes
  endif
endif

ifneq ($(NO_GPU), yes)
  HAS_GPU = yes
endif

# GPU architecture

# CUDA architecture / code version
# Fermi: -gencode=arch=compute_10,code=sm_10 not supported
# Tesla (default): -gencode=arch=compute_20,code=sm_20
# Geforce GT 650m: -gencode=arch=compute_30,code=sm_30
# Kepler (cuda5) : -gencode=arch=compute_35,code=sm_35
GENCODE_20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE_30 = -gencode=arch=compute_30,code=\"sm_30,compute_30\"
GENCODE_35 = -gencode=arch=compute_35,code=\"sm_35,compute_35\"

# CUDA version 5.x
##GENCODE = $(GENCODE_35)
# CUDA version 4.x
#GENCODE = $(GENCODE_20)

# CUDA flags and linking
#NVCC_FLAGS_BASE = $(CUDA_FLAGS) $(CUDA_INC) $(CUDA_DEBUG) $(MPI_INCLUDES)
##NVCC_FLAGS = $(NVCC_FLAGS_BASE) -dc $(GENCODE)
#NVCC_FLAGS = $(NVCC_FLAGS_BASE) -DUSE_OLDER_CUDA4_GPU $(GENCODE)

##NVCCLINK_BASE = $(NVCC) $(CUDA_INC) $(MPI_INCLUDES)
##NVCCLINK = $(NVCCLINK_BASE) -dlink $(GENCODE)
#NVCCLINK = $(NVCCLINK_BASE) -DUSE_OLDER_CUDA4_GPU $(GENCODE)

NVCC_FLAGS =
NVCCLINK = $(NVCC) $(NVCC_FLAGS)

#######################################
####
#### VTK
#### with configure: ./configure --enable-vtk ..
####
#######################################

#VTK = yes
VTK = no

CPPFLAGS += 
LDFLAGS += 
MPILIBS += 


#######################################
####
#### ADIOS
#### with configure: ./configure --with-adios ADIOS_CONFIG=..
####
#######################################

#ADIOS = yes
ADIOS = no

FCFLAGS += 
MPILIBS += 


#######################################
####
#### FORCE_VECTORIZATION
#### with configure: ./configure --with-vec ..
####
#######################################
FORCE_VECTORIZATION = yes
#FORCE_VECTORIZATION = no


#######################################
####
#### CEM
#### with configure: ./configure --with-cem CEM_LIBS=.. CEM_FCFLAGS=..
####
#######################################

#CEM = yes
CEM = no

FCFLAGS += 
MPILIBS += 


#######################################
####
#### directories
####
#######################################

## compilation directories
# B : build directory
B = .
# E : executables directory
E = $B/bin
# O : objects directory
O = $B/obj
# S_TOP : source file root directory
S_TOP = .
# setup file directory
SETUP = $B/setup
# output file directory
OUTPUT = $B/OUTPUT_FILES


#######################################
####
#### targets
####
#######################################

# code subdirectories
SUBDIRS = \
	shared \
	create_header_file \
	meshfem3D \
	specfem3D \
	auxiliaries \
	postprocess \
	tomography \
	$(EMPTY_MACRO)

ifeq ($(HAS_GPU),yes)
  SUBDIRS := gpu $(SUBDIRS)
endif

# default targets
DEFAULT = \
	xcreate_header_file \
	xmeshfem3D \
	xspecfem3D \
	xcombine_AVS_DX \
	xcombine_surf_data \
	xcombine_vol_data \
	xcombine_vol_data_vtk \
	xconvolve_source_timefunction \
	xcreate_movie_AVS_DX \
	xcreate_movie_GMT_global \
	$(EMPTY_MACRO)

ifeq ($(ADIOS), yes)
DEFAULT += 	\
	xcombine_vol_data_adios \
	xcombine_vol_data_vtk_adios \
	$(EMPTY_MACRO)
endif

all: default aux movies postprocess tomography

default: $(DEFAULT)

ifdef CLEAN
clean:
	@echo "cleaning by CLEAN"
	-rm -f $(foreach dir, $(CLEAN), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_SHARED_OBJECTS) $($(dir)_TARGETS))
else
clean:
	@echo "cleaning all"
	-rm -f $(foreach dir, $(SUBDIRS), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_TARGETS)) $O/*
endif

realclean: clean
mrproper: clean

help:
	@echo "usage: make [executable]"
	@echo ""
	@echo "supported main executables:"
	@echo "    xmeshfem3D"
	@echo "    xspecfem3D"
	@echo ""
	@echo "defaults:"
	@echo "    xcreate_header_file"
	@echo "    xmeshfem3D"
	@echo "    xspecfem3D"
	@echo "    xcombine_AVS_DX"
	@echo "    xcombine_surf_data"
	@echo "    xcombine_vol_data"
	@echo "    xcombine_vol_data_vtk"
	@echo "    xconvolve_source_timefunction"
	@echo "    xcreate_movie_AVS_DX"
	@echo "    xcreate_movie_GMT_global"
	@echo ""
	@echo "additional executables:"
	@echo "- auxiliary executables: [make aux]"
	@echo "    xcombine_vol_data"
	@echo "    xcombine_vol_data_vtk"
ifeq ($(ADIOS), yes)
	@echo "    xcombine_vol_data_adios"
	@echo "    xcombine_vol_data_vtk_adios"
endif
	@echo "    xcombine_surf_data"
	@echo "    xcombine_AVS_DX"
	@echo "    xconvolve_source_timefunction"
	@echo ""
	@echo "- movie executables: [make movies]"
	@echo "    xcreate_movie_AVS_DX"
	@echo "    xcreate_movie_GMT_global"
	@echo "    xcombine_paraview_strain_data"
	@echo ""
	@echo "- postprocessing tools: [make postprocess]"
	@echo "    xaddition_sem"
	@echo "    xclip_sem"
	@echo "    xcombine_sem"
	@echo "    xdifference_sem"
	@echo "    xinterpolate_model"
	@echo "    xsmooth_sem"
ifeq ($(ADIOS), yes)
	@echo "    xconvert_model_file_adios"
endif
	@echo ""
	@echo "- tomography tools: [make tomography]"
	@echo "    xadd_model_iso"
	@echo "    xadd_model_tiso"
	@echo "    xadd_model_tiso_cg"
	@echo "    xadd_model_tiso_iso"
	@echo "    xsum_kernels"
	@echo "    xsum_preconditioned_kernels"
	@echo ""

.PHONY: all default clean help

#######################################


# Get dependencies and rules for building stuff
include $(patsubst %, ${S_TOP}/src/%/rules.mk, $(SUBDIRS))


#######################################

##
## Shortcuts
##

# Shortcut for: <prog>/<xprog> -> bin/<xprog>
define target_shortcut
$(patsubst $E/%, %, $(1)): $(1)
.PHONY: $(patsubst $E/%, %, $(1))
$(patsubst $E/x%, %, $(1)): $(1)
.PHONY: $(patsubst $E/x%, %, $(1))
endef

# Shortcut for: dir -> src/dir/<targets in here>
define shortcut
$(1): $($(1)_TARGETS)
.PHONY: $(1)
$$(foreach target, $$(filter $E/%,$$($(1)_TARGETS)), $$(eval $$(call target_shortcut,$$(target))))
endef

$(foreach dir, $(SUBDIRS), $(eval $(call shortcut,$(dir))))

# Other old shortcuts
mesh: $E/xmeshfem3D
spec: $E/xspecfem3D
.PHONY: mesh spec

