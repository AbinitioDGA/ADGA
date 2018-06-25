# This file is part of the Abinitio Dynamical Vertex Approximation (ADGA)
# package. It is an electronic structure code which allows the inclusion of
# non-local correlations beyond DMFT and the calculation of momentum-dependent
# susceptibilities.
#
# The public repository can be found at
# https://github.com/AbinitioDGA/ADGA
#
# The arXiv publication can be found at
# https://arxiv.org/abs/1710.06651
#
# Copyright (C) <2017, 2018>
# <Anna Galler*, Patrick Thunström, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak, Karsten Held>
# * Corresponding author. E-mail address: galler.anna@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

SHELL=/bin/bash

PROG = abinitiodga
PREPROG = setupvertex

include ../make_config

INSTALL_DIR = ../bin

MAIN_OBJECTS := main.o parameters_module.o interaction_module.o kq_tools.o lapack_module.o one_particle_quant_module.o \
		    susc_module.o eom_module.o hdf5_module.o aux.o mpi_org.o lookup_module.o config_module.o
VERTEX_OBJECTS := vertex_chann_sym.o hdf5_module.o parameters_module.o aux.o kq_tools.o

.SUFFIXES:
.SUFFIXES: .o .f90 .F90

%.o: %.f90
	$(F90) $(FFLAGS) $(FINCLUDE) -c $< -o $@

%.o: %.F90
	$(F90) $(FFLAGS) $(FPPFLAGS) $(FINCLUDE) -c $< -o $@

.PHONY: all
all: $(INSTALL_DIR) $(INSTALL_DIR)/$(PROG) $(INSTALL_DIR)/$(PREPROG)

$(INSTALL_DIR):
	if [ ! -d $(INSTALL_DIR) ]; then mkdir -p $(INSTALL_DIR) ; fi

$(INSTALL_DIR)/$(PROG): $(MAIN_OBJECTS)
	$(LD) $^ -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

$(INSTALL_DIR)/$(PREPROG): $(VERTEX_OBJECTS)
	$(LD) $^ -o $@ $(FFLAGS) $(LDINCLUDE) $(LDFLAGS)

main.o: parameters_module.o aux.o kq_tools.o lapack_module.o one_particle_quant_module.o susc_module.o \
				eom_module.o interaction_module.o hdf5_module.o mpi_org.o config_module.o
vertex_chann_sym.o: parameters_module.o hdf5_module.o aux.o kq_tools.o

eom_module.o: parameters_module.o one_particle_quant_module.o
susc_module.o: parameters_module.o hdf5_module.o
one_particle_quant_module.o: aux.o lapack_module.o parameters_module.o mpi_org.o
kq_tools.o: parameters_module.o
interaction_module.o: parameters_module.o hdf5_module.o aux.o
hdf5_module.o: parameters_module.o kq_tools.o aux.o
mpi_org.o: parameters_module.o
lookup_module.o: parameters_module.o
config_module.o: parameters_module.o lookup_module.o

.PHONY: clean pristine
clean:
	rm -f *.o *.mod

pristine:
	rm -f *.o *.mod
	rm -f $(INSTALL_DIR)/$(PROG)
	rm -f $(INSTALL_DIR)/$(PREPROG)
