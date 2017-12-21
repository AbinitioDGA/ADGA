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
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

#!/bin/bash
#Check if we already have a RSPTmake.inc
if [ -f ../make_config ]
then
   exit
fi

#host
myhost=`hostname |awk -F"[0-9.]" '{print $1}'`

#cluster (Add your favourite cluster in the list below!)
cluster="default"
if [[ "$OSTYPE" =~ "darwin" ]]
then
   # Catch Macs
   cluster="mac"
elif [ "$myhost" == "n" ] # hclm is running on node n101
then
   cluster="hclm"
elif [ "$myhost" == "l" ] # vsc3 is running on login nodes l31 to l35
then
   cluster="vsc3"
else
   # Fall back
   cluster="$myhost"
fi

# Check if we already have a working template
template=""
if [ -n "$cluster" ] 
then
   template=`ls -1 make_config_*|grep -i "make_config_$cluster" | head -1`
fi

# Generate the RSPTmake.inc file
if [ -n "$template" ]
then
   cp $template ../make_config
else
   echo "make_config_auto.sh failed to find a suitable make_config. Please look at the available templates in make_configs."
fi

