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

