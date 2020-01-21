#!/bin/bash

######################################################
# input parameters
######################################################
run=$1
dataset=$2

######################################################
# variable parameters
# Change these as needed
######################################################
data_path="/Volumes/xen01/analysis/xurich/DPv1.1/Run"
ml_path="/Users/marrodan/scripts"
mr_path="/Users/marrodan/Xurich/MatROOT_convert"

######################################################
# Output parameters
######################################################

if [ $2 ]
then

${ml_path}/mlConvert ${data_path}${run}/${dataset}.mat 
${mr_path}/exeROOTreadBin091008 ${data_path}${run}/${dataset}.dat ${data_path}${run}/${dataset}.root
rm -f ${data_path}${run}/${dataset}.dat

######################################################
# Input Error Message
######################################################

else
    echo "Script requires 2 parameters:"
    echo "    parameter1: Run (two digits)"
    echo "    parameter2: Dataset (without extension)"
fi
