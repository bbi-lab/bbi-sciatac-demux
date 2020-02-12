#!/bin/bash

#
# Notes:
#   o  run in sciatac_pipeline directory
#   o  when running on cpuid_level=22 nodes in a qlogin session, use the
#      qlogin parameters -l mfree=16G -l cpuid_level=22
#   o  the pypy virtual environment is used only to run 'barcode_correct_sciatac.py'
#
module purge
module load modules modules-init modules-gs 
export DIR=$(dirname `readlink -f $0`)
source $DIR/load_pypy_env_reqs.sh
module load git/latest

echo 'Cleaning cache directories...'
rm -rf ./__pycache__/* src/__pycache__/*  ~/.cache/*

if [ -d $DIR/src/pypy_env ]; then
        echo 'Removing existing pypy virtualenv...'
        rm -rf $DIR/src/pypy_env
fi

# First, the pypy virtualenv
echo 'Building pypy virtualenv...'
export PYTHONPATH=''
virtualenv-pypy $DIR/src/pypy_env

source $DIR/src/pypy_env/bin/activate

pip install -r $DIR/pypy_requirements.txt

# git clone https://github.com/andrewhill157/barcodeutils.git

pushd barcodeutils
pypy setup.py install
popd

deactivate

# Second, the python virtualenv
module purge
module load modules modules-init modules-gs
export DIR=$(dirname `readlink -f $0`)
source $DIR/load_python_env_reqs.sh
module load virtualenv/16.0.0
# 
echo 'Cleaning cache directory...'
rm -rf ./__pycache__/* src/__pycache__/* ~/.cache/*
# 
if [ -d $DIR/src/python_env ]; then
        echo 'Removing existing virtualenv...'
	rm -rf $DIR/src/python_env
fi
# 
echo 'Bulding python virtualenv...'
virtualenv $DIR/src/python_env
export PYTHONPATH=''
source $DIR/src/python_env/bin/activate
# 
pip install -r $DIR/python_requirements.txt

pushd barcodeutils
python setup.py install
popd

deactivate
