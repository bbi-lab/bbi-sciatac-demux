#!/bin/bash

#
# Notes:
#   o  run in sciatac_pipeline directory
#   o  the Shendure cluster has nodes with Intel cpuid_level=11, cpuid_level=13,
#      and cpuid_level=22.
#      The cpuid_level 22 nodes have instructions (vectorized) that are not
#      part of the cpuid_level 11 architecture. So certain software built on
#      cpuid_level 22 nodes may not run on cpuid_level 11 nodes (I am guessing that
#      the cpuid_level 22 instructions are a superset of those on the cpuid_level
#      11 nodes.) This seems to affect at least the numpy module, if I remember
#      correctly. I suggest installing the modules using a cpuid_level=11 node
#      so that this pipeline runs on all Shendure cluster nodes.
#      The Trapnell cluster nodes are all AMD so this note is irrelevant for
#      this pipeline when installed and run on Trapnell cluster nodes.
#      In order to install this pipeline on a cpuid_level 11 node, use a shell
#      obtained with the command
#
#        qlogin -l mfree=16G -l cpuid_level=11
#
#   o  the pypy virtual environment is used only to run 'barcode_correct_sciatac.py'
#

#
# Prepare pypy virtual environment.
#
module purge
module load modules modules-init modules-gs
export DIR=$(dirname `readlink -f $0`)
source $DIR/load_pypy_env_reqs.sh
module load git/2.18.0

echo 'Cleaning cache directories...'
rm -rf ~/.cache/pip

#
# Remove existing pypy virtual environment in the background, if it exists.
#
if [ -d $DIR/src/pypy_env ]; then
        echo 'Removing existing pypy virtualenv...'
        mv $DIR/src/pypy_env $DIR/src/pypy_env.tmp
        rm -rf $DIR/src/pypy_env.tmp &
fi

# First, the pypy virtualenv
echo 'Building pypy virtualenv...'
# export PYTHONPATH=''
virtualenv -p /net/gs/vol3/software/modules-sw-python/3.6.5/pypy/7.3.1/Linux/CentOS7/x86_64/bin/pypy3 $DIR/src/pypy_env

source $DIR/src/pypy_env/bin/activate

pypy3 -m ensurepip
pip install -r $DIR/pypy_requirements.txt

#
# Clone the two repositories once only.
#
git clone https://github.com/andrewhill157/barcodeutils.git

pushd barcodeutils
pypy setup.py install
popd

deactivate

