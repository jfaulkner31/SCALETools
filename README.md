# SCALETools
Tools for SCALE

To setup SCALEDepleter and automatically download the SCALEDepleter into a folder:
First cd into some project directory:

cd project

Then make a directory for your calculation:

mkdir myCalculation

Then run the setup script:

bash setupSCALEDepleter.sh

Then if needed, add SCALEDepleter to your python path:

export PYTHONPATH=<project>/<myCalculation>/SCALEDepleter:$PYTHONPATH

# To install on a system

### CD into your default directory
cd $HOME

### Clone this repo:
git clone https://github.com/jfaulkner31/SCALETools.git

### Create a new conda env
conda create --name scale

conda activate scale

### Install pip for conda
conda install pip

### Install SCALEDepleter
cd SCALETools

pip install -e .

### If you want a modifiable version of SCALEDepleter such that changes reflect at runtime:
pip install -e .

### ... and if you want to not modify changes to the module:
pip install .

### Test that it works:
cd ../

python

import SCALEDepleter

### Exit python with ctrl+D
ctrl+D

### To copy the input file for CELI into a base directory:
cd "yourWorkingProjectFolder"

python

import SCALEDepleter.depletion_python_scripts.CELI as CELI

CELI.makeInput()


