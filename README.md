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

#### CD into your default directory
cd $HOME

#### Clone this repo:
git clone https://github.com/jfaulkner31/SCALETools.git

#### Create a new conda env
conda create --name scale
conda activate scale

#### Install SCALEDepleter
cd SCALETools
pip install .

#### If you want a modifiable version of SCALEDepleter such that changes reflect at runtime:
pip install -e .
