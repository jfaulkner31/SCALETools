#!/bin/bash

# pulls latest version of SCALEDepleter as a folder in current repo.

# copy this folder to one youd like to setup a calculation/case in.
# then run this folder:
# bash setupSCALEDepleter.sh
# may need to also export SCALEDepleter into the python path when you are done.

# WARNING - deleters entire SCALEDepleter folder - only use if absolutely necessary or if starting out for the first time.
rm -fr SCALETools
rm -fr SCALEDepleter
rm -fr celiCalculationInput.py

git clone https://github.com/jfaulkner31/SCALETools.git
cd SCALETools
git sparse-checkout init --cone
git sparse-checkout set SCALEDepleter
git pull origin main
cd ../
cp -r SCALETools/SCALEDepleter SCALEDepleter
rm -fr SCALETools

echo ' '
echo 'Now setting up python path...'
echo ' '
echo 'You may need to add the following line manually in the bashrc, or rerun the export line for every terminal you open'
echo 'export PYTHONPATH=~/projects/semiImplicitBurnup/SCALEDepleter:$PYTHONPATH'
export PYTHONPATH=~/projects/semiImplicitBurnup/SCALEDepleter:$PYTHONPATH

echo ' '
echo 'Now copying input file for SCALEDepleter as celiCalculationInput.py'
cp SCALEDepleter/runCELI.py celiCalculationInput.py

echo ' '
echo 'All done!'
echo ' '
