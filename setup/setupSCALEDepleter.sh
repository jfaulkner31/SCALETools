#!/bin/bash

# pulls latest version of SCALEDepleter as a folder in current repo.

# copy this folder to one youd like to setup a calculation/case in.
# then run this folder:
# bash setupSCALEDepleter.sh

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
cp -r SCALEDepleter/* .
rm -fr SCALETools
rm -fr SCALEDepleter

echo ' '
echo 'Now copying input file for SCALEDepleter as _INPUT_FILE_.py'
cp runCELI.py _INPUT_FILE_.py

echo ' '
echo 'All done! Now work through _INPUT_FILE_.py to begin your case'
echo ' '
