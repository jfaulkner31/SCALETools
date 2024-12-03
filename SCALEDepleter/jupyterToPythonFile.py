import sys
import numpy as np
import subprocess
import time
import os

import getComps
import runAndKillScale
import copyMatAndF33Files
import makeAndRunOrigen
# some notes:
# FMA=fission matrix
#
#
#
##################################################################
########################### USER INPUT ###########################
##################################################################
"""
  Step 0 - ask yourself:
  Is first step length ok in the input file?
  Is the FMA parameter set to yes?
  Are the NPG and other neutron parameters correct?
"""

# fissionable regions - used for origen later
fissionable_mats = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116]
fissionable_mats_vols = [6.02834870915574000000E+04]*16

# add addnux dict location
addnuxdict = 'addnuxDicts/addnux0Dict.dict'

# TODO add volumes, get fission source shape after monte carlo from MC output, use fission source shape to adjust poewrs in each node.

# ORIGEN information
specific_power = [200]
steplength_days = [5]
origen_predictor_divs = 100 # number of time divs for predictor
origen_corrector_divs = 100 # number of time divs for corrector

# file handles
base_triton = 'triton_base.inp'
origen_base = 'baseOrigenFile.inp'

# nprocs
Nprocs = 48


##################################################################
############################ SETUP ###############################
##################################################################

# get temp dir for transport calculations
tmpdir = sys.argv[2]

# get machine file
machinefile = sys.argv[1]

# scale run handles
print("NPROCS", Nprocs)
print("MACHINEFILE",machinefile)
print("BASE_TRITON",base_triton)
scale_run_line = ['scalerte', '-N', str(Nprocs), '-M',machinefile, base_triton, '-T', tmpdir, '-m'] # for cluster.
# scale_run_line = ['scalerte', base_triton, '-T', tmpdir, '-m'] # for pc

# get materials
time_lib = getComps.time_dependent_material_lib()

# start iteration
step_num = 0 # temporary thing

### FOR LOOP STARTS HERE

# set origen temp dir

################################################
################## PREDICTOR ###################
################################################

# STEP 1 - run transport with base file and then copy files from temp dir
runAndKillScale.runAndKillScale(scale_run_line, base_triton)
std_cmp_label = 'PREDICTOR_step'+str(step_num)
f33_label = 'PREDICTOR_step'+str(step_num)
f33_files, std_cmp_files = copyMatAndF33Files.copy_files_from_temp(tmpdir=tmpdir, fissionable_mats=fissionable_mats, std_cmp_label=std_cmp_label, f33_label=f33_label)
# deleteTempDir(tmpdir=tmpdir) # delete temp dir

# Predictor step of origen - first write origen file for this step.
origen_file_list = []
origen_tmpdirs = []
material_lib_at_this_step = getComps.material_lib()
for fiss_mat_id in fissionable_mats:
  # make origen file, append origen file to list of origen files, append Beginning of step mat composition to the material library.
  mat, file_handle, origen_tmpdir = makeAndRunOrigen.makeOrigenFile(origen_base, fiss_mat_id, f33_label, step_num, steplength_days, origen_predictor_divs, specific_power, predictor_corrector_string='PREDICTOR')
  origen_file_list.append(file_handle)

  if step_num == 0:
    # on step 0, we append the input materials to origen - need a baseline nuclide vector to start depletion.
    # for future steps, use predictor corroector to determine what nuclide vector should be
    material_lib_at_this_step.append_mat_to_lib(mat)

  origen_tmpdirs.append(origen_tmpdir)

# append material lib to time lib once its fully made.
time_lib.append_lib(material_lib_at_this_step, time=steplength_days[step_num], step=step_num, PC_flag='C') # labelled C since the BOS comp will come from end of corrector step.

# if step num == 0 - we need to add initial nuclides to addnuxdict
if step_num == 0:
  addnuxdict = getComps.makeNewAddnuxDict(zeromatdict=time_lib.mats_by_steps[0].material_dict, tmpdir=tmpdir, addnuxdict=addnuxdict)


# run origen now, unpack isotopics, then append isotopics to time_lib
for idx, fiss_mat_id in enumerate(fissionable_mats):
  # get temp directory for this origen directory
  temp = origen_tmpdirs[idx]
  # get origen input
  origen_inp = origen_file_list[idx]
  # run origen
  origen_output_loc = makeAndRunOrigen.runOrigenFile(tmpdir=temp, origen_file=origen_inp, material_id=fiss_mat_id)
  # now unpack isotopes from output using ouptput location


# remove some random leftover files from origen:
makeAndRunOrigen.removePattern('F33*.f33')
makeAndRunOrigen.removePattern('ORIGEN_*.f33')
makeAndRunOrigen.removePattern('*.f71')
makeAndRunOrigen.removePattern('ORIGEN_*.msg')
