import numpy as np
import getComps
import subprocess
import time
import glob

def removePattern(pattern):
  files = glob.glob(pattern)
  for file in files:
    move = subprocess.run(['rm -fr '+file], shell=True) # enable shell to allow shell globbing

def runOrigenFile(origen_file, tmpdir, material_id, skipRunning):
  msg = origen_file[0:-4]+'.msg'
  out = origen_file[0:-4]+'.out'
  # run origen
  if not skipRunning:
    origenrun = subprocess.run(['scalerte', origen_file, '-T', tmpdir, '-m'])

  # move origen files for this into its temp directory and then delete the random stuff polluting main dir
  move1 = subprocess.run(['cp', origen_file, tmpdir+'/'+origen_file])
  move2 = subprocess.run(['cp', msg, tmpdir+'/'+msg])
  move3 = subprocess.run(['cp', out, tmpdir+'/'+out])

  # now run some subprocesses to remove files -
  # i could not get the * to work in the subprocess line so we first find files that have the pattern and then kill them file by fil
  pattern1 = origen_file[0:-4]+'*'
  files = glob.glob(pattern1)
  for file in files:
    move = subprocess.run(['rm -fr '+file], shell=True) # enable shell to allow shell globbing

  pattern1 = 'F33*'+str(material_id)+'*.f33'
  files = glob.glob(pattern1)
  for file in files:
    move = subprocess.run(['rm -fr '+file], shell=True) # enable shell to allow shell globbing

  pattern1 = 'STD_CMP_*'
  files = glob.glob(pattern1)
  for file in files:
    move = subprocess.run(['rm -fr '+file], shell=True) # enable shell to allow shell globbing

  move = subprocess.run(['rm -fr '+msg], shell=True) # remove message file specifically.

  origen_output_loc = tmpdir+'/'+out
  return origen_output_loc

def makeOrigenFile(origen_base, fiss_mat_id, f33_label, step_num, steplength_days, origen_predictor_divs, specific_power, predictor_corrector_string):

  # open base file
  with open(origen_base, 'r') as file:
    lines = file.readlines()
  file_handle = 'ORIGEN_'+predictor_corrector_string+'_step'+str(step_num)+'_mat'+str(fiss_mat_id)
  F33_FILE = 'F33_'+str(fiss_mat_id)+'_'+f33_label+'.f33'
  this_file = open(file_handle+'.inp', 'w', encoding="utf-8")

  # first make copy-paste shell section
  origen_tmpdir = 'tmp_origen_'+predictor_corrector_string+'_step'+str(step_num)
  line = '=shell'
  this_file.write(line+'\n')
  line = 'cp $INPDIR/'+F33_FILE+' $TMPDIR/'+F33_FILE
  this_file.write(line+'\n')
  line = 'rm $INPDIR/'+F33_FILE
  this_file.write(line+'\n')
  line = 'rm $INPDIR/'+'STD_CMP_'+str(fiss_mat_id)+'_'+predictor_corrector_string+'_step'+str(step_num)
  this_file.write(line+'\n')
  line = 'end'
  this_file.write(line+'\n\n\n')

  # now do main body of origen:
  for line in lines:
    if 'CASE_TITLE' in line:
      CASE_TITLE = file_handle
      line = '  title="' + CASE_TITLE + '"'
    if 'F33_FILE' in line:
      line = '    file="'+F33_FILE+'"'
    if 'F33_POS_INT' in line:
      F33_POS_INT=1
      line = '    pos='+str(F33_POS_INT)
    if 'TIME_VECTOR' in line:
      TIME_VECTOR = np.linspace(0, steplength_days[step_num], origen_predictor_divs)
      TIME_VECTOR = TIME_VECTOR[1::]
      line = '    t=' + str(TIME_VECTOR)
    if 'POWER VECTOR' in line:
      POWER_VECTOR = specific_power * np.ones(len(TIME_VECTOR))
      line = '  power=' + str(POWER_VECTOR) + ' %MW '
    if 'ISOTOPES_FROM_MATS' in line: # this one is only done for Step 0 -> gets comps from the triton stdcmp file. at tstep > 0, we get comps from previous origen EOS f71 file
      mats_file = 'STD_CMP_'+str(fiss_mat_id)+'_'+predictor_corrector_string+'_step'+str(step_num)
      mat = getComps.get_comps_from_std_mix_file(mats_file)
      line = mat.make_origen_materials()
    if "EOS_OUTPUT_F71_FILE" in line:
      EOS_OUTPUT_F71_FILE = predictor_corrector_string+'_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'.f71'
      line = '    file="'+EOS_OUTPUT_F71_FILE+'"'
    this_file.write(line+"\n")

  this_file.close()
  origen_input = file_handle+'.inp'
  return mat, origen_input, origen_tmpdir
