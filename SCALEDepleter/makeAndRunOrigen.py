import numpy as np
import getComps
import subprocess
import time
import glob
from getComps import material_lib
from getComps import material_normal
import os

def removePattern(pattern):
  files = glob.glob(pattern)
  for file in files:
    move = subprocess.run(['rm -fr '+file], shell=True) # enable shell to allow shell globbing

def runOrigenFile(origen_file, tmpdir, material_id, skipRunning):
  msg = origen_file[0:-4]+'.msg'
  out = origen_file[0:-4]+'.out'
  # run origen
  if not skipRunning:
    origenrun = subprocess.run(['scalerte', tmpdir+'/'+origen_file, '-T', tmpdir, '-m'])

  origen_output_loc = tmpdir+'/'+out
  return origen_output_loc

def makeOrigenFile(origen_base: str, fiss_mat_id: int, f33_files: dict, origenResults_F71dir: str,
                   step_num: int, steplength_days: float, origen_predictor_divs: int,
                   specific_power: float, predictor_corrector_string: str,
                   bos_cmp: material_normal):

  # first make directory for running origen in
  origen_tmpdir = 'tmp_origen_'+predictor_corrector_string+'_step'+str(step_num)
  if origen_tmpdir not in [entry.name for entry in os.scandir('.') if entry.is_file()]:
    mkdir = subprocess.run(['mkdir', origen_tmpdir])

  # open base file
  with open(origen_base, 'r') as file:
    lines = file.readlines()
  file_handle = 'ORIGEN_'+predictor_corrector_string+'_step'+str(step_num)+'_mat'+str(fiss_mat_id)
  F33_FILE = f33_files[fiss_mat_id]
  this_file = open(origen_tmpdir+'/'+file_handle+'.inp', 'w', encoding="utf-8")


  line = '=shell'
  this_file.write(line+'\n')

  line = 'cp $INPDIR/../'+F33_FILE+' this_f33.f33'
  this_file.write(line+'\n')

  # line = 'rm $INPDIR/'+F33_FILE
  # this_file.write(line+'\n')

  # line = 'rm $INPDIR/'+'STD_CMP_'+str(fiss_mat_id)+'_'+predictor_corrector_string+'_step'+str(step_num)
  # this_file.write(line+'\n')

  line = 'end'
  this_file.write(line+'\n\n\n')

  # now do main body of origen:
  for line in lines:
    if 'CASE_TITLE' in line:
      CASE_TITLE = file_handle
      line = '  title="' + CASE_TITLE + '"'
    if 'F33_FILE' in line:
      line = '    file="this_f33.f33"'
    if 'F33_POS_INT' in line:
      F33_POS_INT=1
      line = '    pos='+str(F33_POS_INT)
    if 'TIME_VECTOR' in line:
      TIME_VECTOR = np.linspace(0, steplength_days, origen_predictor_divs)
      TIME_VECTOR = TIME_VECTOR[1::]
      line = '    t=' + str(TIME_VECTOR)
    if 'POWER VECTOR' in line:
      POWER_VECTOR = [specific_power] * np.ones(len(TIME_VECTOR))
      line = '  power=' + str(POWER_VECTOR) + ' %MW '
    if 'ISOTOPES_FROM_MATS' in line:
      line = bos_cmp.make_origen_materials()
    if "EOS_OUTPUT_F71_FILE" in line:
      EOS_OUTPUT_F71_FILE =  '../'+origenResults_F71dir+'/'+predictor_corrector_string+'_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'.f71'
      line = '    file="'+EOS_OUTPUT_F71_FILE+'"'
    this_file.write(line)

  this_file.close()
  origen_input = file_handle+'.inp'
  return origen_input, origen_tmpdir

def origenBlend(origen_f71_results_dir: str, step_num: int, mat_id: int, blended_basefilename: str):

  new_input_name = blended_basefilename.split('.')[0] + '_step'+str(step_num)+'_mat'+str(mat_id)+'.inp'

  blender_dir = 'tmp_blender'
  input_path = blender_dir+'/'+new_input_name
  if blender_dir not in [entry.name for entry in os.scandir('.') if entry.is_file()]:
    mkdir = subprocess.run(['mkdir', blender_dir])
  cp = subprocess.run(['cp', blended_basefilename, blender_dir+'/'+new_input_name])

  corrector_f71_string = '$INPDIR/../'+origen_f71_results_dir+'/CORRECTOR_EOS_step'+str(step_num)+'_mat'+str(mat_id)+'.f71'
  predictor_f71_string = '$INPDIR/../'+origen_f71_results_dir+'/PREDICTOR_EOS_step'+str(step_num)+'_mat'+str(mat_id)+'.f71'
  blended_f71_string ='$INPDIR/../'+origen_f71_results_dir+'/BLENDED_EOS_step'+str(step_num)+'_mat'+str(mat_id)+'.f71'

  # '../../OrigenResults_F71dir/CORRECTOR_EOS_step0_mat101.f71'
  # now make the file
  with open(input_path, 'r') as file:
    lines = file.readlines()

  for idx, line in enumerate(lines):
    if 'PREDICTOR_STRING_HERE' in line:
      line = '  cp '+predictor_f71_string+' p.f71\n'
    if 'CORRECTOR_STRING_HERE' in line:
      line = '  cp '+corrector_f71_string+' c.f71\n'
    if 'BLENDER_SAVE' in line:
      line = 'save{ file="step'+str(step_num)+'_mat'+str(mat_id)+'.f71'+'"    steps=[LAST] }\n'
    if 'BLENDED_STRING_HERE' in line:
      line = '  cp $INPDIR/step'+str(step_num)+'_mat'+str(mat_id)+'.f71'+' '+blended_f71_string+'\n'
    lines[idx] = line

  with open(input_path, 'w') as f:
    for line in lines:
        f.write(f"{line}")

  # returns 1 - temp for blender 2 - input file name inside of that temp, 3 - path to f71
  return blender_dir, new_input_name, origen_f71_results_dir+'/BLENDED_EOS_step'+str(step_num)+'_mat'+str(mat_id)+'.f71'



