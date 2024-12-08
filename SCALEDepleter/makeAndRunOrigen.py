import numpy as np
import getComps
import subprocess
import time
import glob
from getComps import material_normal
from getComps import material_lib
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
                   bos_cmp: material_normal, volume: float,
                   no_blended_name: bool):

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

  # get blended_eos_step_mat.f71 file path - used for cases where step_num > 0 or if we are using initial conditions from followon calcs
  # if we are not doing the CEBM scheme, blended name is just CORRECTOR_.... since we arent doing blending unless we are doing CEBM
  if no_blended_name: # defaults to corrector_....
    blended_filename = 'CORRECTOR_EOS_step'+str(step_num-1)+'_mat'+str(fiss_mat_id)+'.f71'
    blended_filepath = blended_filepath = origenResults_F71dir+'/'+blended_filename
  else: # for the CEBM scheme where we use blending to average
    blended_filename = 'BLENDED_EOS_step'+str(step_num-1)+'_mat'+str(fiss_mat_id)+'.f71'
    blended_filepath = origenResults_F71dir+'/'+blended_filename


  line = '=shell'
  this_file.write(line+'\n')

  this_f33_name = 'mat'+str(fiss_mat_id)+'.f33'
  line = 'cp $INPDIR/../'+F33_FILE+' '+this_f33_name
  this_file.write(line+'\n')

  if step_num > 0:
    line = 'cp $INPDIR/../'+blended_filepath+' '+blended_filename
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
      line = '    file="'+this_f33_name+'"'
    if 'F33_POS_INT' in line:
      line = '    pos=1'
    if 'TIME_VECTOR' in line:
      TIME_VECTOR = np.linspace(0, steplength_days, origen_predictor_divs)
      TIME_VECTOR = TIME_VECTOR[1::]
      line = '    t=' + str(TIME_VECTOR)
    if 'POWER VECTOR' in line:
      POWER_VECTOR = [specific_power] * np.ones(len(TIME_VECTOR))
      line = '  power=' + str(POWER_VECTOR) + ' %MW '
    if 'ISOTOPES_FROM_MATS' in line:
      if step_num == 0:
        line = bos_cmp.make_origen_materials()
      else:
        line = '    load{ file="' + blended_filename + '" pos=1 }'
    if 'VOLUME_HERE' in line:
      if step_num == 0:
        line = '    volume=' + str(volume)+'\n'
      else:
        line = '\n' #  skip printing volumes if we are getting f71 material def
    if ('units=ATOMS-PER-BARN-CM' in line) & (step_num > 0):
      line = '\n' # do not print units= part of line -
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

def f33Interpolate(filepath: str, bos_file: str, eos_file: str, times: list, start_time: float, end_time: float):
  """
  f33_files_bos (dict): f33 files for each mat id at the beginning of this step.
  f33_files_eos (dict): f33 files for each mat id at the end of this step.
  times (list): list of floats
  start_time (float): the start time of the current step.
  end_time (float): the end time of the current step.
  """
  print("Now interpolating f33 files:", bos_file, '('+str(start_time)+') -> '+ eos_file +'('+str(end_time)+')')
  mkdir = os.makedirs(filepath)
  f33_substep_filepath_list = []
  for idx, time in enumerate(times):
    print("\tNow doing:", time)
    f33name = 'substep'+str(idx)+'.f33'
    f33path = filepath+'/'+f33name
    f33_substep_filepath_list.append(f33path)
    p = subprocess.run(['bash', 'interp2files.sh', bos_file, eos_file, str(start_time), str(end_time), str(time), f33path],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)
  return f33_substep_filepath_list

def makeOrigenCELIFile(fiss_mat_id: int,
                       step_num: int,
                       predictor_corrector_string: str,
                       f33_substep_filepath_list: list,
                       origenResults_F71dir: str,
                       LI_starts: list,
                       LI_ends: list,
                       dt: float,
                       del_t: float,
                       origen_steps_per_div: int,
                       specific_power: float,
                       volume: float,
                       bos_cmp: material_normal):
  """
  Makes origen files for CE/LI scheme corrector step.
  For the corrector step, we use linear interpolated f33 files for a set number of steps (dt long)
  Inside each of those steps are substeps that are del_t long.
  """
  # make temp directory for running origen for this step
  origen_tmpdir = 'tmp_origen_'+predictor_corrector_string+'_step'+str(step_num)
  if origen_tmpdir not in [entry.name for entry in os.scandir('.') if entry.is_file()]:
    mkdir = subprocess.run(['mkdir', origen_tmpdir])

  # filename of origen input file and making this_file
  file_handle = 'ORIGEN_'+predictor_corrector_string+'_step'+str(step_num)+'_mat'+str(fiss_mat_id)
  this_file = open(origen_tmpdir+'/'+file_handle+'.inp', 'w', encoding="utf-8")

  # name of starting composition f71 file
  f71filename_start = predictor_corrector_string+'_EOS_step'+str(step_num-1)+'_mat'+str(fiss_mat_id)+'.f71'
  f71filepath_start = '../'+origenResults_F71dir+'/'+f71filename_start

  # copy f33 files
  this_file.write('=shell\n')
  for idx, start_time in enumerate(LI_starts):
    f33name = 'mat'+str(fiss_mat_id)+'_substep'+str(idx)+'.f33'
    this_file.write('  cp '+'$INPDIR/../'+f33_substep_filepath_list[idx]+' '+f33name+'\n')
  this_file.write('end\n')

  # now write origen file
  this_file.write('=origen\nsolver{type=CRAM}\n')
  for idx, start_time in enumerate(LI_starts):
    if idx == len(LI_starts)-1:
      last_idx = True
    else:
      last_idx = False
    CASETITLE = 'mat'+str(fiss_mat_id)+'_step'+str(step_num)+'_substep'+str(idx)
    F33_FILEPATH = 'mat'+str(fiss_mat_id)+'_substep'+str(idx)+'.f33'

    TIME_VECTOR = [0]
    dt_vec = [del_t]*origen_steps_per_div # [del_t] -> [delt delt delt ...]
    for this in dt_vec:
      TIME_VECTOR.append(TIME_VECTOR[-1]+this)
    TIME_VECTOR = TIME_VECTOR[1:] #  delete the leading zero

    POWER_VECTOR = [specific_power]*origen_steps_per_div
    VOLUME = volume

    this_file.write("case{\n")
    this_file.write('  title="'+CASETITLE+'"\n')
    this_file.write('  lib{\n' )
    this_file.write('    file="'+F33_FILEPATH+'"\n')
    this_file.write('    pos=1\n')
    this_file.write('  }\n')
    this_file.write('  time{\n')
    this_file.write('    units=DAYS\n')
    this_file.write('    start=0\n')
    this_file.write('    t='+str(TIME_VECTOR)+'\n')
    this_file.write('  }\n')
    this_file.write('  power='+str(POWER_VECTOR)+' % MW\n')
    if (step_num == 0) & (idx == 0): # step_num == 0 we use from stdcmp library. idx != 0 we use from end of previous origen case.
      ISOTOPES_FROM_MATS = bos_cmp.make_origen_materials()
      this_file.write('  mat{\n')
      this_file.write('    units=ATOMS-PER-BARN-CM\n    volume='+str(VOLUME)+'\n'+ISOTOPES_FROM_MATS+'\n')
      this_file.write('  }\n')
    elif idx == 0: # get material def from previous EOS comp
      this_file.write('  mat{\n')
      this_file.write('    load{ file="'+f71filepath_start+'" pos=1 }\n')
      this_file.write('  }\n')
    else: # get material def from previous origen case
      this_file.write('  mat{\n')
      this_file.write('    previous=LAST\n')
      this_file.write('  }\n')
    if last_idx:
      EOS_OUTPUT_F71_FILE =  '../'+origenResults_F71dir+'/'+predictor_corrector_string+'_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'.f71'
      this_file.write('  save{ file="'+EOS_OUTPUT_F71_FILE+'"\n    steps=[LAST]\n  }\n')
    this_file.write('}\n')


  this_file.write('end\n')
  this_file.close()

  origen_input = file_handle
  return origen_input, origen_tmpdir

def makeOrigenCEPEFile(fission_mat_ids: list,
                       substep_power: dict,
                       substep_idx: int,
                       origen_steps_per_div: int,
                       del_t: float,
                       specific_power_this_step: float,
                       volumes: list,
                       step_num: int,
                       f33_substep_filepath_dict: dict,
                       bos_cmp_lib: material_lib,
                       origenResults_F71dir: str,
                       LI_starts: list):
  """
  Depletes N materials for a single substep
  """
  # make temp directory for running origen for this step
  origen_tmpdir = 'tmp_origen_CORRECTOR_step'+str(step_num)
  if origen_tmpdir not in [entry.name for entry in os.scandir('.') if entry.is_file()]:
    mkdir = subprocess.run(['mkdir', origen_tmpdir])

  # filename of origen input file and making this_file
  file_handle = 'ORIGEN_CORRECTOR_step'+str(step_num)+'_substep'+str(substep_idx)
  this_file = open(origen_tmpdir+'/'+file_handle+'.inp', 'w', encoding="utf-8")

  # get names of starting f71 files from beginning of substep or BOS as well as f33 files for this substep
  bos_f71_by_mat_id = {} # f71 files by mat_id
  eos_f71_by_mat_id = {}
  f33_file_dict = {} # f33 files by mat_id
  f71_filenames = {} # f71 output filenames dict where key = mateerial id
  for mat_id in fission_mat_ids:
    f33_file_dict[mat_id] = f33_substep_filepath_dict[mat_id][substep_idx]

    if substep_idx == 0:
      bos_f71_by_mat_id[mat_id] = 'CORRECTOR_EOS_step'+str(step_num-1)+'_mat'+str(mat_id)+'.f71'
    else:
      bos_f71_by_mat_id[mat_id] = 'CORRECTOR_substep'+str(substep_idx-1)+'_mat'+str(mat_id)+'.f71'

    if substep_idx == (len(LI_starts)-1): #  last step  - save eos file as special
      eos_f71_by_mat_id[mat_id] = 'CORRECTOR_EOS_step'+str(step_num)+'_mat'+str(mat_id)+'.f71'
    else:
      eos_f71_by_mat_id[mat_id] = 'CORRECTOR_substep'+str(substep_idx)+'_mat'+str(mat_id)+'.f71'
    f71_filenames[mat_id] = eos_f71_by_mat_id[mat_id]

  # cp f33 files
  this_file.write('=shell\n')
  for fiss_mat_id in fission_mat_ids:
    f33name = 'mat'+str(fiss_mat_id)+'_substep'+str(substep_idx)+'.f33'
    this_file.write('  cp '+'$INPDIR/../'+f33_file_dict[fiss_mat_id]+' '+f33name+'\n')
  this_file.write('end\n')

  # now write origen file
  this_file.write('=origen\nsolver{type=CRAM}\n')
  for idx, fiss_mat_id in enumerate(fission_mat_ids):
    CASETITLE = 'mat'+str(fiss_mat_id)+'_step'+str(step_num)+'_substep'+str(substep_idx)
    F33_FILEPATH = 'mat'+str(fiss_mat_id)+'_substep'+str(substep_idx)+'.f33'
    TIME_VECTOR = [0]
    dt_vec = [del_t]*origen_steps_per_div # [del_t] -> [delt delt delt ...]
    for this in dt_vec:
      TIME_VECTOR.append(TIME_VECTOR[-1]+this)
    TIME_VECTOR = TIME_VECTOR[1:] #  delete the leading zero
    specific_power = specific_power_this_step * substep_power[fiss_mat_id]
    POWER_VECTOR = [specific_power]*origen_steps_per_div
    VOLUME = volumes[idx]
    this_file.write("case{\n")
    this_file.write('  title="'+CASETITLE+'"\n')
    this_file.write('  lib{\n' )
    this_file.write('    file="'+F33_FILEPATH+'"\n')
    this_file.write('    pos=1\n')
    this_file.write('  }\n')
    this_file.write('  time{\n')
    this_file.write('    units=DAYS\n')
    this_file.write('    start=0\n')
    this_file.write('    t='+str(TIME_VECTOR)+'\n')
    this_file.write('  }\n')
    this_file.write('  power='+str(POWER_VECTOR)+' % MW\n')
    if (step_num == 0) & (substep_idx == 0): # step_num == 0 we use from stdcmp library. idx != 0 we use from end of previous origen case.
      bos_cmp = bos_cmp_lib.material_dict[fiss_mat_id]
      ISOTOPES_FROM_MATS = bos_cmp.make_origen_materials()
      this_file.write('  mat{\n')
      this_file.write('    units=ATOMS-PER-BARN-CM\n    volume='+str(VOLUME)+'\n'+ISOTOPES_FROM_MATS+'\n')
      this_file.write('  }\n')
    else: # get material def from previous EOS comp
      this_file.write('  mat{\n')
      # fix ??? possibly fix the use of step_num in the lines below for pos= but i think we are ok here....
      if substep_idx == 0:
        this_file.write('    load{ file="'+'../'+origenResults_F71dir+'/'+bos_f71_by_mat_id[fiss_mat_id]+'" pos=1 }\n')
      else:
        this_file.write('    load{ file="'+'../'+origenResults_F71dir+'/'+bos_f71_by_mat_id[fiss_mat_id]+'" pos='+str(step_num+1)+' }\n')
      this_file.write('  }\n')
    EOS_OUTPUT_F71_FILE = '../'+origenResults_F71dir+'/'+eos_f71_by_mat_id[fiss_mat_id]
    this_file.write('  save{ file="'+EOS_OUTPUT_F71_FILE+'"\n    steps=[LAST]\n  }\n')
    this_file.write('}\n')
  this_file.write('end\n')
  this_file.close()

  return file_handle, origen_tmpdir, f71_filenames


