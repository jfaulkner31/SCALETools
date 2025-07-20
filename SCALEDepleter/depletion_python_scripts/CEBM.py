import sys
import numpy as np
import subprocess
import time
import os

import getComps
import makeStdCmp
import SCALEDepleter.depletion_python_scripts.runAndKillScale as runAndKillScale
import makeTritonFile
import copyMatAndF33Files
import makeAndRunOrigen
import removeAndMakeDir
import powerFromOutput
import pickledData

"""
Runs a CEBM predictor/corrector scheme.
Step 0 - Run Transport at T0.
Step 1 - Deplete with T0 f33 file until T1 - get T1_predictor.f71.
Step 2 - Run Transport with T1.f71 nuclide composition - get T1.f33
Step 3 - Depletion from T0 -> T1 with T1.f33 - get T1_corrector.f71
Step 4 - Take average of T1_corrector.f71 and T1_predictor.f71 - get BLENDED.f71

Step 5 - Restart at Step 1 using T1.f33 as BOS f33 file and BLENDED.f71 as BOS composition
         ....
"""


def CEBM(fissionable_mats: list,
          fissionable_mats_vols: list,
          residual_number_density: float,
          include_non_fission_material_power: bool,
          print_transport_powers: bool,
          system_IHM_mass_grams: float,
          specific_power: list,
          steplength_days: list,
          origen_predictor_divs: int,
          origen_corrector_divs: int,
          addnuxdictbase: str,
          base_triton: str,
          origen_base: str,
          blender_base: str,
          origenResults_F71dir: str,
          MonteCarloResults_F33dir: str,
          Nprocs: int,
          machinefile: str,
          tmpdir: str,
          is_parallel: bool):

  ##################################################################
  ####################### SANITY CHECKS ############################
  ##################################################################
  num_steps = int
  if len(specific_power) != len(steplength_days):
    raise Exception("Specific power and steplength lengths do not match.")
  else:
    num_steps = len(specific_power)
  if len(fissionable_mats_vols) != len(fissionable_mats):
    raise Exception("Fissionable mat ids and fissionable mat volume lengths do not match.")
  if residual_number_density <= 1e-24:
    raise Exception("Residual number density cannot be less than 1e-24")


  ##################################################################
  ############################ SETUP ###############################
  ##################################################################

  # get temp dir for transport calculations
  tmpdir = tmpdir

  # get machine file
  machinefile = machinefile # defined from def

  # scale run handles
  print("NPROCS", Nprocs)
  print("MACHINEFILE",machinefile)
  print("BASE_TRITON",base_triton)

  # Get materials at time 0 and append
  time_lib = getComps.time_dependent_material_lib()
  fh = open(base_triton, 'r')
  initial_mats = getComps.get_comps(fh)
  fh.close()
  time_lib.append_lib(initial_mats, time=0.0, step=int(0), PC_flag='C') # C flag since it is the "true" solution

  # completely remove and remake some directories - since we are makign a new case
  removeAndMakeDir.removeAndMakeDir(dirct=tmpdir, make=True)
  removeAndMakeDir.removeAndMakeDir(dirct=MonteCarloResults_F33dir, make=True)
  removeAndMakeDir.removeAndMakeDir(dirct=origenResults_F71dir, make=True)
  removeAndMakeDir.removeAndMakeDir(dirct='tmp_blender', make=False)
  for idx in range(num_steps):
    removeAndMakeDir.removeAndMakeDir(dirct='tmp_origen_PREDICTOR_step'+str(idx), make=False)
    removeAndMakeDir.removeAndMakeDir(dirct='tmp_origen_CORRECTOR_step'+str(idx), make=False)

  # we need to add initial nuclides to addnuxdict - addnux dict holds all the addnux nuclides + nuclides defined in the initial material defs
  addnuxdict, new_isotopes_from_addnux = getComps.makeNewAddnuxDict(time_lib.mats_by_steps[0]['C'].material_dict, tmpdir, addnuxdictbase)

  # any isotopes that are in the addnux dict but not the initial comp now needs to be added as residual isotopes
  for iso in new_isotopes_from_addnux:
    mat_lib = time_lib.mats_by_steps[0]['C']
    mat_lib.add_residual_isotope(isotope=iso, numdens=residual_number_density)

  # make a base triton scale input with flags for stdcmp - make for predicted and corrected stdcmps
  corrected_triton_input = makeTritonFile.makeTritonFile(base_triton, fissionable_mats, stdcmp_tag='corrected_cmp') # uses "corrcted" stdcmp material defs
  predicted_triton_input = makeTritonFile.makeTritonFile(base_triton, fissionable_mats, stdcmp_tag='predicted_cmp') # uses "predicted" stdcmp material defs

  # make a stdcmp file based on initial materials - automatically move stdcmps for triton_stdcmp to tmp dir
  for material_index in fissionable_mats:
    thisfilename = 'StdCmpMix'+str(material_index)+'_corrected_cmp'
    makeStdCmp.makeStdCmpFromMatLib(outputFilename=thisfilename, material_lib=initial_mats, material_index=material_index, tmpdir=tmpdir)


  # make scale run line based on serial or parallel
  if is_parallel:
    scale_run_line_c = ['scalerte', '-N', str(Nprocs), '-M',machinefile, corrected_triton_input, '-T', tmpdir, '-m'] # for cluster.
    scale_run_line_p = ['scalerte', '-N', str(Nprocs), '-M',machinefile, predicted_triton_input, '-T', tmpdir, '-m'] # for cluster.
  else:
    scale_run_line_c = ['scalerte', corrected_triton_input, '-T', tmpdir, '-m'] # for pc
    scale_run_line_p = ['scalerte', predicted_triton_input, '-T', tmpdir, '-m'] # for pc

  # initializing output data
  power_by_step = {}
  keff_lines = []

  # start iteration
  for step_num in range(num_steps):
    steplength_days_this_step = steplength_days[step_num]
    specific_power_this_step = specific_power[step_num] * system_IHM_mass_grams / 1e6 # multiply specific power by the amount of tons in system
    step_start_time = sum(steplength_days[0:step_num])
    step_end_time = sum(steplength_days[0:step_num+1])

    ##############################################################
    ################## PREDICTOR - MONTE CARLO ###################
    ##############################################################

    # only run this on step 0 since we already have a transport solution for T0 if step_num > 0
    if step_num == 0:
      # STEP 1 - run transport with base file and then copy files from temp dir
      keff_line = runAndKillScale.runAndKillScale(scale_run_line_c, corrected_triton_input)
      keff_lines.append(keff_line)

      # make a folder with all f33's from the tmpdir moved into the new folder - label folder with step number - f33_files is a dict[material_id]
      f33_files = copyMatAndF33Files.copy_files_from_temp(tmpdir=tmpdir, step_num=step_num, mcf33dir=MonteCarloResults_F33dir)

      # get and append power
      power_by_step[step_num] = powerFromOutput.getPower(printP=print_transport_powers, fission_mat_ids=fissionable_mats,
                                                        include_non_fission_material_power=include_non_fission_material_power,
                                                        filename=corrected_triton_input, total_power_python=specific_power_this_step)
    else:
      # if we are not doing very first step, use f33 files from T1 from previous step
      f33_files = f33_files_next
    ##############################################################
    ################## PREDICTOR - ORIGEN ########################
    ##############################################################

    # Predictor step of origen - first write origen file for this step.
    origen_file_list = []
    origen_tmpdirs = []


    # make origen file, append origen file to list of origen files, append Beginning of step mat composition to the material library.
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      bos_mat_lib = time_lib.material_lib_from_step(requested_step=step_num, PC_flag='C')
      bos_mat_lib_this_mat = bos_mat_lib.material_dict[fiss_mat_id]
      origen_file_handle, origen_tmpdir = makeAndRunOrigen.makeOrigenFile(origen_base=origen_base, fiss_mat_id=fiss_mat_id, f33_files=f33_files,
                                                                          origenResults_F71dir=origenResults_F71dir, step_num=step_num,
                                                                          steplength_days=steplength_days_this_step, origen_predictor_divs=origen_predictor_divs,
                                                                          specific_power=specific_power_this_step*power_by_step[step_num][fiss_mat_id],
                                                                          predictor_corrector_string='PREDICTOR',
                                                                          bos_cmp=bos_mat_lib_this_mat, volume=fissionable_mats_vols[idx],
                                                                          no_blended_name=False)
      origen_file_list.append(origen_file_handle)
      origen_tmpdirs.append(origen_tmpdir)


    # run origen now, unpack isotopics, then append isotopics to time_lib
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      # get temp directory for this origen directory
      temp = origen_tmpdirs[idx]
      # get origen input
      origen_inp = origen_file_list[idx]
      # run origen
      origen_output_loc = makeAndRunOrigen.runOrigenFile(tmpdir=temp, origen_file=origen_inp, material_id=fiss_mat_id, skipRunning=False)
      origen_f71_loc = origenResults_F71dir+'/PREDICTOR_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'.f71'

      # now unpack isotopes from output using output location
      stdCmpFilename = 'StdCmpMix'+str(fiss_mat_id)+'_predicted_cmp'
      temperature = time_lib.material_lib_from_step(requested_step=step_num, PC_flag='C').material_dict[fiss_mat_id].temp # temperature from beginning of simulation for this matid

      # make a predicted stdcmp from the f71 file at the end of the predictor step
      _ = makeStdCmp.makeStdCmpFromF71(materialNumber=fiss_mat_id, temperature=temperature, filename=origen_f71_loc, dictionaryFilename=addnuxdict, outputFilename=stdCmpFilename, tmpdir=tmpdir)


    ##############################################################
    ################## CORRECTOR - MONTE CARLO ###################
    ##############################################################

    # now, we have to run eigenvalue calc after depleting - get f33 at T1
    keff_line = runAndKillScale.runAndKillScale(scale_run_line_p, predicted_triton_input) # run scale for the 'predicted' triton input
    keff_lines.append(keff_line)

    # now copy over the new f33 files
    f33_files = copyMatAndF33Files.copy_files_from_temp(tmpdir=tmpdir, step_num=step_num+1, mcf33dir=MonteCarloResults_F33dir)
    f33_files_next = f33_files # set f33 files for the BOS for the next step....

    # append power
    power_by_step[step_num+1] = powerFromOutput.getPower(printP=print_transport_powers, fission_mat_ids=fissionable_mats,
                                                      include_non_fission_material_power=include_non_fission_material_power,
                                                      filename=predicted_triton_input, total_power_python=specific_power_this_step)
    ##############################################################
    ################## CORRECTOR - ORIGEN ########################
    ##############################################################

    # now deplete using origen with the new f33 files from step_num+1

    origen_file_list = []
    origen_tmpdirs = []

    # make origen file, append origen file to list of origen files, append Beginning of step mat composition to the material library -
    # use f33 files from the corrector monte carlo - f33_files from above
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      bos_mat_lib = time_lib.material_lib_from_step(requested_step=step_num, PC_flag='C')
      bos_mat_lib_this_mat = bos_mat_lib.material_dict[fiss_mat_id]
      # todo - power_by_step[step_num] or step_num+1 for use in makeOrigenFile for the corrector? it was originalkly step_num+1 but i think i fixed it?
      origen_file_handle, origen_tmpdir = makeAndRunOrigen.makeOrigenFile(origen_base=origen_base, fiss_mat_id=fiss_mat_id, f33_files=f33_files,
                                                                          origenResults_F71dir=origenResults_F71dir, step_num=step_num,
                                                                          steplength_days=steplength_days_this_step, origen_predictor_divs=origen_corrector_divs,
                                                                          specific_power=specific_power_this_step*power_by_step[step_num][fiss_mat_id],
                                                                          predictor_corrector_string='CORRECTOR',
                                                                          bos_cmp=bos_mat_lib_this_mat, volume=fissionable_mats_vols[idx],
                                                                          no_blended_name=False)
      origen_file_list.append(origen_file_handle)
      origen_tmpdirs.append(origen_tmpdir)


    # run origen now, unpack isotopics, then append isotopics to time_lib
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      # get temp directory for this origen directory
      temp = origen_tmpdirs[idx]
      # get origen input
      origen_inp = origen_file_list[idx]
      # run origen
      origen_output_loc = makeAndRunOrigen.runOrigenFile(tmpdir=temp, origen_file=origen_inp, material_id=fiss_mat_id, skipRunning=False)
      origen_f71_loc = origenResults_F71dir+'/CORRECTOR_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'.f71'




    # now blend files together and calculate the average of the predictor and corrector nuclide vector using Origen - predictor and corrector blended - 0.5
    corrected_mat_lib = getComps.material_lib()
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      blender_tmp, new_input_name, path_to_f71 = makeAndRunOrigen.origenBlend(origen_f71_results_dir=origenResults_F71dir, step_num=step_num, mat_id=fiss_mat_id, blended_basefilename=blender_base)
      makeAndRunOrigen.runOrigenFile(tmpdir=blender_tmp, origen_file=new_input_name, material_id=-1, skipRunning=False)

      # make a standard cmp from the f71 - these are the corrected isotopes after averaging
      stdCmpFilename = 'StdCmpMix'+str(fiss_mat_id)+'_corrected_cmp'
      temperature = time_lib.material_lib_from_step(requested_step=0, PC_flag='C').material_dict[fiss_mat_id].temp # temperature from beginning of simulation for this matid
      origen_f71_loc = path_to_f71

      # make corrector stdcmp and place them in tmpdir
      _ = makeStdCmp.makeStdCmpFromF71(materialNumber=fiss_mat_id, temperature=temperature, filename=origen_f71_loc, dictionaryFilename=addnuxdict, outputFilename=stdCmpFilename, tmpdir=tmpdir)

      # get material
      corrected_mat_lib.append_mat_to_lib(getComps.get_comps_from_std_mix_file(tmpdir+'/'+stdCmpFilename))

    # append corrected mat lib to time lib
    time_lib.append_lib(time=step_end_time, lib=corrected_mat_lib, step=step_num+1, PC_flag='C')

  for line in keff_lines:
    print(line)


  # Save the data to a pickle file
  import pickle

  with open('keff.pkl', 'wb') as file:
    pickle.dump(keff_lines, file)
  with open('power.pkl', 'wb') as file:
    pickle.dump(power_by_step, file)
  with open('nuclides.pkl', 'wb') as file:
    pickle.dump(time_lib, file)

  pickledData.makeOutput(keff_lines=keff_lines, power_by_step=power_by_step,
                         specific_power=specific_power, steplength_days=steplength_days,
                         time_lib=time_lib, pkl_filename='output.pkl')
