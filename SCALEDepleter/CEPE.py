import sys
import numpy as np
import subprocess
import time
import os

import getComps
import makeStdCmp
import runAndKillScale
import makeTritonFile
import copyMatAndF33Files
import makeAndRunOrigen
import removeAndMakeDir
import powerFromOutput
import pickledData

def normalizePower(bos_power: dict, eos_power: dict, previous_substep_power: dict,
                   denomVariableT0: dict, denomVariableT1: dict,
                   substepX0: dict, substepX1: dict):

  """
  bos_power (dict): beginning of macro step power
  eos_power (dict): end of macro step power
  previous_substep_power (dict): power used in previous substep
  denomVariableT0: usually xenon conc at beginning of macrostep
  denomVariableT1: usually xenon conc at beginning of macr0step
  """

  # get normalization factor to normalize power vector with
  bos_total = 0.0
  for key in bos_power.keys():
    bos_total += bos_power[key]

  # make dict for each material id of ratio deltaP/dX
  ratioDict = {}
  new_power_dict_unnormalized = {}
  for key in bos_power.keys():
    deltaV = denomVariableT1[key] - denomVariableT0[key]
    deltaP = eos_power[key] - bos_power[key]
    ratioDict[key] = deltaP/deltaV
    new_power_dict_unnormalized[key] = ratioDict[key] * (substepX1[key]-substepX0[key]) + previous_substep_power[key]

  # now renormalize
  new_total = 0.0
  for key in new_power_dict_unnormalized.keys():
    new_total += new_power_dict_unnormalized[key]

  # ...
  new_power_dict = {}
  for key in new_power_dict_unnormalized.keys():
    new_power_dict[key] = new_power_dict_unnormalized[key] / new_total * bos_total


  return new_power_dict

def print_substepping_arrays(xenon: dict, power: dict, step_num: int,
                             MC_power_eos: dict, MC_xenon_eos: dict):
  print('\n\n================================================================')
  print("Substepping now complete for Step Number", step_num)
  for substep_key in xenon[step_num]:
    print("================================================")
    print("Substep Number:", substep_key)
    print("MatID - Xenon BOSubstep | Power BOSubstep")
    print("Last Xenon entry is the EOS of the last substep - this no assoicated power vector.")
    for mat_id in xenon[step_num][substep_key]:
      try:
        p = power[step_num][substep_key][mat_id]
        x = xenon[step_num][substep_key][mat_id]
        print(x, '\t|\t', p)
      except:
        x = xenon[step_num][substep_key][mat_id]
        print(mat_id, '-', x, '\t\t|\t\t', 'N\A')
  last_xe_substep = max([x for x in xenon[step_num].keys()])
  print("================================================")
  print("Xenon Concentration by step num - BOSubstep concentrations")
  for mat_id in xenon[step_num][0]:
    line = str(mat_id)+' '
    for substep in xenon[step_num].keys():
      line += str(xenon[step_num][substep][mat_id])+' '
    print(line)
  print("================================================")
  print("Power by step num - BOS value - constant from BOSubstep to EOSubstep")
  for mat_id in power[step_num][0]:
    line = str(mat_id)+' '
    for substep in power[step_num].keys():
      line += str(power[step_num][substep][mat_id])+' '
    print(line)
  print("================================================")
  print("Predictor Xenon at T1 | Xenon from Corrector at T1 | Predictor Power from Transport | Corrected power in last substep")
  for mat_id in xenon[step_num][last_xe_substep]:
    print(mat_id, '-', MC_xenon_eos[mat_id], xenon[step_num][last_xe_substep][mat_id], MC_power_eos[mat_id], power[step_num][last_xe_substep-1][mat_id])
  print('\n\n================================================================')


def CEPE(fissionable_mats: list,
          fissionable_mats_vols: list,
          residual_number_density: float,
          include_non_fission_material_power: bool,
          print_transport_powers: bool,
          system_IHM_mass_grams: float,
          specific_power: list,
          steplength_days: list,
          origen_predictor_divs: int,
          addnuxdictbase: str,
          base_triton: str,
          origen_base: str,
          origenResults_F71dir: str,
          MonteCarloResults_F33dir: str,
          Nprocs: int,
          machinefile: str,
          tmpdir: str,
          is_parallel: bool,
          origen_LI_divs: int,
          origen_steps_per_div: int,
          print_substepping: bool):

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
  removeAndMakeDir.removeAndMakeDir(dirct='tmp_interpolated', make=False)
  for idx in range(num_steps):
    removeAndMakeDir.removeAndMakeDir(dirct='tmp_origen_PREDICTOR_step'+str(idx), make=False)
    removeAndMakeDir.removeAndMakeDir(dirct='tmp_origen_CORRECTOR_step'+str(idx), make=False)

  # we need to add initial nuclides to addnuxdict - addnux dict holds all the addnux nuclides + nuclides defined in the initial material defs
  addnuxdict, new_isotopes_from_addnux = getComps.makeNewAddnuxDict(time_lib.mats_by_steps[0]['C'].material_dict, tmpdir, addnuxdictbase)

  # any isotopes that are in the addnux dict but not the initial comp now needs to be added as residual isotopes
  mat_lib = time_lib.mats_by_steps[0]['C']
  for iso in new_isotopes_from_addnux:
    mat_lib.add_residual_isotope(isotope=iso, numdens=residual_number_density)

  # make a base triton scale input with flags for stdcmp - make for predicted and corrected stdcmps
  corrected_triton_input = makeTritonFile.makeTritonFile(base_triton, fissionable_mats, stdcmp_tag='corrected_cmp') # uses "corrcted" stdcmp material defs
  predicted_triton_input = makeTritonFile.makeTritonFile(base_triton, fissionable_mats, stdcmp_tag='predicted_cmp') # uses "predicted" stdcmp material defs

  xenon_start = {}

  # make a stdcmp file based on initial materials - automatically move stdcmps for triton_stdcmp to tmp dir
  for material_index in fissionable_mats:
    thisfilename = 'StdCmpMix'+str(material_index)+'_corrected_cmp'
    makeStdCmp.makeStdCmpFromMatLib(outputFilename=thisfilename, material_lib=initial_mats, material_index=material_index, tmpdir=tmpdir)

    # append xenon to current xenon vector
    xenon_start[material_index] = float(mat_lib.material_dict[material_index].return_iso_atom_dens(iso='xe-135'))

  # make scale run line based on serial or parallel
  if is_parallel:
    scale_run_line_c = ['scalerte', '-N', str(Nprocs), '-M',machinefile, corrected_triton_input, '-T', tmpdir, '-m'] # for cluster.
    scale_run_line_p = ['scalerte', '-N', str(Nprocs), '-M',machinefile, predicted_triton_input, '-T', tmpdir, '-m'] # for cluster.
  else:
    scale_run_line_c = ['scalerte', corrected_triton_input, '-T', tmpdir, '-m'] # for pc
    scale_run_line_p = ['scalerte', predicted_triton_input, '-T', tmpdir, '-m'] # for pc

  # initializing output data
  power_by_step = {}
  xenon_by_step = {} # xenon by step corresponds to power_by_step statepoints - it is the same xenon used in the MC calcs to generate power vector
  xenon_by_step[0] = xenon_start
  keff_lines = []

  # xenon and power tracking for substepping
  xenon_by_substeps = {} # xenon concentration for each substep
  power_by_substeps = {}

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
      f33_files_bos = f33_files # set bos f33 files

      # get and append power
      power_by_step[step_num] = powerFromOutput.getPower(printP=print_transport_powers, fission_mat_ids=fissionable_mats,
                                                        include_non_fission_material_power=include_non_fission_material_power,
                                                        filename=corrected_triton_input, total_power_python=specific_power_this_step)
    else:
      # if we are not doing very first step, use f33 files from T1 from previous step
      f33_files = f33_files_next
      f33_files_bos = f33_files # bos f33_files
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
                                                                          no_blended_name=True)
      origen_file_list.append(origen_file_handle)
      origen_tmpdirs.append(origen_tmpdir)


    # run origen now, unpack isotopics, then append isotopics to time_lib
    xenon_predicted = {}
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

      # update xenon_predicted
      pred_mat = getComps.get_comps_from_std_mix_file(tmpdir+'/'+stdCmpFilename) # gets comp from this material
      xenon_predicted[fiss_mat_id] = pred_mat.return_iso_atom_dens(iso='xe-135')

    # update xenon -
    xenon_by_step[step_num+1] = xenon_predicted

    ##############################################################
    ################## CORRECTOR - MONTE CARLO ###################
    ##############################################################

    # now, we have to run eigenvalue calc after depleting - get f33 at T1
    keff_line = runAndKillScale.runAndKillScale(scale_run_line_p, predicted_triton_input) # run scale for the 'predicted' triton input
    keff_lines.append(keff_line)

    # now copy over the new f33 files
    f33_files = copyMatAndF33Files.copy_files_from_temp(tmpdir=tmpdir, step_num=step_num+1, mcf33dir=MonteCarloResults_F33dir)
    f33_files_next = f33_files # set f33 files for the BOS for the next step....
    f33_files_eos = f33_files # eos f33 files

    # append power
    power_by_step[step_num+1] = powerFromOutput.getPower(printP=print_transport_powers, fission_mat_ids=fissionable_mats,
                                                      include_non_fission_material_power=include_non_fission_material_power,
                                                      filename=predicted_triton_input, total_power_python=specific_power_this_step)
    ##############################################################
    ################## CORRECTOR - ORIGEN ########################
    ##############################################################

    # we now have both f33 files. need to perform interpolations based on num substeps
    # origen_LI_divs # number of f33 interpolations
    # origen_steps_per_div # number of origen timesteps per interpolation
    dt = steplength_days_this_step / origen_LI_divs # how long averaging step is
    del_t = dt  / origen_steps_per_div  # how long each microstep is

    # get midpoints of the LI steps 0 is BOS and 5 is EOS - get 0.5,1.5,etc.
    # these are the times we interpolate f33 files to.
    #   0.5   1.5   2.5   3.5   4.5
    # |     |     |     |     |     |
    # 0     1     2     3     4     5
    LI_ends = np.linspace(0,steplength_days_this_step,origen_LI_divs+1)[1:]
    LI_starts = np.linspace(0,steplength_days_this_step,origen_LI_divs+1)[0:-1]
    LI_times = ( LI_ends + LI_starts ) /2.0
    start_time = 0.0
    end_time = start_time + steplength_days_this_step

    # first make f33 libraries that we use to interpolate
    f33_substep_filepath_dict = {} # keys are going to be fission material ids. values are lists of filepaths for f33 files for that mateiral
                                   # {101: [...substep0.f33, ....substep2.f33 ....]}
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      eos_file = f33_files_eos[fiss_mat_id]
      bos_file = f33_files_bos[fiss_mat_id]
      filepath = "tmp_interpolated/step"+str(step_num)+"/mat"+str(fiss_mat_id) # folder that holds the interpolated f33 files
      f33_substep_filepath_list = makeAndRunOrigen.f33Interpolate(filepath=filepath, bos_file=bos_file, eos_file=eos_file, times=LI_times, start_time=start_time, end_time=end_time)
      f33_substep_filepath_dict[fiss_mat_id] = f33_substep_filepath_list

    # 1. Make origen file - 1 big file containing depletion in that substep for ALL materials.
    # 2. Get power, XeStart, and XeEnd from previous step -> normalize power for this step
    # 3. Deplete with origen. Get Xe at end of each substep
    # 3. Update powers and xenons for next step. P0,X0 ->> P1,X1
    # 4. Go back to 1.
    xenon_by_substeps[step_num] = {}
    power_by_substeps[step_num] = {}

    for substep_idx, LI_time in enumerate(LI_times):
      # now recalculate power
      bos_power = power_by_step[step_num]
      eos_power = power_by_step[step_num+1]
      denomVariableT1 = xenon_by_step[step_num+1] #  dict w/ matid as keys
      denomVariableT0 = xenon_by_step[step_num] # dict w/ matid as keys
      if substep_idx == 0:
        # first substep so use P = P0
        substep_power = bos_power
        xStart = denomVariableT0
        if step_num == 0:
          xStart = denomVariableT0 # xstart should not be the MC-used xenon value - it should be the xenon value corresponding to the 'corrector' f71 file we did last step
        else:
          xStart = xEnd # sets to last value from previous step number that we got. this way xenon IC is equiv to f71 file in origen - otherwise negative number madness!
        xenon_by_substeps[0][substep_idx] = xStart # ste to xStart since this dict contains BOS information
      else:
        substep_power = normalizePower(bos_power=bos_power, eos_power=eos_power,
                                       previous_substep_power=previous_substep_power,
                                       denomVariableT0=denomVariableT0, denomVariableT1=denomVariableT1,
                                       substepX0=xStart, substepX1=xEnd)
      # save some data and for next step..
      power_by_substeps[step_num][substep_idx] = substep_power
      previous_substep_power = substep_power # set for next step




      bos_mat_lib = time_lib.material_lib_from_step(requested_step=step_num, PC_flag='C')
      substep_time_start = LI_starts[substep_idx]
      substep_time_end = LI_ends[substep_idx]
      # filehandle = input file name, tmpdir  = tmpdir for running origen, f71 names = list of names of output f71 files at end of
      origen_file_handle, origen_tmpdir, f71_names = makeAndRunOrigen.makeOrigenCEPEFile(fission_mat_ids=fissionable_mats, substep_power=substep_power,
                                                                                         substep_idx=substep_idx, origen_steps_per_div=origen_steps_per_div,
                                                                                         del_t=del_t, specific_power_this_step=specific_power_this_step,
                                                                                         volumes=fissionable_mats_vols, step_num=step_num,
                                                                                         f33_substep_filepath_dict=f33_substep_filepath_dict,
                                                                                         bos_cmp_lib=bos_mat_lib, origenResults_F71dir=origenResults_F71dir,
                                                                                         LI_starts=LI_starts)
      origen_output_loc = makeAndRunOrigen.runOrigenFile(tmpdir=origen_tmpdir, origen_file=origen_file_handle, material_id=-1, skipRunning=False)

      # now get Xenon at end of the origen run -> sets for next step
      if substep_idx > 0: # only reset xStart to preivous step xEnd if substep != 0
        xStart = xEnd # set xStart to previous value of xEnd for next step
      xEnd = {}
      for fiss_mat_id in fissionable_mats:
        # todo f71 files in tmp_origen_substepN_matN.f71 -> these have more than 1 pos - can keep it this way or
        xEnd[fiss_mat_id] = makeStdCmp.grabNuclideFromF71(filename=origenResults_F71dir+'/'+f71_names[fiss_mat_id],
                                                          nuclide='xe-135', precision=12)
      xenon_by_substeps[step_num][substep_idx+1] = xEnd # set to xEnd

    # update last index of xenon vector

    if print_substepping:
      print_substepping_arrays(xenon=xenon_by_substeps, power=power_by_substeps, step_num=step_num,
                               MC_power_eos=eos_power, MC_xenon_eos=denomVariableT1)


    # depletion now finished based on substepping approach - time to wrap things up since this is now EOS
    corrected_mat_lib = getComps.material_lib()
    for idx, fiss_mat_id in enumerate(fissionable_mats):
      # get stdcmp + tmeperature
      stdCmpFilename = 'StdCmpMix'+str(fiss_mat_id)+'_corrected_cmp'
      temperature = time_lib.material_lib_from_step(requested_step=0, PC_flag='C').material_dict[fiss_mat_id].temp # temperature from beginning of simulation for this matid
      # make a standard cmp from the f71
      origen_f71_loc = origenResults_F71dir+'/'+f71_names[fiss_mat_id]
      _ = makeStdCmp.makeStdCmpFromF71(materialNumber=fiss_mat_id, temperature=temperature, filename=origen_f71_loc, dictionaryFilename=addnuxdict, outputFilename=stdCmpFilename, tmpdir=tmpdir)
      corrected_mat_lib.append_mat_to_lib(getComps.get_comps_from_std_mix_file(tmpdir+'/'+stdCmpFilename))

    # append coirrected lib to time_lib
    time_lib.append_lib(time=step_end_time, lib=corrected_mat_lib, step=step_num+1, PC_flag='C')

    print('ok')


    # append coirrected lib to time_lib
    # time_lib.append_lib(time=step_end_time, lib=corrected_mat_lib, step=step_num+1, PC_flag='C')


  ##############################################################
  ###################### DATA OUTPUT ###########################
  ##############################################################


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
