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
import shutil
import results_data

"""
Runs a CELI predictor/corrector scheme.
Step 0 - Run Transport at T0.
Step 1 - Deplete with T0 f33 file until T1 - get T1_predictor.f71.
Step 2 - Run Transport with T1.f71 nuclide composition - get T1.f33
Step 3 - Depletion from T0 -> T1 with T1.f33 - get T1_corrector.f71
Step 4 - Take average of T1_corrector.f71 and T1_predictor.f71 - get BLENDED.f71

Step 5 - Restart at Step 1 using T1.f33 as BOS f33 file and BLENDED.f71 as BOS composition
         ....
"""


def CELI(fissionable_mats: list,
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
          corrector_iterations: int,
          relaxation_factor: float,
          case_name: str,
          include_predictor_in_blender: bool):

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
  if corrector_iterations <= 0:
    raise Exception("Corrector iterations must be greater than 1")
  if (relaxation_factor > 1.0) | (relaxation_factor <=0.0):
    raise Exception("Relaxation factor is bounded (0,1]")
  if (origen_predictor_divs < 2):
    raise Exception("Origen predictor divs must be 2 or more - 2 reccomended since CRAM methodology is being used and CRAM accuracy is not strongly related to timestep size.")

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
  initial_mats = getComps.get_comps(fh) # returns a matlib
  fh.close()
  time_lib.append_lib(initial_mats, time=0.0, step=int(0), PC_flag='C') # C flag since it is the "true" solution

  # make a results lib specifically made for CELI results
  celi_results_lib = results_data.ResultsCELI(
              fissionable_mats=fissionable_mats,
              fissionable_mats_vols=fissionable_mats_vols,
              residual_number_density=residual_number_density,
              include_non_fission_material_power=include_non_fission_material_power,
              print_transport_powers=print_transport_powers,
              system_IHM_mass_grams=system_IHM_mass_grams,
              specific_power=specific_power,
              steplength_days=steplength_days,
              origen_predictor_divs=origen_predictor_divs,
              addnuxdictbase=addnuxdictbase,
              base_triton=base_triton,
              origen_base=origen_base,
              origenResults_F71dir=origenResults_F71dir,
              MonteCarloResults_F33dir=MonteCarloResults_F33dir,
              Nprocs=Nprocs,
              machinefile=machinefile,
              tmpdir=tmpdir,
              is_parallel=is_parallel,
              origen_LI_divs=origen_LI_divs,
              origen_steps_per_div=origen_steps_per_div,
              corrector_iterations=corrector_iterations,
              relaxation_factor=relaxation_factor,
              case_name=case_name,
              include_predictor_in_blender=include_predictor_in_blender)


  # completely remove and remake some directories - since we are makign a new case
  removeAndMakeDir.removeAndMakeDir(dirct=tmpdir, make=True)
  removeAndMakeDir.removeAndMakeDir(dirct=MonteCarloResults_F33dir, make=True)
  removeAndMakeDir.removeAndMakeDir(dirct=origenResults_F71dir, make=True)
  removeAndMakeDir.removeAndMakeDir(dirct='tmp_blender', make=False)
  removeAndMakeDir.removeAndMakeDir(dirct='tmp_interpolated', make=False)
  removeAndMakeDir.removeAndMakeDir(dirct='blended_CELI_tmp', make=False)
  removeAndMakeDir.removeAndMakeDir(dirct='CASE_'+case_name, make=True)

  for idx in range(num_steps):
    removeAndMakeDir.removeAndMakeDir(dirct='tmp_origen_PREDICTOR_step'+str(idx), make=False)
    for ci in range(corrector_iterations):
      removeAndMakeDir.removeAndMakeDir(dirct='tmp_origen_CORRECTOR_step'+str(idx)+'_corrIter'+str(ci), make=False)



  # we need to add initial nuclides to addnuxdict - addnux dict holds all the addnux nuclides + nuclides defined in the initial material defs
  addnuxdict, new_isotopes_from_addnux = getComps.makeNewAddnuxDict(time_lib.mats_by_steps[0]['C'].material_dict, tmpdir, addnuxdictbase)

  # any isotopes that are in the addnux dict but not the initial comp now needs to be added as residual isotopes
  for iso in new_isotopes_from_addnux:
    mat_lib = time_lib.mats_by_steps[0]['C']
    mat_lib.add_residual_isotope(isotope=iso, numdens=residual_number_density)

  # append IC results for CELI datakeeping
  celi_results_lib.append_isotopics_data(step_num=0, iteration_num=-1, matlib=mat_lib)

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

  # start iteration for each burnup step
  for step_num in range(num_steps):
    steplength_days_this_step = steplength_days[step_num]
    specific_power_this_step = specific_power[step_num] * system_IHM_mass_grams / 1e6 # multiply specific power by the amount of tons in system
    step_start_time = sum(steplength_days[0:step_num])
    step_end_time = sum(steplength_days[0:step_num+1])

    ##############################################################
    ################## PREDICTOR - MONTE CARLO ###################
    ##############################################################

    # STEP 1 - run transport with base file and then copy files from temp dir
    # this part runs MC at BOS - the correct MC solution for the beginning of this step
    # since it comes either from IC nuclide density or a nuclide density converged
    # at the corrector at the previous steps EOS convergence procedure
    keff_line = runAndKillScale.runAndKillScale(scale_run_line_c, corrected_triton_input)
    keff_lines.append(keff_line)

    # append results to celi mat lib
    celi_results_lib.append_keff_data(step_num=step_num, iteration_num=-1, keff_line=keff_line)


    # make a folder with all f33's from the tmpdir moved into the new folder - label folder with step number - f33_files is a dict[material_id]
    f33_files = copyMatAndF33Files.copy_files_from_temp(tmpdir=tmpdir, step_num=step_num, mcf33dir=MonteCarloResults_F33dir, appendThis='_pred')
    f33_files_bos = f33_files # set bos f33 files

    # get and append power - this overwrites values from the corrector iteration steps - correct power at T0 for this step
    # power_by_step has fractional powers. (all powers in all materials will sum to 1.0,
    # if not including nonfissionable materials, they will sum to slightly less than 1.0)

    _, _, mass_HM_this_step = powerFromOutput.get_triton_depletion_table(filename=corrected_triton_input)

    power_by_step[step_num] = powerFromOutput.getPower(printP=print_transport_powers, fission_mat_ids=fissionable_mats,
                                                      include_non_fission_material_power=include_non_fission_material_power,
                                                      filename=corrected_triton_input, total_power_python=specific_power_this_step)
    # for celi results data storage
    celi_results_lib.append_power_data(step_num=step_num, iteration_num=-1, power_data_block=power_by_step[step_num])

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

    origen_f71_locs_all = {} # make empty dict -> origen_f71_locs_all[iteration_num][fiss_mat_id]
    for ci in range(corrector_iterations):
      mat_lib_this_step = getComps.material_lib()

      # now, we have to run eigenvalue calc after depleting - get f33 at T1
      keff_line = runAndKillScale.runAndKillScale(scale_run_line_p, predicted_triton_input) # run scale for the 'predicted' triton input
      keff_lines.append(keff_line)

      # append celi results for datakeeping
      celi_results_lib.append_keff_data(step_num=step_num, iteration_num=ci, keff_line=keff_line)

      # now copy over the new f33 files
      f33_files = copyMatAndF33Files.copy_files_from_temp(tmpdir=tmpdir, step_num=step_num+1, mcf33dir=MonteCarloResults_F33dir, appendThis='_corrIter'+str(ci))
      f33_files_eos = f33_files # eos f33 files

      # append power - progressively updates power as we iterate - OK - does not record power during Stochastic iterations since it overwrites.
      power_by_step[step_num+1] = powerFromOutput.getPower(printP=print_transport_powers, fission_mat_ids=fissionable_mats,
                                                        include_non_fission_material_power=include_non_fission_material_power,
                                                        filename=predicted_triton_input, total_power_python=specific_power_this_step)

      celi_results_lib.append_power_data(step_num=step_num, iteration_num=ci, power_data_block=power_by_step[step_num+1])

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
      origen_file_list = []
      origen_tmpdirs = []

      # interpolate power beteen T0 and T1 @ starts of substeps (LI_starts): interpd_power_dict[mat_id][idx of timestep]
      power_t0 = power_by_step[step_num] # power at t0
      power_t1 = power_by_step[step_num+1] # power at t1
      interpd_power_dict = powerFromOutput.interpolatePower(power_by_step_t0_pre=power_t0, power_by_step_t1_pre=power_t1,
                                        times=LI_starts, start_time=0.0, end_time=steplength_days_this_step,
                                        specific_power_this_step=specific_power_this_step,
                                        printP=False)


      # now doing origen
      for idx, fiss_mat_id in enumerate(fissionable_mats):
        # first make interpolated f33 files.
        eos_file = f33_files_eos[fiss_mat_id]
        bos_file = f33_files_bos[fiss_mat_id]
        filepath = "tmp_interpolated/step"+str(step_num)+"/mat"+str(fiss_mat_id) # folder that holds the interpolated f33 files
        f33_substep_filepath_list = makeAndRunOrigen.f33Interpolate(filepath=filepath, bos_file=bos_file, eos_file=eos_file, times=LI_times, start_time=start_time, end_time=end_time, appendThis='_corrIter'+str(ci))
        bos_mat_lib = time_lib.material_lib_from_step(requested_step=step_num, PC_flag='C')
        bos_mat_lib_this_mat = bos_mat_lib.material_dict[fiss_mat_id]
        # interpolated files in "tmp_interpolated/step"+str(step_num)+"/mat"+str(fiss_mat_id)+'/'+substep+str(sumstepInt)+'.f33' ---> tmp_interpolated/step0/mat101/substep0.f33

        # now make origen files that we are going to run
        # todo - power_by_step[step_num] or step_num+1 for use in makeOrigenFile for the corrector? it was originalkly step_num+1 but i think i fixed it?
        # todo --- resolved ---> what i need to do is interpolate the power using the CELI scheme.
        # todo ---- see Flux renormalization in constant power burnup calculations by Isotalo for more information on the above.
        # todo ---- what we want to basically do is add power interpolation between T0 and T1,
        # todo ---- good sens. analysis can also be to test if using all T0, all T1, or a 0.7T0 + 0.3T1 might be better? reasonable question to ask what power to use when depleting.
        # TODO:: |V
        # my predictor corrector scheme is free to sya the power is anything - whether it is T0, T1, or T1/2 or interp(T) power or a weighted crank nicolson style power.
        # just note that the fundamental way origen does normalization during each origen substep is very unique and should be noted and considered in any paper.

        origen_file_handle, origen_tmpdir = makeAndRunOrigen.makeOrigenCELIFile(fiss_mat_id=fiss_mat_id, step_num=step_num, predictor_corrector_string='CORRECTOR',
                                                                                f33_substep_filepath_list=f33_substep_filepath_list, origenResults_F71dir=origenResults_F71dir,
                                                                                LI_starts=LI_starts, LI_ends=LI_ends, dt=dt, del_t=del_t, origen_steps_per_div=origen_steps_per_div,
                                                                                specific_power=specific_power_this_step*power_by_step[step_num][fiss_mat_id],
                                                                                volume=fissionable_mats_vols[idx], bos_cmp=bos_mat_lib_this_mat,
                                                                                appendThis='_corrIter'+str(ci),
                                                                                interpd_power_dict=interpd_power_dict)
        origen_file_list.append(origen_file_handle)
        origen_tmpdirs.append(origen_tmpdir)


      # now run origen, unpack isotopics, and then append isotopics to the time_lib
      corrected_mat_lib = getComps.material_lib()

      origen_f71_locs_all[ci] = {} # make empty dict for storing the location of f71 files for this iteration

      for idx, fiss_mat_id in enumerate(fissionable_mats):
        # get temp directory for this origen directory
        temp = origen_tmpdirs[idx]
        # get origen input
        origen_inp = origen_file_list[idx]
        # run origen
        origen_output_loc = makeAndRunOrigen.runOrigenFile(tmpdir=temp, origen_file=origen_inp, material_id=fiss_mat_id, skipRunning=False)
        origen_f71_loc = origenResults_F71dir+'/CORRECTOR_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'_corrIter'+str(ci)+'.f71' # OrigenResults_F71dir/corrector_eos_step0_mat101_corrIter2.f71
        origen_f71_locs_all[ci][fiss_mat_id] = origen_f71_loc

        ##### right here we have to update N for next MC iteration based on all previous iterations - use origen blend function and origen_f71_loc

        # f71 name: '../'+origenResults_F71dir+'/BLENDED'+'_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+appendThis+'.f71'
        blended_f71_filename = makeAndRunOrigen.blendCELIOrigenFiles(corrector_iteration=ci, step_num=step_num, fiss_mat_id=fiss_mat_id,
                                                                     origenResults_F71dir=origenResults_F71dir, predictor_corrector_string='CORRECTOR',
                                                                     appendThis='_corrIter'+str(ci),
                                                                     relaxation_factor=relaxation_factor,
                                                                     origen_f71_locs_all=origen_f71_locs_all,
                                                                     include_predictor_in_blender=include_predictor_in_blender)

        # get std cmp filename for next MC run
        if ci == corrector_iterations - 1: # (last iteration no we have need to setup next corrector iteration)
          stdCmpFilename = 'StdCmpMix'+str(fiss_mat_id)+'_corrected_cmp' # use a corrected comp for predictor iteration
        else:
          stdCmpFilename = 'StdCmpMix'+str(fiss_mat_id)+'_predicted_cmp' # use predicted comp for corrector iteration

        # get temperature for this material id
        temperature = time_lib.material_lib_from_step(requested_step=0, PC_flag='C').material_dict[fiss_mat_id].temp # temperature from beginning of simulation for this matid

        # make the standard cmp from the f71 - whether next is P or C we are running a calc at T1.
        _ = makeStdCmp.makeStdCmpFromF71(materialNumber=fiss_mat_id, temperature=temperature, filename=blended_f71_filename, dictionaryFilename=addnuxdict, outputFilename=stdCmpFilename, tmpdir=tmpdir)

        # append the blended material lib for this iteration and this maetrial
        mat_lib_this_step.append_mat_to_lib(getComps.get_comps_from_std_mix_file(tmpdir+'/'+stdCmpFilename))

        # add to mat lib if this is last iteration
        # CHANGED HERE !!! Original -> ci == corrector_iterations - 1 --> if True
        # changed so that this always runs and all substep isotopic libraries are always kept track of and updated.
        if True:
          corrected_mat_lib.append_mat_to_lib(getComps.get_comps_from_std_mix_file(tmpdir+'/'+stdCmpFilename))

        # setup f71 to use as IC in BOS predictor calculation for next step since we are done iterating and supposedly converged
        if ci == corrector_iterations - 1:
          shutil.copyfile(blended_f71_filename, origenResults_F71dir+'/'+'CORRECTOR_EOS_step'+str(step_num)+'_mat'+str(fiss_mat_id)+'.f71')

      # once all materials are done - append for record keeping for celi data at this iteration
      celi_results_lib.append_isotopics_data(iteration_num=ci, step_num=step_num, matlib=corrected_mat_lib)

    # append coirrected lib to time_lib now that we are all done iterating
    time_lib.append_lib(time=step_end_time, lib=corrected_mat_lib, step=step_num+1, PC_flag='C') # time lib step_num=1 is T1 at step 0

    # append final lib at BOS for next step for celi data storage
    celi_results_lib.append_isotopics_data(iteration_num=-1, step_num=step_num+1, matlib=corrected_mat_lib)

    # write all data in this state to pkl file before EOS
    celi_results_lib.write_state_to_pkl(step_num, is_final=False, folder='CASE_'+case_name)

    # TODO now append converged lib to time_lib



  ##############################################################
  ###################### DATA OUTPUT ###########################
  ##############################################################

  # writes final output as pkl file
  celi_results_lib.write_state_to_pkl(step_num='asdasd', is_final=True, folder='CASE_'+case_name)


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
