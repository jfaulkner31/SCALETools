"""
This folder makes data structures for outputting and manipulating CELI results in python.
"""
from getComps import material_lib
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt

def open_pkl(filename):
  # opens the CELI pkl file and returns a ResultsCELI object from a completed run
  with open(filename, "rb") as openfile:
    return pickle.load(openfile)


class ResultsCELI:
  def __init__(self,
              fissionable_mats=list,
              fissionable_mats_vols=list,
              residual_number_density=float,
              include_non_fission_material_power=bool,
              print_transport_powers=bool,
              system_IHM_mass_grams=float,
              specific_power=list,
              steplength_days=list,
              origen_predictor_divs=int,
              addnuxdictbase=str,
              base_triton=str,
              origen_base=str,
              origenResults_F71dir=str,
              MonteCarloResults_F33dir=str,
              Nprocs=int,
              machinefile=str,
              tmpdir=str,
              is_parallel=bool,
              origen_LI_divs=int,
              origen_steps_per_div=int,
              corrector_iterations=int,
              relaxation_factor=float,
              case_name=str):

    # if we are making the class for the first time then we must make a variable case settings and attach everything
    self.case_settings = {}
    self.case_settings['fissionable_mats'] = fissionable_mats
    self.case_settings['fissionable_mats_vols'] = fissionable_mats_vols
    self.case_settings['residual_number_density'] = residual_number_density
    self.case_settings['include_non_fission_material_power'] = include_non_fission_material_power
    self.case_settings['print_transport_powers'] = print_transport_powers
    self.case_settings['system_IHM_mass_grams'] = system_IHM_mass_grams
    self.case_settings['specific_power'] = specific_power
    self.case_settings['steplength_days'] = steplength_days
    self.case_settings['origen_predictor_divs'] = origen_predictor_divs
    self.case_settings['addnuxdictbase'] = addnuxdictbase
    self.case_settings['base_triton'] = base_triton
    self.case_settings['origen_base'] = origen_base
    self.case_settings['origenResults_F71dir'] = origenResults_F71dir
    self.case_settings['MonteCarloResults_F33dir'] = MonteCarloResults_F33dir
    self.case_settings['Nprocs'] = Nprocs
    self.case_settings['machinefile'] = machinefile
    self.case_settings['tmpdir'] = tmpdir
    self.case_settings['is_parallel'] = is_parallel
    self.case_settings['origen_LI_divs'] = origen_LI_divs
    self.case_settings['origen_steps_per_div'] = origen_steps_per_div
    self.case_settings['corrector_iterations'] = corrector_iterations
    self.case_settings['relaxation_factor'] = relaxation_factor
    self.case_settings['case_name'] = case_name

    # now dealing with the actual data and results
    self.keff_all = {} # keff by steps and then by another internal dict by predictor iterations : keff_all[step][iteration]
    self.BOS_keff_lines = {} # correct keff at BOS for each real step
    self.isotopics_all = {} # isotopics by steps and also by corrector iterations
    self.BOS_isotopics = {} # correct isotopics at BOS for each real step
    self.power_all = {}
    self.BOS_power = {}
  def append_isotopics_data(self, step_num, iteration_num, matlib):
    try:
      self.isotopics_all[step_num]
    except:
      self.isotopics_all[step_num] = {}
    self.isotopics_all[step_num][iteration_num] = matlib
    if iteration_num == -1:
      self.BOS_isotopics[step_num] = matlib

  def append_keff_data(self, keff_line, step_num, iteration_num):
    try:
      self.keff_all[step_num] # try to access for this step number and if not we make it a blank dict.
    except:
      self.keff_all[step_num] = {}
    self.keff_all[step_num][iteration_num] = keff_line
    # if predictor iteration
    if iteration_num == -1:
      self.BOS_keff_lines[step_num] = keff_line

  def append_power_data(self, power_data_block, step_num, iteration_num):
    try:
      self.power_all[step_num]
    except:
      self.power_all[step_num] = {}
    self.power_all[step_num][iteration_num] = power_data_block
    if iteration_num == -1:
      self.BOS_power[step_num] = power_data_block

  def write_state_to_pkl(self, step_num, is_final, folder):
    # removes old savestate and writes new one
    if is_final:
      print("\n\nWriting pkl state for final results...")
      name = self.case_settings['case_name']
      with open(folder+'/'+name+'_FINAL_.pkl', 'wb') as file:
        pickle.dump(self, file)


    else:
      print("\n\nWriting pkl state for step_num....", step_num)
      name = self.case_settings['case_name']
      with open(folder+'/'+name+'_upTo_stepNum_'+str(step_num)+'.pkl', 'wb') as file:
        pickle.dump(self, file)

      print("Removing pkl state for step_num", step_num-1)
      if os.path.exists(folder+'/'+name+'_upTo_stepNum_'+str(step_num-1)+'.pkl'):
        os.remove(folder+'/'+name+'_upTo_stepNum_'+str(step_num-1)+'.pkl')
      else:
        print("Tried removing", folder+'/'+name+'_upTo_stepNum_'+str(step_num-1)+'.pkl', "but could not as it does not exist.")

  def print_case_settings(self):
    for key in self.case_settings.keys():
      print(key+':', self.case_settings[key])
    print("\n\nAvailable Variables: ")
    print("\n- BOS values are the 'true' values at the T0 part of the predictor corrector iteration.")
    print("- 'all' values are the true values in addition to the corrector values in the corrector iterations.")
    print("- In the 'all' values, values with a key of -1 refer to predictor calcs at BOS.")
    print("- This means any key of -1 is the BOS values at T0 of that step.")

    print("\nkeff_all - Every keff value including predictor + corrector iterations")
    print("isotopics_all - Every isotopic library including predictor + corrector iterations")
    print("power_all - Every power library including predictor + corrector iterations")

    print("BOS_keff_lines - BOS keff lines from the BOS predictor calculation")
    print("BOS_isotopics - BOS isotopics for each step.")
    print("BOS_power - BOS power from the BOS power predictor calculations.")

  def print_methods(self):
    print("available methods for postprocessing results:")
    print("out.print_case_settings()")
    print("out.get_BOS_keffs()")
    print("out.get_BOS_isotope()")
    print("out.get_BOS_power()")
    print("out.get_corrector_keffs()")
    print("out.get_BOS_AO")
    print("\n\nneed to make the following still")
    print("out.get_corrector_power()")
    print("out.get_corrector_isotope()")
    print("out.get_corrector_AO()")

  def get_BU_vectors(self):
    """
    Outputs:
      time (np.array): Time in days
      bu (np.array): Burnup in GWD/MTIHM (GWdays per metric ton initial heavy metal)
    """
    time = np.array([0])
    bu = np.array([0])
    for idx, d in enumerate(self.case_settings['steplength_days']):
      time = np.append(time, time[-1]+d) # DAYS
      bu = np.append(bu, bu[-1]+d*self.case_settings['specific_power'][idx] / 1000) # GWD/MTIHM

    return time, bu

  def get_BOS_keffs(self):
    full_time, full_bu = self.get_BU_vectors()
    BOS_keffs = np.array([])
    BOS_sigmas = np.array([])
    BOS_BU = np.array([])
    BOS_time = np.array([])
    for idx, key in enumerate(self.BOS_keff_lines.keys()):
      line = self.BOS_keff_lines[key].split()
      BOS_keffs = np.append(BOS_keffs, float(line[4]))
      BOS_sigmas = np.append(BOS_sigmas, float(line[-1]))
      BOS_BU = np.append(BOS_BU, full_bu[idx])
      BOS_time = np.append(BOS_time, full_time[idx])
    return BOS_keffs, BOS_sigmas, BOS_BU, BOS_time

  def get_BOS_isotope(self, material_id, isotope):
    full_time, full_bu = self.get_BU_vectors()
    BOS_BU = np.array([])
    BOS_time = np.array([])
    iso_by_time = np.array([])
    for idx, step in enumerate(self.BOS_isotopics):
      try:
        val = self.BOS_isotopics[step].material_dict[material_id].return_iso_atom_dens(iso=isotope)
        iso_by_time = np.append(iso_by_time, val)
        BOS_BU = np.append(BOS_BU, full_bu[idx])
        BOS_time = np.append(BOS_time, full_time[idx])
      except:
        raise Exception("Error occured - likely due to incorrect isotope being requested.")
    return iso_by_time, BOS_BU, BOS_time

  def get_BOS_power(self):
    full_time, full_bu = self.get_BU_vectors()
    BOS_BU = np.array([])
    BOS_time = np.array([])
    for idx, step in enumerate(self.BOS_power.keys()):
      BOS_BU = np.append(BOS_BU, full_bu[idx])
      BOS_time = np.append(BOS_time, full_time[idx])
    return BOS_BU, BOS_time, self.BOS_power

  def get_BOS_AO(self):
    full_time, full_bu = self.get_BU_vectors()
    BOS_BU = np.array([])
    BOS_time = np.array([])
    ao = np.array([])
    for idx, key in enumerate(self.BOS_power.keys()):
      BOS_BU = np.append(BOS_BU, full_bu[idx])
      BOS_time = np.append(BOS_time, full_time[idx])

      this_powers = self.BOS_power[key]
      mat_ids = this_powers.keys()
      numkeys = len(mat_ids)
      bottom = 0.0
      top = 0.0
      for z, id in enumerate(mat_ids):
        this_power = self.BOS_power[key][id]
        if z < numkeys/2:
          bottom += this_power
        else:
          top += this_power
      this_ao = (top - bottom) / (top+bottom)
      ao = np.append(ao, this_ao)


    return BOS_BU, BOS_time, ao

  def get_corrector_keffs(self, step_num):
    """
    Returns the corrector-iterated keffs for the T0 value for this step.
    """
    final_keff = float
    final_sigma = float
    keffs = np.array([])
    sigmas = np.array([])
    time = float
    bu = float
    self.keff_all

    step_num -= 1
    try:
      for it in self.keff_all[step_num].keys():
        if it >= 0:
          this_k = float(self.keff_all[step_num][it].split()[4])
          this_sigma = float(self.keff_all[step_num][it].split()[-1])
          keffs = np.append(keffs, this_k)
          sigmas = np.append(sigmas, this_sigma)
    except:
      raise Exception("keff_all[step_num] does not exist for step_num"+str(step_num+1))

    try:
      final_keff = float(self.keff_all[step_num+1][-1].split()[4])
      final_sigma = float(self.keff_all[step_num+1][-1].split()[-1])
    except:
      final_keff = None
      final_sigma = None
    keffs = np.append(keffs, final_keff)
    sigmas = np.append(sigmas, final_sigma)

    return keffs, sigmas

  def plot_BOS_power(self, figsize: tuple, normalize: bool, fontsize=14, fontname='Cambria'):
    BOS_BU, BOS_time, power_dict = self.get_BOS_power()

    # materials
    mats = power_dict[0].keys()

    # first make array for power:
    parr = np.array([])
    for materialKey in power_dict[0].keys():
      mat_vs_time = []
      for timestepKey in power_dict.keys():
        p = power_dict[timestepKey][materialKey]
        mat_vs_time.append(p)
      try:
        parr = np.vstack((parr, mat_vs_time))
      except:
        parr = mat_vs_time

    # normalize power if requested
    if normalize:
      # power is either normalied to 1.0 from this or the power fraction as returned by scale.
      # power fraction not always 1.0 since gamma heating accounted for in nonfissile materials as well.
      # see runCELI.py input option - include_non_fission_material_power
      parr = parr / np.sum(parr, axis=0)


    # now make pcolormesh
    fig, ax = plt.subplots(figsize=figsize)
    c = ax.pcolormesh(BOS_time, mats, parr)

    # colorbar
    if normalize:
      colorbar = plt.colorbar(c, ax=ax)
      colorbar.set_label('Power (normalized to unity)',
                    fontdict={'fontsize': fontsize,
                              'fontname': fontname})
    else:
      colorbar = plt.colorbar(c, ax=ax)
      colorbar.set_label('Power fraction',
                    fontdict={'fontsize': fontsize,
                              'fontname': fontname})

    # other
    ax.set_xlabel('Burnup step number',
                  fontdict={'fontsize': fontsize,
                            'fontname': fontname})

    ax.set_ylabel('Material ID',
                  fontdict={'fontsize': fontsize,
                            'fontname': fontname})
    return parr






