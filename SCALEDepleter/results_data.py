"""
This folder makes data structures for outputting CELI results in python.
"""
from getComps import material_lib
import pickle
import os

class ResultsCELI:
  def __init__(self,
              fissionable_mats: list,
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
              case_name: str):
    # make a variable case settings and attach everything
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

  def write_state_to_pkl(self, step_num, is_final):
    # removes old savestate and writes new one
    if is_final:
      print("\n\nWriting pkl state for final results...")
      name = self.case_settings['case_name']
      with open(name+'_FINAL_.pkl', 'wb') as file:
        pickle.dump(self, file)


    else:
      print("\n\nWriting pkl state for step_num....", step_num)
      name = self.case_settings['case_name']
      with open(name+'_upTo_stepNum_'+str(step_num)+'.pkl', 'wb') as file:
        pickle.dump(self, file)

      print("Removing pkl state for step_num", step_num-1)
      if os.path.exists(name+'_upTo_stepNum_'+str(step_num-1)+'.pkl'):
        os.remove(name+'_upTo_stepNum_'+str(step_num-1)+'.pkl')
      else:
        print("Tried removing", name+'_upTo_stepNum_'+str(step_num-1)+'.pkl', "but could not as it does not exist.")

