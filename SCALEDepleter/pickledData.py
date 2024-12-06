import pickle
from getComps import time_dependent_material_lib

def makeOutput(keff_lines: list, power_by_step: dict, specific_power: list, steplength_days: list, time_lib: time_dependent_material_lib, pkl_filename: str):

  output = {}
  output['keff_lines'] = keff_lines
  output['power_by_step'] = power_by_step
  output['specific_power'] = specific_power
  output['steplength_days'] = steplength_days
  output['time_lib'] = time_lib


  # make keff and sigma.
  keff_vec = []
  sigma_vec = []
  for line in keff_lines.split():
    keff_vec.append(float(line[4]))
    sigma_vec.append(float(line[8]))

  output['keff'] = keff_vec
  output['sigma'] = sigma_vec

  # make time
  time_vec = [0]
  for idx, this in enumerate(steplength_days[1:]):
    time_vec.append(this+time_vec[idx-1])
  output['time'] = time_vec

  # dump to a pickle file
  with open(pkl_filename, 'wb') as file:
    pickle.dump(output, file)

  # with open('keff.pkl', 'wb') as file:
  #   pickle.dump(keff_lines, file)
  # with open('power.pkl', 'wb') as file:
  #   pickle.dump(power_by_step, file)
  # with open('nuclides.pkl', 'wb') as file:
  #   pickle.dump(time_lib, file)

def getOutput(pkl_filename: str):
  with open(pkl_filename, 'rb') as file:
    output = pickle.load(file)
  return output
