import numpy as np

def getPower(filename: str, include_non_fission_material_power: bool, fission_mat_ids, printP: bool, total_power_python: float):
  filename = filename.split('.')[0]+'.out'

  target_string_631 = "Number (MW/MTIHM)    (---)   (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)  n/(cm^2*sec)"
  target_string_632 = "Number     (MW/MTIHM)    (---)       (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)     n/(cm^2*sec)"
  found_line = False
  power_dict = {}
  fissionable_mat_power = 0.0
  nonfissionable_mat_power = 0.0
  with open(filename, 'r') as file:
    lines = file.readlines()
  for line_number, line in enumerate(lines):  # Use enumerate for line numbers
    if (target_string_631 in line) | (target_string_632 in line):
      found_line = True
    elif found_line:
      powerline = line.split()
      if 'Total' in powerline[0]:
        break
      else:
        normalized_power = float(powerline[1])
        mix = int(powerline[0])
        if mix in fission_mat_ids:
          power_dict[mix] = normalized_power
          fissionable_mat_power += normalized_power
        else:
          nonfissionable_mat_power += normalized_power
          continue
  # now we have:
  # output_total_power = total poewr in the output of ALL materials
  # fissionable_mat_power = total power across all materials that are fissionable - aka all materials in our depletion input
  # power_dict[mat_id] = dict with powers for each material.
  fissionable_to_total_ratio = fissionable_mat_power / (fissionable_mat_power + nonfissionable_mat_power)

  # if we are including non fissionable materials under the power definition - we need to normalize power_dict to output_total_power
  if include_non_fission_material_power:
    for key in power_dict.keys():
      fractional_power_from_fissionable_materials = power_dict[key] / fissionable_mat_power
      power_dict[key] = fractional_power_from_fissionable_materials * fissionable_to_total_ratio

  else: # if we are NOT including non fissionable materials in the power definition
    for key in power_dict.keys():
      power_dict[key] = power_dict[key]/fissionable_mat_power

  if printP:
    print("NOW PRINTING POWERS")
    for key in power_dict.keys():
      print("\t",key,"|",power_dict[key],"|",power_dict[key]*total_power_python)

  return power_dict

def interpolatePower(power_by_step_t0: dict, power_by_step_t1: dict , times, start_time, end_time):
  """
  Takes in power_by_step for t0 and t1 and interpolates based on time value

  Times - midpoints of substep times.

  Start time - t0 (Beginning of big step - when MC is first ran)

  End time - t1 (End of big step - when MC is ran second)
  """
  mat_ids = power_by_step_t0.keys()
  norm_t0 = 0.0
  norm_t1 = 0.0

  # get norms (aka total power) at t0 and t1
  for mat_id in mat_ids:
    norm_t0 += power_by_step_t0[mat_id]
    norm_t1 += power_by_step_t1[mat_id]

  # interpolate norms to get a norm for every substep.
  interp_norms = np.interp(times, [start_time, end_time], [norm_t0, norm_t1])

  # make a dict with keys of mat_id and each entry is a vector of interpolated times.
  power_vs_time_dict = {}
  for mat_id in mat_ids:
    p0 = power_by_step_t0[mat_id]
    p1 = power_by_step_t1[mat_id]
    interp_power_this_mat = np.interp(times, [start_time, end_time], [p0, p1])
    power_vs_time_dict[mat_id] = interp_power_this_mat

  # now get sums / norms of each timestep for the interpolated power dict
  for idx, time in enumerate(power_vs_time_dict[   int( list(mat_ids)[0] )   ]):
    this_timesteps_sum = 0.0
    for mat_id in mat_ids:
      this_timesteps_sum += power_vs_time_dict[mat_id][idx]

    # normalize all the materials in this timestep
    for mat_id in mat_ids:
      power_vs_time_dict[mat_id][idx] *= interp_norms[idx] / this_timesteps_sum


  # we now have power_vs_time_dict which is a dict with mat_ids where each entry is a vector for the power vs time for that timestep.
  return power_vs_time_dict
