import numpy as np

def get_triton_depletion_table(filename: str):
  """
  Gets power for each material from scale output and
  also gets normalization factors for HM mass

  Parameters
  ----------
  filename : str
    name of scale output file from TRITON

  Returns
  -------
  table : dict
    table of power and flux information from TRITON
  normalization_factor : float
    total mass of HM

  """
  table = {}

  found_line_631_632 = False
  found_line_624 = False

  target_string_631 = "Number (MW/MTIHM)    (---)   (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)  n/(cm^2*sec)"
  target_string_632 = "Number     (MW/MTIHM)    (---)       (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)     n/(cm^2*sec)"
  target_string_624 = "Number (MW/MTIHM)    (---)   (MW/MTIHM)  n/(cm^2*sec)  n/(cm^2*sec)"


  filename = filename.split('.')[0]+'.out'
  with open(filename, 'r') as file:
    lines = file.readlines()

  for idx , line in enumerate(lines):
    if "Masses will be normalized by a factor of" in line:
      # this is the normalization line
      normline = line.split()
      for this_idx, word in enumerate(normline):
        if word == "of":
          idx_of_factor = this_idx + 1
          break
      # the normalization factor
      norm_factor = float(normline[idx_of_factor])
    if (target_string_631 in line) | (target_string_632 in line):
      found_line_631_632 = True
    if (target_string_624 in line):
      found_line_624 = True
    ### FOUND LINE 6.3.1 or 6.3.2 ###
    elif found_line_631_632:
      powerline = line.split()
      # print(powerline)
      if 'Total' in powerline:
        break
      else:
        mixture_number = int(powerline[0])
        total_power = float(powerline[1])
        frac_power = float(powerline[2])
        if "N/A" in powerline[3]:
          mixture_power = 0.0
        else:
          mixture_power = float(powerline[3])
        if "N/A" in powerline[4]:
          mixture_burnup = 0.0
        else:
          mixture_burnup = float(powerline[4])
        thermal_flux = float(powerline[5])
        total_flux = float(powerline[6])
        item = {"total_power": total_power,
                "frac_power": frac_power,
                "mixture_power": mixture_power,
                "mixture_burnup": mixture_burnup,
                "thermal_flux": thermal_flux,
                "total_flux": total_flux}
        table[mixture_number] = item
    ### FOUND LINE 6.2.4 ###
    elif found_line_624:
      powerline = line.split()
      # print(powerline)
      if 'Total' in powerline:
        break
      else:
        mixture_number = int(powerline[0])
        total_power = float(powerline[1])
        frac_power = float(powerline[2])
        if "N/A" in powerline[3]:
          mixture_power = 0.0
        else:
          mixture_power = float(powerline[3])
        mixture_burnup = 0.0 # in 6.2.4 always just set bu to 0.0 since it doesnt return this info
        thermal_flux = float(powerline[4])
        total_flux = float(powerline[5])
        item = {"total_power": total_power,
                "frac_power": frac_power,
                "mixture_power": mixture_power,
                "mixture_burnup": mixture_burnup,
                "thermal_flux": thermal_flux,
                "total_flux": total_flux}
        table[mixture_number] = item


  # total system mass will then be = 1e6 / norm_factor
  # since norm_factor = 1e6 / TOTAL_MASS

  #      table  norm_factor   total_mass ()
  return table, norm_factor, 1e6/norm_factor


def getPower(filename: str, include_non_fission_material_power: bool, fission_mat_ids, printP: bool, total_power_python: float):
  filename = filename.split('.')[0]+'.out'

  target_string_631 = "Number (MW/MTIHM)    (---)   (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)  n/(cm^2*sec)"
  target_string_632 = "Number     (MW/MTIHM)    (---)       (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)     n/(cm^2*sec)"
  target_string_624 = "Number (MW/MTIHM)    (---)   (MW/MTIHM)  n/(cm^2*sec)  n/(cm^2*sec)"
  found_line = False
  power_dict = {}
  fissionable_mat_power = 0.0
  nonfissionable_mat_power = 0.0
  with open(filename, 'r') as file:
    lines = file.readlines()
  for line_number, line in enumerate(lines):  # Use enumerate for line numbers
    if (target_string_631 in line) | (target_string_632 in line) | (target_string_624 in line):
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
    print("NOW PRINTING POWERS", "")
    for key in power_dict.keys():
      print("\t",key,"|",power_dict[key],"|",power_dict[key]*total_power_python)

  return power_dict

def interpolatePower(power_by_step_t0_pre: dict, power_by_step_t1_pre: dict , times, start_time, end_time,
                     specific_power_this_step: float, printP: bool):
  """
  Takes in power_by_step for t0 and t1 and interpolates based on time value

  Times - midpoints of substep times.

  Start time - t0 (Beginning of big step - when MC is first ran)

  End time - t1 (End of big step - when MC is ran second)
  """
  mat_ids = power_by_step_t0_pre.keys()
  norm_t0 = 0.0
  norm_t1 = 0.0

  # first multiply each element by specific power to get true power that we want to work with
  power_by_step_t0 = {}
  power_by_step_t1 = {}
  for mat_id in mat_ids:
    power_by_step_t0[mat_id] = power_by_step_t0_pre[mat_id] * specific_power_this_step
    power_by_step_t1[mat_id] = power_by_step_t1_pre[mat_id] * specific_power_this_step

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
