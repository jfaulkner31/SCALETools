def getPower(filename: str, include_non_fission_material_power: bool, fission_mat_ids, printP: bool, total_power_python: float):
  filename = filename.split('.')[0]+'.out'

  target_string = "Number (MW/MTIHM)    (---)   (MW/MTIHM)   (GWd/MTIHM)   n/(cm^2*sec)  n/(cm^2*sec)"
  found_line = False
  power_dict = {}
  fissionable_mat_power = 0.0
  nonfissionable_mat_power = 0.0
  with open(filename, 'r') as file:
    lines = file.readlines()
  for line_number, line in enumerate(lines):  # Use enumerate for line numbers
    if target_string in line:
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




