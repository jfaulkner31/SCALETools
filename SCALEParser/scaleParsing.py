import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['figure.figsize'] = [10, 4]

def burnupPowerPlot(powers):
  first_key = next(iter(powers))
  grid_powers = np.zeros(len(powers[first_key]))
  for key in powers.keys():
    grid_powers = np.vstack([grid_powers, powers[key]])
  grid_powers = grid_powers[1:]
  plt.pcolor(grid_powers)
  return grid_powers

# a variety of definitions that are used for parsing scale message (.msg) files.
def open_message_file(file):
  with open(file, 'r') as fh:
    lines = fh.readlines()
  return lines

def find_index_of_string(lines, string):
  index = int(len(lines)-1)
  for idx,this in enumerate(lines):
    if string in this:
      index = idx
      break
  return index

def get_keff_generations(lines, run_number=1):
  # run number = keff run number - aka kenovi run number
  # lines = lines fro open_message_file
  if run_number < 1:
      raise Exception("get_keff_generation run_number must be 1 or larger.")
  if run_number ==1:
    start_counting_index = find_index_of_string(lines, "generation     average      avg k-eff      generation    elapsed time")
    stop_counting_index = find_index_of_string(lines[start_counting_index:], "best estimate system k-eff")
  else:
    raise Exception ("run_number = 1 is only current supported number of ones currently...")

  # data arrays
  gen = np.array([])
  genkeff = np.array([])
  avgkeff = np.array([])
  avgdev = np.array([])
  entropy = np.array([])
  time = np.array([])
  for idx,this in enumerate(lines[start_counting_index:stop_counting_index]):
    if "E+0" in this:
      values = this.split()
      gen = np.append(gen, float(values[0]))
      genkeff = np.append(genkeff, float(values[1]))
      avgkeff = np.append(avgkeff, float(values[2]))
      avgdev = np.append(avgdev, float(values[3]))
      entropy = np.append(entropy, float(values[4]))
      time = np.append(time, float(values[5]))
  linedata = {
    "Generation": gen,
    "Generation keff": genkeff,
    "Average keff": avgkeff,
    "Average Deviation": avgdev,
    "Shannon Entropy": entropy,
    "Runtime": time
  }

  return linedata

def get_keff_depletion(lines):
  keff_list = np.array([])
  sigma_list = np.array([])
  for idx,this in enumerate(lines):
    if "best estimate system k-eff" in this:
      words = this.split()
      keff_list = np.append(keff_list, float(words[4]))
      sigma_list = np.append(sigma_list, float(words[8]))

  return keff_list, sigma_list

def get_mixture_powers(lines, is_msg_file=True):
  mixture_powers = {}
  for idx,this in enumerate(lines):
    if " multi-zone depletion..." in this:
      break
    if "mixture=" in this:
      words = this.split()
      mixture = words[1]
      mixnum = int(mixture[8:])
      mixture_powers[mixnum] = np.array([])

  for idx,this in enumerate(lines):
    if "mixture=" in this:
      words = this.split()
      mixture = words[1]
      power = float(words[3][6:-1])
      mixnum = int(mixture[8:])
      mixture_powers[mixnum] = np.append(mixture_powers[mixnum], power)

  keys = list(mixture_powers.keys())
  # output returns substep information as well so we correct to only get powers after transport step is done
  if is_msg_file: # do for .msg files
    for key in keys:
      mixture_powers[key] = mixture_powers[key][::2]
      mixture_powers[key] = np.append(mixture_powers[key], mixture_powers[key][-1])
  else:
    for key in keys:
      mixture_powers[key] = mixture_powers[key][::1]
      mixture_powers[key] = np.append(mixture_powers[key], mixture_powers[key][-1])

  return mixture_powers

def get_axial_offset(mixture_powers, bottom_half, top_half):
  axial_offset = np.array([])
  for timestep, garbage in enumerate(mixture_powers[bottom_half[0]]):
    bottom_sum = 0
    top_sum = 0
    for index in bottom_half:
      bottom_sum += mixture_powers[index][timestep]
    for index in top_half:
      top_sum += mixture_powers[index][timestep]
    val = (top_sum - bottom_sum) / (top_sum + bottom_sum)
    axial_offset = np.append(axial_offset, val)

  return axial_offset

def make_data_dict(file, top_half, bottom_half, csvname="outputCSV.csv", do_csv = False):

  lines = open_message_file(file) # gets line by line from file
  keff_list, sigma_list = get_keff_depletion(lines) # gets keff list and sigma list for depletion
  mixture_powers = get_mixture_powers(lines)  # gets dict with mixture powers by index
  axial_offset = get_axial_offset(mixture_powers, top_half, bottom_half)

  mixture_powers["keff"] = keff_list
  mixture_powers["sigma"] = sigma_list
  mixture_powers["axial_offset"] = axial_offset
  df = pd.DataFrame(mixture_powers)
  if do_csv:
    df.to_csv(csvname, index=False)
    print("CSV of name", csvname, "created with pandas!")

  return df

def combine_csv_files(csv_files):
  # lists all keffs, sigmas, axial offsets. then lists all power densities as averages. then combines this into a csv file
  # first open all csv files and list all keffs and their sigmas
  num_files = len(csv_files)
  data = {}
  keys = []
  # first list all keffs together
  for idx, file in enumerate(csv_files):
    csv_base = file[0:-4]
    dict_data = pd.read_csv(file)
    num_steps = len(dict_data["keff"].values)
    data["stepNum"] = np.linspace(1, num_steps, num_steps)
    data["keff_"+csv_base] = dict_data["keff"].values
    keys = list(dict_data.keys())


  # next do all sigmas
  for idx, file in enumerate(csv_files):
    csv_base = file[0:-4]
    dict_data = pd.read_csv(file)
    data["sigma_"+csv_base] = dict_data["sigma"].values
    keys = list(dict_data.keys())

  # next do all axial offsets together
  for idx, file in enumerate(csv_files):
    csv_base = file[0:-4]
    dict_data = pd.read_csv(file)
    data["axOff_"+csv_base] = dict_data["axial_offset"].values
    keys = list(dict_data.keys())


  # next we want to average all power density. we keep the keys from prior
  for key in keys:
    if ("keff" in key) | ("sigma" in key) | ("axial" in key):
      a = 1 # do nothing

    else: # else it is a power density key we are woirking with
      if key not in list(data.keys()): # if key isnt in data.keys, add it
        data[key] = None
      for idx, file in enumerate(csv_files): # start iterating over the files since we are on a correct key
        dict_data = pd.read_csv(file) # read file
        if idx == 0:
          data[key] = dict_data[key].values/num_files
        else:
          data[key] = data[key]  + dict_data[key].values/num_files

  df = pd.DataFrame(data)
  df.to_csv("csv_combined.csv", index=False)
  return "csv_combined.csv"



def csv_set_to_book(csv_files, output_excel_file):

  # write everything to a book.xlsx
  with pd.ExcelWriter(output_excel_file, engine="openpyxl") as writer:
    for csv_file in csv_files:
      # Read each CSV file into a DataFrame
      df = pd.read_csv(csv_file)

      # Use the file name (without extension) as the sheet name
      sheet_name = csv_file.split(".")[0]

      # Write the DataFrame to a sheet
      df.to_excel(writer, sheet_name=sheet_name, index=False)

  print(f"Excel workbook '{output_excel_file}' created successfully!")

