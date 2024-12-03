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

