# a variety of definitions that are used for parsing scale message (.msg) files.

def open_message_file(file):
  with open(file, 'r') as fh:
    lines = fh.readlines()
  return lines

def find_index_of_string(lines, string):
  index = int(-2)
  for idx,this in enumerate(lines):
    if string in this:
      index = idx
      break
  return index

def get_keff_generations(lines):
  # first, use lines to get the index of the first generation.
  start_counting_index = find_index_of_string(lines, "generation     average      avg k-eff      generation    elapsed time")
  for idx,this in enumerate(lines[start_counting_index:]):
    print(this)
  return 1
