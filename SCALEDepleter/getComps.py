import numpy as np
import re

# definitions for getting comps and stuffs.
class material_normal:
  def __init__(self):
    self.isotope_list = np.array([])
    self.mat_id = int(-1)
    self.biasid = int(-1)
    self.atom_dens = []
    self.temp = -1.0
  def make_material_lines(self):
    lines = []
    for idx, iso in enumerate(self.isotope_list):
      line = ''
      line+=iso+"   "
      line+=str(self.mat_id)+"   "
      line+="0   "
      line+=str(self.atom_dens[idx])+"   "+str(self.temp)+"   end"
      lines+=[line]
    return lines

  def append_mat_info(self, iso: str, atom_density: float):
    if iso == 'C' or iso == 'c':
      iso = 'c12' # converts iso of C from stdcomp to carbon 12
    elif len(iso) == 1:
      raise Exception("len(iso) == 1 for iso with label"+iso+" - clearly an edge case that hasnt been considered and properly converted to zaid from stdcmp")
    self.isotope_list = np.append(self.isotope_list, iso)
    self.atom_dens = np.append(self.atom_dens, atom_density)

  def make_origen_materials(self):
    line = '    iso=['
    for idx, iso in enumerate(self.isotope_list):
      line += ' '+iso+'='+str(self.atom_dens[idx])
      if (idx+1) % 10 == 0:
        line+="\n"
    line+=']'
    return line

  def return_iso_atom_dens(self, iso: str):
    try:
      lower_iso = iso.lower()
      lower_list = [item.lower() for item in self.isotope_list]
      idx = lower_list.index(lower_iso)
      return self.atom_dens[idx]
    except:
      raise Exception("Requested isotope "+iso+" not found!")







class material_lib:
  def __init__(self):
    self.material_dict = {}
    self.time = float
    self.step = int

  def check_if_id_exists(self, test_id):
    if test_id in self.material_dict.keys():
      return True
    else:
      return False

  def append_mat_to_lib(self, material: material_normal):
    self.material_dict[material.mat_id] = material

  def merge_mat(self, material_new):
      self.material_dict[material_new.mat_id].isotope_list = np.append(self.material_dict[material_new.mat_id].isotope_list, material_new.isotope_list)
      self.material_dict[material_new.mat_id].atom_dens = np.append(self.material_dict[material_new.mat_id].atom_dens, material_new.atom_dens)

  def add_residual_isotope(self, numdens: float, isotope: str):
    """
    adds a residual isotope amount to every single material in the material lib
    """
    for key in self.material_dict.keys():
      self.material_dict[key].append_mat_info(iso=isotope, atom_density=numdens)



class time_dependent_material_lib:
  def __init__(self):
    self.mats_by_steps = {}
    self.steps = []
    self.times = []

  def append_lib(self, lib: material_lib, step: int, time: float, PC_flag: str):
    lib.time = time
    self.times.append(time)

    lib.step = step
    self.steps.append(step)

    self.mats_by_steps[step] = {PC_flag: lib}

  def dump_info(self):
    print("Timelib:"+time_dependent_material_lib)
    for step in self.steps:
      print("\t", self.mats_by_steps[step])

  def material_lib_from_time(self, requested_time: float, PC_flag: 'str'):
    for idx, t in enumerate(self.times):
      if t <= requested_time+1e-12 | t>= requested_time-1e-12:
        # found the time and step
        return self.mats_by_steps[idx][PC_flag]
    raise Exception("Could not find material library for requested time - "+str(time)+" with requested PC_Flag - " + PC_flag)

  def material_lib_from_step(self, requested_step: int, PC_flag: 'str'):
    try:
      returned_lib = self.mats_by_steps[requested_step][PC_flag]
      return returned_lib
    except:
      raise Exception("Could not find a material library for requested step - "+str(requested_step)+" and with reques ted PC_Flag - " + PC_flag)




def get_comps(file):
  # returns a material_lib object with all materials and their associated data based on the input file
  matlib = material_lib()
  do_comp = False
  for wholeline in file.readlines():
    line = wholeline.split() # -> line = ['this', 'is', 'how', 'line', 'split', 'looks']

    # starts and stops reading of compositions
    if 'read' in line and 'comp' in line:
      if line.index('read') == line.index('comp') - 1:
        print("read comp found! Now starting comp directory.")
        do_comp = True # -> activates composition mode
    if 'end' in line and 'comp' in line:
      if line.index('end') == line.index('comp') - 1:
        print("end comp found! Now finalizing comp directory.")
        break # -> break the for loop


    if (do_comp & ('end' in line)):
      if (line[2] == 'end'):
        # this case is custom scale standard comp
        isotope = line[0]
        mat_id = line[1]
        dummy = ''
        dens = ''
        temp = ''
      else:
        # get values in current line for this material of type [iso 1 0 dens temp end]
        isotope = line[0]
        mat_id = int(line[1])
        dummy = line[2]
        # print(line)
        dens = line[3]
        temp = line[4]

      # make new material def
      new = material_normal()
      new.isotope_list = isotope
      new.mat_id = mat_id
      new.biasid = dummy
      new.atom_dens = dens
      new.temp = temp
      if (matlib.check_if_id_exists(mat_id) == False):
        # make new material and append to library if material doesnt already exist
        # print('Append material ', mat_id, "to library!")
        matlib.append_mat_to_lib(new)
      else:
        # matid already exists so we need to add it to the library - merges with material with same id.
        # print('Merging another instance of material ', mat_id, "to library!")
        matlib.merge_mat(new)

  return matlib

def get_comps_from_std_mix_file(filename):
  mat = material_normal()
  with open(filename, 'r') as file:
    lines = file.readlines()
  for idx, line in enumerate(lines):
    if ('concentrations in mixture' in line) | ('material number' in line): # skip these lines -
      continue
    else:
      parts = line.strip().split()
      nuclide = parts[0]
      dash = nuclide.find('-')
      nuclide = nuclide[0:dash].lower()+'-'+nuclide[dash+1::] # probably best to keep the dash for consistency?
      mat_id = int(parts[1])
      zero = int(parts[2])
      atomdens = float(parts[3])
      temp = float(parts[4])
      end = parts[5]
      mat.append_mat_info(iso=nuclide, atom_density=atomdens)
  mat.biasid = 0
  mat.temp = temp
  mat.mat_id = mat_id
  return mat


def split_isotope(isotope):
  match = re.match(r"([a-zA-Z]+)(\d+)", isotope)
  if '-' in isotope:
    return -1, -1, isotope
  elif match:
    letters = match.group(1)  # Extract the letters
    numbers = match.group(2)  # Extract the numbers
    return letters, numbers, letters+'-'+numbers
  else:
    letters = match.group(1)  # Extract the letters
    numbers = match.group(2)  # Extract the numbers
    print(letters)
    print(numbers)
    print(isotope)
    raise ValueError(f"Invalid isotope format: {isotope} {letters} {numbers}")

def makeNewAddnuxDict(zeromatdict, tmpdir, addnuxdict):
  extraAddnuxIsotopes = []
  in_dict_but_not_input = []
  zero_mat_dict = zeromatdict
  for key in zero_mat_dict.keys():
    iso_list = zero_mat_dict[key].isotope_list
    for isotope in iso_list:
      _1, _2, new_iso = split_isotope(isotope) # returns proper format as variable new_iso li6 -> li-6
      if new_iso.lower() not in extraAddnuxIsotopes:
        extraAddnuxIsotopes.append(new_iso.lower())

        # now make list of unique isotopes in addnuxdict
  with open(addnuxdict, 'r') as file:
    lines = file.readlines()

  addnuxList = []
  for line in lines:
    line = line.replace("\n", "")
    if (line.lower() not in extraAddnuxIsotopes) & (len(line) > 0):
      extraAddnuxIsotopes.append(line)
      in_dict_but_not_input.append(line)

  # now make new file for addnuxdict and write to it
  userAddNuxDict = open(tmpdir+'/usraddnuxdict.dict', 'w', encoding="utf-8")
  for iso in extraAddnuxIsotopes:
    userAddNuxDict.write(iso+'\n')






  userAddNuxDict.close()
  addnuxdict = tmpdir+'/usraddnuxdict.dict' # set equals for use later - new addnuxdict is now in the temp directory. under tmpdir+'/usraddnuxdict.dict'

  return addnuxdict, in_dict_but_not_input
