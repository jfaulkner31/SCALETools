# definitions for getting comps and stuffs.
class material_normal:
  def __init__(self):
    self.isotope_list = np.array([])
    self.mat_id = int(-1)
    self.biasid = int(-1)
    self.atom_dens = []
    self.temp = -1.0

class material_lib:
  def __init__(self):
    self.material_dict = {}

  def check_if_id_exists(self, test_id):
    if test_id in self.material_dict.keys():
      return True
    else:
      return False

  def append_mat_to_lib(self, material):
    self.material_dict[material.mat_id] = material

  def merge_mat(self, material_new):
      self.material_dict[material_new.mat_id].isotope_list = np.append(self.material_dict[material_new.mat_id].isotope_list, material_new.isotope_list)
      self.material_dict[material_new.mat_id].atom_dens = np.append(self.material_dict[material_new.mat_id].atom_dens, material_new.atom_dens)

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
        print(line)
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
        print('Append material ', mat_id, "to library!")
        matlib.append_mat_to_lib(new)
      else:
        # matid already exists so we need to add it to the library - merges with material with same id.
        print('Merging another instance of material ', mat_id, "to library!")
        matlib.merge_mat(new)

  return matlib

