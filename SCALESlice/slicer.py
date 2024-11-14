



# todo
# 2. add hexagons, xcyls, ycyls.
# 3. do radial subdivisions of a simple cylinder IMPORTANT FIX
# 5. do vol flag for perfect geometries


# imports
import numpy as np
import copy
from geometry_collection import *
import re
import sys
from colorama import Fore


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

# give a unit number and a file and it returns the unit block that is specified with all lines in a dictionary.
def get_unit(file, unit_num):
  # input an input file and a unit number and it returns a dict with line by line of that unit
  do_unit = False
  start_grabbing_lines = False
  linedict = {}
  idx = 0
  for wholeline in file.readlines():
    line = wholeline.split() # -> line = ['this', 'is', 'how', 'line', 'split', 'looks']

    # starts and stops reading of compositions
    if 'read' in line and 'geometry' in line:
      if line.index('read') == line.index('geometry') - 1:
        print("read geometry found! Now starting comp directory.")
        do_unit = True # -> activates composition mode
        # print(line)
    if 'end' in line and 'geometry' in line:
      if line.index('end') == line.index('geometry') - 1:
        print("end geometry found! Now finalizing comp directory.")
        break # -> break the for loop

    # once we access are in the unit block we care about - go line by line and get a dict with all the lines in it
    if (do_unit & (str(unit_num) in line) & ('unit' in line)):
      start_grabbing_lines = True

    if start_grabbing_lines:
      line = wholeline.split()
      linedict[idx] = line
      idx += 1
      if ('boundary' in line):
        break
  return linedict

def get_media_index(unitbank, medialine):
  # give a unit dictionary and it returns which key holds the media_line of the user input
  # inputs:
  # unitbank =  {0: line1, 1: line2 ...}
  # medialine = 'media 1 1 1'
  a = False
  sliced = medialine.split()
  newline = ''
  for word in sliced:
    newline += ' ' + word
  medialine = newline[1::]

  for key in unitbank.keys():
    newline = ''
    for word in unitbank[key]:
      newline += ' ' + word
    newline = newline[1::]
    # print(newline)

    if medialine == newline:
      a = True
      break

  if a == False:
    raise Exception("Medialine == newline -- line was never found as we searched! This is usually due to a wrong 'media' input in the user input! Make sure you copy the media card exactly as it appears in the .inp file. Also make sure that your unit_number is correct in your input")
  return newline, key


def surfaces_from_media_line(medialine):
  # pass in a medialine
  medialine = medialine.split()
  mat_id = medialine[1]
  importance = medialine[2]
  surfaces = medialine[3::]

  inside = np.array([])
  outside = np.array([])
  for this in surfaces:
    if int(this) < 0:
      inside = np.append(inside, int(this))
    elif int(this) > 0:
      outside = np.append(outside, int(this))

  all_surfaces = np.unique(  np.concatenate( (inside,outside) ) )
  return inside, outside, all_surfaces, importance, mat_id

def get_surface_cards(card_array, unitbank):
  # if not returning surfaces check geometry_collection.py
  surface_line_dict = {}
  for surf_float in card_array:
    # start at each surface number we want to get lines for
    surface_number = abs(int(surf_float))
    for key in unitbank.keys():
      thisline = unitbank[key]
      for geomtype in geometry_data.keys():
        if geomtype in thisline:
          # if it is a recognized geomtype -> proceed with search
          if surface_number == int(thisline[1]): # -> if surface number we are looking for is the surface number in the line
            # if true then we got a hit and the surface number we are looking for is this line!
            line_to_write = []
            idx = -1
            do_write_next = True
            for word in thisline[0::]:
              idx +=1
              if ('p' in word) & (idx !=0): # this indicates we are doing 6p or 2p -  avoid edge case where name of the shape might have a 'p' in it
                do_write_next = False
                index_p = word.index('p')
                integer_half = int(word[:index_p])
                dimension = float(thisline[idx+1]) # number to repeat
                for number in range(0,integer_half):
                  line_to_write.append( dimension* ((-1)**number)   )
              else:
                if do_write_next:
                  line_to_write.append(word)
                do_write_next = True

            surface_line_dict[surface_number] = line_to_write

  for surf_float in card_array:
    if int(abs(surf_float)) not in surface_line_dict.keys():
      raise ValueError("Surface number", int(surf_float), "was not found in card_array")
  return surface_line_dict


def get_origin_from_surface_card(card):
  # takes in a surface card. if it has 'origin' in it then it returns origin xyz components.
  new_item_list = []
  keep_going = False
  for item in card:

    if isinstance(item, str):
      if 'origin' in item:
        keep_going = True
      if keep_going:
        new_item_list.append(item)


  if keep_going:
    new_line = ''
    for item in new_item_list[1::]:
      new_line += item

    do_x = False
    do_y = False
    do_z = False
    locs = np.array([])
    if 'x' in new_line:
      do_x = True
      xpos = new_line.find('x')
      locs = np.append(locs, xpos)
    if 'y' in new_line:
      do_y = True
      ypos = new_line.find('y')
      locs = np.append(locs, ypos)
    if 'z' in new_line:
      do_z = True
      zpos = new_line.find('z')
      locs = np.append(locs, zpos)

    order = np.argsort(locs) # -> gives order of sorting
    final_positions = np.empty_like(order)
    final_positions[order] = np.arange(len(locs))
    order = final_positions

    new_line_1 = new_line.split('origin')

    new_line_1 = ''.join(new_line_1)

    new_line_1 = new_line_1.split('x')
    new_line_1 = ''.join(new_line_1)
    new_line_1 = new_line_1.split('y')
    new_line_1 = ''.join(new_line_1)
    new_line_1 = new_line_1.split('z')
    new_line_1 = ''.join(new_line_1)
    new_line_1 = new_line_1.split('=')
    new_line_1 = [string for string in new_line_1 if string != '']

    output = np.array([])
    for item in new_line_1:
      output = np.append(output, float(item))

    if do_x & do_y & do_z:
      x = output[order[0]]
      y = output[order[1]]
      z = output[order[2]]
    elif do_x & do_y:
      x = output[order[0]]
      y = output[order[1]]
      z = 0.0
    elif do_x & do_z:
      x = output[order[0]]
      y = 0.0
      z = output[order[1]]
    elif do_y & do_z:
      x = 0.0
      y = output[order[0]]
      z = output[order[1]]
    elif do_x:
      x = output[order[0]]
      y = 0.0
      z = 0.0
    elif do_y:
      x = 0.0
      y = output[order[0]]
      z = 0.0
    elif do_z:
      x = 0.0
      y = 0.0
      z = output[order[1]]
    return x,y,z
  else:
    return 0.0, 0.0, 0.0

class cuboid:
  def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, shiftx, shifty, shiftz, media_id):
    self._xmin = xmin
    self._xmax = xmax
    self._ymin = ymin
    self._ymax = ymax
    self._zmin = zmin
    self._zmax = zmax
    self._shiftx = shiftx
    self._shifty = shifty
    self._shiftz = shiftz
    self._media_id = media_id

  @property
  def media_id(self):
    return self._media_id

  @property
  def xmin(self):
    return self._xmin+self._shiftx
  @property
  def xmax(self):
    return self._xmax+self._shiftx

  @property
  def ymin(self):
    return self._ymin+self._shifty
  @property
  def ymax(self):
    return self._ymax+self._shifty

  @property
  def zmin(self):
    return self._zmin+self._shiftz
  @property
  def zmax(self):
    return self._zmax+self._shiftz

class cylinder:
  def __init__(self, radius, minh, maxh, orientation, shiftx, shifty, shiftz, media_id):
    self._radius = radius
    self._orientation = orientation
    self._minh = minh
    self._maxh = maxh
    self._shiftx = shiftx
    self._shifty = shifty
    self._shiftz = shiftz
    self._media_id = media_id

  @property
  def media_id(self):
    return self._media_id

  @property
  def orientation(self):
    return self._orientation

  @property
  def radius(self):
    return self._radius

  @property
  def minh(self):
    return self._minh

  @property
  def maxh(self):
    return self._maxh

  @property
  def xmin(self):
    if (self.orientation == 'z' or self.orientation == 'y'):
      return -self.radius+self._shiftx
    elif (self.orientation == 'x'):
      return self.minh+self._shiftx

  @property
  def xmax(self):
    if (self.orientation == 'z' or self.orientation == 'y'):
      return self.radius+self._shiftx
    elif (self.orientation == 'x'):
      return self.maxh+self._shiftx

  @property
  def ymin(self):
    if (self.orientation == 'z' or self.orientation == 'x'):
      return -self.radius+self._shifty
    elif (self.orientation == 'y'):
      return self.minh+self._shifty

  @property
  def ymax(self):
    if (self.orientation == 'z' or self.orientation == 'x'):
      return self.radius+self._shifty
    elif (self.orientation == 'y'):
      return self.maxh+self._shifty

  @property
  def zmin(self):
    if (self.orientation == 'y' or (self.orientation == 'x')):
      return -self.radius+self._shiftz
    elif (self.orientation == 'z'):
      return self.minh+self._shiftz

  @property
  def zmax(self):
    if (self.orientation == 'x' or self.orientation == 'y'):
      return -self.radius+self._shiftz
    elif (self.orientation == 'z'):
      return self.maxh+self._shiftz

def objects_from_surface_cards(surface_cards, all_surfaces):
  # returns array of objects from surface cards
  #   # surface_cards example: {4: ['cylinder', '4', '65.1175778839795', 65.1175778839795, -65.1175778839795],
  #                             99999: ['cuboid', '99999', 75.1175778839795, -75.1175778839795, 75.1175778839795, -75.1175778839795, 75.1175778839795, -75.1175778839795]}
  objectarray = np.array([])
  for idx, key in enumerate(surface_cards.keys()):
    object_type = surface_cards[key][0]
    this_object = surface_cards[key]
    shiftx,shifty,shiftz = get_origin_from_surface_card(this_object)

    if object_type == 'cylinder':
      name = this_object[0]
      surf_id = key
      r = float(this_object[2])
      minh = min([float(this_object[3]), float(this_object[4])])
      maxh = max([float(this_object[3]), float(this_object[4])])
      objectarray = np.append(objectarray, cylinder(r, minh, maxh, 'z', shiftx, shifty, shiftz, int(all_surfaces[idx]) ))

    if object_type == 'cuboid':
      name = this_object[0]
      surf_id = key
      xmin = min([float(this_object[2]), float(this_object[3])])
      xmax = max([float(this_object[2]), float(this_object[3])])
      ymin = min([float(this_object[4]), float(this_object[5])])
      ymax = max([float(this_object[4]), float(this_object[5])])
      zmin = min([float(this_object[6]), float(this_object[7])])
      zmax = max([float(this_object[6]), float(this_object[7])])

      objectarray = np.append(objectarray, cuboid(xmin, xmax, ymin, ymax, zmin, zmax, shiftx, shifty, shiftz, int(all_surfaces[idx] ) ))

  return objectarray

def CartesianGridSplit(surface_object_collection, medialine, ui_num_x_divs, ui_num_y_divs, ui_num_z_divs, plane_id_start, media_id_offset):
  # medialine example: 'media 10 1 -1 2'
  # surface_cards example: {4: ['cylinder', '4', '65.1175778839795', 65.1175778839795, -65.1175778839795],
                            # 99999: ['cuboid', '99999', 75.1175778839795, -75.1175778839795, 75.1175778839795, -75.1175778839795, 75.1175778839795, -75.1175778839795]}

  # for all media cards, we will have this list of 'senses' to start
  surface_list = ''
  # first determine the xyz
  for item in surface_object_collection:
    xmax = -9e20
    xmin = 9e20
    ymax = -9e20
    ymin = 9e20
    zmax = -9e20
    zmin = 9e20
    # find a bounding box around the domain with positive sense
    if item.media_id > 0:
      # print(item.media_id)
      xmax = max(xmax,item.xmax)
      ymax = max(ymax,item.ymax)
      zmax = max(zmax,item.zmax)
      xmin = min(xmin,item.xmin)
      ymin = min(ymin,item.ymin)
      zmin = min(zmin,item.zmin)

    surface_list += str(int(item.media_id)) + ' ' # append so we get a list of the senses / surfaces for our media cards

  # make planes for the split - 3 zones means 4 divisions so we do +1
  xplane_locs = np.linspace(xmin, xmax, ui_num_x_divs+1)
  yplane_locs = np.linspace(ymin, ymax, ui_num_y_divs+1)
  zplane_locs = np.linspace(zmin, zmax, ui_num_z_divs+1)

  plane_idx = int(plane_id_start)
  plane_input = {}
  xplane_list = {}
  yplane_list = {}
  zplane_list = {}

  # Ax + By + Cz + con = 0 ---> x = -con; A = 1
  for idx, item in enumerate(xplane_locs):
    xplane_list[plane_idx] = 'plane ' + str(plane_idx) + '  xpl=1.0 ypl=0.0 zpl=0.0 con=' + str(-item)
    plane_idx += 1
  for idx, item in enumerate(yplane_locs):
    yplane_list[plane_idx] = 'plane ' + str(plane_idx) + '  xpl=0.0 ypl=1.0 zpl=0.0 con=' + str(-item)
    plane_idx += 1
  for idx, item in enumerate(zplane_locs):
    zplane_list[plane_idx] = 'plane ' + str(plane_idx) + '  xpl=0.0 ypl=0.0 zpl=1.0 con=' + str(-item)
    plane_idx += 1

  # make a list of all of our planes.
  plane_input.update(xplane_list)
  plane_input.update(yplane_list)
  plane_input.update(zplane_list)

  # media_id_offset - how much to offset from the original media id. keep at 0 if media is 0 already
  media = int(medialine.split()[1])
  importance = int(medialine.split()[2])
  medialist = {}
  media_id_list = np.array([])

  for (zkey1, zp), (zkey2, zn) in zip(list(zplane_list.items())[:-1], list(zplane_list.items())[1:]):
    zplane_id_posi = zp.split()[1]
    zplane_id_nega = zn.split()[1]
    for (zkey1, yp), (zkey2, yn) in  zip(list(yplane_list.items())[:-1], list(yplane_list.items())[1:]):
      yplane_id_posi = yp.split()[1]
      yplane_id_nega = yn.split()[1]
      for (zkey1, xp), (zkey2, xn) in  zip(list(xplane_list.items())[:-1], list(xplane_list.items())[1:]):
        xplane_id_posi = xp.split()[1]
        xplane_id_nega = xn.split()[1]
        if media == 0:
          media_input = 0
        else:
          media_input = media_id_offset
        medialist[media_id_offset] = 'media ' + str(media_input) + ' ' + str(importance) + ' '  + surface_list +  ' +'+str(xplane_id_posi) + ' -'+str(xplane_id_nega) +  ' +'+str(yplane_id_posi) + ' -'+str(yplane_id_nega) +  ' +'+str(zplane_id_posi) + ' -'+str(zplane_id_nega)
        media_id_list = np.append(media_id_list, media_input)
        media_id_offset += 1


  media_id_offset += 1
  return plane_input, medialist, media_id_list


def make_and_write_media_id_list(material, media_id_list, outputfilename):
  output = {}
  for label in media_id_list:
    line = ''
    for idx, iso in enumerate(material.isotope_list):
      line += '\n' + iso + ' ' + str(int(label)) + ' ' + str(material.biasid) + ' ' + str(material.atom_dens[idx]) + ' ' + str(material.temp) + ' end'
      output[label] = line


  with open(outputfilename, "w") as file:
    if media_id_list[0] != 0:
      for this in output.keys():
        file.write(output[this])
    else:
      file.write("Media was found to be zero aka void - nothing to write.")

  return output

def delete_old_media_and_replace_with_new_then_write(media_key, unitbank, media_input, surface_input, filename):
  # replace media key with header
  old_medialine = ''
  for word in unitbank[media_key]:
    old_medialine+= word + ' '

  newfill = "' =====================================================================\n' THIS SECTION GENERATED BY SCALESLICE\n' =====================================================================\n" + "'OLD LINE = "+old_medialine
  for this in surface_input.keys():
    item = surface_input[this]
    newfill += '\n'+item
  newfill += '\n'
  newfill += "\n'NOW TIME FOR MEDIA CARDS"
  newfill += '\n'
  for this in media_input.keys():
    item = media_input[this]
    newfill += '\n'+item

  newfill+= '\n'+"' =====================================================================\n' END OF SCALESLICE GENERATED SECTION\n' =====================================================================\n"

  unitbank[media_key] = newfill


  with open(filename, "w") as file:
    for this in unitbank.keys():
      if isinstance(unitbank[this], list):
        to_write = ' '.join(unitbank[this]) + '\n'
      else:
        to_write = unitbank[this] + '\n'
      file.write(to_write)


  return unitbank

def parse_input_file(file_path):
  parameters = {}
  with open(file_path, 'r') as file:
    current_label = None
    for line in file:
      if line.strip().startswith("[") and line.strip().endswith("]"):
        current_label = line.strip()[1:-1]
        parameters[current_label] = {}
      elif line.strip() and not line.strip().startswith("#"):
        key, value = re.findall(r'\b(\w+)\s*=\s*(.+)\b', line)[0]
        parameters[current_label][key] = value.strip("'")
  return parameters

def process_parameters(parameters):
  print(f"\nProcessing parameters....\n")
  print("Parameters:")
  for box in parameters.keys():
    if len(box) > 0:
      print("========================================")
      print("Parameter object: ", box.split(']')[0])
      for item in parameters[box].keys():
        print(item, " = ", parameters[box][item])
  print("========================================")
  print("Calling your code with these parameters...")
  print()

def return_list_from_box(parameters, box):
  para = parameters[box]
  splittype = para['type']
  filename = para['filename']
  medialine_ui =para['media_line']
  unit_int = int(para['unit_number'])
  nx = int(para['nx'])
  ny = int(para['ny'])
  nz = int(para['nz'])
  surfaces_offset = int(para['surfaces_offset'])
  media_offset = int(para['media_offset'])
  do_vols = para['do_vols']
  geometry_output_name = filename[0:-4] + '_COMPOSITION_.txt'
  unit_output_name = filename[0:-4] + '_UNIT_.txt'
  if filename[-4::] != '.inp':
    raise Exception("File end MUST be .inp. Did you make a spelling mistake?", filename[-4::])
  if 'true' or 'True' in do_vols:
    do_vols = True
  else:
    do_vols = False
  return splittype, filename, medialine_ui, unit_int, nx,ny,nz, surfaces_offset,media_offset, do_vols, geometry_output_name, unit_output_name

def print_intro_garbage(inp):
  n = len(inp)
  print("\n\n========================================")
  print("Total arguments passed:", n)
  print("\nName of Python script:", inp[0])
  if n > 0:
    print("\nName of Input file", inp[1])
  else:
    raise Exception("Need input file to run this script!")

#############################################
# USER INPUT
#############################################
# EXAMPLE INPUT PARAMETERS FOR REST OF SCRIPT
# filename = 'Tank_HC_ActivationDose_Monaco.inp'
# medialine_ui = 'media 50000 1 82 -50 -401 -411 -81'
# unit_int = int(1000)
# nx = 20
# ny = 20
# nz = 10
# surfaces_offset = 100001
# media_offset = 200001
# do_vols = True
# geometry_output_name = filename[0:-4] + '_COMPOSITION_.txt'
# unit_output_name = filename[0:-4] + '_UNIT_.txt'
#############################################

#############################################
# RUNTIME
#############################################
# steps to run:
# 1. open file
# 2. get material database
# 3. get a dict of the unit line by line where dict = {0: 'line1here', 1: 'line2hre'} ...
# 4. get media line / card of the media input from the user.
# 5. grab media number.
# 6. get all surfaces from the media card -
# 7. get all surface cards from the unit dictionary.
# 8. make objects for each surface card (currently cylinders and cuboids)
# 9. run cartesian grid split.
# 10. delete old media and replace it with the new media and surfaces
# 11. write to file
#############################################


print_intro_garbage(sys.argv)
parameters = parse_input_file(sys.argv[1])
process_parameters(parameters)
mybox = list(parameters.keys())[0]

splittype, filename, medialine_ui, unit_int, nx,ny,nz, surfaces_offset,media_offset, do_vols, geometry_output_name, unit_output_name = return_list_from_box(parameters, mybox)

file = open(filename)
material_database = get_comps(file=file)
file = open(filename)

# get a dict full of all the lines from the user defined unit
unitbank = get_unit(file, unit_int)

# gets the line and key of unitbank dict for the user specified media card that we are replacing
medialine, media_key = get_media_index(unitbank, medialine_ui)
unitbank[media_key] # this returns the medialine line from the unit bank - media key is the key from unitbank the media line occurs on.
media_number = medialine.split()[1]

# gets the surfaces from the media line that we are using - labels inside, outside, and all surfaces. also gets importances.
inside, outside, all_surfaces, importance, medias_mat_id = surfaces_from_media_line(medialine)

# gets the individual lines defining the surfaces that we are interested in.
surface_cards = get_surface_cards(all_surfaces, unitbank) # list of surfaces that appear in the media card -- both positive and negative values here.

# make a collection of surface objects for every surface in the relevant media card.
surface_object_collection = objects_from_surface_cards(surface_cards, all_surfaces)




############################################################################
# THIS IS WHERE WE STORE SLICING METHODS
# slicing method should return: media_id_list - list of media id's we need to make new compositions for.
# surface input - dict of surfaces for input into the new input file
# media input dict of media card for input into the new input file

# now we need to use the old media line, surface cards, and then make new surfaces for xyz planes
# use maxes and mins from surface cards to determine boundaries of the split. Only use surface cards with negative ids to decide mins and maxes

if splittype == 'CartesianGridSplit':
  surface_input, media_input, media_id_list = CartesianGridSplit(surface_object_collection,
                                                                  medialine, nx, ny, nz, surfaces_offset, media_offset)

elif splittype = 'CylinderRZSplit':
  surface_input, media_inpt, media_id_list = CylinderRZSplit(surface_object_collection,
                                                                  medialine, nx, ny, nz, surfaces_offset, media_offset)
else:
  raise Exception("Split type unknown - current valid split types are: CartesianGridSplit\nDid you make a spelling or case mistake in the 'type' parameter?")



############################################################################

mat_dummy = material_normal
mat_dummy.mat_id = 0
mat_dummy.isotope_list = np.array([])
mat_dummy.biasid = int(-1)
mat_dummy.atom_dens = []
mat_dummy.temp = -1.0
material_database.append_mat_to_lib(mat_dummy)
media_replacement = make_and_write_media_id_list(material_database.material_dict[int(media_number)], media_id_list, geometry_output_name)

# ok now delete unitbank[mediakey] and replace with new medias and surfaces.
unitbank_updated = delete_old_media_and_replace_with_new_then_write(media_key, unitbank, media_input, surface_input, unit_output_name)
RESET = '\033[0m'
RED = '\033[31m'
BRIGHT_CYAN = '\033[96m'

print("========================================")
print(RED + "Finished" + RESET)
print("Outputs are in: ")
print("Composition: " +BRIGHT_CYAN + geometry_output_name + RESET)
print("Replacement for input unit: " + BRIGHT_CYAN + unit_output_name + RESET)
