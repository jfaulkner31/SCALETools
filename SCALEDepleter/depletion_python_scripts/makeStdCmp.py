import subprocess
from getComps import material_lib
import shutil
import os


def makeStdCmpFromF71(materialNumber: int, temperature: float, filename: str, dictionaryFilename: str, outputFilename: str, tmpdir: str, scale_version: int):
  """
  Makes a standard cmp file using OBIWAN script.
  materialNumber = material id
  temperature = material temperature
  filename = name of f71 file / path to f71 file e.g. tmp/thisfile.f71
  dictionaryFilename = dictionary of all the nuclides to include in the stdcmp input
  outputFilename = output name of the created stdcmp library
  scale_version = scale version number from original input file.
  """

  # Get the path to the f71_to_comp.sh script
  this_dir = os.path.dirname(__file__)
  bash_f71_path = this_dir + "/../depletion_bash_scripts/f71_to_comp.sh"
  input_filename_path = os.getcwd()+'/'+filename # INPUTDIR/OrigenResults_F71dir/PREDICTOR_EOS_STEPNNN_MATMMM.f71
  dict_filename_path = os.getcwd()+'/'+dictionaryFilename

  # scaleDepleterPath = os.path.dirname(os.path.realpath(__file__))

  print("Now making standard comp file for material", materialNumber, " | output=", outputFilename)
  # rm = subprocess.run(['rm', outputFilename])
  origenrun = subprocess.run(["bash", bash_f71_path, str(materialNumber), str(temperature), input_filename_path, dict_filename_path, outputFilename, str(scale_version)],
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.STDOUT)
  print("Stdcmp file with name", outputFilename, 'was successfully made!')

  # cp = subprocess.run(['cp', outputFilename, tmpdir+'/'+outputFilename])
  # rm = subprocess.run(['rm', outputFilename])

  shutil.move(outputFilename, tmpdir+'/'+outputFilename)
  return

def makeStdCmpFromMatLib(outputFilename: str, material_lib: material_lib, material_index: int, tmpdir: str):
  """
  Makes a stdcmp file from a material library.
  """
  print("Now making standard comp file for material", material_index, " | output=", outputFilename)
  fh = open(outputFilename, 'w')
  this_lib = material_lib.material_dict[material_index]
  fh.write("' material number "+str(material_index)+" in the stdcmp format with filename "+outputFilename+'\n')
  for idx, iso in enumerate(this_lib.isotope_list):
    atmdens = this_lib.atom_dens[idx]
    temp = this_lib.temp
    fh.write(iso+" "+str(material_index)+" 0 "+atmdens+" "+str(temp)+" end\n")
  fh.close()
  cp = subprocess.run(['cp', outputFilename, tmpdir+'/'+outputFilename])
  rm = subprocess.run(['rm', outputFilename])
  # {nuclide} $materialNumber 0 ${atomdensity} $temperature end
  # initial_mats.material_dict[101].isotope_list
  # initial_mats.material_dict[101].atom_dens
  # initial_mats.material_dict[101].temp
  return

def grabNuclideFromF71(filename: str, nuclide: str, precision: int):
  """"""
  # obiwan view -format=csv -units=atom -prec=8 -idform='{:ee}-{:A}{:m}' PREDICTOR_EOS_step0_mat101.f71 | grep "xe-135,"
  if "," not in nuclide:
    nuclide += ','
  if '-' not in nuclide:
    raise Exception("Nuclide not in correct format - needs to be in format such as: xe-135")
  p = subprocess.run(['bash', 'f71_nuclide_dump.sh', filename, str(precision), nuclide],
                      capture_output = True,
                      text = True)
  out = p.stdout
  nd = float(out.split(',')[-1]) # get number density - works for f71 files that have more than 1 pos - e.g. out = '              xe-135, 1.225176631824E-10, 2.081296201948E-10'
  return nd
