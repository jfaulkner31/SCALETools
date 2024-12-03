import subprocess
from getComps import material_lib

def makeStdCmpFromF71(materialNumber: int, temperature: float, filename: str, dictionaryFilename: str, outputFilename: str, tmpdir: str, pctag: str):
  """
  Makes a standard cmp file using OBIWAN script.
  materialNumber = material id
  temperature = material temperature
  filename = name of f71 file / path to f71 file e.g. tmp/thisfile.f71
  dictionaryFilename = dictionary of all the nuclides to include in the stdcmp input
  outputFilename = output name of the created stdcmp library
  """
  updatedoutputfilename = outputFilename+pctag

  print("Now making standard comp file for material", materialNumber, " | output=", updatedoutputfilename)
  rm = subprocess.run(['rm', updatedoutputfilename])
  origenrun = subprocess.run(["bash", "f71_to_comp.sh", str(materialNumber), str(temperature), filename, dictionaryFilename, updatedoutputfilename])

  cp = subprocess.run(['cp', updatedoutputfilename, tmpdir+'/'+updatedoutputfilename])
  rm = subprocess.run(['rm', updatedoutputfilename])
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
