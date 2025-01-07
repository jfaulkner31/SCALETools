import os
import shutil
import subprocess
def copy_files_from_temp(tmpdir: str, step_num: int, mcf33dir: str, appendThis=''):
  # get all files in directory
  all_files = [entry.name for entry in os.scandir(tmpdir) if entry.is_file()]
  f33_files = {}

  # if monte carlo results doesnt exist, make a directory
  if mcf33dir not in [entry.name for entry in os.scandir('.') if entry.is_file()]:
    mkdir = subprocess.run(['mkdir', mcf33dir])

  # copy f33 files
  for filename in all_files:
    if ("f33" in filename) & ('mix' in filename):
      index_dot = filename.find('.')
      index_x = filename.find('x')
      number = int(filename[index_x+1 : index_dot])
      new_filename = str(number)+'.f33'

      new_fileloc = mcf33dir+'/'+'step'+str(step_num)+appendThis+'/'+new_filename # location of the new file from the running directory
      if 'step'+str(step_num)+appendThis not in [entry.name for entry in os.scandir(mcf33dir) if entry.is_file()]:
        mkdir = subprocess.run(['mkdir', mcf33dir+'/step'+str(step_num)+appendThis])

      shutil.move(tmpdir+'/'+filename, new_fileloc)
      print(f"Moved: {filename} from temporary directory as {new_filename}")
      f33_files[number] = new_fileloc

  return f33_files







