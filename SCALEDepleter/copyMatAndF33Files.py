import os
import shutil

def copy_files_from_temp(tmpdir: str, fissionable_mats: list, std_cmp_label: str, f33_label: str):
  # working_dir
  working_dir = tmpdir+'/../'
  # get all files in directory
  all_files = [entry.name for entry in os.scandir(tmpdir) if entry.is_file()]
  std_cmp_files = {}
  f33_files = {}

  # first do standard comp files
  for filename in all_files:
    if "StdCmpMix" in filename:
      index_u = filename.find('_')
      index_x = filename.find('x')
      number = int(filename[index_x+1 : index_u])
      new_filename = 'STD_CMP_'+str(number)+'_'+std_cmp_label
      shutil.copy(tmpdir+'/'+filename, new_filename)
      print(f"Copied: {filename} from temporary directory as {new_filename}")
      std_cmp_files[number] = new_filename
    if ("f33" in filename) & ('mix' in filename):
      index_dot = filename.find('.')
      index_x = filename.find('x')
      number = int(filename[index_x+1 : index_dot])
      new_filename = 'F33_'+str(number)+'_'+f33_label+'.f33'
      shutil.copy(tmpdir+'/'+filename, new_filename)
      print(f"Copied: {filename} from temporary directory as {new_filename}")
      f33_files[number] = new_filename
  return f33_files, std_cmp_files







