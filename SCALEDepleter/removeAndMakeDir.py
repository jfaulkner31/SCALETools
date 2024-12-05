import subprocess
import shutil
import time
def removeAndMakeDir(dirct: str):
  """
  Completely deletes tmp dir and remakes it as an empty folder.
  """
  try:
    rm = shutil.rmtree(dirct)
  except:
    mkdir = subprocess.run(["mkdir", dirct])
  time.sleep(0.1)
