import subprocess
import shutil
import time
def removeAndMakeDir(dirct: str, make=True):
  """
  Completely deletes given directory and remakes it as an empty folder if make=True.
  Removes but does not make if make=False
  """
  try:
    rm = shutil.rmtree(dirct)
  except:
    if make:
      mkdir = subprocess.run(["mkdir", dirct])
  time.sleep(0.3)
