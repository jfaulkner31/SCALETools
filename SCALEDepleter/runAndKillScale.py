import sys
import numpy as np
import subprocess
import time
import os
import signal

import getComps

def runAndKillScale(scale_input_line, scale_input_file_name):

  dot = scale_input_file_name.find('.')
  msg_file_name = scale_input_file_name[0:dot]+'.msg'

  process = subprocess.Popen(scale_input_line , preexec_fn=os.setsid)

  print(f"Started process with PID: {process.pid}", flush=True)

  print("Now sleeping python for 5 seconds", flush=True)

  time.sleep(5)

  try:
    # Step 2: Monitor the criteria
    start_time = time.time()
    terminate_process = False
    while not terminate_process:
      # print("Checking",msg_file_name,'if step 1 is finished yet....')
      # Example criteria: Stop after 100 seconds
      # if time.time() - start_time > 100:
      #   print("Criteria met. Terminating the process.")
      #   break
      # time.sleep(1)  # Simulate monitoring delay
      # criteria  stop after origen finishes running.
      with open(msg_file_name, 'r') as file:
        lines = file.readlines()
      for line in lines[::-1]:
        if "best estimate" in line:
          keff_line = line
        if "ORIGEN multi-zone depletion global step 1" in line:
          terminate_process = True
          raise_exception = False
          print("Termination of process detected", flush=True)
        if "Output is stored in" in line:
          terminate_process = True
          raise_exception = True
          print("Termination of process detected", flush=True)

      time.sleep(1) # sleep for 1 seconds to


  except KeyboardInterrupt:
    print("Manual interruption received.")

  # Step 3: Kill the process
  time.sleep(6) # sleep for 20 seconds to let things finish writing/copying
  os.killpg(os.getpgid(process.pid), signal.SIGTERM)
  process.wait() # garbage collection for subprocess to finish killing
  print("Process terminated after", time.time()-start_time, "seconds.", flush=True)

  if raise_exception:
    raise Exception("Output terminated - exception raised due to SCALE error.")

  return keff_line
