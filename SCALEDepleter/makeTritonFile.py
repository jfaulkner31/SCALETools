def makeTritonFile(baseFilename, fissionable_mats):
  """
  Makes a triton file using standard composition files instead of listing nuclides out line by line.
  """
  file = open(baseFilename, 'r')
  newFilename = baseFilename.split('.')[0] + '_stdcmp.inp'
  newfile = open(newFilename, 'w')

  do_comp = False
  for wholeline in file.readlines():
    line = wholeline.split() # -> line = ['this', 'is', 'how', 'line', 'split', 'looks']
    additional_lines = ''

    # starts and stops reading of compositions
    if 'read' in line and 'comp' in line:
      if line.index('read') == line.index('comp') - 1:
        print("read comp found! Now starting comp directory.")
        do_comp = True # -> activates composition mode
        for mat in fissionable_mats:
          additional_lines += '<StdCmpMix'+str(mat)+'\n'

    if 'end' in line and 'comp' in line:
      if line.index('end') == line.index('comp') - 1:
        print("end comp found! Now finalizing comp directory.")
        do_comp = False

    # now check each line...
    anyline = any(str(x) in line for x in fissionable_mats)
    if (do_comp) & ("end" in line) & (anyline): # if this is a line with a depletable material - do nothing
      continue
    elif len(line) == 0:
      continue
    elif line[0][0] == "'":
      continue # commented line
    else:
      newfile.write(wholeline+additional_lines)

  newfile.close()
  file.close()
  return newFilename
