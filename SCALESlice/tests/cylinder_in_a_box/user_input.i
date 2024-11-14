[media_label_here]]
  type = CartesianGridSplit
  filename = 'tests/cylinder_in_a_box/cylinder_in_a_box.inp'
  media_line = 'media 0 1 99999 -4'
  unit_number = 1
  nx = 20
  ny = 20
  nz = 10
  surfaces_offset = 100001
  media_offset = 200001
  do_vols = True
  composition_suffix = '_COMPOSITION_.txt'
  unit_suffix = '_UNIT_.txt'
[]

