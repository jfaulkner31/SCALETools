Geometry slicer for scale.
Can slice geometry rectangles/boxes, cylinders, and more!
# SCALESlice
A geometry subdivider for SCALE geometries


Steps to use:

In python:
python slicer.py tests/cylinder_in_a_box/user_input.i

Then look at new .txt outputs.
Replace geometry unit with the newly generated unit.
Then add the materials block into your composition block.
We do not regenerate the materials block due to potential innaccuries with interpreting many different scale material types.
Verify that the materials in the replicated grid are the same definition as what you would expect them to be.

See different slicing options below with examples:
1. CartesianGridSplit - splits a material into a cartesian mesh
[media_label_here]
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
