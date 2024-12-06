
##################################################################
########################### USER INPUT ###########################
##################################################################
"""
  Step 0 - ask yourself:
  Is first step length ok in the input file?
  Is the FMA parameter set to yes?
  Are the NPG and other neutron parameters correct?
  Is addnux set to 0 in the triton base file?

  To run via python command line: python ...py machinefile tmpdir numProcsTransport
"""

# fissionable regions - used for origen later
fissionable_mats = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116] # list of fissionable materials we are depleting -
fissionable_mats_vols = [6.02834870915574000000E+04]*16 # list of vols for each fissionable material
residual_number_density = 1e-20 # residual number density for trace nuclides in initial conditions

# ORIGEN information
include_non_fission_material_power = True # include power effects from non-fissionable materials?
print_transport_powers = True  # print powers after every transport step?
system_IHM_mass_grams = 7.213356e+04 # heavy metal mass in grams i the initial system - the ENTIRE system.
specific_power = [200, 200, 200, 200, 200, 200] # MW/TIHM
steplength_days = [5, 5, 5, 5, 5, 5] # length of each step in days
origen_predictor_divs = 100 # number of time divs for predictor - use for CEBM model only
origen_LI_divs = 10 # how many times during origen calculation the f33 is divided - value of 10 means we use 10 interpolated f33 files per material - use for CELI scheme
origen_steps_per_div = 10 # how many origen substeps are used for each division of f33 files

# File Handles - keep as is - only consider changing the addnux dictionary
addnuxdictbase = 'addnuxDicts/addnux0Dict.dict'
base_triton = 'triton_base.inp'
origen_base = 'baseOrigenFile.inp'
origenResults_F71dir = 'OrigenResults_F71dir'
MonteCarloResults_F33dir = 'MonteCarloResults_F33'


##################################################################
########################### EXECUTION ############################
##################################################################

# command line args
import sys
machinefile = sys.argv[1]
tmpdir = sys.argv[2]
Nprocs = int(sys.argv[3])
# machinefile = 'asdasdasd'
# tmpdir = 'tmp'
# Nprocs = 1


if Nprocs <= 0:
  raise Exception("Nprocs for transport cannot be less than or equal to zero!")
elif Nprocs == 1:
  is_parallel = False
else:
  is_parallel = True

import CELI
CELI.CELI(fissionable_mats=fissionable_mats,
          fissionable_mats_vols=fissionable_mats_vols,
          residual_number_density=residual_number_density,
          include_non_fission_material_power=include_non_fission_material_power,
          print_transport_powers=print_transport_powers,
          system_IHM_mass_grams=system_IHM_mass_grams,
          specific_power=specific_power,
          steplength_days=steplength_days,
          origen_predictor_divs=origen_predictor_divs,
          addnuxdictbase=addnuxdictbase,
          base_triton=base_triton,
          origen_base=origen_base,
          origenResults_F71dir=origenResults_F71dir,
          MonteCarloResults_F33dir=MonteCarloResults_F33dir,
          Nprocs=Nprocs,
          machinefile=machinefile,
          tmpdir=tmpdir,
          is_parallel=is_parallel,
          origen_LI_divs=origen_LI_divs,
          origen_steps_per_div=origen_steps_per_div)
