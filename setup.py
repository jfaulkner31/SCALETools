from setuptools import setup, find_packages

setup(
  name='SCALEDepleter',
  version='0.1',
  packages=find_packages(),
  install_requires=["numpy>=1.26.4", # add any dependencies here
                    "pandas>=2.2.2",
                    "mpi4py>=4.0.0",
                    "matplotlib>=3.10.0"],
  package_data={
    'SCALEDepleter': ['addnuxDicts/*', 'depletion_bash_scripts/*', 'origen_base_files/*']
  }
)
