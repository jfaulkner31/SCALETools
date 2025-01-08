"""
This folder makes data structures for outputting CELI results in python.
"""
class ResultsCELI:
  def __init__(self):
    self.steps = [] # list of num_steps
    self.keff_all = {} # keff by steps and then by another dict by substeps
    self.isotopics_all = {} # isotopics by steps and also by substeps
  def append_keff_line(self, keff_line, step_num, substep):
    try:
      self.keff_all[step_num]
      # now append substep
    except:
      self.keff

