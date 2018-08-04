# -*- coding: UTF-8 no BOM -*-

from .solver import Solver
import damask
import subprocess,re

class Abaqus(Solver):

  def __init__(self,version=''):                                                                    # example version string: 2017
    self.solver='Abaqus'
    if version =='':
      version = damask.Environment().options['ABAQUS_VERSION']
    else:
      self.version = version

  def return_run_command(self,model):
    env=damask.Environment()
    shortVersion = re.sub('[\.,-]', '',self.version)
    try:
      cmd='abq'+shortVersion
      subprocess.check_output(['abq'+shortVersion,'information=release'])
    except OSError:                                                                                 # link to abqXXX not existing
      cmd='abaqus'
      process = subprocess.Popen(['abaqus','information=release'],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      detectedVersion = process.stdout.readlines()[1].split()[1]
      if self.version != detectedVersion:
        raise Exception('found Abaqus version %s, but requested %s'%(detectedVersion,self.version))
    return '%s -job %s -user %s/src/DAMASK_abaqus interactive'%(cmd,model,env.rootDir())
