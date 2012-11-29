#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,numpy,skfmm
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP
from scipy import ndimage

# -----------------------------
class extendedOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
    
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            lvalue = value.split(",")
            values.ensure_value(dest, []).extend(lvalue)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)



def outStdout(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    print cmd
  else:
    print cmd
  return

def outFile(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    locals['filepointer'].write(cmd+'\n')
  else:
    locals['filepointer'].write(cmd+'\n')
  return


def output(cmds,locals,dest):
  for cmd in cmds:
    if isinstance(cmd,list):
      output(cmd,locals,dest)
    else:
      {\
      'File': outFile,\
      'Stdout': outStdout,\
      }[dest](str(cmd),locals)
  return


# +++++++++++++++++++++++++++++++++++++++++++++++++++
def vtk_writeASCII_mesh(dim,res,origin,data):
# +++++++++++++++++++++++++++++++++++++++++++++++++++
  """ function writes data array defined on a rectilinear grid """
  N  = res[0]*res[1]*res[2]
  
  cmds = [\
          '# vtk DataFile Version 3.1',
          string.replace('powered by $Id: spectral_geomCheck.py 1575 2012-06-26 18:07:38Z MPIE\p.eisenlohr $','\n','\\n'),
          'ASCII',
          'DATASET RECTILINEAR_GRID',
          'DIMENSIONS %i %i %i'%(res[0]+1,res[1]+1,res[2]+1),
          'X_COORDINATES %i float'%(res[0]+1),
          ' '.join(map(str,[i*dim[0]/res[0]+origin[0] for i in range(res[0]+1)])),
          'Y_COORDINATES %i float'%(res[1]+1),
          ' '.join(map(str,[i*dim[1]/res[1]+origin[1] for i in range(res[1]+1)])),
          'Z_COORDINATES %i float'%(res[2]+1),
          ' '.join(map(str,[i*dim[2]/res[2]+origin[2] for i in range(res[2]+1)])),
          'CELL_DATA %i'%N,
         ]
  
  for datatype in data:
    for item in data[datatype]:
      cmds += [\
               '%s %s float'%(datatype.upper()+{True:'',False:'S'}[datatype.lower().endswith('s')],item),
               'LOOKUP_TABLE default',
               [[['\t'.join(map(str,data[datatype][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])]
              ]

  return cmds


# ----------------------- MAIN -------------------------------

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
        'origin':     ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
        'origin':     lambda x: float(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog [geomfile[s]]', description = """
Produce Euclidean distance map from geom description

""" + string.replace('$Id: spectral_geomCheck.py 1575 2012-06-26 18:07:38Z MPIE\p.eisenlohr $','\n','\\n')
)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name)})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

  info = {'resolution': [0,0,0],
          'dimension':  [0.0,0.0,0.0],
          'origin':     [0.0,0.0,0.0],
         }
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] in identifiers.keys():
      for i in xrange(len(identifiers[headitems[0]])):
        info[headitems[0]][i] = mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])

  if info['resolution'] == [0,0,0]:
    print 'no resolution info found.'
    sys.exit(1)
  if info['dimension'] == [0.0,0.0,0.0]:
    print 'no dimension info found.'
    sys.exit(1)

  print 'resolution: %s'%(' x '.join(map(str,info['resolution'])))
  print 'dimension:  %s'%(' x '.join(map(str,info['dimension'])))
  print 'origin:     %s'%(' : '.join(map(str,info['origin'])))

  dx = info['dimension'][0]/info['resolution'][0]
  
  data = {'scalar':{'perimeter':numpy.zeros(info['resolution'],'i'),
              		'distance':numpy.zeros(info['resolution'],'i')}}
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      data['scalar']['perimeter'][i%info['resolution'][0],(i/info['resolution'][0])%info['resolution'][1],i/info['resolution'][0]/info['resolution'][1]] = item
      i += 1

#  data['scalar']['perimeter'] = numpy.where(ndimage.morphology.grey_dilation(data['scalar']['perimeter'],size=(3,3,3))-data['scalar']['perimeter']>0,0,1)
  
  FDstencil_x = numpy.zeros([3,3,3])
  FDstencil_x[:,1,1] = [-1,0,1]
  FDstencil_y = numpy.zeros([3,3,3])
  FDstencil_y[1,:,1] = [-1,0,1]
  FDstencil_z = numpy.zeros([3,3,3])
  FDstencil_z[1,1,:] = [-1,0,1]
  data['scalar']['perimeter'] = numpy.where(numpy.abs(ndimage.convolve(data['scalar']['perimeter'], FDstencil_x)) + numpy.abs(ndimage.convolve(data['scalar']['perimeter'], FDstencil_y)) + numpy.abs(ndimage.convolve(data['scalar']['perimeter'], FDstencil_z))>0,0,1)
  data['scalar']['distance'] = skfmm.distance(data['scalar']['perimeter'], dx=dx)
  
  out = {}
  out['mesh'] = vtk_writeASCII_mesh(info['dimension'],info['resolution'],info['origin'],data)
  
  for what in out.keys():
    if file['name'] == 'STDIN':
      output(out[what],{},'Stdout')
    else:
      (head,tail) = os.path.split(file['name'])
      vtk = open(os.path.join(head,what+'_'+os.path.splitext(tail)[0]+'.vtk'), 'w')
      output(out[what],{'filepointer':vtk},'File')
      vtk.close()
