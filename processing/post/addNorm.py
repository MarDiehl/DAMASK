#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,math,string,damask
from optparse import OptionParser, Option

# -----------------------------
class extendableOption(Option):
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



# definition of element-wise p-norms for matrices

# p = 1
def normAbs(object):
  return sum(map(abs, object))

# p = 2
def normFrobenius(object):
  return math.sqrt(sum([x*x for x in object]))

# p = infinity
def normMax(object):
  return max(map(abs, object))


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing norm of requested column(s) being either vectors or tensors.

""" + string.replace('$Id$','\n','\\n')
)

normChoices = ['abs','frobenius','max']

parser.add_option('-n','--norm',        dest='norm', action='store', type='choice', choices=normChoices, \
                                        help='type of element-wise p-norm (%s) [2]'%(','.join(map(str,normChoices))))
parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', \
                                        help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', \
                                        help='heading of columns containing tensor field values')
parser.add_option('-s','--special',     dest='special', action='extend', type='string', \
                                        help='heading of columns containing field values of special dimension')
parser.add_option('-d','--dimension',   dest='N', action='store', type='int', \
                                        help='dimension of special field values [%default]')

parser.set_defaults(norm = 'frobenius')
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(special = [])
parser.set_defaults(N = 12)

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) + len(options.special)== 0:
  parser.error('no data column specified...')

datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
             'special':    {'len':options.N,
                            'label':[]},
           }


if options.vector  != None:    datainfo['vector']['label']  += options.vector
if options.tensor  != None:    datainfo['tensor']['label']  += options.tensor
if options.special != None:    datainfo['special']['label'] += options.special


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

# --------------- figure out columns to process
  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        sys.stderr.write('column %s not found...\n'%key)
      else:
        if datatype not in active: active[datatype] = []
        if datatype not in column: column[datatype] = {}
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)                   # remember columns of requested data
        table.labels_append('norm%s(%s)'%(options.norm.capitalize(),label)) # extend ASCII header with new labels

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table
    
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested norms
        eval("table.data_append(norm%s(map(float,table.data[column[datatype][label]:column[datatype][label]+datainfo[datatype]['len']])))"%options.norm.capitalize())
    
    table.data_write()                                                      # output processed line

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
