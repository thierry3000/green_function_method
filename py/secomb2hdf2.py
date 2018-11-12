# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 17:18:15 2015

@author:  Thierry Fredrich

@brief:   This file was used to convert secombs .dat file to a .h5 file used by II
          
@details: Welter M, Fredrich T, Rinneberg H, Rieger H (2016) 
          Computational Model for Tumor Oxygenation Applied to Clinical Data 
          on Breast Tumor Hemoglobin Concentrations Suggests Vascular Dilatation and Compression.
          PLOS ONE 11(8): e0161267. 
          https://doi.org/10.1371/journal.pone.0161267
          (I)
          Fredrich, T., Welter, M. & Rieger, H. 
          Tumorcode
          Eur. Phys. J. E (2018) 41: 55. 
          https://doi.org/10.1140/epje/i2018-11659-x
          (II)
          Secomb, T.W., Hsu, R., Park, E.Y.H. and Dewhirst, M.W. 
          Green's function methods for analysis of oxygen delivery 
          to tissue by microvascular networks. 
          Annals of Biomedical Engineering, 32: 1519-1529 (2004)
          (III)
          
          At time of publication (I), the software in this repository was already present (III). 
          In (I) the oxygen distribution within artificial blood vessel networks was calculated. 
          At time of publication the software in this repository was available and 
          we compared run time and accuracy to that. 
          The present file was used to transform the topological network data provided by secomb 
          to a .h5 file used by tumorcode (II).

"""

import krebsutils as ku
import h5py
import numpy as np
import sys
import os

def correct_bc_type_from_secomb_to_MW(anArray):
  for i,v in enumerate(anArray):
    if v == 0:
      anArray[i] = 1
  return anArray

def get_network_from_file():
  with open('Network.dat','r') as network_file:
    data = network_file.readlines()
  nice_name = str(data[0])
  #edges start at 8
  # need to find end
  for (i,line) in enumerate(data):
    if 'total number of segments' in line:
      beginOfEdges = i
    if 'total number of nodes' in line or '\ttotal\tnumber\tof\tnodes' in line:
      endOfEdges = i
      
  cropped_edges = data[beginOfEdges+2:endOfEdges]

  label=[]
  mw_vessel_flag= []
  node_a_index=[]
  node_b_index=[]
  radii = []
  length = []
  hema = []
  flow = []
  
  #to be calculated
  flow = []

  for line in cropped_edges:
    vesselId = line.split()[0]
    vesselType = line.split()[1]
    vesselA = line.split()[2]
    vesselB = line.split()[3]
    vesselDiameter = line.split()[4]
    vesselFlow = line.split()[5]
    vesselHema = line.split()[6]
    
    label.append(int(vesselId)-1)
    mw_vessel_flag.append(ku.CAPILLARY)
    node_a_index.append(int(vesselA)-1)
    node_b_index.append(int(vesselB)-1)
    radii.append(float(vesselDiameter)/2)
    flow.append(np.fabs(float(vesselFlow))*1e6/60.)
    hema.append(float(vesselHema))

  #edge stuff
  N_edges= len(cropped_edges)
  ''' there are double entries, which maybe spoil the computation '''
  
  edgegrp.attrs.create('COUNT',N_edges)
  ds_nodeA = edgegrp.create_dataset('node_a_index', data= node_a_index)
  ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
  ds_radius = edgegrp.create_dataset('radius', data=radii)
  ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
  ds_flow = edgegrp.create_dataset('flow', data=flow)
  ds_flags = edgegrp.create_dataset('flags', data= mw_vessel_flag)
  ds_vessel_label = edgegrp.create_dataset('vessel_label', data= label)
  return nice_name


def get_roots_from_file():
  with open('Network.dat','r') as network_file:
    data = network_file.readlines()
  
  for (i,line) in enumerate(data):
    if 'total number of boundary nodes' in line or 'total\tnumber\tof\tboundary\tnodes' in line:
      beginOfRoots = i
      
  cropped_roots = data[beginOfRoots+2:] 
  
  indeces_of_roots = []
  bctyp_of_roots = []
  value_of_bc = []
  roots_hema = []
  roots_PO2 = []
  for line in cropped_roots:
    indeces_of_roots.append(int(line.split()[0])-1)
    type_secomb = int(line.split()[1])
    bctyp_of_roots.append(type_secomb)
    value = float(line.split()[2])
    if(type_secomb == 2):
      value_of_bc.append(-1*value*1e6/60.)
    if(type_secomb == 0):
      value_of_bc.append(np.fabs(value)*0.13)
  print(value_of_bc)
    
  ds_value_of_bc = nodegrp.create_dataset('bc_value', data = value_of_bc)        
  ds_bctyp_of_roots = nodegrp.create_dataset('bc_type', data= correct_bc_type_from_secomb_to_MW(bctyp_of_roots))
  ds_bc_node_index = nodegrp.create_dataset('bc_node_index', data = indeces_of_roots)
  ds_bc_conductivity_value = nodegrp.create_dataset('bc_conductivity_value', data=np.zeros_like(value_of_bc)) 
  ds_roots = nodegrp.create_dataset('roots', data = indeces_of_roots)
        
def get_nodes_from_file():
  with open('Network.dat','r') as network_file:
    data = network_file.readlines()
  nice_name = str(data[0])
  # need to find end
  for (i,line) in enumerate(data):
    if 'total number of nodes' in line or '\ttotal\tnumber\tof\tnodes' in line:
      beginOfNodes = i
    if 'total number of boundary nodes' in line or 'total\tnumber\tof\tboundary\tnodes' in line:
      endOfNodes = i
      
  
    
  cropped_nodes = data[beginOfNodes+2:endOfNodes]


  node_index = []
  #node_label2index = dict()
  #positions_of_nodes = np.zeros([ 972, 3])
  positions_of_nodes = []
  max_x=0
  max_y=0
  max_z=0
    
  for line in cropped_nodes:
    node_index.append(int(line.split()[0])-1)
    x_pos =float(line.split()[1])
    y_pos =float(line.split()[2])
    z_pos =float(line.split()[3])
    if x_pos>max_x:
      max_x=x_pos
    if y_pos>max_y:
      max_y=y_pos
    if z_pos>max_z:
      max_z=z_pos
    coordinate=[x_pos,y_pos,z_pos]
    positions_of_nodes.append(coordinate)


  #node stuff
  N_nodes = len(node_index)
  nodegrp.attrs.create('COUNT', N_nodes)
  ds_world_pos = nodegrp.create_dataset('world_pos', data = np.array(positions_of_nodes))

def correct_vessel_indeces(node_label2index, vesselgrp):
    va = vesselgrp['edges/node_a_index']
    vb = vesselgrp['edges/node_b_index']
    del vesselgrp['edges/node_a_index']
    del vesselgrp['edges/node_b_index']
    va_new = []
    vb_new = []
    
    for old_index in va:
        va_new.append(node_label2index[old_index])
    for old_index in vb:
        vb_new.append(node_label2index[old_index])
        
    vesselgrp.create_dataset('edges/node_a_index', data= va_new)
    vesselgrp.create_dataset('edges/node_b_index', data= vb_new)
    
    print("done... renumbering")

''' I found that there are some edges twice, 
    delete them
    '''
def delete_double_edges(vesselgrp):
  #176 deleted for convergenz issues
  by_hand_identified_doubles = [279,283,284,294,312,325,328]
  va = np.asarray(vesselgrp['edges/node_a_index'])
  vb = np.asarray(vesselgrp['edges/node_b_index'])
  flags = np.asarray(vesselgrp['edges/flags'])
  radii = np.asarray(vesselgrp['edges/radius'])
  hema = np.asarray(vesselgrp['edges/hematocrit'])
  no_edges = len(va)
  good_indeces = np.ones(no_edges)
  for i in by_hand_identified_doubles:
    good_indeces[i] = 0
  good_indeces = good_indeces.astype(bool)
  del vesselgrp['edges/node_a_index']
  del vesselgrp['edges/node_b_index']
  del vesselgrp['edges/flags']
  del vesselgrp['edges/radius']
  del vesselgrp['edges/hematocrit']
#  matrix = np.asarray([va,vb])
#  unique_matrix = np.unique(matrix, axis=1)
#  va_new = unique_matrix[0,:]
#  vb_new = unique_matrix[1,:]
  va_new = va[good_indeces]
  vb_new = vb[good_indeces]
  vesselgrp.create_dataset('edges/node_a_index', data= va_new)
  vesselgrp.create_dataset('edges/node_b_index', data= vb_new)
  N_edges= len(va_new)
  edgegrp = vesselgrp['edges']
  edgegrp.attrs['COUNT']= N_edges
  flags_new = flags[good_indeces]
  vesselgrp.create_dataset('edges/flags', data= flags_new)
  radii_new = radii[good_indeces]
  vesselgrp.create_dataset('edges/radius', data= radii_new)
  hema_new = hema[good_indeces]
  vesselgrp.create_dataset('edges/hematocrit', data= hema_new)
  
  print("done... deleting double entries")

if __name__ == '__main__':
  if not os.path.isfile("Network.dat"):
    print("Network.dat not found in current folder");
    sys.exit(1);
  else:
    print("reading data from Network.dat");
  
  
  #write the data to a h5 file
  fn = 'output.h5'
  with h5py.File(fn,'w') as f3:
    vesselgrp = f3.create_group('vessels')
    vesselgrp.attrs['CLASS']= 'REALWORLD'
    nodegrp = f3['vessels'].create_group("nodes")
    edgegrp = vesselgrp.create_group("edges")
    
    get_network_from_file()
    get_nodes_from_file()
    get_roots_from_file()
    #delete_double_edges(vesselgrp)
    #correct_vessel_indeces(node_label2index, vesselgrp)
    

  '''***** tumorcode stuff *****'''
  
  ''' see 
      py/krebsjobs/parameters/
  '''
  calcflow_param_dict_for_tumorcode = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect = 1,)
  ##CALCULATE!!!!
  dd = ku.calc_vessel_hydrodynamics_Ccode(fn, 'vessels', True, calcflow_param_dict_for_tumorcode, False, True)
  #print(dd)