#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 13:44:32 2018

@author: Thierry Fredrich

@brief: a quick script to recreate a vessel file from previous tumorcode versions to current standards.
"""

import h5py
import krebsutils as ku

old_file_name = 'vessels-q2d-8mm-P10-typeE-7x2L130-sample01.h5'
new_file_name = 'recreated_' + old_file_name

with h5py.File(old_file_name, 'r') as f:
  print f.keys()
  with h5py.File(new_file_name, 'w') as f_out:
    f_out.create_group('vessels')
    f_out['vessels'].attrs['CLASS'] = 'GRAPH'
    f_out['vessels'].create_group('edges')
    f_out['vessels/edges'].attrs['COUNT'] = f['vessels/edges'].attrs['COUNT']
    f_out['vessels/edges'].create_dataset('flags', data=f['vessels/edges/flags'])
    f_out['vessels/edges'].create_dataset('node_a_index', data=f['vessels/edges/node_a_index'])
    f_out['vessels/edges'].create_dataset('node_b_index', data=f['vessels/edges/node_b_index'])
    f_out['vessels/edges'].create_dataset('radius', data=f['vessels/edges/radius'])
    
    f_out['vessels'].create_group('lattice')
    print(f['vessels/lattice'].attrs.keys())
    for aKey in f['vessels/lattice'].attrs.keys():
      f_out['vessels/lattice'].attrs[aKey] = f['vessels/lattice'].attrs[aKey]
      
    f_out['vessels'].create_group('nodes')
    f_out['vessels/nodes'].attrs['COUNT'] = f['vessels/nodes'].attrs['COUNT']
    f_out['vessels/nodes'].create_dataset('lattice_pos', data=f['vessels/nodes/lattice_pos'])
    f_out['vessels/nodes'].create_dataset('roots', data=f['vessels/nodes/roots'])
    f_out['vessels/nodes'].create_dataset('pressure', data=f['vessels/nodes/pressure'])
    
    f_out['vessels/nodes'].create_dataset('bc_value', data = [3.59805,5.23394])
    f_out['vessels/nodes'].create_dataset('bc_type', data = [1,1])
    f_out['vessels/nodes'].create_dataset('bc_node_index', data = [0,1])
    f_out['vessels/nodes'].create_dataset('bc_conductivity_value', data = [0,0])

calcflow_param_dict_for_tumorcode = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = 0,)

dd = ku.calc_vessel_hydrodynamics_Ccode(new_file_name, 'vessels', True, calcflow_param_dict_for_tumorcode, False, True)
print(dd)