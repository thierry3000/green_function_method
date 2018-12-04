# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:56:26 2014

@author:  Thierry Fredrich

@brief:   This file was used to convert "vessels-q2d-8mm-P10-typeE-7x2L130-sample01_secombComparison.h5"
          into a format usable by the secomb software.
          
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
          
          At time of publication (I), the software in this repository was already present. 
          In (I) the oxygen distribution within artificial blood vessel networks was calculated. 
          At time of publication the software in this repository was available and 
          we compared run time and accuracy to that. 
          The present file was used to transform the simulation result of tumorcode (II) 
          into a format that could be used by the present software.
"""

import h5py
import krebsutils
import sys
import numpy as np

percentage_of_boundary_security = 1. #0 is no security, 1 means one time the systemsize

if __name__ == '__main__':
    fn = sys.argv[1]
    group = sys.argv[2]
    fout = open('Network.dat','w+')
      
    ''' the single capillary case was used to test the software '''
    single_capillary_case = fn == 'vessel-single-all.h5'
    with h5py.File(fn,'r') as f:
      if single_capillary_case:
          node_grp = f['thierry_case2/vessels/nodes']
          vessel_grp = f['thierry_case2/vessels/']
      else:
          node_grp = f[group + '/nodes']
          vessel_grp = f[group]
      roots_list = np.asarray(node_grp['roots'])
      segment_list,real_world_positions, pressure_at_node,radius_data, flow_data, hematocrit_data, flag_data_edges, a, b= krebsutils.read_vessels_from_hdf(vessel_grp,['position','pressure','radius','flow','hematocrit','flags','node_a_index','node_b_index'])
    
    
    pressure_at_node = pressure_at_node[:,0]
    radius_data =radius_data[:,0]
    flow_data = flow_data[:,0]
    hematocrit_data = hematocrit_data[:,0]
    flag_data_edges = flag_data_edges[:,0]
    a = a[:,0]
    b = b[:,0]
    ''' *************** Unit change **********************'''
    ''' secomb uses flow measured in nl as input
        see flowfac in his implementation
    '''
    flow_data  = flow_data * 60/1e6
    ''' secomb measures pressure in mmHg
        therefore 1/7.5= 0.1333333
    '''
    pressure_at_node = pressure_at_node * 0.133333333333333333333333333333
    
    
    print ('Creating intro...')
    
    #data needed for intro
    xmin = min(real_world_positions[:,0])
    xmax = max(real_world_positions[:,0])
    ymin = min(real_world_positions[:,1])
    ymax = max(real_world_positions[:,1])
    zmin = min(real_world_positions[:,2])
    zmax = max(real_world_positions[:,2])
    ex_x = xmax-xmin
    ex_y = ymax-ymin
    ex_z = zmax-zmin
    #percentage_of_boundary_security = 10. #0 is no security, 1 means one time the systemsize
    for (i,apoint) in enumerate(real_world_positions):
      x = apoint[0]+ 1000. #percentage_of_boundary_security*ex_x
      y = apoint[1]+ 1000. #percentage_of_boundary_security*ex_y
      z = apoint[2]+ 1050. #percentage_of_boundary_security*ex_z
      real_world_positions[i,:]=[x,y,z]
    
    xmax_box = max(real_world_positions[0,:])+percentage_of_boundary_security*ex_x
    ymax_box = max(real_world_positions[1,:])+percentage_of_boundary_security*ex_y
    zmax_box = max(real_world_positions[2,:])+percentage_of_boundary_security*ex_z
    
    print('shifted x: %f, extension x: %f, max box x: %f'% (percentage_of_boundary_security*ex_x, ex_x, xmax_box))
    print('shifted x: %f, extension x: %f, max box y: %f'% (percentage_of_boundary_security*ex_y, ex_y, ymax_box))
    print('shifted x: %f, extension x: %f, max box z: %f'% (percentage_of_boundary_security*ex_z, ex_z, zmax_box))    
    print('suggest plane z=%f for projection!' %(ex_z/2+percentage_of_boundary_security*ex_z))
    
    if single_capillary_case:
        fout.write('Imported from %s hdf Network\n' % fn[-40:] )
        fout.write('  %f %f %f box dimensions in microns\n' % (4000, 400, 400))
        fout.write('  %i %i %i number of tissue points in x,y,z directions\n' % (50, 10, 10))
        fout.write('  %f outer bound distance (100. in test case)\n' %100. )
        fout.write('  %f max. segment length.  (5. in test case)\n'% 15.)
        fout.write('  %i nodsegm, max. allowed number of segments per node\n' %3)
    else:
        fout.write('Imported from %s hdf Network\n' % fn[-40:] )
        fout.write('  %f %f %f box dimensions in microns\n' % (xmax_box,ymax_box,zmax_box))
        fout.write('  %i %i %i number of tissue points in x,y,z directions\n' % (30,30,15))
        fout.write('  %f outer bound distance (100. in test case)\n' %1000. )
        fout.write('  %f max. segment length.  (5. in test case)\n'% 55.)
        fout.write('  %i nodsegm, max. allowed number of segments per node\n' %3)
    print('Finished writting intro')
    
    
    #get MW stuff
    is_set = lambda flags_,flag: np.asarray(np.bitwise_and(flags_, flag), np.bool)
    
    #finding bad segments
    circulated_segments = is_set(flag_data_edges,krebsutils.CIRCULATED)
    connected_segments = is_set(flag_data_edges,krebsutils.CONNECTED)
    boundary_data = is_set(flag_data_edges,krebsutils.BOUNDARY)
    good_segments=circulated_segments*connected_segments
    
    segment_list=segment_list[good_segments]
    radius_data=radius_data[good_segments]
    flow_data=flow_data[good_segments]
    hematocrit_data=hematocrit_data[good_segments]
    
    print('Deleted %i uncirculated segments.'%(len(good_segments)-sum(good_segments)))
    
    fout.write('  %i total number of segments\n'%len(segment_list))
    
    fout.write('name	type	from	to	diam		flow		HD\n')
    how_many_sign_switch=0
    for (i, segment) in enumerate(segment_list):
        #check for direction
        if(pressure_at_node[segment[0]]>pressure_at_node[segment[1]]):
            fout.write('%i 5 %i %i %f %5.10f %f\n'%(i+1, segment[0]+1,segment[1]+1, 2*radius_data[i], flow_data[i],hematocrit_data[i]))        
        else:
           how_many_sign_switch=how_many_sign_switch+1
           fout.write('%i 5 %i %i %f %5.10f %f\n'%(i+1, segment[0]+1,segment[1]+1, 2*radius_data[i], -1*flow_data[i],hematocrit_data[i]))        
    print('We switche %i times the sign of the flow!'%how_many_sign_switch)
    
    fout.write('  %i total number of nodes\n'%real_world_positions[:,:].shape[0])
    fout.write('name x y z\n')
    for i in range(0,real_world_positions[:,:].shape[0]):
        fout.write('%i %f %f %f\n'%(i+1,real_world_positions[i,0], real_world_positions[i,1], real_world_positions[i,2]))
        
    
    ''' create boundary conditions '''
    fout.write('%i total number of boundary nodes\n' % len(roots_list))
    fout.write('node	bctyp	press/flow	HD	PO2\n')
    for i in range(0,len(roots_list)):
        #fout.write('%i 0 %f 0.37 55. 1\n'%(roots_list2[i]+1,pressure_at_nodes[roots_list2[i]]))
        fout.write('%i 2 %f 0.37 55. 1\n'%(roots_list[i]+1,pressure_at_node[roots_list[i]]))
        fout.write('Secomb flow is 2, Sebomb pressure is 0, pressure: %f, flow: %f\n' %(pressure_at_node[roots_list[i]],flow_data[roots_list[i]]))
    print('Minimal flow in system: %f'% (min(abs(flow_data))))
    fout.close()