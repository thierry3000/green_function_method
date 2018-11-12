# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 17:18:15 2015

@author:  Thierry Fredrich
"""
import sys
import os
import string
import h5py
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import krebsutils


def get_o2_field(fn):
  with open('TissueLevels.out','r') as f:
    read_data = f.read()
  f.closed
  print("we read the tissue oxygen levels")
  niceList=[]
  start = False
  for (linenumber,line) in enumerate(read_data.split(os.linesep)):          
      if start:
          print line
          for wert in line.split():
              niceList.append(wert)
      if line =='Solute 1':
          start = True
      if not string.find(line, 'pmean,'):
          start = False
  '''on this stage secombs data is in array nice List
  now we get Michaels data...
  '''
  niceList.pop()
  niceList.pop()
  niceList.pop()
  niceList.pop()
  niceList.pop()
  niceList.pop()
  niceList.pop()
  niceList.pop()
  
  print("we load the tissue data")
  with open('TissueSources.out','r') as f:
    read_data = f.read()
  f.closed
  discretization_points=[]
  X=[]
  Y=[]
  Z=[]

  myindex=0
  for (linenumber,line) in enumerate(read_data.split(os.linesep)):          
      if linenumber==0:
          print line
          for (pos,wert) in enumerate(line.split()):
            if pos<3:
              discretization_points.append(int(wert))
      if linenumber>1 and (not 'Tissue point xyz' in line):
          print line
          for (pos,wert) in enumerate(line.split()):
            if myindex<discretization_points[0]:
              X.append(float(wert))
              myindex=myindex+1
              continue
            if myindex>discretization_points[0]-1 and myindex<(discretization_points[0]+discretization_points[1]):
              Y.append(float(wert))
              myindex=myindex+1
              continue
            if myindex>(discretization_points[0]+discretization_points[1]-1) and myindex<(discretization_points[0]+discretization_points[1]+discretization_points[2]):
              Z.append(float(wert))
              myindex=myindex+1
              continue
  
  
  discretization_points=np.asarray(discretization_points)

  niceList2 = []
  start = False
  for (linenumber,line) in enumerate(read_data.split(os.linesep)):
        if 'source strengths for solute' in line:
          start = False
            
        if start:
            print line
            for wert in line.split():
                niceList2.append(int(wert))
        if line =='Tissue point xyz indices':
            start = True
  print('Found %i datapoint, got %i discretization points' % (len(niceList),discretization_points[0]*discretization_points[1]*discretization_points[2]))  
  
  niceList2=np.asarray(niceList2)
  dim=0
  number=0
  coords=[]
  coords2=np.zeros([3,len(niceList2)])
  
  for aint in niceList2:
    if dim==0:
      x=X[aint-1]
      coords2[dim,number]=x
      dim=dim+1
      continue
    if dim==1:
      y=Y[aint-1]
      coords2[dim,number]=y
      dim=dim+1
      continue
    if dim==2:
      z=Z[aint-1]
      coords2[dim,number]=z
      dim=dim+1
    if dim==3:
      coords.append([x,y,z])
      number=number+1
      dim=0
      
    
  with h5py.File(fn, 'r') as f:
    #need that data from hdf
    po2_field_mw = np.asarray(f['/po2/vessels/po2field'])
    ld = krebsutils.read_lattice_data_from_hdf_by_filename(fn,'/po2/vessels/field_ld')
    po2_field_mw_at_secomb_positons = krebsutils.sample_field(np.asarray(coords2,dtype=np.float32).transpose(),po2_field_mw,ld)
  
  return coords2, niceList2, po2_field_mw_at_secomb_positons


if __name__=='__main__':
  
  fn = sys.argv[1]
  pos, o2_secomb, o2_mw = get_o2_field(fn)
  
  fig1=plt.figure()
  ax = fig1.add_subplot(111, projection='3d')
  
  ax.scatter(pos[0,:],pos[1,:],pos[2,:], c=o2_secomb)
  
  plt.show()
  
  