#!//usr/bin/python

import MDAnalysis
import numpy
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
import numpy as np
import numpy.linalg 
from numpy.linalg import norm
from collections import OrderedDict
from pprint import pprint

import MDAnalysis.analysis.distances
from  numpy import matrix

print ' Size and Shape Characteristics of Polystyrene and Poly(ethylene oxide) Star Polymer Melts Studied By Atomistic Simulations '

print 'Asphericity a fluctuates between the values of 0 and 1 and for a perfectly spherical shape a = 0  '

print 'Prolateness, on the other hand, reaches values between -1 and 1. For -1 the object is fully oblate and 1 implies a perfectly prolate shape '


print 'need cluster index , to find out largest cluster and will only print the largest clusterr  result '

print 'must consider all atoms, other wise we will have a lot of calvity '


structurefile    = sys.argv[1]

cluster_id      =  sys.argv[2]

size            =  int(sys.argv[3])

##########################

print('read psf, list all M1')

ifile1 = open(structurefile,'r')


res_index    =  OrderedDict()


counter = 1

for line in ifile1:

    columns = line.split()
     
    if len(columns) > 8 and columns[3] == 'GLM1' : #
           
            
       if columns[4] == 'BB'  :  # 1 st atom in GLM1 , record +2 +3 +4 4 
          res_index[counter] = OrderedDict()
          res_index[counter]['index'] = [int(columns[0])+2, int(columns[0])+3, int(columns[0])+4 ]
          res_index[counter]['resid'] = columns[2]
          res_index[counter]['segid'] = columns[1]
 

          counter += 1 

    if len(columns) > 8 and columns[1][0] == 'W':     
       break

if counter != 0:
   print('total ', counter  , 'M1 found ')
else:

   print('no peptide found, something wrong')
   quit()       
##########################
##########################

print('now read cluster index ,  it has time, follow by cluster id of each peptide,  ')

cluster_index    =  OrderedDict()

ifile_cluster = open(cluster_id,'r')

for line in ifile_cluster:

    columns = line.split()

    if len(columns) > 5 : #

       time = int(float(columns[0]))   # we will not have ps unit time 
       cluster_index[time]  = []

       # now find largest cluster 
       
       i = 0 
       temp_index = {}
       current_max_cluster_number = 0
 
       for current_id in columns :
         
           if i == 0 :  # i will be the M1 index from 1
              pass

           else:

              current_id = int(float(current_id)) 
              if current_id not in temp_index :
                 temp_index[current_id] = {}
                 temp_index[current_id]['count'] = 1
                 temp_index[current_id]['list']  = []
                 temp_index[current_id]['list'].append(i) 
                 temp_index[current_id]['index'] = []
                 temp_index[current_id]['index'] += res_index[i]['index']   # SC2 SC3 SC4

              else :
                 temp_index[current_id]['count'] += 1
                 temp_index[current_id]['list'].append(i)

                 temp_index[current_id]['index'] += res_index[i]['index'] 
              #  now  compare with current max 

             # if len(temp_index[current_id]['list'])  > current_max_cluster_number :

             #    current_max_cluster_number = len(temp_index[current_id]['list'])
             #    current_max_cluster_id     = current_id  
           i += 1     
           # after 1st time column,  i will equal to the resid 

       # now record result of this frame 
       #cluster_index[time] = temp_index[current_max_cluster_id]

       for current_cluster in  temp_index:

        if temp_index[current_cluster]['count'] > size:
 
           print  current_cluster, ' contain ', temp_index[current_cluster]['count'],  ' M1'

           # convert to vmd formart:
           temp = 'index '
           for index in  temp_index[current_cluster]['index'] :
               temp += str (' ' + str(index-1) ) 

           print temp 





















