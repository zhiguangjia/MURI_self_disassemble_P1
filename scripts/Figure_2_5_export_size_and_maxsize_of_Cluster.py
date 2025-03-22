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


cluster_id      =  sys.argv[1]


outfilename   = sys.argv[2]

try:
    initial_time = int(sys.argv[3])

    try:
       final_time = int(sys.argv[4])

       try:
            skip  =  int(sys.argv[5])
       except:
            skip  =  1

    except:
       final_time = 9999999999
       skip  =  1

except:
    initial_time = 0
    final_time = 9999999999
    skip  =  1



##########################


   





##########################
##########################

print('now read cluster index ,  it has time, follow by cluster id of each peptide,  ')

cluster_index    =  OrderedDict()

ifile_cluster = open(cluster_id,'r')
ofile  = open(str(outfilename + '_clustersize_list.dat'),'w') # open file for writing
ofile_max  = open(str(outfilename + '_clustersize-max_list.dat'),'w') # open file for writing
ofile_average  = open(str(outfilename + '_clustersize-average_list.dat'),'w') # open file for writing



for line in ifile_cluster:

    columns = line.split()

    if len(columns) > 5 : #

       current_time = int(float(columns[0]))   # we will not have ps unit time 


       if current_time < initial_time or float(current_time)%skip != 0 :
           pass
       elif current_time > final_time :

           print current_time, 'exceed ', final_time, 'fnishing and output final result' 
           break
       else :
           print current_time, 'ns'


           cluster_index[current_time]  = []

          # now find all cluster 
       
           i = 0 
           temp_index = {}
           current_max_cluster_number = 0
 
           for current_id in columns :
         
               if i == 0 :  # i will be the M1 index from 1  , 1st is time
                   pass

               else:

                   current_id = int(float(current_id))

                   if current_id == -1 :  # noise
                      pass
                   else :  
                      if current_id not in temp_index :    
                         temp_index[current_id] = {}
                         temp_index[current_id]['count'] = 1
                         temp_index[current_id]['list']  = []
                         temp_index[current_id]['list'].append(i) 

   
                      else :
                         temp_index[current_id]['count'] += 1
                         temp_index[current_id]['list'].append(i)


                     #  now  compare with current max 

                      if len(temp_index[current_id]['list'])  > current_max_cluster_number :

                         current_max_cluster_number = len(temp_index[current_id]['list'])
                         current_max_cluster_id     = current_id  
               i += 1     
           # after 1st time column,  i will equal to the resid 

       # now record result of this frame 
       #cluster_index[time] = temp_index[current_max_cluster_id]
 

           ofile.write(str(current_time,) + ' ' )
           ofile_max.write(str(current_time,) + ' ' + str(current_max_cluster_number) + "\n")
           ofile_average.write(str(current_time,) + ' ')
 
           average = 0.0
           count   = 0.0
           for current_id in  temp_index:
               ofile.write( str(temp_index[current_id]['count']) +' '  ) 
               average +=  temp_index[current_id]['count']
               count   += 1.0

           if count == 0.0 :
               ofile_average.write(str(0) + "\n")
           else :
               ofile_average.write(str(average/count) + "\n")
           ofile.write("\n")



            

















