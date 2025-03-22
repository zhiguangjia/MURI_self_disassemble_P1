#!//usr/bin/python

import MDAnalysis
import numpy
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
import numpy.linalg
from numpy.linalg import norm
from collections import OrderedDict
from pprint import pprint

import MDAnalysis.analysis.distances
from  numpy import matrix



structurefile    = sys.argv[1]
targettrajectory = sys.argv[2]

#skip  =  int(sys.argv[4])
outfilename   = sys.argv[3]

try:
    initial_time = int(sys.argv[4])

    try:
       final_time = int(sys.argv[5])

       try:
            skip  =  int(sys.argv[6])
       except:
            skip  =  1

    except:
       final_time = 9999999999
       skip  =  1

except:
    initial_time = 0
    final_time = 9999999999
    skip  =  1

print '  '

print '  '



# unit in A!


#######################################
#  relations :
  ############  whole calculation will be made ,VSD and RCK just means the contact method between two chain (who give VSD and who give RCK)

##########################

#  constants
# kB = 3.2976268E-24  cal/K
kB = 1.3806488E-23    # j/K
An = 6.02214179E23
T  = float(300)
R  = 1.987  # (cal/mol.degree)
#frame = int(0)
frame = float(0)
#  later it will be divided, so float may be better??

##########################

print('read psf, list all M1')

ifile1 = open(structurefile,'r')


res_index    =  OrderedDict()

index = 0
counter = 0

for line in ifile1:

    columns = line.split()
     
    if len(columns) > 8 and columns[1][0] == 'P': # because all GLM has in P section
           
       if columns[3] == 'GLM1' :

        #  if columns[4] in ['SC2', 'SC3', 'SC4']:   # no need, simple BB +2 +3 +4
        #     res_index[]            
             
          if columns[4] == 'BB' :
             res_index[counter] = OrderedDict()
             res_index[counter]['BBindex'] = index
             res_index[counter]['resid']   = int(columns[2])
             res_index[counter]['segname'] = columns[1]
             counter += 1 

       index += 1

    if len(columns) > 8 and columns[1][0] == 'W':     
       break

if counter != 0:
   print('total ', counter , 'M1 found ')
else:

   print('no M1 found, something wrong')
   quit()       
##########################


targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)

print  ' total frame number ',  targetprotein.trajectory.n_frames
#half_frame = int(targetprotein.trajectory.n_frames/2)

##  write file title  , set up average array 

ofile_COM_list  = open(str(outfilename + '_M1_min_list.dat'),'w') # open file for writing

ofile_COM_list.write('info ' + ' number  ' + str(counter)  + '  ' + "\n")

#ofile_COM_distribute  = open(str(outfilename + '_COM_distribute.dat'),'w') # open file for writing


segids_list = []
M1_list     = []

inital_frame_flag = 0
for ts in targetprotein.trajectory:
    current_time = ts.time/1000 # to ns
    #print current_time

    if current_time < initial_time or float(ts.frame)%skip != 0 :
         pass
    elif current_time > final_time :

         print current_time, 'exceed ', final_time, 'fnishing and output final result' 
         break
    else :
         print current_time , 'ns'
         ofile_COM_list.write('time' + ' ' + str(current_time) + "\n")         
         ###################################################################################################
         ####################  a  long test see if it is correct  


      #   if inital_frame_flag == 0 :

      #      inital_frame_flag = 1

      #      print 'test selection'
      #      test = targetprotein.select_atoms('segid P*  and resname GLM1')
      #      print test
          #  pprint(vars(test))  <== print alll atom list
      #      print vars(test)

      #      print dir(test)

         #   for att in dir(test) :
         #       print (att, getattr(test,att))   #  this cuould not done as some attribe in it is unfinished 

         #   for segids in test.segids:
         #       if segids not in segids_list:
         #           print segids  
         #           segids_list.append(segids)   #   this give segids PPP1 to end, each has a single line

          #  print test.n_residues           # single number , for 80 mer we have 480 . each has 6 M1 so correct

       #     counter = 0
       #     for residue in test.segids:
       #         print residue, test.resnames[counter], test.resnums[counter],test.resids[counter]
       #         counter += 1                              

          # print test.n_residues           # single number
       #     for segments in test.segments:
       #         print segments
          #  for resnums in  test.resnums:   
          #      print resnums              # still atom based

          #  for fragments in  test.fragments  :
          #      print fragments            #  strange, Cter is cut off so you have 2 n  list ?


             
          #    print 'now loop M1 SC2 SC3 SC4'


          #     Current_sel   = targetprotein.select_atoms(" (resid %d:%d and name CA Ca C N  ) and protein "%(resid_1, resid_2) )

         #box = MDAnalysis.coordinates.base.Timestep.dimensions
         #box = MDAnalysis.coordinates.timesteps.Timestep.dimensions
         box = ts.dimensions
         #print  (box)
         counter = 0
         for current_BB in res_index :        

             #print dir(Current_sel)
             #print Current_sel.positions   


             counter_2 = 0
             for current_BB_2 in res_index :
             
                     
                 if counter < counter_2 :
                     #Current_sel   = targetprotein.select_atoms(" (resid %d:%d and name CA Ca C N  ) and protein "%(resid_1, resid_2) )
                     Current_sel   = targetprotein.select_atoms(" (bynum %d:%d) and resname GLM1 "%(res_index[counter]['BBindex']+3, res_index[counter]['BBindex']+5) )
                                                      # bynum is 1 based

                     Current_sel_2 = targetprotein.select_atoms(" (bynum %d:%d) and resname GLM1 "%(res_index[counter_2]['BBindex']+3, res_index[counter_2]['BBindex']+5) ) 

                     #print dir(Current_sel_2) 
                     x = MDAnalysis.analysis.distances.distance_array(Current_sel.positions,Current_sel_2.positions,box)      
                     new_dist =  round(  x.min(), 3)
                        
                     #new_dist = "{:.3f}".formart(new_dist)
                     ofile_COM_list.write(' ' + str(new_dist) + ' ' ) 

                     res_index[counter][counter_2] = str(new_dist)

                 elif counter == counter_2 :

                     ofile_COM_list.write(' ' + str(0.0) + ' ' )
                 
                 else :

                     ofile_COM_list.write(' ' + str(res_index[counter_2][counter]) + ' ' )

                 counter_2 += 1 
             
          #   print(Current_sel) 
          #   print new_dist

             counter += 1

             ofile_COM_list.write("\n") 



       #  x_array_.append(Current_sel.center_of_mass()[0])
       #  y_array_.append(Current_sel.center_of_mass()[1])
       #  z_array_.append(Current_sel.center_of_mass()[2] - membrane.center_of_mass()[2] )

print('if not consider box reverge')

print('average, minus thi value to PMF y')
#print numpy.average(z_array_)











