#!/home/zgjia/Software/Python3/miniconda3/bin/python3
import sys,string
import numpy as np
import math
#import MDAnalysis
import math
#from  MDAnalysis import analysis
from sklearn.cluster import DBSCAN

from collections import OrderedDict


try:
     distfile      = sys.argv[1]

     current_eps         = float(sys.argv[2])
     current_min_samples = int(sys.argv[3])
     outfilename   = sys.argv[4]

except:
     print("Usage:",sys.argv[0], "pdb-for-coordinate pdb-for-formart  outfile")
     sys.exit(1)


print('......')



##### Variable Initializations ##########

ofile_cluster_noise = open(str(outfilename + '_noise.dat'),'w') # open file for writing
ofile_cluster_number = open(str(outfilename + '_cluster_number.dat'),'w') # open file for writing
ofile_cluster_index  = open(str(outfilename + '_cluster_index.dat'),'w') # open file for writing
ofile_cluster_number_of_AA_in_cluster = open(str(outfilename + '_number_of_AA_in_cluster.dat'),'w') # open file for writing




#############################################################
###  molecule part [ moleculetype ]   [ atoms ]


ifile1 = open(distfile,'r') # open file for reading


ini_flag = 1

##  atom part, genergae bond list

for line in ifile1:

    columns = line.split()
    if len(columns) == 0:
       pass

   # elif columns[0] == 'CRYST1':
   #    ofile.write(line)

    else:

       if columns[0] == 'info' :
          total_number = int(columns[2])  # either M1 or peptide number
          
       elif columns[0] != 'time'  : 


          if len(columns) != total_number   :
             print('number of res not match unit number', len(columns),total_number )
             print(line)
             quit()
          current_column = 0
          for i in columns:

             i = float(i) 
             x[current_row][current_column] = i  # current_row is inislzed later in time part

             if current_column == current_row and i != 0 :

                print('crossterm non 0 at time, current_column current_row ', current_time,current_column,current_row) 
                quit()
             current_column += 1

          current_row += 1

          ####
       if columns[0] == 'time' :

          current_time = columns[1]

          if ini_flag == 1 :

             # first frame, skip

             ini_flag = 0

             # initialize or clean the array
             x = np.zeros(shape=(total_number,total_number))
             current_row = 0

          else:


             ofile_cluster_number.write(str(current_time) + ' ')
             ofile_cluster_index.write(str(current_time) + ' ')
             ofile_cluster_noise.write(str(current_time) + ' ')
             ofile_cluster_number_of_AA_in_cluster.write(str(current_time) + ' ')


             # treat last frame data
          #   print(x)
             db = DBSCAN(eps=current_eps, min_samples=current_min_samples, metric='precomputed'  ).fit(x)
             
          #   print (db)
          #   print (dir(db))

             labels = db.labels_
             n_cluster = len(set(labels)) - (1 if -1 in labels else 0)   # set remove duplicate element
             n_noise  = list(labels).count(-1)
 
          #   print(labels)
          #   print('n_cluster',n_cluster) 
          #   print('n_noise',n_noise) 

          #   print ('done time', current_time)
             
             for j in labels :
                 ofile_cluster_index.write(' ' + str(j))
             ofile_cluster_index.write("\n")


             ofile_cluster_number.write(' ' + str(n_cluster) + ' ' + "\n")
             ofile_cluster_noise.write(' ' + str(n_noise) + ' ' + "\n")


             #  now check how many cluster 
             ofile_cluster_number_of_AA_in_cluster.write(' ' + str(total_number - n_noise) + ' ' + "\n") 
             

              
          #   ofile_cluster_number.write("\n")
     

             # initialize or clean the array
             x = np.zeros(shape=(total_number,total_number))
             current_row = 0


          #  ** for BSA, all coil became blank stranvge











