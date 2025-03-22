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

rectange_range = float(sys.argv[4])
axis_length    = float(sys.argv[5])

bins           = int(sys.argv[6])

blocksize        = int(sys.argv[7])

try:
    initial_time = int(sys.argv[8])

    try:
       final_time = int(sys.argv[9])

       try:
            skip  =  int(sys.argv[10])
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

#res_index    =  OrderedDict()

##########################

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)

print  ' total frame number ',  targetprotein.trajectory.n_frames
#half_frame = int(targetprotein.trajectory.n_frames/2)

##  write file title  , set up average array 

ofile_COM_list  = open(str(outfilename + '_hist.dat'),'w') # open file for writing


ofile_COM_list.write('# probability  BB   SC_M1  SC_M2   SC_M3   W  DII_head DII_tail DII  ' + ' ' +  "\n")

#ofile_COM_distribute  = open(str(outfilename + '_COM_distribute.dat'),'w') # open file for writing

x_BB     = []
x_SC_M1     = []
x_SC_M2     = []
x_SC_M3     = []

x_W     = []
x_DII_head     = []
x_DII_tail     = []
x_DII     = []


inital_frame_flag = 0
for ts in targetprotein.trajectory:
    current_time = ts.time/1000 # to ns

    if current_time < initial_time or float(ts.frame)%skip != 0 :
         pass
    elif current_time > final_time :

         print ts.time, 'exceed ', final_time, 'fnishing and output final result' 
         break
    else :
         print ts.time/1000, 'ns'
         #ofile_COM_list.write('time' + ' ' + str(current_time) + "\n")         
         ###################################################################################################
         ####################  a  long test see if it is correct  


         #box = MDAnalysis.coordinates.base.Timestep.dimensions
         #box = MDAnalysis.coordinates.timesteps.Timestep.dimensions
         box = ts.dimensions
         #print  (box)
         counter = 0
   
         current_center = [ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2]

         #print current_center, float(current_center[0])
         ##  laver our x is along the x axis, so only x is needed

         ##  BB

         xmin_current = current_center[0] - axis_length 
         xmax_current = current_center[0] + axis_length
         ymin_current = current_center[1] - rectange_range
         ymax_current = current_center[1] + rectange_range
         zmin_current = current_center[2] - rectange_range
         zmax_current = current_center[2] + rectange_range



         Current_sel   = targetprotein.select_atoms(" name BB CTE and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

 
         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_BB.append( position[0] - current_center[0] )

       #  W = MDAnalysis.Writer("BB.pdb")
       #  W.write(Current_sel)

         Current_sel   = targetprotein.select_atoms(" (resname GLM1 and not name BB CTE) and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )


         for position  in  Current_sel.positions :
             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :
                 x_SC_M1.append( position[0] - current_center[0] )

         Current_sel   = targetprotein.select_atoms(" (resname GLM2 and not name BB CTE) and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

      #   W = MDAnalysis.Writer("M2.pdb")
      #   W.write(Current_sel)
      #   quit()
         
         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_SC_M2.append( position[0] - current_center[0] )

         Current_sel   = targetprotein.select_atoms(" (resname GLM3 and not name BB CTE) and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

       #  W = MDAnalysis.Writer("M3.pdb")
       #  W.write(Current_sel)

         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :
    
                 x_SC_M3.append( position[0] - current_center[0] )
 
         Current_sel   = targetprotein.select_atoms(" name W and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

       #  W = MDAnalysis.Writer("Protein.pdb")
       #  W.write(Current_sel)
       #  quit()

         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_W.append( position[0] - current_center[0] )



         Current_sel   = targetprotein.select_atoms(" name W and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

       #  W = MDAnalysis.Writer("Protein.pdb")
       #  W.write(Current_sel)
       #  quit()

         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_W.append( position[0] - current_center[0] )

         Current_sel   = targetprotein.select_atoms(" (resname DII and name C1 SC1A SC2A SC3A SC4A SC1B SC2B SC3B SC4B ) and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_DII_head.append( position[0] - current_center[0] )

         Current_sel   = targetprotein.select_atoms(" (resname DII and name C1A  C1B  C2A C3A C4A C5A C2B C3B C4B C5B ) and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_DII_tail.append( position[0] - current_center[0] )

         Current_sel   = targetprotein.select_atoms(" resname DII  and ( prop x > %d and prop x < %d and prop y > %d and prop y < %d and prop z > %d and prop z < %d )  " %(xmin_current,xmax_current,ymin_current,ymax_current,zmin_current,zmax_current  )  )

         for position  in  Current_sel.positions :

             if (position[1]-current_center[1])*(position[1]-current_center[1]) + (position[2]-current_center[2])*(position[2]-current_center[2]) <= rectange_range*rectange_range :

                 x_DII.append( position[0] - current_center[0] )


       #  z_array_.append(Current_sel.center_of_mass()[2] - membrane.center_of_mass()[2] )
#print x_BB
histBB, bin_edgesBB    = numpy.histogram(x_BB,bins,(-axis_length, axis_length))
#print histBB
#print bin_edgesBB
histSC_M1, bin_edgesSC_M1    = numpy.histogram(x_SC_M1,bins,(-axis_length, axis_length))
histSC_M2, bin_edgesSC_M2    = numpy.histogram(x_SC_M2,bins,(-axis_length, axis_length))
histSC_M3, bin_edgesSC_M3    = numpy.histogram(x_SC_M3,bins,(-axis_length, axis_length))

histW, bin_edgesW    = numpy.histogram(x_W,bins,(-axis_length, axis_length))
histDII_head, bin_edgesDII_head    = numpy.histogram(x_DII_head,bins,(-axis_length, axis_length))
histDII_tail, bin_edgesDII_tail    = numpy.histogram(x_DII_tail,bins,(-axis_length, axis_length))
histDII, bin_edgesDII    = numpy.histogram(x_DII,bins,(-axis_length, axis_length))

###  now  consider 

# for some group, it may not ecist, so make it 1
BBmax = sum(histBB)
SC_M1max = sum(histSC_M1)
SC_M2max = sum(histSC_M2)
SC_M3max = sum(histSC_M3)
Wmax = sum(histW)
DII_headmax = max(sum(histDII_head), 1)
DII_tailmax = max(sum(histDII_tail), 1)
DIImax = max(sum(histDII),1)


#BBmax = max(histBB)
#SC_M1max = max(histSC_M1)
#SC_M2max = max(histSC_M2)
#SC_M3max = max(histSC_M3)
#Wmax = max(histW)


for x in range(bins):

  if x > blocksize and x < bins-blocksize-1 :

     # order  BB   SC_M1  SC_M2   SC_M3   W   DII_head DII_tail DII 
     ofile_COM_list.write( str( round(bin_edgesBB[x],2)  ) + " " )

     ofile_COM_list.write(str(round(float( sum (histBB[x-blocksize:x+blocksize+1])   )/BBmax,3)) + " ")
                                           #  upper raneg need +1, python sum role
                                         
     ofile_COM_list.write(str(round(float( sum (histSC_M1[x-blocksize:x+blocksize+1]) )/SC_M1max,3)) + " " )
     ofile_COM_list.write(str(round(float(sum (histSC_M2[x-blocksize:x+blocksize+1]) )/SC_M2max,3)) + " " )
     ofile_COM_list.write(str(round(float(sum (histSC_M3[x-blocksize:x+blocksize+1]) )/SC_M3max,3)) + " " )
     ofile_COM_list.write( str(round(float(sum (histW[x-blocksize:x+blocksize+1])  )/Wmax,3)) + " "  +  str(round(float(sum ( histDII_head[x-blocksize:x+blocksize+1]) )/DII_headmax,3)) + " " +  str(round(float(sum (histDII_tail[x-blocksize:x+blocksize+1]) )/DII_tailmax,3)) + " " +  str(round(float(sum (histDII[x-blocksize:x+blocksize+1]))/DIImax,3)) + " " +"\n")















