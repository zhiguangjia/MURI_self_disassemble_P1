#!/bin/bash

################################################################

         script formart

**  .py is python script
    .inp is charmm input script

     groamcs has no specific input file but may need same index (.ndx) topology file (.itp .top)

################################################################

          data formart

 ** .agr will be a xmgrace save file, which has raw data inside it


 

################################################################





################################################################

   inital structure ----  AA part

################################################################

 1, add new residue type M1 M2 M3 in AA force field, build AA peptide

  working dir:  /home/zgjia/Project/MURI/2015ACS/AA/monomer_11_re
 
          ** in  Readme.txt there is discuss about parameters and peptide length/sequence

          ** output  random_*.pdb  will be continued in nex step:  /home/zgjia/Project/MURI/2015ACS/CG_monomer/random_11mer_new_V2_re 


  in scripts folder, we have example charmm input file:  chain_random1.inp


################################################################

   inital structure ----  CG part

################################################################

 1, add new residue type M1 M2 M3 in CG force field, build CG peptide

  working dir:  /home/zgjia/Project/MURI/2015ACS/CG_monomer/random_11mer_new_V2_re

  output is random_CG_*.pdb and correspoding itp

  in topology_inital_structure , we have an example itp:   random_CG_1.itp


   ** some analysis require psf formart as it is a clean topology , but this can easy transferemed from itp+gro/pdb

################################################################

   self aggregation part (Figure 2, also part of Figure 5)

################################################################

****  In this part, the stragedy are same for aggregation of system contains 10 to 320 monomers.

****  all system has three replica in /home/zgjia/Project/MURI/2015ACS/CG_monomer/random_11mer_new_V3_re/

****  all system has simular constract:

      5*n peptide is random placed in box, only limitation is the minimal distance between peptides (5 different P1 is added in a 1:1:1:1:1 ratio, e.g. for 40 mer system, n=4).

      The system is solved in CG water and run 1-10 us as descript is SI table 1

***************************************************************
 
      analysis

*****  Figure 2

in Figure_2
    

   ** for the DBSCAN analysis
 
     ==> first, we need get minus distance between peptidews (or M1) , there are two script

            cp  ../../CG_simulation_11mer_re/Analysis_script/M1_mindist_V2.py scripts/Figure_2_5_M1_mindist_V2.py 
            cp  ../../CG_simulation_11mer_re/Analysis_script/pep_mindist_V2.py scripts/Figure_2_5_peptide_mindist_V2.py  
 
            ** optaion step, analysis the histogram of the distribution, set cut-off after 1st pea

     ==> now, use the min distance, we perform cluster

            cp ../../CG_simulation_11mer_re//Analysis_script/DBSCAN_mindist.py scripts/Figure_2_5_DBSCAN_mindist.py

     ==> In current figure, we focus on the largest cluster
 
            cp ../../CG_simulation_11mer_re//Analysis_script/simple_export_size_and_maxsize_of_Cluster.py scripts/Figure_2_5_export_size_and_maxsize_of_Cluster.py 

     ==> raw data

            cp  ../Fig_cluster_number/320mer/sim1_8-5_V2.agr  Figure_2/Figure_2_320mer_cluster.agr


   ** the density profile:
    
     ==> script cp ../../CG_simulation_11mer_re/Analysis_script/profile.py scripts/Figure_1_density_profile.py
   
     ==> data   cp ../Fig_2_profile/320mer/320mer_distribution_v2.agr  Figure_2/
                cp ../Fig_2_profile/10mer/distribution_V2.agr          Figure_2/10mer_distribution_v2.agr      

**** 



################################################################

   M1 cluster diffuse (Figure 3)

################################################################

*** motivation is tae the DBSCAN clustered based on M1, then on vmd color M1 by current id, load snap shot a few time later, see how thos mixed

*** step 1, get a index file which sorted by id , and which atom are included in this id 

    ****  1  origin DBSCAN file

           ** raw data: cluster id file

          cp /home/zgjia/Project/MURI/2015ACS/CG_simulation_11mer_re/V2_320_rep1/analysis/DBSCAN/25ns_0-10000ns_M1_DBSCAN_from_min_8-5_cluster_index.dat Figure_3/320_sim1_DBSCAN_cluster_index.dat
 
          ** simple copy the 5800 ns data line to a new file: 5800_8-5_id  

          cp /home/zgjia/Project/MURI/2015ACS/Manuscript_11mer/Fig_M1_cluster_by_id/5800_8-5_id  Figure_3/320_sim1_DBSCAN_cluster_index_5800ns.dat

 
          ** script convert this to index file

          cp  ../../CG_simulation_11mer_re/Analysis_script/simple_export_index_in_M1_Cluster_with_sizelimit.py   scripts/Figure_3_clust_id_to_atom_index.py




                 




################################################################

   DiI inbedding (Figure 4)

################################################################

***   we take aggregated polymer  from last step, put DiI in random position

  ***  raw date:

       sas

       cp /home/zgjia/Project/MURI/2015ACS/CG_simulation_11mer_re/V2_80_DII_5_rep1/sas.agr  Figure_4/DII_SAS.agr

       profile
 
       cp /home/zgjia/Project/MURI/2015ACS/Manuscript_11mer/Fig_2_profile/DII_new.agr         Figure_4/DII_profile.agr


  ** no new script, profile calculated same as figure 2

  ** sas by gmx tool:   gmx_mpi sasa


################################################################

   Dissemble (Figure 5)

################################################################

  ** no new script, cluster calculated same as figure 2 ; sasa calculated as figure 4

  **  raw data

     ************   80 mer

      for both bCA and control 

           **  cluster

           cp /home/zgjia/Project/MURI/2015ACS/Manuscript_11mer/Fig_dissemble_80mer/disassemble_cluster_average_80_running_average.agr Figure_5/80mer_cluster.agr

           ** SASA
 
           cp /home/zgjia/Project/MURI/2015ACS/Manuscript_11mer/Fig_dissemble_80mer/disassemble_sasa_average_80_running_average.agr  Figure_5/10mer_SASA.agr



     ************   10 mer

      for both bCA and control

          cp /home/zgjia/Project/MURI/2015ACS/Manuscript_11mer/Fig_dissemble_10mer/disassemble_cluster_average_80_running_average.agr Figure_5/10mer_cluster.agr

          cp /home/zgjia/Project/MURI/2015ACS/Manuscript_11mer/Fig_dissemble_80mer/disassemble_sasa_average_80_running_average.agr    Figure_5/10mer_SASA.agr 



################################################################

   peptide sequence correction

################################################################

 1  as the peptide 1 mistakely use peptide 2 sequence, after copy sequence here, delete peptide 1 , and shift the number


***************************************************


git hub part


**** git

git init
git add  readme.txt
git commit -m "MURI_self_disassemble_P1"

 git add *
git commit -m "V3"


git remote add origin git@github.com:zhiguangjia/MURI_self_disassemble_P1.git
git branch -M main
git push -u origin main

***  for later change:

git add .

git commit -m "revision_1"

git push origin main


username zhiguangjia
 token see laptoop code dataase

#### 6.21  change to version 2

git add .

git commit -m "revision_1_V2"

git push origin main --force  # because delete some file in github online 



git add .

git commit -m "revision_1_V3"

git push origin main













