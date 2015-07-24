# centering_clusters

Implementation of the mis-centering corrections outlined in the document 
'How to use the centering probabilities of galaxy clusters for cosmological investigations' by E. Rozo and M. Becker. 
This writeup can be found in 'documents/miscentering_writeup_v6.pdf'. In particular, the equations implemented are
72-78, for two centers with probability p1 and p2, with p1 + p2 =1. The input catalogs are the cross-correlation signal of stacked clusters about the most likely and second most likely centers, as a function of radial distnace R (Mpc/h). The catalogs used were provided by E. Sheldon and are located in 'catalogs_erin/'.
The program 'code/centering.py' has the relevant functions, and it is called by 
the program 'code/test_centering_module.py'. 
