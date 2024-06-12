# Hypervariability of Ramachandran maps

This repo contains the scripts for generating Bond geometry-specific Ramachandran maps from MD trajectory of a given protein. The trajectory should be PDB format. The main script for the map generation is *rotation_parallel.py* . I have also included the GROMACS commands and inputs I used for calculating the energy of every two-linked peptide unit for a given residue (*get_potential_from_gromacs.py*). 

To understand why I did this, please read **Ravikumar, A., & Srinivasan, N. (2021). Hypervariability of accessible and inaccessible conformational space of proteins. Current research in structural biology, 3, 229-238.**
