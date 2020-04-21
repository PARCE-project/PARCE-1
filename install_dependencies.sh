#!/bin/bash


#PARCE: Protocol for Amino acid Refinement through Computational Evolution

#From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
#Computer Physics Communications
#Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio
#Year: 2020


########################################################################################
# Ubuntu 16.04 tested commands
########################################################################################

sudo apt-get update
sudo apt-get install pdb2pqr
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-yaml
sudo apt-get install python3-tk
pip3 install GromacsWrapper

########################################################################################
# Ubuntu 16.04 script to test additional functionalities
########################################################################################

python3 test.py
