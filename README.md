---------------------------------
BASIC VASCULAR BRANCHING ANALYSIS
---------------------------------

This example takes a segmentation of a vascular branching network (a 3D binary image stack) and demonstrates how to convert this segmentation into a branching network of nodes and elements, and how to analyse structure and function in that network.

This example works for data "F8", which due to size is not available on Github. Please contact us for more information. 

Python requirements are listed in requirements.txt, please note that this code does not work exactly the same way with skan version 0.9 as it does with version 0.10, due to imprpovements in skan's algorithms. This code is not checked with skan <0.10

*segmentation2skeleton.py* needs to be run first. It creates a  skeleton in the data/sample-skeleton directory, initially a complete tree, and then a series of edited trees which assign radii to vessels and remove loops.

Next run *branch_analysis.py* which creates spreadsheets in the data/sample-branch directory and allows you to see statistics about the nature of vascular branching.

To come.... conversion to a computational model



