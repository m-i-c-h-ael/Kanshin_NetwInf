# Network inference from high-resolution time course data
This repo contains code, input and output files for a network inference project. Input is phosphoproteomics data from Kanshin, 2015. These are mapped to a STRINGdb network which is expanded by one shell byond significant proteins. For each time-slice (two adjacent timepoints), a Steiner tree is created, which are then stitched together to a temporal network. Paths are extracted. Evaluation of the network via edge confidence scores and kinase-significant protein edges is performed.

The release from 28/4/23 is associated with the BioRxiv submission "Network inference from temporal phosphoproteomics informed by protein-protein interactions"
