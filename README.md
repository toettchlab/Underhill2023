# Underhill2023
Movies and quantification code for paper "Control of gastruloid patterning and morphogenesis by the Erk and Akt signaling pathways," Development (2023)

This code is used for quantifying the distribution of a fluorescent readout along the major axis of a gastruloid. 
The Matlab script takes as input a .tif stack with individual gastruloid images, along with .txt files containing the coordinates
of the edges and major axis for each gastruloid. 
The code outputs a plot of the fluorescence profile for each individual gastruloid, along with an aggregated plot
containing the mean and standard deviation of the gastruloids in a batch.
