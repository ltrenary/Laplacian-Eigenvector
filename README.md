# Laplacian-Eigenvector
In geosciences, we often want/need to reduce the dimensionality of the data. One way of doing this is to perform an EOF analysis
on the field of interest. Unfortunately, EOF analysis can yield different results for different datasets, which makes direct comparison 
between climate model and observational records challenging. One way to circumvent this is to project the data onto a common basis set. 
The code in this directory projects a 2-dimensional field onto the eigenvectors of the Laplacian operator. Codes are written in MATLAB. Example of file also provided.
