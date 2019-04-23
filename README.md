# The-GMD-biplot
This repository contains two files: "GMD_source.R" contains main functions for creating the GMD-biplot and the GMD scree plot; "examples.R" gives a simple example on how to apply the main functions based on a simulated data set. Usages of the main functions in the "GMD_source.R" are described below. 

1. GMD(X, H, Q, K): 
Input: 
(1) X: an nXp matrix; (2) H: an nXn positive semi-definite matrix; (3) Q: a pXp positive semi-definite matrix; (4) K: the dimension of the GMD
Output: a GMD object. 
(1) U: the nXK left GMD vectors; (2) V: the pXK right GMD vectors; (3) D: the KX1 vector of GMD values.

2. biplot.GMD(fit, index, names, sample.col, sample.pch):
Input: 
(1) fit: a GMD object.
(2) index: optional; a vector containing the index of the variables that are plotted in the GMD-biplot. The default is to plot all variables.
(3) names: optional; a vector containing the names of the variables that are plotted in the GMD-biplot. The default names for the j-th variable is "Vj".
(4) sample.col: optional; the color of the sample points in the GMD-biplot.
(5) sample.pch: optional; the symbol of the sample points in the GMD-biplot.
Output: the GMD-biplot

3. screeplot.gmd(fit, K)
Input: 
(1) fit: a GMD object.
(2) K: the number of components plotted in the scee plot.
Output: the GMD scree plot




