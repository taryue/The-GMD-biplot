# An example of constructing the GMD-biplot and the corresponding scree plot
source('GMD_source.R')

library(MASS)
# simulate data matrix
X.org = matrix(rnorm(1000), 20, 50)
autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}
H.sqrt.n = autocorr.mat(20, 0.6)
X = H.sqrt.n%*%X.org
H = ginv(H.sqrt.n%*%H.sqrt.n)

labels = (rowSums(X) > 0) # simulate labels

#--------------------------------------------
# step1: GMD
GMD.fit = GMD(X, H, diag(rep(1,50)), 10) # Q is identity; we only calculate top 10 components
# step 2: biplot
par(mar = c(4.1, 4.1, 2.1, 2.1)) # set your graphic parameters
plot.index = 1:5 # choose the variables that you want to display
plot.names = paste0("V",plot.index) # add names for the variables
biplot.gmd(fit = GMD.fit, index = plot.index, names = plot.names, sample.col = labels + 1, sample.pch = 19)

# step 3: screeplot
screeplot.gmd(GMD.fit, 10) # plot top 10 components

