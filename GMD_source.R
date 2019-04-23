#---------------------------------------------------------------------------------------------------------
# This piece of code includes functions for GMD, the GMD-biplot,  the GMD-screeplot
# Yue Wang: 20190422
#---------------------------------------------------------------------------------------------------------

get_uv = function(X, H, Q, u_0, v_0){
  
  u = X%*%Q%*%v_0/as.numeric(sqrt(t(v_0)%*%t(Q)%*%t(X)%*%H%*%X%*%Q%*%v_0))
  v = t(X)%*%H%*%u/as.numeric(sqrt(t(u)%*%t(H)%*%X%*%Q%*%t(X)%*%H%*%u))

  return(list(u = u, v = v))
    
}


GMD =function(X, H, Q, K){
  
  n = dim(X)[1]
  p = dim(X)[2]
  # output matrix/vec
  U = matrix(0, n, K)
  V = matrix(0, p, K)
  D = rep(0,K)
  
  X_0 = X
  u_0 = c(1, rep(0,n-1))
  v_0 = c(1, rep(0,p-1))
  
  for(iter in 1:K){
    
    
    error = 1
 
    while(error > 1e-5){
      
      temp.uv = get_uv(X_0, H, Q, u_0, v_0)
      
      u = temp.uv$u
      v = temp.uv$v
      
      error = norm(u - u_0, "2") + norm(v - v_0, "2")
      #print(error)
      
      u_0 = u
      v_0 = v
   
    }
    
    U[,iter] = u
    V[,iter] = v
    d = t(u)%*%H%*%X_0%*%Q%*%v
    D[iter] =  d
    
    X_0 = X_0 - u%*%d%*%t(v)
  
  }
 
  return(list(U = U, V = V, D = D, H = H, Q = Q, X = X))
  
}



#-----------------------------------------------------
# screeplots for GMD
# fit: a GMD object
# K: the number of component that is plotted in the scree plot; the default of K is min(10,p).
#-----------------------------------------------------
screeplot.gmd = function(fit, K){
  
  p = dim(fit$V.s)[1]
  
  if(missing(K)){K = min(10,p)}
  
  D = fit$D
  X = fit$X
  H = fit$H
  Q = fit$Q
  
  D.plot = D^2/sum(diag(t(X)%*%H%*%X%*%Q))
  
  plot(1:K, D.plot[1:K], xaxt = 'n', ylab = 'Percentage of variance explained', type = 'b', xlab = ' ', ylim = c(0,ceiling(D.plot[1]*10)/10), cex.axis = 1.2, cex.lab = 1.2)
  axis(1,1:K, paste0("PC", 1:K))
  
}


#---------------------------------------------------
# GMD biplots (unsupervised)
# fit: a GMD object
# index (optional): a vector indicating which variable you want to plot in the biplots; the default is all variables.
# names (optional): a vector of variable names which has the same length as the index; if not provided, the jth variable will be labeled as "Vj". 
# sample.col (optional): sample colors ;
# sample.pch (optional): sample symbols (optional);
#---------------------------------------------------
biplot.gmd = function(fit, index, names, sample.col, sample.pch){
  
  U = fit$U
  D = fit$D
  V = fit$V
  
  U = U[,order(D, decreasing = T)]
  V = V[,order(D, decreasing = T)]
  
  k1 = order(D, decreasing = T)[1]
  k2 = order(D, decreasing = T)[2]
  
  D = sort(D, decreasing = T)
  
  eta = U%*%diag(D)
  
  
  max.xlab = max(abs(eta[,1]))
  max.ylab = max(abs(eta[,2]))  
  
  if(missing(sample.col)){sample.col = 'grey50'; arrow.col = 'lightgreen'; legend.col = 'blue'}
  if(!missing(sample.col)){arrow.col = 'grey50'; legend.col = 'black'}
  if(missing(sample.pch)){sample.pch = 19}
  
  plot(eta[,1], eta[,2], xlab = paste0('PC',k1), ylab = paste0('PC',k2), pch = sample.pch, xlim = c(-1.1*max.xlab, 1.1*max.xlab ), ylim =  c(-1.1*max.ylab, 1.1*max.ylab) , col = sample.col, cex.axis = 1.2, cex.lab = 1.2 )
  xaxp = axTicks(1)
  yaxp = axTicks(2)

  
  if(missing(index)){index = 1:dim(V)[1]}
  # only plot these user-specified variables
  
  #calculate coordinates
  Q = fit$Q
  V.plot = Q%*%V
  arrow.x = V.plot[,1]
  arrow.y = V.plot[,2]
  
  # original code (equivalent!)
  #arrow.x = diag(rep(1, dim(V)[1]))%*%Q%*%V[,1] 
  #arrow.y = diag(rep(1, dim(V)[1]))%*%Q%*%V[,2] 
  
  if(missing(names)){names = paste0("V", index)}
  iter = 1
  
  max.xarrow = max(abs(arrow.x))
  max.yarrow = max(abs(arrow.y))
  xratio = max.xarrow/max.xlab
  yratio = max.yarrow/max.ylab

  
  xsci = as.numeric(unlist(strsplit(formatC(xratio, format = 'e'),"e")))
  xlab.arrow = round(xaxp*xsci[1]*10^(xsci[2]), digits = 2)

  
  ysci = as.numeric(unlist(strsplit(formatC(yratio, format = 'e'),"e")))
  ylab.arrow = round(yaxp*ysci[1]*10^(ysci[2]), digits = 2)
  
  for(i in index){
    
   # if(arrow.x[i]^2 + arrow.y[i]^2 >= 0.1^2){
      
      arrows(x0 = 0,y0 = 0,x1 = arrow.x[i]/xratio, y1 = arrow.y[i]/yratio, length = 0.05, col = arrow.col)
      if(!missing(names)){
        text(arrow.x[i]/xratio, arrow.y[i]/yratio*1.1, names[iter], cex = 1, col = legend.col)
      }
      
    #}
    
    iter = iter + 1
  }
  
  # add new axis
  axis(3, at = xaxp, labels = as.character(xlab.arrow), cex.axis = 1.2)
  axis(4, at = yaxp, labels = as.character(ylab.arrow), cex.axis = 1.2)
  
  points(eta[,1], eta[,2], col = sample.col, pch = sample.pch)
}

