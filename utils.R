# localization operator
localization.op <- function(evalues, evectors, g, i=NULL){ # support of T_i g is centered at node i
  if(is.null(i)){
    res <- evectors %*% (t(Conj(evectors))*g(evalues)) # ith column: T_i g
  } else{
    res <- as.vector(evectors %*% (Conj(evectors)[i,]*g(evalues))) # T_i g vector when i is specified
  }
  return(res)
}

# graph windowed Fourier transform
graph.window.FT <- function(x, S, g, M){
  # eigenres <- eigen(S)
  # evalues <- eigenres$values
  # evectors <- eigenres$vectors
  # lmax <- max(evalues)
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  lmax <- max(evalues)
  
  C <- NULL
  normconst <- c()
  for(m in 1:M){
    gm <- function(lambda){
      return(g(lambda, sigma.sq=lmax*(M+1)/M^2, m, tau=lmax*(M+1)/M^2))
    }
    Tgm <- localization.op(evalues, evectors, gm)
    C <- cbind(C, t(Conj(Tgm))%*%x) # where C[i,m] = <x, Ti gm>
    normconst <- c(normconst, norm(Tgm, type="F")^2)
  }
  return(list(C=C, normconst=normconst))
}

# cpsd estimate
cpsd.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # x, y should be R realizations (N x R matrix)
    # eigenres <- eigen(S)
    # x.tilde <- t(Conj(eigenres$vectors))%*%x
    # y.tilde <- t(Conj(eigenres$vectors))%*%y
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    x.tilde <- t(Conj(evectors))%*%x
    y.tilde <- t(Conj(evectors))%*%y
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  } else if(method=="window"){
    # x, y should be one realization (N x 1 vector)
    C1 <- graph.window.FT(x, S, g, M)
    C2 <- graph.window.FT(y, S, g, M)
    cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
  } else if(method=="random"){
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    x.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(x))
    y.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(y))
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  }
  return(cpsd)
}

cpsd.graph.fast <- function(x, y, evalues, evectors, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # x, y should be R realizations (N x R matrix)
    # eigenres <- eigen(S)
    # x.tilde <- t(Conj(eigenres$vectors))%*%x
    # y.tilde <- t(Conj(eigenres$vectors))%*%y
    x.tilde <- t(Conj(evectors))%*%x
    y.tilde <- t(Conj(evectors))%*%y
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  } else if(method=="window"){
    # x, y should be one realization (N x 1 vector)
    C1 <- graph.window.FT(x, S, g, M)
    C2 <- graph.window.FT(y, S, g, M)
    cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
  } else if(method=="random"){
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    x.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(x))
    y.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(y))
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  }
  return(cpsd)
}


# psd estimate
psd.graph <- function(x, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  return(cpsd.graph(x=x, y=x, S=S, g=g, M=M, sigma=sigma, method=method, seed=seed))
}

# coherence estimate
coherence.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, seed=seed)
  } else if(method=="window"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, method=method, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, method=method, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, method=method, seed=seed)
  } else if(method=="random"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, M=M, sigma=sigma, method=method, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, M=M, sigma=sigma, method=method, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, M=M, sigma=sigma, method=method, seed=seed)
  }
  return(cpsd*Conj(cpsd) / psd.x / psd.y)
}

# cross spectrum analysis All in one
cross.spectrum.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  return(list(cpsd=cpsd, psd.x = psd.x, psd.y = psd.y, coherence=cpsd*Conj(cpsd) / psd.x / psd.y))
}

graph.stationary.level <- function(cov, S){
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  gamma <- t(Conj(evectors))%*%cov%*%evectors
  return(norm(diag(gamma), type=2) / norm(gamma, type="F"))
}

# robust cpsd estimate 
huberloss <- function(t, c){
  res <- c()
  for(i in t){
    if(abs(i) <= c){
      res <- c(res, i^2)
    } else{
      res <- c(res, 2*c*(abs(i)-c/2))
    }
  }
  return(res)
}

# random window generation
windowbank.random <- function(N, M, V, sigma, seed=1){
  res <- matrix(0, nrow=M, ncol=N)
  set.seed(seed)
  for(i in 1:M){
    W.tilde <- diag(N) + rnorm(N^2, 0, sigma)
    res[i,] <- diag(V %*% W.tilde %*% t(Conj(V)))
  }
  return(res)
}


# cpsd estimate
# robust.cpsd.graph.window <- function(x, y, S, g, M){
#   C1 <- graph.window.FT(x, S, g, M)
#   C2 <- graph.window.FT(y, S, g, M)
#   cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
#   return(cpsd)
# }



plot_graph_custom <- function (z, size = 0.75, edge_color, vertex_color=NULL) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$color <- factor(edge_color, levels = unique(edge_color))
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, colour = color), 
                                              linewidth = 2, data = df2) +
    scale_color_manual(values=color.cand, labels = paste("Line", 1:8, sep=""),
                       name = "Line number") + 
    geom_point(aes(fill=vertex_color), size = size, shape=21) + 
    scale_fill_gradient(low="white", high="black", na.value = "yellow", name = "People") +
    theme_void() +
    theme(legend.margin = margin(10,10,10,10))
  print(p1)
}

plot_graph_custom3 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  # if(signal==FALSE){
  #   p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
  #                                                   xend = y1, yend = y2),  linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
  #     geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
  #   print(p1)
  # }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w),  linewidth=e.size, data = df2) + ggtitle(title) +
      scale_color_gradient(low="lightgray", high="black", name="Edge Weight", guide="none") + 
      geom_point(size = v.size, fill="white", colour="black", shape=21, stroke=1.2) + theme_void()+
      geom_label_repel(aes(label = w), size = 3, max.overlaps = Inf, box.padding = 0.3, segment.color = "gray50")+
      theme(plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5),aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom4 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom5 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2, show.legend=FALSE) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(color=guide_colourbar(order=1)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}



gFreqReg <- function(X, Y, S, evalues, evectors, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1, lambda=0){
  if(is.null(method)){
    # X should be n x p x R array (R realizations)
    # Y should be n x R array
    # p-dimensional graph signal
    if(length(dim(X))!=3){
      stop("X should be 3-dimensional n x p x R array!")
    }
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    R <- dim(X)[3]
    
    if(n!=nrow(S)){
      stop("The number of rows of X[,,r] should be equal to the number of rows of GSO!")
    }
    if(nrow(Y)!=nrow(S)){
      stop("The number of rows of Y should be equal to the number of rows of GSO!")
    }
    if(R!=ncol(Y)){
      stop("The number of columns of X[,,r] should be equal to the number of columns of Y!")
    }
    
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    out_X <- X # n x p x R
    out_Y <- array(0, c(n, 1, R)) # n x 1 x R
    out_Y[,1,] <- Y 
    
    
    print("calculating P_array_X...")
    P_array_X <- array(NA, dim = c(n,p,p))
    
    # apply evectors to all slices simultaneously
    out_X.tilde <- t(Conj(evectors)) %*% matrix(out_X, nrow=n, ncol=p*R)
    
    # reshape back into n x p x m
    out_X.tilde <- array(out_X.tilde, dim = c(n, p, R))
    
    print("calculating P_array_Y...")
    P_array_Y <- array(NA, dim = c(n,1,1))
    # apply evectors to all slices simultaneously
    out_Y.tilde <- t(Conj(evectors)) %*% matrix(out_Y, nrow=n, ncol=1*R)
    
    # reshape back into n x p x m
    out_Y.tilde <- array(out_Y.tilde, dim = c(n, 1, R))
    
    print("calculating P_array_XY...")
    P_array_XY <- array(NA, dim = c(n,p,1))
    
    rho.hat <- vector(length=n)
    mu.tilde.hat <- vector(length=n)
    beta.hat <- matrix(0,nrow=n, ncol=p)
    y.tilde.hat <- array(0, c(n, R)) # n x 1 x R
    for (k in 1:n) {
      if(p==1){
        P_array_X[k,,] <- sum(out_X.tilde[k,,] * Conj(t(out_X.tilde[k,,]))) / R
      }
      else{
        P_array_X[k,,] <- out_X.tilde[k,,] %*% Conj(t(out_X.tilde[k,,])) / R
      }
      
      P_array_Y[k,,] <- sum(out_Y.tilde[k,,] * Conj(out_Y.tilde[k,,])) / R
      P_array_XY[k,,] <- out_X.tilde[k,,] %*% Conj(out_Y.tilde[k,,]) / R
      rho.hat[k] <- as.numeric(t(P_array_XY[k,,]) %*%  ginv(P_array_X[k,,]+diag(lambda,p)) %*% P_array_XY[k,,] / P_array_Y[k,,])
      beta.hat[k,] <- ginv(P_array_X[k,,]+diag(lambda,p)) %*% Conj(P_array_XY[k,,])
      mu.tilde.hat[k] <- mean(out_Y.tilde[k,1,]) - sum(beta.hat[k,]*out_X.tilde[k,,]) / R
      
      if(p==1){
        y.tilde.hat[k,] <- mu.tilde.hat[k] + beta.hat[k,]*out_X.tilde[k,,]
      } else{
        y.tilde.hat[k,] <- mu.tilde.hat[k] + t(out_X.tilde[k,,]) %*% beta.hat[k,]
      }
    }
    
    Fstat <- rho.hat * (R-p-1) / (p*(1-rho.hat))
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$mu.tilde.hat <- mu.tilde.hat ; res$beta.hat <- beta.hat ; res$rho.hat <- rho.hat
    res$Px <- P_array_X ; res$Py <- P_array_Y ; res$Pxy <- P_array_XY
    res$Y.tilde.hat <- y.tilde.hat ; res$Fstat <- Fstat
    return(res)
  } 
  
  else if(method=="random"){
    # X should be n x p array (1 realization)
    # Y should be n-dimensional vector
    # p-dimensional graph signal
    if(length(dim(X))!=2){
      stop("X should be 2-dimensional n x p array!")
    }
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    if(n!=nrow(S)){
      stop("The number of rows of X should be equal to the number of rows of GSO!")
    }
    if(length(Y)!=nrow(S)){
      stop("The length of Y should be equal to the number of rows of GSO!")
    }
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    
    
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    
    m <- nrow(WB)
    out_X <- array(0, c(nrow(X), ncol(X), m)) # n x p x m
    out_Y <- array(0, c(length(Y), 1, m)) # n x 1 x m
    
    
    for (i in 1:m) {
      out_X[,,i] <- sweep(X, 1, WB[i, ], `*`)
      out_Y[,,i] <- sweep(as.matrix(Y), 1, WB[i, ], `*`)
    }
    
    
    print("calculating P_array_X...")
    P_array_X <- array(NA, dim = c(n,p,p))
    
    
    # apply evectors to all slices simultaneously
    out_X.tilde <- t(Conj(evectors)) %*% matrix(out_X, nrow=n, ncol=p*m)
    
    # reshape back into n x p x m
    out_X.tilde <- array(out_X.tilde, dim = c(n, p, m))
    
    print("calculating P_array_Y...")
    P_array_Y <- array(NA, dim = c(n,1,1))
    # apply evectors to all slices simultaneously
    out_Y.tilde <- t(Conj(evectors)) %*% matrix(out_Y, nrow=n, ncol=1*m)
    
    # reshape back into n x p x m
    out_Y.tilde <- array(out_Y.tilde, dim = c(n, 1, m))
    
    print("calculating P_array_XY...")
    P_array_XY <- array(NA, dim = c(n,p,1))
    
    rho.hat <- vector(length=n)
    mu.tilde.hat <- vector(length=n)
    beta.hat <- matrix(0,nrow=n, ncol=p)
    y.tilde.hat <- vector(length=n) 
    
    X.tilde <- t(Conj(evectors)) %*% X
    for (k in 1:n) {
      if(p==1){
        P_array_X[k,,] <- sum(out_X.tilde[k,,] * Conj(t(out_X.tilde[k,,]))) / m
      }
      else{
        P_array_X[k,,] <- out_X.tilde[k,,] %*% Conj(t(out_X.tilde[k,,])) / m
      }
      
      P_array_Y[k,,] <- sum(out_Y.tilde[k,,] * Conj(out_Y.tilde[k,,])) / m
      P_array_XY[k,,] <- out_X.tilde[k,,] %*% Conj(out_Y.tilde[k,,]) / m
      rho.hat[k] <- as.numeric(t(P_array_XY[k,,]) %*%  ginv(P_array_X[k,,]+diag(lambda,p)) %*% P_array_XY[k,,] / P_array_Y[k,,])
      beta.hat[k,] <- ginv(P_array_X[k,,]+diag(lambda,p)) %*% Conj(P_array_XY[k,,])
      mu.tilde.hat[k] <- mean(out_Y.tilde[k,1,]) - sum(beta.hat[k,]*out_X.tilde[k,,]) / m
      y.tilde.hat[k] <- mu.tilde.hat[k] + sum(X.tilde[k,] * beta.hat[k,])
    }
    
    Fstat <- rho.hat * (m-p-1) / (p*(1-rho.hat))
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$mu.tilde.hat <- mu.tilde.hat ; res$beta.hat <- beta.hat ; res$rho.hat <- rho.hat
    res$Px <- P_array_X ; res$Py <- P_array_Y ; res$Pxy <- P_array_XY
    res$Y.tilde.hat <- y.tilde.hat ; res$Fstat <- Fstat
    return(res)
  }
}


