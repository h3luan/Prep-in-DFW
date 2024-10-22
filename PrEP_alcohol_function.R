## Create spatial neighbors: account for islands

addNB <- function(shp){
  shp_neighbs <- poly2nb(shp,queen=T)
  count <- card(shp_neighbs)
  if(!any(count==0)){
    return(shp_neighbs)
  }
  # nnbs <- knearneigh(coordinates(shp))$nn
  nnbs <- knearneigh(st_centroid(st_geometry(shp)))$nn
  islands <- which(count==0)
  for(i in islands){
    shp_neighbs[[i]] <- nnbs[i]
    shp_neighbs[[nnbs[i]]] <- c(shp_neighbs[[nnbs[i]]],i)
  }
  return(shp_neighbs)
}


## Create the matrix in INLA for Leroux prior

leroux_matrix <- function(g.path){
  g <- inla.read.graph(g.path)
  num <- g$n
  R <- matrix(0,num,num)
  for(i in 1:num){
    R[i,i] <- g$nnbs[[i]]
    R[i,g$nbs[[i]]] <- -1
  }
  R.Leroux <- diag(length(g$nbs))-R
  return(R.Leroux)
}

## Correlation analysis
corr.matrix <- function(x){
  num <- dim(x)[2]
  corr.m <- matrix(NA,num,num)
  for(i in 1:num){
    for(j in 1:num){
      corr.m[i,j] <- rcorr(x[,i],x[,j],"pearson")$r[1,2]
    }
  }
  return(corr.m)
}

## Posterior predictive probability checking
ppp_check <- function(obs,model){
  ppp <- vector(mode = "numeric", length = length(obs))
  for(i in 1:length(ppp)){
    ppp[i] <- inla.pmarginal(q = obs[i],marginal = model$marginals.fitted.values[[i]])
  }
  return(ppp)
}


# Testing zero-inflation --------------------------------------------------

ZI.test <- function(fits){
  # obs <- fits$y
  obs <- fits$model[,1]
  N <- length(obs)
  mean.obs <- mean(obs)
  fitted <- fits$fitted.values
  p.zero.fitted <- exp(-fitted)
  a <- sum((as.numeric(obs==0)-p.zero.fitted)/p.zero.fitted)^2
  b <- sum((1-p.zero.fitted)/p.zero.fitted)-N*mean.obs
  test.stat <- a/b
  p.value <- 1-pchisq(test.stat,1)
  res <- c(test.stat,p.value)
  names(res) <- c("Test Statistic","p-value")
  return(res)
}

## Leroux cumulative
cum_fun <- function(num){
  cum <- vector()
  cum[1] <- 0
  for(i in 1:length(num)){
    cum[i+1] <- sum(num[1:i])
  }
  return (cum)
}


# Plot figures ------------------------------------------------------------

plot_prep <- function(layer,field,class,breaks,color,title,labels){
  plot <- tm_shape(layer)+
    tm_fill(col=field,palette=color,title="",breaks=breaks,labels = labels)+
    tm_borders(lwd = 0.3)+
    tm_scale_bar(position = c("right","bottom"),width = 0.2)+
    tm_compass(type = "arrow",position = c("left","top"))+
    tm_layout(main.title = title,main.title.position = "center",frame = F,legend.position = c(0,0.5),legend.text.size = 0.4,legend.height = 0.4,legend.title.size = 0.5,main.title.size = 0.6,main.title.fontface = "bold")
  return(plot)
}
