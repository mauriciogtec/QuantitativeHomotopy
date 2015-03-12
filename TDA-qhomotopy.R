# LET'S GENERATE A CIRCLE AND ITS SUBDIVISIONS
circle.pts <- function(n.sub, size=2^n.sub){
  # INPUT:
  # n.sub = number of subdivisions
  # size = sample size
  # OUTPUT:
  # "size" points chosen from "2^n.sub" homogenuous points of the circle
  theta <- seq(0, 2*pi, length.out=2^n.sub+1)[-(2^n.sub+1)]
  coords <- data.frame(x=cos(theta), y=sin(theta))
  return(coords[sample(2^n.sub, size) ,])
}
plot(circle.pts(5))

# Now points from the torus!
torus.pts <- function(n.sub, size=2^n.sub, noise=0){
  # INPUT:
  # n.sub = number of subdivisions
  # size = sample size
  # noise = white gaussian noise to points
  # OUTPUT:
  # "size" points chosen from "2^n.sub" homogenuous points of the circle
  N <- 2^n.sub
  theta <- runif(N,0,2*pi)
  phi <- runif(N,0,2*pi)
  x = (2 + cos(theta))*cos(phi) + rnorm(N, 0, noise)
  y = (2 + cos(theta))*sin(phi) + rnorm(N, 0, noise)
  z = sin(theta)  + rnorm(N, 0, noise)
  coords <- data.frame(x, y, z)
  return(coords[sample(2^n.sub, size) ,])
}
library(rgl)
library(car)
coords <- torus.pts(11, noise=.1)
scatter3d(coords[ ,1], coords[ ,2], coords[ ,3],surface=FALSE,axis.scales=TRUE, point.col = "black")


genustwo.pts <- function(n.sub, noise=0){
  # INPUT: see torus.pts
  # OUTPUT:
  # "size" points chosen from "2^(n.sub+1)" homogenuous points of the circle
  N <- 2^n.sub
  theta <- runif(N,0,2*pi)
  phi <- runif(N,0,2*pi)
  x = (2 + cos(theta))*cos(phi) + rnorm(N, 0, noise)
  y = (2 + cos(theta))*sin(phi) + rnorm(N, 0, noise)
  z = 2*sin(theta)  + rnorm(N, 0, noise)
  coords <- coords2 <- data.frame(x, y, z)
  coords$x <- coords$x - 2
  coords2$x <- coords2$x + 2
  coords <- coords[coords$x < 0, ]
  coords2 <- coords2[coords2$x > 0, ]
  coords <- rbind(coords, coords2)
  return(coords)
}
coords <- genustwo.pts(9, noise=0)
library(scatterplot3d)
scatterplot3d(coords[ ,1], coords[ ,3], coords[ ,2], angle=80,
              highlight.3d = TRUE, col.axis = "blue",              
              col.grid = "lightblue", main = "", pch =20)



d.test <- function(i1, i2, functions, land){
  f1 <- functions[[i1]]
  f2 <- functions[[i2]]
  return(max(abs(land[f2, ]-land[f1, ])))
}

distance.matrix <- function(test, land, nsim=0, max.funs=100){
    require(plyr)
    # 1) create nsim functions
    N.test <- dim(test)[1]
    N.land <- dim(land)[1]
    # Gotta be carefull with the HUGE number of possible functions!!!!!!
    if (nsim==0) nsim <- N.land^N.test
    nsim <- min(nsim, max.funs)
    # each element of this list will be a function
    functions <- lapply(1:nsim, function(x) sample.int(n=N.land, size=N.test, replace=TRUE))
    combs <- combn(nsim, 2)
    d.matrix <- matrix(0, nrow=nsim, ncol=nsim)
    distances <- apply(combs, 2, function(x)  d.test(x[1], x[2], functions, land))
    for(x in 1:ncol(combs)){
      d.matrix[combs[1, x], combs[2, x]] <- distances[x] 
      d.matrix[combs[2, x], combs[1, x]] <- distances[x]
    }
    return(d.matrix)
}


# Now let's see the persistent homology of this object!

d1 <- distance.matrix(test=circle.pts(5), land=genustwo.pts(5, noise=0), max.funs=500)
library(phom)
max_dim <- 2
max_f <- 9
intervals <- pHom(d1, max_dim, max_f, metric="distance_matrix")
plotPersistenceDiagram(intervals, max_dim, max_f)
