#-----------------------------
# Stat 322
# Name: David Cao
# Instructor: Cristina Anton
# Markov Project
#-----------------------------
install.packages("markovchain","diagram","diagrammeR")

require(markovchain)
library(expm)
library(DiagrammeR)
library(markovchain)
library(diagram)
library(matlab)

# Find Matrix T
getRT <- function(M,type="T"){
  if(length(absorbingStates(M)) == 0) stop("Not Absorbing Matrix")
  tm <- M@transitionMatrix
  d <- diag(tm)
  m <- max(which(d == 1))
  n <- length(d)
  ifelse(type=="T",
         A <- tm[(m+1):n,(m+1):n],
         A <- tm[(m+1):n,1:m])
  return(A)
}

byRow <- TRUE
States = c('0','1','2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'10','W' ,'L' ,'A6','A7','A8','A9','A10')
data = c(1/13,1/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,  0 ,  0 ,  0 ,  0 ,
           0 ,1/13,1/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,  0 ,  0 ,  0 ,
           0 ,  0 ,1/13,1/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,  0 ,  0 ,
           0 ,  0 ,  0 ,1/13,1/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,  0 ,
           0 ,  0 ,  0 ,  0 ,1/13,1/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,
           0 ,  0 ,  0 ,  0 ,  0 ,1/13,1/13,2/13,2/13,2/13,4/13,1/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,2/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,2/13,2/13,2/13,2/13,4/13,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,2/13,2/13,2/13,6/13,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,2/13,2/13,8/13,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,1/13,2/13,10/13, 0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,4/13,  0 ,1/13,2/13,2/13,2/13,2/13,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,4/13,  0 ,  0 ,  0 ,2/13,  0 ,  0 ,1/13,2/13,2/13,2/13,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,2/13,4/13,  0 ,  0 ,2/13,  0 ,  0 ,  0 ,1/13,2/13,2/13,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,2/13,2/13,4/13,  0 ,2/13,  0 ,  0 ,  0 ,  0 ,1/13,2/13,
           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,2/13,2/13,2/13,4/13,2/13,  0 ,  0 ,  0 ,  0 ,  0 ,1/13)
Matrix <- matrix(data , byrow = byRow, nrow = 18,dimnames = list(States, States))

mc <- new("markovchain", states = States, byrow = byRow, transitionMatrix = Matrix, name = "BlackJack")
show(mc)

# Determine transient states
transientStates(mc)

# determine absorbing states
absorbingStates(mc)  

#Put the transition matrix into Canonical Form
P <- canonicForm(mc)
P
T <- getRT(P)
T

# Find Fundamental Matrix
I <- diag(dim(T)[2])
I
N <- solve(I - T)
N
write.table(N, 'clipboard', sep='\t')

#Calculate time to absorption
e <- rep(1,dim(N)[2])
m <- N %*% e
m
R <- getRT(P,"R")
R
A <- N %*% R
A

#Draw a graphical representation of the transition matrix
require(igraph)
mcIgraph <- as(mc, "igraph")
plot.igraph(mcIgraph, vertex.color="green" )
mcIgraph

#Run multiple simulation of the markov chain and plot the sequence states until the absorption state
MCts <- rmarkovchain(n=10, object=mc)
MCts
MCtsDf <- as.data.frame(MCts,stringsAsFactors = TRUE)
MCtsDf
MCtsDf$index <- 1:10
MCtsDf$MCts <- as.numeric(MCtsDf$MCts)
require(ggplot2)

library(ggplot2)
p <- ggplot(MCtsDf,aes(index,MCts))
p  +geom_point(colour="dark red")+
  xlab("time") +
  ylab("state") +
  ggtitle("Markov Chain")
