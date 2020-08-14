# create Z matrix
p <- apply(M, 2, mean)/2
P <- matrix(rep(p*2, nrow(M)), ncol=ncol(M), nrow=nrow(M), byrow=TRUE)
Z <- M - P
# sum 2pq to scale G to the A matrix
q <- 1 - p
sum2pq <- 2*sum(p*q)  # note you can pull out the 2 (redundant)
print(sum2pq)
# calculate G!
G <- (Z %*% t(Z)) / sum2pq
print(round(G, 3))
