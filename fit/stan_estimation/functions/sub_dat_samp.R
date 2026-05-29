sample_indices <- function(n, S, P) {

  xx = matrix(1:n, ncol = S)
  nc = nrow(xx)
  
  y = as.vector(apply(xx, 2, sample, size = P*nc))
  
  return(sort(y))
}

# n = d_stan$n; S = d_stan$N; P = .5
