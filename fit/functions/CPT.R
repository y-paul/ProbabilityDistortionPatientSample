# value fun
v <- function(x, l, a = 1) {
  L = ifelse(x < 0, l, 1) # set loss aversion
  vx = sign(x) * L * abs(x)^a 
  return(vx)
}

# pwf
pwf_prtk <- function(p, g) {
  (p^g) / (  (  p^g  + ((1-p)^g)  )^(g^(-1))   ) 
}

# choice rule
softmax <- function(va, vb, t) { 1 / ( 1 + exp(-t * (va - vb))) }

cpt_mod_pr <- function(XA, PA, XB, PB, l, g, t,
                       out = 'pa') {
  
  # subjective values
  vXA = v(XA, l = l); vXB = v(XB, l = l)
  
  # decision weights
  wPA = pwf_prtk(PA, g = g); wPB = pwf_pr(PB, g = g)

  # subjective values
  va = sum(vXA * wPA)
  vb = sum(vXB * wPB)
  
  if(out == 'ev') return( cbind(va = va, vb = vb))
  
  
  # choice A probability
  pa <- softmax(va, vb, t = t)
  
  if(out == 'pa') return(pa)
  
  
  # predicted choice
  co <- rbinom(length(pa), 1, pa)
  
  if(out == 'co') return(co)
  
  
}
