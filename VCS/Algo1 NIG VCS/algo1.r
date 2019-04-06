
Uniform <- function() {
  #set.seed(seed);
    U = runif(N,0,1)
    zt <- rep(0,N)
    for (i in 1:N) {
    if (U[i] < dt/(dt+gamma*zeta[i])) {
      zt[i] = zeta[i] }
    else  {
      zt[i] = ((delta^2)*(t^2)) / (gamma^2*zeta[i]) }}
    return(zt)} 