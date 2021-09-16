def.quadprog <- function(qp_pars){
  qp_pars
}

update_quadprog <- function(Dmat, dvec, Cmat, bvec, qp_pars){
  res <- quadprog::solve.QP(Dmat, dvec, t(Cmat), bvec = bvec)
  active <- rep(FALSE, length(res$solution))
  active[res$iact] <- TRUE
  list(alpha = res$solution, active = active)
}