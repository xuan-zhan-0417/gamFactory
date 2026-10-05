# Cross blocks of the Hessian between two effects (o1 horizontal, o2 vertical): crossprod(J2, llk$d2 * J1),
# where J is d eta / d(parameters) of the effect. Called by DllkDbeta.linpreds as .Hess.<class o1>_<class o2>.

.Hess.stand_stand <- function(o1, o2, llk){
  crossprod(o2$store$X, llk$d2 * o1$store$X)
}

# J of a nested effect: (df/dz_k) dz_k/dalpha_k for each margin k with parameters, then X0
.jac_nest <- function(o) cbind(do.call("cbind", Map("*", o$store$f1, o$store$g1)), o$store$X0)

.Hess.nest_stand <- function(o1, o2, llk){
  crossprod(o2$store$X, llk$d2 * .jac_nest(o1))
}
.Hess.stand_nest <- function(o1, o2, llk){
  t( .Hess.nest_stand(o1 = o2, o2 = o1, llk) )
}

.Hess.nest_nest <- function(o1, o2, llk){
  crossprod(.jac_nest(o2), llk$d2 * .jac_nest(o1))
}
