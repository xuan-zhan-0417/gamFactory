############################
# Returns list of effects: standard (eff_stand) or nested (eff_nest)
#
.build_effects <- function(X, info, outer){
  lapply(seq_along(info$type), function(ii){
    if( info$type[[ii]][1] == "stand" ){
      eff_stand(X[ , info$iec[[ii]], drop = FALSE])
    } else {
      eff_nest(margin = info$extra[[ii]]$si$margin, basis = info$extra[[ii]]$basis)
    }
  })
}
