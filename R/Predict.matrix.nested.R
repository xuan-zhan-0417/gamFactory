#'
#' Predict using nested effects
#' 
#' @name Predict.matrix.nested
#' @rdname Predict.matrix.nested
#' @export
#'
#'
Predict.matrix.nested <- function(object, data, ...){
  
  if(class(object)[1] == "si"){
    return( .predict.matrix.si(object, data, ...) ) 
  }
  if(class(object)[1] == "nexpsm"){
    return( .predict.matrix.nexpsm(object, data, ...) ) 
  }
  if(class(object)[1] == "mgks"){
    return( .predict.matrix.mgks(object, data, ...) ) 
  }
  if(class(object)[1] == "si_nexpsm"){
    return( .predict.matrix.si_nexpsm(object, data, ...) )
  }
  if(class(object)[1] == "inter_linear"){
    return( .predict.matrix.inter_linear(object, data, ...) )
  }
  if(class(object)[1] == "inter_nexp"){
    return( .predict.matrix.inter_nexp(object, data, ...) )
  }
  if(class(object)[1] == "inter_mgks"){
    return( .predict.matrix.inter_mgks(object, data, ...) )
  }
  
  stop("Predict.matrix.nested --- I do not know this effect type")
  
}