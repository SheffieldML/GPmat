modelOut <- function(model, X, Phi.return=FALSE, ...) {
  fhandle <- get(paste(model$type, "Out", sep=""), mode="function")
  #fhandle = str2func(paste(model$type,"Out",sep=""))
  
  if (Phi.return) {
    temp = fhandle(model, X, ...) #[Y, Phi] = fhandle(model, X, varargin{:})
    Y=temp$Y; Phi=temp$Phi
  } else
    Y = fhandle(model, X, ...)
  
  if ("indexOut" %in% names(model) && length(model$indexOut)>0)
    Y[,setdiff(c(1:dim(Y)[2]),model$indexOut)] = NaN

  if (Phi.return)
    return(list(Y=Y, Phi=Phi))
  else
    return(Y)
}
