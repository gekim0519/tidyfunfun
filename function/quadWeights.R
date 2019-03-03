quadWeights = function (argvals, method = "trapezoidal") 
{
  ret <- switch(method, trapezoidal = {
    D <- length(argvals)
    1/2 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 
                                                                 2)], argvals[D] - argvals[D - 1])
  }, midpoint = c(0, diff(argvals)), stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
  return(ret)
}