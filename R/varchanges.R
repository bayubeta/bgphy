# change of variable methods for constrained parameters
varchange <- function(f, ...){
  UseMethod("varchange")
}

# bounded under (a,b)
varchange.prior <- function(f){
  # retrieve bounds
  bounds <- attr(f, "bounds")

  if (is.null(bounds)){
    return(f)
  }

  else{
    bound_type <- !is.infinite(bounds)

    # lower and upper bound
    if (bound_type[1] & bound_type[2]){
      a <- bounds[1]
      b <- bounds[2]

      fy <- function(y){
        invlogit_y <- invlogit(y)
        f(a + (b-a)*invlogit_y) + (b-a) + log(invlogit_y) + log(1- invlogit_y)
      }

      attributes(fy) <- attributes(f)[-1]
      class(fy) <- c(class(f), "transformed")
      attr(fy, "untransformed") <- attr(f, "srcref")

      return(fy)
    }


    # lower bound
    else if (bound_type[1]){
      a <- bounds[1]
      fy <- function(y){
        f(exp(y) + a) + y #log(f(exp(y) + a) * exp(y))
      }

      attributes(fy) <- attributes(f)[-1]
      class(fy) <- c(class(f), "transformed")
      attr(fy, "untransformed") <- attr(f, "srcref")

      return(fy)
    }

    # upper bound
    else{
      b <- bounds[2]
      fy <- function(y){
        f(b - exp(y)) + y #log(f(exp(y) + a) * exp(y))
      }

      attributes(fy) <- attributes(f)[-1]
      class(fy) <- c(class(f), "transformed")
      attr(fy, "untransformed") <- attr(f, "srcref")

      return(fy)
    }
  }
}











