# change of variable methods for constrained parameters

#' @export
varchange <- function(prior){
  UseMethod("varchange")
}

# bounded under (a,b)
#' @export
varchange.priorpdf <- function(prior){
  # retrieve bounds
  bounds <- attr(prior, "bounds")

  if (is.null(bounds)){
    return(prior)
  }

  else{
    bound_type <- !is.infinite(bounds)

    # lower and upper bound
    if (bound_type[1] & bound_type[2]){
      a <- bounds[1]
      b <- bounds[2]

      fy <- function(y){
        invlogit_y <- invlogit(y)
        prior(a + (b-a)*invlogit_y) + log(b-a) + log(invlogit_y) + log(1- invlogit_y)
      }

      attributes(fy) <- attributes(prior)[-1]
      class(fy) <- c("transformed", class(prior))
      attr(fy, "btype") <- "lowup"

      return(fy)
    }


    # lower bound
    else if (bound_type[1]){
      a <- bounds[1]
      fy <- function(y){
        prior(exp(y) + a) + y #log(f(exp(y) + a) * exp(y))
      }

      attributes(fy) <- attributes(prior)[-1]
      class(fy) <- c("transformed", class(prior))
      attr(fy, "btype") <- "low"

      return(fy)
    }

    # upper bound
    else{
      b <- bounds[2]
      fy <- function(y){
        prior(b - exp(y)) + y #log(f(exp(y) + a) * exp(y))
      }

      attributes(fy) <- attributes(prior)[-1]
      class(fy) <- c("transformed", class(prior))
      attr(fy, "btype") <- "up"

      return(fy)
    }
  }
}


#' @export
varchange.mgpm_prior <- function(priors){
  priors_tr <- lapply(priors, varchange)
  attr(priors_tr, "class") <- "mgpm_prior"

  return(priors_tr)
}


#' @export
print.transformed <- function(transformed, unit = TRUE){

  if (unit){
    cat("  ~ ")
  }

  disttype <- class(transformed)[3]
  btype <- attr(transformed, "btype")
  bounds <- attr(transformed, "bounds")

  if (btype == "lowup"){
    a <- bounds[1]
    b <- bounds[2]
    cat("logit((")
    print.priorpdf(transformed, unit = FALSE)
    cat(" - a)/(", "b - a)),", "\n", sep = "")
    cat("    a = ", a, ", b = ", b, sep = "")
  }

  else if (btype == "low"){
    a <- bounds[1]
    cat("log(")
    print.priorpdf(transformed, unit = FALSE)
    cat(" - a),", "\n")
    cat("    a =", a, "\n")
  }

  else{
    b <- bounds[2]
    cat("log(b - ", sep = "")
    print.priorpdf(transformed, unit = FALSE)
    cat(" - a),", "\n")
    cat("    a =", a, "\n")
  }
}



#' @export
trfunc <- function(priors_tr){
  # function that returns a function
  # to transform a parameter to/from the unbounded space
  # using the information of supports from priors

  # priors: priors_tr


  # parameter names and order in the vector p
  parnames <- names(priors_tr)

  # number of parameters
  npars <- length(parnames)

  # get information about the bounds from priors
  # bound types
  btypes <- lapply(parnames, function(x){attr(priors_tr[[x]], "btype")})
  pbounds <- lapply(parnames, function(x){attr(priors_tr[[x]], "bounds")})


  # placeholder list for functions
  # f: transform to unbounded
  # g: transform from unbounded
  f <- stats::setNames(vector("list", npars), parnames)
  g <- stats::setNames(vector("list", npars), parnames)


  for (i in 1:npars){
    # assign transformation functions for different types of bounds
    b_type <- btypes[[i]]

    if (is.null(b_type)){
      f[[i]] <- function(x){
        return(x)
      }

      g[[i]] <- function(y){
        return(y)
      }
    }

    else if (b_type == "low"){
      a <- pbounds[[i]][1]
      force(a)

      f[[i]] <- f_low(a)
      g[[i]] <- g_low(a)
    }

    else if (b_type == "up"){
      b <- pbounds[[i]][2]
      force(b)

      f[[i]] <- f_up(b)
      g[[i]] <- g_up(b)
    }

    else if (b_type == "lowup"){
      a <- pbounds[[i]][1]
      b <- pbounds[[i]][2]

      f[[i]] <- f_lowup(a, b)
      g[[i]] <- g_lowup(a, b)
    }
  }


  # combine transformation functions for each parameter into one
  t <- list()

  t$f <- function(p){
    for (i in 1:npars){
      p[i] <- f[[i]](p[i])
    }
    return(p)
  }

  t$g <- function(p){
    for (i in 1:npars){
      p[i] <- g[[i]](p[i])
    }
    return(p)
  }

  return(t)
}




# transformation functions for the data
# X to Y = f(X) or Y to X = f(Y)
f_lowup <- function(a, b){
  force(a)
  force(b)
  function(x){
    logit((x - a) / (b - a))
  }
}


g_lowup <- function(a, b){
  force(a)
  force(b)
  function(y){
    a + (b-a)*invlogit(y)
  }
}


f_low <- function(a){
  force(a)
  function(x){
    log(x - a)
  }
}


g_low <- function(a){
  force(a)
  function(y){
    exp(y) + a
  }
}


f_up <- function(b){
  force(b)
  function(x){
    log(b - x)
  }
}


g_up <- function(b){
  force(b)
  function(y){
    b - exp(y)
  }
}







