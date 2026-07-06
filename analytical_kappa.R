lambda_D <- 0.5  
lambda_R <- 0.75
tau <- 2 

F_D <- function(x, lambda_D) { 1 - exp(-lambda_D * x) }
F_R <- function(x, lambda_R) { 1 - exp(-lambda_R * x) }

f_DR <- function(x, y, lambda_D, lambda_R) {
  if (length(x) != length(y)) stop("x and y must be of same length")
  
  Fd <- F_D(x, lambda_D)
  Fr <- F_R(y, lambda_R)
  
  result <- rep(0, length(x))
  valid <- (Fd > 1e-10) & (Fr > 1e-10)
  
  if (any(valid)) {
    A <- 1 / Fd[valid]
    B <- 1 / Fr[valid]
    A_prime <- -lambda_D * exp(-lambda_D * x[valid]) / (Fd[valid]^2)
    B_prime <- -lambda_R * exp(-lambda_R * y[valid]) / (Fr[valid]^2)
    denom <- (A + B - 1)^3
    
    result[valid] <- ifelse(denom <= 1e-10, 0, 2 * A_prime * B_prime / denom)
  }
  
  return(result)
}

f_DR_vec <- Vectorize(f_DR, vectorize.args = c("x", "y"))

P_D_greater_u_less_tau_R_less_u_numerator <- function(u, lambda_D, lambda_R, tau) {
  integrand_outer <- function(d) {
    sapply(d, function(d_scalar) {  # VECTORIZE HERE
      upper_r <- min(d_scalar, u)
      integrand_inner <- function(r) {
        f_DR_vec(d_scalar, r, lambda_D, lambda_R)
      }
      integrate(integrand_inner, lower = 0, upper = upper_r, rel.tol = 1e-6)$value
    })
  }
  integrate(integrand_outer, lower = u, upper = tau, rel.tol = 1e-6)$value
}

P_D_greater_u_less_tau_R_less_u_numerator(1, lambda_D, lambda_R, 2)

P_R_leq_D <- function(lambda_D, lambda_R, upper = 10) {
  integrand_outer <- function(d) {
    sapply(d, function(d_scalar) {
      integrand_inner <- function(r) {
        f_DR_vec(d_scalar, r, lambda_D, lambda_R) 
      }
      integrate(integrand_inner, lower = 0, upper = d_scalar, rel.tol = 1e-6)$value
    })
  }
  integrate(integrand_outer, lower = 0, upper = upper, rel.tol = 1e-6)$value
}

P_R_leq_D(lambda_D, lambda_R)

# P(u < D <= tau, R <= u | R <= D)
P_D_greater_u_less_tau_R_less_u <- function(u, lambda_D, lambda_R, tau) {
  sapply(u, function(u_scalar) {
    numerator <- P_D_greater_u_less_tau_R_less_u_numerator(u_scalar, lambda_D, lambda_R, tau)
    denominator <- P_R_leq_D(lambda_D, lambda_R)
    numerator / denominator
  })
}

P_D_greater_u_less_tau_R_less_u (1, lambda_D, lambda_R, 2)

# calculate P(C < D <= tau, R <= C, C <= tau)
integral_1 <- integrate(function(u) {
  P_D_greater_u_less_tau_R_less_u(u, lambda_D, lambda_R, tau)* 1/3
}, lower = 0, upper = tau)$value


P_D_greater_u_R_less_u_numerator  <- function(u, lambda_D, lambda_R) {
  integrand_outer <- function(d) {
    sapply(d, function(d_scalar){
      upper_r <- min(d_scalar, u)
      integrand_inner <- function(r) {
        f_DR_vec(d_scalar, r, lambda_D, lambda_R)
      }
      integrate(integrand_inner, lower = 0, upper = upper_r, rel.tol = 1e-6)$value
    })
  }
  integrate(integrand_outer, lower = u, upper = 20, rel.tol = 1e-6)$value
}

P_D_greater_u_R_less_u_numerator(1, lambda_D, lambda_R)

# P(u < D, R <= u | R <= D)
P_D_greater_u_R_less_u <- function(u, lambda_D, lambda_R) {
  sapply(u, function(u_scalar) {
    numerator <- P_D_greater_u_R_less_u_numerator(u_scalar, lambda_D, lambda_R)
    denominator <- P_R_leq_D(lambda_D, lambda_R)
    numerator / denominator
  })
}

# calculate P(D > C, R <= C, C <= tau)
integral_2 <- integrate(function(u) {
    P_D_greater_u_R_less_u(u, lambda_D, lambda_R) * 1/3
}, lower = 0, upper = tau)$value

# calculate P(D > 2, R <= 2| R <= D) P(C > tau)
term_1 <- P_D_greater_u_R_less_u(tau,lambda_D, lambda_R) /3

integral_1/(integral_2 + term_1)



P_D_greater_u_less_tau_R_greater_u_numerator <- function(u, lambda_D, lambda_R, tau) {
  integrand_outer <- function(d) {
    sapply(d, function(d_scalar) {  # VECTORIZE HERE
      integrand_inner <- function(r) {
        f_DR_vec(d_scalar, r, lambda_D, lambda_R)
      }
      integrate(integrand_inner, lower = u, upper = d_scalar, rel.tol = 1e-6)$value
    })
  }
  integrate(integrand_outer, lower = u, upper = tau, rel.tol = 1e-6)$value
}

# P(u < D <= tau, R > u | R <= D)
P_D_greater_u_less_tau_R_greater_u <- function(u, lambda_D, lambda_R, tau) {
  sapply(u, function(u_scalar) {
    numerator <- P_D_greater_u_less_tau_R_greater_u_numerator(u_scalar, lambda_D, lambda_R, tau)
    denominator <- P_R_leq_D(lambda_D, lambda_R)
    numerator / denominator
  })
}

integral_3 <- integrate(function(u) {
  P_D_greater_u_less_tau_R_greater_u(u, lambda_D, lambda_R, tau)* 1/3
}, lower = 0, upper = tau)$value

P_D_greater_u_R_greater_u_numerator  <- function(u, lambda_D, lambda_R) {
  integrand_outer <- function(d) {
    sapply(d, function(d_scalar){
      integrand_inner <- function(r) {
        f_DR_vec(d_scalar, r, lambda_D, lambda_R)
      }
      integrate(integrand_inner, lower = u, upper = d_scalar, rel.tol = 1e-6)$value
    })
  }
  integrate(integrand_outer, lower = u, upper = 20, rel.tol = 1e-6)$value
}

# P(u < D, R > u | R <= D)
P_D_greater_u_R_greater_u <- function(u, lambda_D, lambda_R) {
  sapply(u, function(u_scalar) {
    numerator <- P_D_greater_u_R_greater_u_numerator(u_scalar, lambda_D, lambda_R)
    denominator <- P_R_leq_D(lambda_D, lambda_R)
    numerator / denominator
  })
}

# calculate P(D > C, R > C, C <= tau)
integral_4 <- integrate(function(u) {
  P_D_greater_u_R_greater_u(u, lambda_D, lambda_R) * 1/3
}, lower = 0, upper = tau)$value

term_2 <- P_D_greater_u_R_greater_u(tau,lambda_D, lambda_R) /3

integral_3/(integral_4 + term_2)
