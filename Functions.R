get_occupancy <- function(time, data, duration = FALSE){
  markers <- levels(data$StreetMarker)
  occupancy <- d0 <- rep(NA, length(markers))
  for (i in 1:length(markers)) {
    data.marker <- filter(data, StreetMarker == markers[i])
    ind <- which(time %within% data.marker$interval)
    if (length(ind) == 1){
      occupancy[i] <- data.marker$State[ind]
      d0[i] <- as.numeric(difftime(time, int_start(data.marker$interval[ind]), units = "mins"))
    } else if (length(ind) == 0){
      occupancy[i] <- -1
    } else {
      occupancy[i] <- -length(ind)
    }
  }
  out <- list()
  out$occupancy <- occupancy
  if (duration) out$duration <- d0
  out
}

get_occupancy2 <- function(time, data, ind, duration = FALSE){
  occupancy <- d0 <- rep(NA, length(ind))
  markers <- levels(data$StreetMarker)[ind]
  data2 <- filter(data, StreetMarker %in% markers,
                  difftime(time, ArrivalTime, units = "days") >= 0,
                  difftime(time, ArrivalTime, units = "days") < 1)
  for (i in 1:length(ind)) {
    data.marker <- filter(data2, StreetMarker == markers[i])
    ind.marker <- which(time %within% data.marker$interval)
    if (length(ind.marker) == 1){
      occupancy[i] <- data.marker$State[ind.marker]
      d0[i] <- as.numeric(difftime(time, int_start(data.marker$interval[ind.marker]), units = "mins"))
    } else if (length(ind.marker) == 0){
      occupancy[i] <- -1
    } else {
      occupancy[i] <- -length(ind.marker)
    }
  }
  out <- tibble(ind = ind, marker = markers, occupancy = occupancy)
  if (duration) out$duration <- d0
  out
}

get_par <- function(time, model, data, occupancy.pre, state, distance){
  markers <- levels(data$StreetMarker)
  ind.markers <- match(unique(data$StreetMarker), markers)
  X <- model.matrix(model)
  X <- X[, -grep("frailty", colnames(X))]
  frails <- -model$frail/model$scale
  alpha <- 1/model$scale
  
  weekday <- factor(weekdays(time, abbreviate = TRUE))
  h <- paste0(")", hour(time))
  ind <- match(markers, data$StreetMarker)
  ind <- ind[which(!is.na(ind))]
  sos <- paste0("Code", data$SideOfStreetCode[ind])
  par <- list()
  par$alpha <- alpha 
  
  beta <- -model$coefficients/model$scale
  par$gamma <- rep(NA, length(markers))
  k <- 0
  for (i in ind.markers) {
    k <- k + 1
    occ <- occupancy.pre
    ind <- which(distance[i, ] <= 50)
    occ[ind] <- NA
    if (state == 0){
      frac <- length(which(occ == 1))/length(markers)
    } else {
      frac <- length(which(occ == 0))/length(markers)
    }
    beta_frac <- as.numeric(tail(beta, 1))
    ind <- c(1, grep(weekday, colnames(X)), grep(h, colnames(X)), grep(sos[k], colnames(X)))
    par$gamma[i] <- exp(sum(beta[ind]) + frac*beta_frac + 
                          frails[k])
  }
  par
}

get_par2 <- function(time, model, data, pre, state, distance){
  
  par <- list()
  par$alpha <- 1/model$scale
  
  X <- model.matrix(model)
  X <- X[, -grep("frailty", colnames(X))]
  frails <- -model$frail/model$scale
  weekday <- factor(weekdays(time, abbreviate = TRUE))
  h <- paste0(")", hour(time))
  
  markers <- levels(data$StreetMarker)
  sos <- paste0("Code", data$SideOfStreetCode[match(markers, data$StreetMarker)])

  beta <- -model$coefficients/model$scale
  par$gamma <- rep(NA, length(markers))
  for (i in 1:length(markers)) {
    ind.marker <- pre$ind[which(pre$marker %in% markers[i])]
    ind.nearby <- which(distance[ind.marker, ] <= 50)
    if (state == 0){
      frac <- length(which(pre$occupancy[pre$ind %in% ind.nearby] == 1))/length(ind.nearby)
    } else {
      frac <- length(which(pre$occupancy[pre$ind %in% ind.nearby] == 0))/length(ind.nearby)
    }
    if (length(ind.nearby) == 0) frac <- 0
    beta_frac <- as.numeric(tail(beta, 1))
    ind <- c(1, grep(weekday, colnames(X)), grep(h, colnames(X)), grep(sos[i], colnames(X)))
    par$gamma[i] <- exp(sum(beta[ind]) + frac*beta_frac + 
                          frails[i])
  }
  par
}

f_star <- function(x, par, d0){
  return(par$b*par$k*(d0 + x)^(par$k - 1)*exp(-par$b*((d0 + x)^par$k - d0^par$k)))
}
f <- function(x, par, d0){
  return(par$b*par$k*(x)^(par$k - 1)*exp(-par$b*x^par$k))
}
F_star <- function(x, par, d0){
  return(1 - exp(-par$b*((d0 + x)^par$k - d0^par$k)))
}
Fj <- function(x, par, d0){
  return(1 - exp(-par$b*x^par$k))
}


P_0star_0 <- function(s, par, d0, ...){
  val <- LaPlace(s, f_star, par = list(b = par$b1, k = par$k1), d0 = d0)*LaPlace(s, f, par = list(b = par$b2, k = par$k2), d0 = d0)*
    (1/s - LaPlace(s, Fj, par = list(b = par$b1, k = par$k1), d0 = d0))/
    (1 - LaPlace(s, f, par = list(b = par$b1, k = par$k1), d0 = d0)*LaPlace(s, f, par = list(b = par$b2, k = par$k2), d0 = d0))
  return(val)
}
P_1star_0 <- function(s, par, d0, ...){
  val <- LaPlace(s, f_star, par = list(b = par$b2, k = par$k2), d0 = d0)*
    (1/s - LaPlace(s, Fj, par = list(b = par$b1, k = par$k1), d0 = d0))/
    (1 - LaPlace(s, f, par = list(b = par$b1, k = par$k1), d0 = d0)*LaPlace(s, f, par = list(b = par$b2, k = par$k2), d0 = d0))
  return(val)
}
P_0star_1 <- function(s, par, d0, ...){
  val <- LaPlace(s, f_star, par = list(b = par$b1, k = par$k1), d0 = d0)*
    (1/s - LaPlace(s, Fj, par = list(b = par$b2, k = par$k2), d0 = d0))/
    (1 - LaPlace(s, f, par = list(b = par$b1, k = par$k1), d0 = d0)*LaPlace(s, f, par = list(b = par$b2, k = par$k2), d0 = d0))
  return(val)
}
P_1star_1 <- function(s, par, d0, ...){
  val <- LaPlace(s, f_star, par = list(b = par$b2, k = par$k2), d0 = d0)*LaPlace(s, f, par = list(b = par$b1, k = par$k1), d0 = d0)*
    (1/s - LaPlace(s, Fj, par = list(b = par$b2, k = par$k2), d0 = d0))/
    (1 - LaPlace(s, f, par = list(b = par$b1, k = par$k1), d0 = d0)*LaPlace(s, f, par = list(b = par$b2, k = par$k2), d0 = d0))
  return(val)
}

invlap <- function (Fs, t1, t2, nnt, a = 6, ns = 20, nd = 19, ...) 
{
  stopifnot(is.numeric(t1), length(t1) == 1, is.numeric(t2), 
            length(t2) == 1, is.numeric(nnt), length(nnt) == 1)
  Fs <- match.fun(Fs)
  radt <- linspace(t1, t2, nnt)
  if (t1 == 0) {
    radt <- radt[2:nnt]
    nnt <- nnt - 1
  }
  alfa <- beta <- numeric(ns + 1 + nd)
  for (n in 1:(ns + 1 + nd)) {
    alfa[n] <- a + (n - 1) * pi * (0+1i)
    beta[n] <- -exp(a) * (-1)^n
  }
  n <- 1:nd
  bdif <- rev(cumsum(gamma(nd + 1)/gamma(nd + 2 - n)/gamma(n)))/2^nd
  beta[(ns + 2):(ns + 1 + nd)] <- beta[(ns + 2):(ns + 1 + nd)] * 
    bdif
  beta[1] = beta[1]/2
  ft <- numeric(nnt)
  for (kt in 1:nnt) {
    tt <- radt[kt]
    s <- alfa/tt
    bt <- beta/tt
    btF <- bt * Fs(s, ...)
    ft[kt] <- sum(Re(btF))
  }
  return(list(x = radt, y = ft))
}

LaPlace <- function(s, g, par, d0){
  f.re <- function(x, s, g, par, d0){ 
    c <- Re(s)
    w <- Im(s)
    val <- exp(-c*x)*g(x, par, d0)*cos(w*x) 
    return(val)
  }
  f.im <- function(x, s, g, par, d0){ 
    c <- Re(s)
    w <- Im(s)
    val <- exp(-c*x)*g(x, par, d0)*sin(w*x) 
    return(val)
  }
  re <- im <- rep(0, length(s))
  for (j in 1:length(s)) {
    re[j] <- tryCatch(
      integrate(f.re, 0, Inf, s = s[j], g = g, d0 = d0, par = par, rel.tol = .Machine$double.eps^0.5)$value,
      error = function(cond){
        return(Inf)
      }
    )
    im[j] <- tryCatch(
      -integrate(f.im, 0, Inf, s = s[j], g = g, d0 = d0, par = par, rel.tol = .Machine$double.eps^0.5)$value,
      error = function(cond){
        return(Inf)
      }
    )
  }
  return(re + im*1i)
}

get_prediction_exp <- function(occupancy.pre, par.0, par.1, dur){
  
  prediction <- rep(NA, length(occupancy.pre))
  
  for (i in 1:length(prediction)) {
    if (is.element(occupancy.pre[i], c(0, 1))){ 
      P <- matrix(0, 2, 2)
      lambda.0 <- par.0$gamma[i]
      lambda.1 <- par.1$gamma[i]
      P[1, 1] <- lambda.1 + lambda.0*exp(-dur*(lambda.0 + lambda.1))
      P[1, 2] <- lambda.0*(1 - exp(-dur*(lambda.0 + lambda.1)))
      P[2, 1] <- lambda.1*(1 - exp(-dur*(lambda.0 + lambda.1)))
      P[2, 2] <- lambda.0 + lambda.1*exp(-dur*(lambda.0 + lambda.1))
      P <- P/(lambda.0 + lambda.1) 
      prediction[i] <- P[occupancy.pre[i] + 1, 1]
    } 
  }
  prediction
}

get_prediction_weibull <- function(occupancy.pre, d0, par.0, par.1, dur, star = TRUE) {
  print(star)
  prediction <- rep(NA, length(occupancy.pre))
  for (i in 1:length(prediction)) {
    if (is.element(occupancy.pre[i], c(0, 1))){ 
      par <- list(b1 = round(par.0$gamma[i], 3),
                  b2 = round(par.1$gamma[i], 3),
                  k1 = round(par.0$alpha, 3),
                  k2 = round(par.1$alpha, 3))
      if (star) {
        P <- rep(0, 4)
        if (occupancy.pre[i] == 0){
          P[1] <- 1 - F_star(dur, par = list(b = par$b1, k = par$k1), d0 = d0[i])
          P[3] <- invlap(P_0star_0, t1 = dur, t2 = dur, nnt = 1, d0 = d0[i], par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y
          P[4] <- invlap(P_0star_1, t1 = dur, t2 = dur, nnt = 1, d0 = d0[i], par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y
        } else {
          P[2] <- 1 - F_star(dur, par = list(b = par$b2, k = par$k2), d0 = d0[i])
          P[3] <- invlap(P_1star_0, t1 = dur, t2 = dur, nnt = 1, d0 = d0[i], par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y
          P[4] <- invlap(P_1star_1, t1 = dur, t2 = dur, nnt = 1, d0 = d0[i], par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y
        }
        if (sum(P) > 1.01 | sum(P) < 0.99) {
          warning("Rowsums are significantly different from 1 occurred.")
          print(i)
        }
        P <- P/sum(P)
        prediction[i] <- sum(P[c(1, 3)])
      } else {
        P <- rep(0, 2)
        if (occupancy.pre[i] == 0){
          P[1] <- invlap(P_0star_0, t1 = dur, t2 = dur, nnt = 1, d0 = 0, 
                         par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y + 
            1 - F_star(dur, par = list(b = par$b1, k = par$k1), d0 = 0)
          P[2] <- invlap(P_0star_1, t1 = dur, t2 = dur, nnt = 1, d0 = 0, par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y
        } else {
          P[1] <- invlap(P_1star_0, t1 = dur, t2 = dur, nnt = 1, d0 = 0, par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y
          P[2] <- invlap(P_1star_1, t1 = dur, t2 = dur, nnt = 1, d0 = 0, 
                         par = par, f = f, f_star = f_star, F_star = F_star, Fj = Fj)$y + 
            1 - F_star(dur, par = list(b = par$b2, k = par$k2), d0 = 0)
        }
        if (sum(P) > 1.01 | sum(P) < 0.99) {
          warning("Rowsums are significantly different from 1 occurred.")
          print(i)
        }
        P <- P/sum(P)
        prediction[i] <- P[1]
      }
    }
  }
  prediction
}

smoothConfidence <- function(theta, V, X, q = 0.05, R = 1000){
  set.seed(1)
  mu_sim <- matrix(0, R, nrow(X))
  for (i in 1:R) {
    theta_sim <- rmvn(1, theta, V)
    mu_sim[i, ] <- as.vector(X%*%theta_sim)
  }
  lower <- upper <- rep(0, ncol(mu_sim))
  for (j in 1:ncol(mu_sim)) {
    lower[j] <- quantile(mu_sim[, j], probs = q/2)
    upper[j] <- quantile(mu_sim[, j], probs = 1 - q/2)
  }
  list(lower = lower, upper = upper)
}

my.roc <- function(observed, predicted) {
  predicted2 <- unique(predicted)
  ord <- order(predicted2)
  thres <- c(-Inf, (predicted2[ord][2:length(predicted2)] + predicted2[ord][-length(predicted2)])/2, Inf)
  TPR <- FPR <- PPV <- rep(NA, length(thres))
  for (i in 1:length(thres)) {
    pred.0 <- which(predicted > thres[i])
    pred.1 <- which(predicted <= thres[i])
    obs.0 <- which(observed == 0)
    obs.1 <- which(observed == 1)
    TPR[i] <- length(intersect(pred.0, obs.0))/length(obs.0)
    TNR <- length(intersect(pred.1, obs.1))/length(obs.1)
    FPR[i] <- 1 - TNR
    PPV[i] <- length(intersect(pred.0, obs.0))/length(pred.0)
  }
  out <- tibble(threshold = thres, TPR = TPR, FPR = FPR, PPV = PPV)
  out
}
