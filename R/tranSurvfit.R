uncondKendall <- function(x, y) {
    n <- length(x)
    res <- vector("double", 1)
    .C("uCondKendall", as.double(x), as.double(y), as.integer(n),
       out = as.double(res), PACKAGE = "tranSurv")$out
}

condKendall <- function(trun, obs, delta = NULL, method = "MB",
                        weights = NULL, a = 0, trans = "linear", ...) {
    methName <- c("MB", "IPW1", "IPW2")
    if (!(method %in% methName)) stop("Invalid method name", call. = FALSE)
    ## Weights arranged by c(trun, obs)
    out <- NULL
    out$Call <- match.call()
    n <- length(trun)
    if (class(trans) == "character") {
        if (trans == "linear") FUN <- function(X, T, a) (T + a * X) / (1 + a)
        if (trans == "log") FUN <- function(X, T, a) exp((log(replace(T, 0, 1)) + a * log(X)) / (1 + a))
        if (trans == "log2") FUN <- function(X, T, a) exp((1 + a) * log(replace(T, 0, 1)) - a * log(X))
        if (trans == "exp") FUN <- function(X, T, a) log((exp(T) + a * exp(X)) / (1 + a))
    } else {
        FUN <- match.fun(trans)
    }
    trun <- mapply(FUN, X = obs, T = trun, a = a)
    if (is.null(delta)) delta <- rep(1, length(trun))
    if (is.null(weights) & method == "MB") weights <- rep(1, 2 * n)
    if (is.null(weights) & method != "MB") {
        sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
        if (length(table(delta)) > 1 & 
            sum(head(sc$n.event[sc$n.event > 0]/sc$n.risk[sc$n.event > 0]) == 1) <= 2) {
            sc$time <- sc$time[sc$n.event > 0]
            sc$surv <- exp(-cumsum(sc$n.event[sc$n.event > 0]/sc$n.risk[sc$n.event > 0]))
        }
    }
    if (length(weights) == length(trun)) weights <- rep(weights, 2)
    if (is.null(weights) & method != "MB") {
        weights <- approx(sc$time, sc$surv, method = "constant",
                          xout = c(trun, obs), yleft = 1, yright = min(sc$surv))$y
    }
    res <- vector("double", 2)
    if (method != "IPW2") {
        tmp <- .C("condKendall", as.double(trun), as.double(obs), as.double(delta),
                  as.integer(n), as.double(weights), as.integer(which(method == methName)), 
                  tmp = as.double(res), PACKAGE = "tranSurv")$tmp
    } else {
        event <- delta == 1
        tmp <- .C("condKendall", as.double(trun[event]), as.double(obs[event]),
                  as.double(delta[event]), as.integer(sum(event)), as.double(weights[rep(event, 2)]),
                  as.integer(which(method == methName)), 
                  tmp = as.double(res), PACKAGE = "tranSurv")$tmp
    }
    out$PE <- tmp[1]
    out$SE <- ifelse(tmp[2] >= 0, sqrt(tmp[2]), NA)
    out$STAT <- abs(out$PE) / out$SE
    out$p.value <- 2 - 2 * pnorm(out$STAT)
    out$trans <- trans
    out$a <- a
    class(out) <- "condKendall"
    out
}

getA <- function(a, trun, obs, delta = NULL, sc = NULL, FUN, test = "CK") {
    if (is.null(delta)) delta <- rep(1, length(trun))
    FUN <- match.fun(FUN)
    ta <- mapply(FUN, X = obs, T = trun, a = a)
    if (is.null(sc)) {
        if (test == "CK") tmp <- condKendall(ta, obs, delta)
        if (test == "PC") tmp <- pmcc(ta[delta == 1], obs[delta == 1])
     } else {
        weights <- approx(sc$time, sc$surv, method = "constant", xout = c(ta, obs),
                          yleft = 1, yright = min(sc$surv))$y
        tmp <- condKendall(ta, obs, delta, method = "IPW2", weights = weights)
    }
    return(list(PE = tmp$PE, p.value = tmp$p.value))
}

tauEm <- function(f0, aij, bij, obs, ta) {
    uij <- t(t(aij) * f0) / colSums(t(aij) * f0)
    vij <- t(t(1 - bij) * f0) / colSums(t(bij) * f0)
    colSums(uij + vij) / sum(uij + vij)
}

tauLik <- function(f0, aij, bij, obs, ta) {
    yi <- unique(obs)
    yi <- yi[order(yi)]
    de <- rev(cumsum(rev(f0 * c(0, diff(yi)))))
    L <- sapply(1:length(f0), function(x) f0[x] / de[min(which(ta[order(obs)][x] <= yi))])
    -sum(log(ifelse(L <= 0, 1, L)))
}

tranSurv.control <- function(interval = c(-1, 50), lower = min(interval), upper = max(interval)) {
    list(lower = lower, upper = upper)
}

tranSurvfit <- function(trun, obs, delta = NULL, trans = "linear", plots = FALSE,
                        control = tranSurv.control(), ...) {
    ## trun = truncation time
    ## obs = observed failure time
    ## delta = censoring indicator
    if (is.null(delta)) delta <- rep(1, length(trun))
    if (class(trans) == "character") {
        if (trans == "linear") FUN <- function(X, T, a) (T + a * X) / (1 + a)
        if (trans == "log") FUN <- function(X, T, a) exp((log(replace(T, 0, 1)) + a * log(X)) / (1 + a))
        if (trans == "log2") FUN <- function(X, T, a) exp((1 + a) * log(replace(T, 0, 1)) - a * log(X))
        if (trans == "exp") FUN <- function(X, T, a) log((exp(T) + a * exp(X)) / (1 + a))
    } else {
        FUN <- match.fun(trans)
    }
    lower <- ifelse(control$lower == -Inf, -.Machine$integer.max, control$lower)
    upper <- ifelse(control$upper == Inf, .Machine$integer.max, control$upper)
    ini <- condKendall(trun, obs, delta)
    ini.ipw <- condKendall(trun, obs, delta, method = "IPW2")    
    sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
    if (plots) S0 <- survfit(Surv(trun, obs, delta) ~ 1) 
    if (length(table(delta)) > 1 & 
        sum(head(sc$n.event[sc$n.event > 0]/sc$n.risk[sc$n.event > 0]) == 1) <= 2) {
        sc$time <- sc$time[sc$n.event > 0]
        sc$surv <- exp(-cumsum(sc$n.event[sc$n.event > 0]/sc$n.risk[sc$n.event > 0]))
    }
    y0 <- sort(obs)
    trun1 <- trun[order(obs)][delta[order(obs)] == 1]
    obs1 <- sort(obs[delta == 1])
    delta1 <- delta[delta == 1]
    yi <- unique(obs1)
    byTau <- byP <- NULL
    p0 <- -as.numeric(coef(lm(trun ~ obs))[2])
    byTau$par <- uniroot.all(f = function(x) sapply(x, function(y)
        getA(y, trun1, obs1, sc = sc, FUN = FUN)$PE),
        interval = c(lower + 1e-5, upper))
    if (length(byTau$par) > 0) {
        byTau$par <- unique(c(uniroot.all(f = function(x)
            sapply(x, function(y) getA(y, trun1, obs1, sc = sc, FUN = FUN)$PE),
            interval = c(lower + 1e-5, byTau$par[1])), byTau$par))
        byTau$par <- sort(byTau$par)
        byTau$obj <- sapply(byTau$par, function(x) getA(x, trun1, obs1, sc = sc, FUN = FUN)$PE)
    } else {
        grids <- seq(lower + 1e-5, upper, length.out = 30)
        tmp <- sapply(1:29, function(y) 
            optimize(f = function(x) abs(getA(x, trun1, obs1, delta1, sc = sc, FUN = FUN)$PE),
                     tol = .01, interval = c(grids[y], grids[y + 1])))
        byTau$par <- as.numeric(tmp[1, which.min(tmp[2,])])
        byTau$obj <- as.numeric(tmp[2, which.min(tmp[2,])])
    }
    suppressWarnings(tmpP <- optimize(f = function(x) getA(x, trun1, obs1, delta1, sc = sc, FUN = FUN)$p.value, interval = c(lower, 2 * byTau$par[1] + 1), maximum = TRUE))
    byP$par <- tmpP$maximum
    byP$obj <- tmpP$objective
    tmp <- getA(byTau$par[1], trun1, obs1, sc = sc, FUN = FUN)$p.value
    byP$par <- ifelse(!is.na(tmp) & tmp > byP$obj, byTau$par[1], byP$par)
    byP$obj <- ifelse(!is.na(tmp) & tmp > byP$obj, tmp, byP$obj)
    if (abs(byTau$obj[1]) > 0.1 || byP$obj < 0.6)
        warning("Best estimate does not achieve tau_c = 0. This model may not be appropriate.", call. = FALSE)
    a <- byTau$par[1]
    ta <- mapply(FUN, X = obs1, T = trun1, a = a)
    wgtT <- approx(sc$time, sc$surv, ta, "constant", yleft = 1, yright = min(sc$surv))$y
    wgtX <- approx(sc$time, sc$surv, obs1, "constant", yleft = 1, yright = min(sc$surv))$y
    ## Turnbull's algorithm
    scX <- approx(sc$time, sc$surv, method = "constant", f = 0,
                  xout = yi, yleft = 1, yright = min(sc$surv))$y
    aij <- sapply(1:length(yi), function(x) 1 * (yi[x] == obs1))
    bij <- sapply(1:length(yi), function(x) 1 * (yi[x] >= ta))
    f0 <- rep(1 / length(yi), length(yi))
    tb <- squarem(par = f0, fixptfn = tauEm, objfn = tauLik, aij = aij, bij = bij,
                  obs = obs1, ta = ta)
    f1 <- tb$par / scX
    f1 <- f1 / sum(f1)
    Sy <- approx(yi, 1 - cumsum(f1), method = "constant", xout = y0, yleft = 1,
                   yright = min(1 - cumsum(f1)))$y
    ## Weighted KM
    ## s0 <- survfit(Surv(ta, yi, rep(1, length(yi))) ~ 1, weights = 1/scX)
    ## if (sum(s0$n.risk == 1) > 0) {
    ## s0 <- survfit(Surv(ta, yi, rep(1, length(yi))) ~ 1, weights = 1/scX,
    ##               subset = -which.min(abs(yi - s0$time[s0$n.risk == 1])))
    ## }
    ## Sy <- approx(s0$time, s0$surv, method = "constant", xout = y0, yleft = 1, yright = 0)$y

    ## Make plots 
    if (plots) {
        op1 <- par(mfrow = c(2,1), oma = c(1,1,1,1) + 0.1, mar = c(3.7,3,1,1) + 0.2)
        plot(trun, obs, cex = .4, main = "", xlab = "", ylab = "", pch = 4)
        mtext("x: uncensored   o: censored", 3, line = .1, cex = .9)
        mtext(expression(bold("Before transformation")), 3, line = .8, cex = 1.2)
        points(trun[delta == 0], obs[delta == 0], cex = .4, pch = 1)
        title(xlab = "Truncation times", ylab = "Failure times", line = 2, cex.lab = 1)
        plot(ta, obs1, cex = .4, main = "", xlab = "", ylab = "", pch = 4)
        mtext("x: uncensored   o: censored", 3, line = .1, cex = .9)
        mtext(expression(bold("After transformation (uncensored objects only)")), 3, line = .8, cex = 1.2)
        points(trun[delta == 0], obs[delta == 0], cex = .4, pch = 1)
        title(xlab = "Truncation times", ylab = "Failure times", line = 2, cex.lab = 1)
        par(op1)

        op1 <- par(mfrow = c(2,1), oma = c(1,1,1,1) + 0.1, mar = c(3.7,3,1,1.5) + 0.2, ask = TRUE)
        if (a > lower + .2)
            run <- c(seq(lower, max(lower, a - .2), length = 30),
                     seq(max(lower, a - .2), a + .2, .01),
                     seq(a + .2, 2 * a + 1, length = 30)) + 1e-5
        if (a <= lower + .2)
            run <- sort(unique(c(seq(lower, 2 * a + 1.01, length = 80) + 1e-5, a)))
        Ytau <- sapply(run, function(x) getA(x, trun1, obs1, delta1, sc, FUN)$PE)
        if (min(abs(Ytau), na.rm = TRUE) < byTau$obj[1]) {
            a <- byTau$par[1] <- run[which.min(abs(Ytau))]
            byTau$obj[1] <- min(abs(Ytau))
        }
        plot(run, Ytau, "l", xlab = "", ylab = "", xlim = c(max(min(run),lower), min(max(run), upper)), 
             ylim = range(Ytau, na.rm = TRUE), 
             main = paste("Estimate a by restricted IPW Kendall's tau", " (a = ", round(a, 3), ")",
                          sep = ""))
        title(xlab = "a", ylab = expression(tau[c]), line = 2, cex.lab = 1)
        abline(v = a, lty = 3)
        abline(h = byTau$obj[1], lty = 3)
        text(x = max(run), y = byTau$obj[1], labels = format(round(byTau$obj[1], 2), nsmall = 2),
             pos = 4, cex = .8, xpd = TRUE, srt = 0, offset = 1.5)
        if (byP$par > lower + .2)
            run <- c(seq(lower, max(lower, byP$par - .2), length = 30),
                     seq(max(lower, byP$par - .2), byP$par + .2, .01),
                     seq(byP$par + .2, 2 * byP$par + 1, length = 30)) + 1e-5
        if (byP$par <= lower + .2)
            run <- sort(unique(c(seq(lower, 2 * byP$par + 1.01, length = 80) + 1e-5, byP$par)))
        Yp <- sapply(run, function(x) getA(x, trun1, obs1, delta1, sc, FUN)$p.value)
        if (max(Yp, na.rm = TRUE) > byP$obj) {
            byP$par <- run[which.max(Yp)]
            byP$obj <- max(Yp, na.rm = TRUE)
        }
        plot(run, Yp, "l", xlab = "", ylab = "", xlim = c(max(min(run),lower), min(max(run), upper)),
             ylim = range(Yp, na.rm = TRUE), 
             main = paste("Estimate a by p-value", " (a = ", round(byP$par, 3), ")", sep = ""))
        title(xlab = "a", ylab = "p-values", line = 2, cex.lab = 1)
        abline(v = byP$par, lty = 3)
        abline(h = byP$obj, lty = 3)
        if (byP$obj >= 0.01)
            text(x = max(run), y = byP$obj, labels = format(round(byP$obj, 2), nsmall = 2),
                 pos = 4, cex = .8, xpd = TRUE, srt = 0, offset = 1.5)
        par(op1)

        op2 <- par(mar = c(3.5, 3.5, 2.5, 2.5), ask = TRUE)
        plot(sort(obs), Sy, xlab = "", ylab = "", "s")
        mtext(expression(bold("Survival estimation")), 3, line = .5, cex = 1.2)
        title(xlab = "Time", ylab = "Survival probability", line = 2, cex.lab = 1)
        lines(S0$time, S0$surv, lty = 2, "s")
        legend("topright", c("Transformation", "Kaplan-Meier"), lty = 1:2, bty = "n")
        par(op2)
    }
    out <- list(Sy = Sy, byTau = byTau, byP = byP,
                qind = data.frame(trun = ta, obs = obs1, wgtT = wgtT, wgtX = wgtX))
    out$Call <- match.call()
    out$iniKendall <- ini$PE
    out$iniKendall.ipw <- ini.ipw$PE
    out$iniP <- ini$p.value
    out$iniP.ipw <- ini.ipw$p.value
    class(out) <- "tranSurvfit"
    out
}

gofPlot <- function(trun, obs, delta = NULL) {
    n <- length(obs)
    if (is.null(delta)) delta <- rep(1, length(trun))
    sc <- survfit(Surv(trun, obs, 1 - delta) ~ 1)
    if (length(table(delta)) > 1 & 
        sum(head(sc$n.event[sc$n.event > 0]/sc$n.risk[sc$n.event > 0]) == 1) <= 2) {
        sc$time <- sc$time[sc$n.event > 0]
        sc$surv <- exp(-cumsum(sc$n.event[sc$n.event > 0]/sc$n.risk[sc$n.event > 0]))
    }
    fit <- tranSurvfit(trun, obs, delta)
    ta <- fit$qind$trun
    xa <- fit$qind$obs
    Fw <- survfit(Surv(-obs, -trun, rep(1, length(obs))) ~ 1, data = fit$qind)
    Fwz <- approx(-Fw$time, Fw$surv, method = "constant", xout = unique(sort(xa)), yleft = 0, yright = 1)$y
    fx <- diff(c(0, 1 - fit$Sy))
    fxz <- approx(sort(obs), fx, method = "constant", yleft = 0, yright = 0, xout = unique(sort(xa)))$y
    Scz <- approx(sc$time, sc$surv, method = "constant", yleft = 0, yright = 1, xout=unique(sort(xa)))$y
    gof <- fxz * Fwz * Scz
    s0 <- survfit(Surv(obs, rep(1, n)) ~ 1)
    plot(unique(sort(xa)), cumsum(gof) / sum(gof), "s", xlab = "Time", main = "Goodness of fit")
    lines(s0$time, 1 - s0$surv, "s", col = 2)
    legend("topleft", c("Transformation GOF", "ECDF"), col = 1:2, lty = 1, bty = "n")
}

pmcc <- function(trun, obs, a = 0, trans = "linear", ...) {
    out <- NULL
    out$Call <- match.call()
    n <- length(trun)
    if (class(trans) == "character") {
        if (trans == "linear") FUN <- function(X, T, a) (T + a * X) / (1 + a)
        if (trans == "log") FUN <- function(X, T, a) exp((log(replace(T, 0, 1)) + a * log(X))/(1 + a))
        if (trans == "log2") FUN <- function(X, T, a) exp((1 + a) * log(replace(T, 0, 1)) - a * log(X))
        if (trans == "exp") FUN <- function(X, T, a) log((exp(T) + a * exp(X)) / (1 + a))
    } else {
        FUN <- match.fun(trans)
    }
    trun <- mapply(FUN, X = obs, T = trun, a = a)
    res <- vector("double", 2)
    pmc <- .C("pmcc", as.double(trun), as.double(obs), as.integer(n),
              tmp = as.double(res), PACKAGE = "tranSurv")$tmp
    out$PE <- pmc[1]
    out$SE <- pmc[2]
    out$STAT <- pmc[1] / sqrt(pmc[2])
    out$p.value <- 2 - 2 * pnorm(abs(pmc[1]) / sqrt(pmc[2]))
    out$trans <- trans
    out$a <- a
    class(out) <- "pmcc"
    return(out)
}
