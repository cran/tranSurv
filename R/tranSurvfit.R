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
    if (trans == "linear") trun <- (trun + a * obs) / (1 + a)
    if (trans == "log" || trans == "log2") trun <- ifelse(trun == 0, 1, trun)
    if (trans == "log") trun <- exp((log(trun) + a * log(obs)) / (1 + a))
    if (trans == "log2") trun <- exp((1 + a) * log(trun) - a * log(obs))
    if (trans == "exp") trun <- log((exp(trun) + a * exp(obs)) / (1 + a))
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
    out$SE <- sqrt(tmp[2])
    out$STAT <- abs(tmp[1]) / sqrt(tmp[2])
    out$p.value <- 2 - 2 * pnorm(out$STAT)
    out$trans <- trans
    out$a <- a
    class(out) <- "condKendall"
    out
}

getA <- function(a, trun, obs, delta = NULL, sc = NULL, trans = "linear", criterion = "PE") {
    if (is.null(delta)) delta <- rep(1, length(trun))
    if (trans == "linear") ta <- (trun + a * obs) / (1 + a)
    if (trans == "log" || trans == "log2") trun <- ifelse(trun == 0, 1, trun)
    if (trans == "log") ta <- exp((log(trun) + a * log(obs)) / (1 + a))
    if (trans == "log2") ta <- exp((1 + a) * log(trun) - a * log(obs))
    if (trans == "exp") ta <- log((exp(trun) + a * exp(obs)) / (1 + a))
    if (is.null(sc)) {
        tmp <- condKendall(ta, obs, delta)
     } else {
        weights <- approx(sc$time, sc$surv, method = "constant", xout = c(ta, obs),
                          yleft = 1, yright = min(sc$surv))$y
        tmp <- condKendall(ta, obs, delta, method = "IPW2", weights = weights)
    }
    if (criterion == "PE") return(tmp$PE)
    if (criterion != "PE") return(tmp$p.value)
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

tranSurvfit <- function(trun, obs, delta = NULL, trans = "linear", plots = FALSE, ...) {
    ## trun = truncation time
    ## obs = observed failure time
    ## delta = censoring indicator
    if (is.null(delta)) delta <- rep(1, length(trun))
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
    ## solve for the transformation parameter, a
    byTau <- byP <- NULL
    iniTau <- upTau <- getA(-.9999, trun1, obs1, delta1, sc, trans)
    for (i in 1:100) {
        upTau <- getA(i / 2 - 1, trun1, obs1, delta1, sc, trans)
        if (sign(iniTau) != sign(upTau) & sign(iniTau) * sign(upTau) != 0) break
    }
    if (i < 100) {
        tmp <- uniroot(f = getA, trun = trun1, obs = obs1, delta = delta1,
                       sc = sc, trans = trans, interval = c(-.9999, i / 2 - 1))
        byTau$par <- tmp$root
        byTau$obj <- as.numeric(tmp$f.root)
    } else {
        ## warning("Optimization over the interval (-1, 50).", call. = FALSE)
        tmp <- optimize(f = function(x) abs(getA(x, trun1, obs1, delta1, sc, trans)),
                        interval = c(-.9999, i / 2 - 1))
        byTau$par <- tmp$minimum
        byTau$obj <- as.numeric(tmp$objective) * sign(upTau)
        if (abs(byTau$obj) > 0.05) 
            warning("Best estimate does not achieve tau_c = 0. This model may not be appropriate.", call. = FALSE)
    }
    tmp <- optimize(f = function(x) getA(x, trun1, obs1, delta1, sc, trans, "P"),
                    interval = c(-.999, i / 2), maximum = TRUE)
    byP$par <- tmp$maximum
    byP$obj <- as.numeric(tmp$objective)
    a <- byTau$par
    if (trans == "linear") ta <- (trun1 + a * obs1) / (1 + a)
    if (trans == "log" || trans == "log2") trun1 <- ifelse(trun1 == 0, 1, trun1)
    if (trans == "log") ta <- exp((log(trun1) + a * log(obs1)) / (1 + a))
    if (trans == "log2") ta <- exp((1 + a) * log(trun1) - a * log(obs1))
    if (trans == "exp") ta <- log((exp(trun1) + a * exp(obs1)) / (1 + a))
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
        plot(trun1, obs1, cex = .4, main = "", xlab = "", ylab = "", pch = 4)
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
        if (byTau$par > -.8)
            run <- c(seq(-1, max(-1, byTau$par - .2), length = 20),
                     seq(max(-1, byTau$par - .2), byTau$par + .2, .015),
                     seq(byTau$par + .2, 2 * byTau$par + 1, length = 20)) + 1e-5
        if (byTau$par <= -.8)
            run <- sort(unique(c(seq(-1, 2 * byTau$par + 1.01, length = 60) + 1e-5, byTau$par)))
        Ytau <- sapply(run, function(x) getA(x, trun1, obs1, delta1, sc, trans))
        plot(run, Ytau, "l", xlab = "", ylab = "",
             main = paste("Estimate a by restricted IPW Kendall's tau", " (a = ", round(byTau$par, 3), ")", sep = ""))
        title(xlab = "a", ylab = expression(tau[c]), line = 2, cex.lab = 1)
        abline(v = byTau$par, lty = 3)
        abline(h = byTau$obj, lty = 3)
        text(x = max(run), y = byTau$obj, labels = format(round(byTau$obj, 2), nsmall = 2),
             pos = 4, cex = .8, xpd = TRUE, srt = 0, offset = 1.5)
        if (byP$par > -.8)
            run <- c(seq(-1, max(-1, byP$par - .2), length = 20),
                     seq(max(-1, byP$par - .2), byP$par + .2, .015),
                     seq(byP$par + .2, 2 * byP$par + 1, length = 20)) + 1e-5
        if (byP$par <= -.8)
            run <- sort(unique(c(seq(-1, 2 * byP$par + 1.01, length = 60) + 1e-5, byP$par)))
        Yp <- sapply(run, function(x) getA(x, trun1, obs1, delta1, sc, trans, "p"))
        plot(run, Yp, "l", xlab = "", ylab = "",
             main = paste("Estimate a by p-value", " (a = ", round(byP$par, 3), ")", sep = ""))
        title(xlab = "a", ylab = "p-values", line = 2, cex.lab = 1)
        abline(v = byP$par, lty = 3)
        abline(h = byP$obj, lty = 3)
        if (byP$obj >= 0.01)
            text(x = max(run), y = byP$obj, labels = format(round(byP$obj, 2), nsmall = 2),
                 pos = 4, cex = .8, xpd = TRUE, srt = 0, offset = 1.5)
        par(op1)

        op2 <- par(mar = c(3.5, 3.5, 2.5, 2.5), ask = TRUE)
        plot(sort(obs), Sy, "l", xlab = "", ylab = "")
        mtext(expression(bold("Survival estimation")), 3, line = .5, cex = 1.2)
        title(xlab = "Time", ylab = "Survival probability", line = 2, cex.lab = 1)
        lines(S0$time, S0$surv, lty = 2)
        legend("topright", c("Transformation", "Kaplan-Meier"), lty = 1:2, bty = "n")
        par(op2)
    }
    out <- list(Sy = Sy, byTau = byTau, byP = byP, qind = data.frame(trun = ta, obs = obs1))
    out$Call <- match.call()
    out$iniKendall <- ini$PE
    out$iniKendall.ipw <- ini.ipw$PE
    out$iniP <- ini$p.value
    out$iniP.ipw <- ini.ipw$p.value
    class(out) <- "tranSurvfit"
    out
}

