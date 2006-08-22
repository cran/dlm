
###### Function to create block diagonal matrices (from R-help, I suppose)
bdiag <- function(...)
{
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)
    n <- length(x)
    if(n==0) return(NULL)
    x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
                stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1,]
    cc <- d[2,]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1,-1] <- rcum[-n]
    ind[2,] <- rcum
    ind[3,-1] <- ccum[-n]
    ind[4,] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                                 (y[3]+1):y[4]], imat=imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
} 

###### From R^p to the AR parameters of a stationary AR(p) (univariate)
###### For the multivariate analog, see Ansley and Kohn, 1986,
###### J. Statist. Comput. Simul. 24. 
ARtransPars <- function(raw) {
    p <- length(raw)
    return(.Call("ARtranspar", p, as.double(raw), PACKAGE="dlm"))
}

###### Constructor for dlm objects
dlm <- function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)
    ## required components
    nm <- c("m0", "C0", "FF", "V", "GG", "W") 
    nmInd <- match(nm, names(x))
    if (any(is.na(nmInd)))
        stop(paste("Component(s)", paste(nm[is.na(nmInd)], collapse=", "),
                   "is (are) missing"))
    x[nmInd[-1]] <- lapply(x[nmInd[-1]], as.matrix)
    if (!is.numeric(x$FF))
        stop("Component FF must be numeric")
    m <- nrow(x$FF)
    p <- ncol(x$FF)
    if (!is.numeric(x$V))
        stop("Component V must be numeric")
    if (!( nrow(x$V) == m && ncol(x$V) == m))
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$GG))
        stop("Component GG must be numeric")
    if (!( nrow(x$GG) == p && ncol(x$GG) == p))
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$W))
        stop("Component W must be numeric")
    if (!( nrow(x$W) == p && ncol(x$W) == p))
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$C0))
        stop("Component C0 must be numeric")
    if (!( nrow(x$C0) == p && ncol(x$C0) == p))
        stop("Incompatible dimensions of matrices")
    if (!( is.numeric(x$m0) && NCOL(x$m0) == 1 && NROW(x$m0) == p ))
        stop(paste("Component m0 must be a numeric vector of length",
                   "\n equal to ncol of component FF, or a matrix with one column and",
                   "\n number of rows equal to ncol of component FF"))
    if (!(all.equal(x$C0, t(x$C0)) && all(eigen(x$C0)$values >= 0)))
        stop("C0 is not a valid variance matrix")
    if (any( c(is.na(x$m0), is.na(x$C0))))
        stop("Missing values are not allowed in components m0 and C0")
    ## extra components for time-varying dlm 
    nm1 <- c("JFF", "JV", "JGG", "JW") 
    nm1Ind <- match(nm1, names(x))
    if (all(is.na(nm1Ind))) {
        if (!(all.equal(x$V, t(x$V)) && all(eigen(x$V)$values >= 0)))
            stop("V is not a valid variance matrix")
        if (!(all.equal(x$W, t(x$W)) && all(eigen(x$W)$values >= 0)))
            stop("W is not a valid variance matrix")
        mod <- x[nmInd]
        class(mod) <- "dlm"
        return(mod)
    }
    x[nm1Ind[!is.na(nm1Ind)]] <- lapply(x[nm1Ind[!is.na(nm1Ind)]], as.matrix)
    if (!is.null(x$JFF)) {
        if (!(is.numeric(x$JFF) &&
              nrow(x$JFF) == m && ncol(x$JFF) == p))
            stop("Invalid component JFF")
        JFF <- round(x$JFF)
        if (all(JFF == 0)) JFF <- NULL
    } else
    JFF <- NULL
    if (!is.null(x$JV)) {
        if (!(is.numeric(x$JV) &&
              nrow(x$JV) == m && ncol(x$JV) == m))
            stop("Invalid component JV")
        JV <- round(x$JV)
        if (all(JV == 0)) JV <- NULL
    } else
    JV <- NULL
    if (!is.null(x$JGG)) {
        if (!(is.numeric(x$JGG) &&
              nrow(x$JGG) == p && ncol(x$JGG) == p))
            stop("Invalid component JGG")
        JGG <- round(x$JGG)
        if (all(JGG == 0)) JGG <- NULL
    } else
    JGG <- NULL
    if (!is.null(x$JW)) {
        if (!(is.numeric(x$JW) &&
              nrow(x$JW) == p && ncol(x$JW) == p))
            stop("Invalid component JW")
        JW <- round(x$JW)
        if (all(JW == 0)) JW <- NULL
    } else
    JW <- NULL
    mx <- max(c(JFF, JV, JGG, JW))
    if ( mx <= 0 ) {
        mod <- x[nmInd]
        class(mod) <- "dlm"
        return(mod)
    }
    if ( is.null(x$X) )
        stop("Component X must be provided for time-varying models")
    if (!(is.numeric(x$X) && is.matrix(x$X) && ncol(x$X) >= mx))
        stop("Invalid component X")
    mod <- c(x[nmInd], list(JFF=JFF, JV=JV, JGG=JGG, JW=JW, X=x$X))
    class(mod) <- "dlm"
    return(mod)
}

is.dlm <- function(obj) inherits(obj, "dlm")

as.dlm <- function(obj)
    if (is.dlm(obj)) obj else dlm(obj)

is.dlmFiltered <- function(obj) inherits(obj, "dlmFiltered")

###### Regression
dlmModReg <- function(X, addInt=TRUE, dV=1, dW=rep(0,NCOL(X)+addInt),
                      m0=rep(0,length(dW)), C0=1e7*diag(nrow=length(dW)))
{
    ## sanity checks
    p <- NCOL(X) + addInt
    if (!( length(dV)==1 && length(dW)==p &&
          length(m0)==p && nrow(C0)==p && ncol(C0)==p))
        stop("Inconsistent dimensions of arguments")
    X <- as.matrix(X)
    JFF <- matrix(1:ncol(X),nr=1)
    if (addInt)
        JFF <- cbind(0,JFF)
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(1,1,p),
                V = as.matrix(dV),
                GG = diag(nrow=p),
                W = diag(x=dW, nrow=p),
                JFF = JFF,
                JV = NULL,
                JGG = NULL,
                JW = NULL,
                X = X)
    class(mod) <- "dlm"
    return(mod)
}

###### Polynomial trend
dlmModPoly <- function(order=2, dV=1,
                       dW=c(rep(0,order-1),1),
                       m0=rep(0,order), C0=1e7*diag(nrow=order))
{
    C0 <- as.matrix(C0)
    ## sanity checks
    if (!( length(dV)==1 && length(dW)==order &&
          length(m0)==order && nrow(C0)==order && ncol(C0)==order))
        stop("Inconsistent dimensions of arguments")
    GG <- diag(order)
    GG[row(GG) == col(GG)-1] <- 1
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(c(1,rep(0,order-1)),nr=1),
                V = as.matrix(dV),
                GG = GG,
                W = diag(dW,nrow=length(dW)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### Seasonal factors
dlmModSeas <- function(frequency, m0=rep(0,frequency-1),
                       C0=1e7*diag(nrow=frequency-1), dV=1,
                       dW=c(1,rep(0,frequency-2))) {
    frequency <- as.integer(frequency)
    p <- frequency - 1
    if (!( length(dV) == 1 && length(dW) == p && length(m0) == p &&
          nrow(C0) == p && ncol(C0) == p ))
        stop("Inconsistent dimensions of arguments")
    GG <- matrix(0,p,p)
    GG[row(GG) == col(GG)+1] <- 1
    GG[1,] <- -1
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(c(1,rep(0,p-1)),nr=1),
                V = as.matrix(dV),
                GG = GG,
                W = diag(dW,nrow=length(dW)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### Fourier representation
dlmModTrig <- function(s, q, om, tau, m0, C0, dV=1, dW=0) {
    if ( hasArg(s) ) {
        if ( hasArg(om) )
            stop("Cannot specify both 's' and 'om'")
        if ( hasArg(tau) )
            stop("Cannot specify both 's' and 'tau'")
        s <- as.integer(s)
        sHalf <- s %/% 2
        if ( hasArg(q) ) {
            q <- as.integer(q)
            if ( q > sHalf )
                stop(paste("Can use",sHalf,"components at most (",q,"were asked for )"))
            if ( q <= 0 ) stop("'q' must be positive")
        } else {
            q <- sHalf }
        om <- 2 * base::pi / s
        evenAll <- (q == sHalf) && !(s %% 2)
        if ( sHalf == 1 && evenAll ) {
            FF <- matrix(1,1,1)
            GG <- matrix(-1,1,1)
        } else {
            H1 <- diag(x=cos(om), nrow=2)
            H1[1,2] <- sin(om)
            H1[2,1] <- - H1[1,2]
            if ( q > 1 ) {
                h <- vector("list",q)
                h[[1]] <- H1
                i <- 2
                while ( i < q ) {
                    h[[i]] <- h[[i-1]] %*% H1
                    i <- i + 1
                }
                h[[q]] <- if ( evenAll ) matrix(-1,1,1) else h[[q-1]] %*% H1
                GG <- bdiag(h)
            } else {
                GG <- H1 }
        }
    } else {
        if ( !hasArg(om) )
            if ( !hasArg(tau) )
                stop("One of 's', 'om' or 'tau' must be specified")
            else om <- 2 * base::pi / tau
        if ( !hasArg(q) )
            stop("When 'om' or 'tau' is specified, 'q' must be supplied as well")
        q <- as.integer(q)
        if ( q <= 0 ) stop("'q' must be positive")
        evenAll <- FALSE
        H1 <- diag(x=cos(om), nrow=2)
        H1[1,2] <- sin(om)
        H1[2,1] <- - H1[1,2]
        if ( q > 1 ) {
            h <- vector("list",q)
            h[[1]] <- H1
            for ( i in 2:q ) 
                h[[i]] <- h[[i-1]] %*% H1
            GG <- bdiag(h)
        } else {
            GG <- H1 }
    }
    if ( hasArg(m0) ) {
        if ( length(m0) != nrow(GG) ) stop("Wrong length of 'm0'")
    } else {
        m0 <- rep(0,nrow(GG)) }
    if ( hasArg(C0) ) {
        if ( (nrow(C0) != nrow(GG)) || (ncol(C0) != nrow(GG)) )
            stop("Wrong dimension of 'C0'")
    } else {
        C0 <- diag(x=1e7, nrow=nrow(GG)) }
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(if (evenAll) c(rep(c(1,0),q-1),1)
                else c(1,0),1,nrow(GG)),
                V = matrix(dV,1,1),
                GG = GG,
                W = diag(x=dW,nrow=nrow(GG)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### ARMA
dlmModARMA <- function(ar=NULL, ma=NULL, sigma2=1, m0, C0, dV)
{
    if (is.matrix(sigma2) && (m <- nrow(sigma2)) > 1)
    { ## multivariate
        if (ncol(sigma2) != m)
            stop("sigma2 must be a square matrix")
        r <- max(p <- length(ar), (q <- length(ma)) + 1)
        k <- m * r
        if (hasArg("dV")) {
            if ( length(dV) != m )
                stop("Incompatible dimensions of arguments")
        }
        else
            dV <- rep(0,m)
        if (hasArg("m0")) {
            if ( length(m0) != k )
                stop("Incompatible dimensions of arguments")
        }
        else 
            m0 <- rep(0,k)
        if (hasArg("C0")) {
            if ( !( nrow(C0) == k && ncol(C0) == k ))
                stop("Incompatible dimensions of arguments")
        }
        else 
            C0 <- 1e7*diag(nrow=k)
        FF <- matrix(0,m,k)
        FF[row(FF) == col(FF)] <- 1
        GG <- matrix(0,k,k)
        for (i in seq(length.out=p))  GG[ (1 + (i-1) * m):(i * m), 1:m ] <- ar[[i]]
        GG[row(GG) == col(GG) - m] <- 1
        R <- matrix(0,k,m)
        R[row(R) == col(R)] <- 1
        for (i in seq(length.out=q)) R[ (1 + i * m):((i+1) * m), ] <- ma[[i]]
        mod <- list(
                    m0 = m0,
                    C0 = C0,
                    FF = FF,
                    V = diag(dV,nrow=m),
                    GG = GG,
                    W = R %*% sigma2 %*% t(R),
                    JFF = NULL,
                    JV = NULL,
                    JGG = NULL,
                    JW = NULL)
    }
    else 
    { ## univariate
        r <- max(p <- length(ar), (q <- length(ma)) + 1)
        if (hasArg("dV")) {
            if ( length(dV) != 1 )
                stop("Incompatible dimensions of arguments")
        }
        else
            dV <- 0
        if (hasArg("m0")) {
            if ( length(m0) != r )
                stop("Incompatible dimensions of arguments")
        }
        else 
            m0 <- rep(0,r)
        if (hasArg("C0")) {
            if ( !( nrow(C0) == r && ncol(C0) == r ))
                stop("Incompatible dimensions of arguments")
        }
        else 
            C0 <- 1e7*diag(nrow=r)
        GG <- matrix(0,r,r)
        if ( p > 0 ) GG[ 1:p, 1 ] <- ar
        GG[row(GG) == col(GG) - 1] <- 1
        R <- rep(0,r)
        R[1] <- 1
        if ( q > 0 ) R[ 1 + 1:q ] <- ma
        mod <- list(
                    m0 = m0,
                    C0 = C0,
                    FF = matrix(c(1,rep(0,r-1)),nr=1),
                    V = matrix(dV),
                    GG = GG,
                    W = sigma2 * crossprod(matrix(R,nr=1)),
                    JFF = NULL,
                    JV = NULL,
                    JGG = NULL,
                    JW = NULL)
    }
    class(mod) <- "dlm"
    return(mod)
}


###### Addition for "dlm" objects
"+.dlm" <- function(mod1, mod2)
{
    if ( (m <- nrow(mod1$FF)) != nrow(mod2$FF) )
        stop("Incompatible models")
    if ( (x1 <- !is.null(mod1$X)) | (x2 <- !is.null(mod2$X)) ) {
        if ( x1 && x2 && ( nrow(mod1$X) != nrow(mod2$X) ) )
            stop("Number of rows of mod1$X and mod2$X must match")
        plus <- function(u,v) ifelse( u==0, 0, u+v )
        p1 <- ncol(mod1$FF)
        p2 <- ncol(mod2$FF)
        r <- if ( x1 ) ncol(mod1$X) else 0
        if ( (one <- !is.null(mod1$JFF)) | (two <- !is.null(mod2$JFF)) ) {
            JFF1 <- if ( one ) mod1$JFF else matrix(0,m,p1)
            JFF2 <- if ( two ) plus(mod2$JFF,r) else matrix(0,m,p2)
            JFF <- cbind(JFF1, JFF2)
        } else
        JFF <- NULL
        if ( (one <- !is.null(mod1$JGG)) | (two <- !is.null(mod2$JGG)) ) {
            JGG1 <- if ( one ) mod1$JGG else matrix(0,p1,p1)
            JGG2 <- if ( two ) plus(mod2$JGG,r) else matrix(0,p2,p2)
            JGG <- bdiag(list(JGG1, JGG2))
        } else
        JGG <- NULL
        if ( (one <- !is.null(mod1$JW)) | (two <- !is.null(mod2$JW)) ) {
            JW1 <- if ( one ) mod1$JW else matrix(0,p1,p1)
            JW2 <- if ( two ) plus(mod2$JW,r) else matrix(0,p2,p2)
            JW <- bdiag(list(JW1, JW2))
        } else
        JW <- NULL
        if ( (one <- !is.null(mod1$JV)) | (two <- !is.null(mod2$JV)) ) {
            if ( one && two )
                stop("mod1$V and mod2$V cannot be both time-varying")
            if ( one ) {
                if ( any(mod2$V != 0) ) {
                    mod2$V[] <- 0
                    warning("the value of mod2$V has been discarded")
                }
                JV <- mod1$JV
            }
            if ( two ) {
                if ( any(mod1$V != 0) ) {
                    mod1$V[] <- 0
                    warning("the value of mod1$V has been discarded")
                }
                JV <- plus(mod2$JV,r)
            }
        } else
        JV <- NULL
        mod <- list(
                    m0 = c(mod1$m0, mod2$m0),
                    C0 = bdiag(list(mod1$C0, mod2$C0)),
                    FF = cbind(mod1$FF, mod2$FF),
                    V = mod1$V + mod2$V,
                    GG = bdiag(list(mod1$GG, mod2$GG)),
                    W = bdiag(list(mod1$W, mod2$W)),
                    JFF = JFF,
                    JV = JV,
                    JGG = JGG,
                    JW = JW,
                    X = switch( x1 + 2 * x2,
                    mod1$X,
                    mod2$X,
                    cbind( mod1$X, mod2$X) ))
    } else
    mod <- list(
                m0 = c(mod1$m0, mod2$m0),
                C0 = bdiag(list(mod1$C0, mod2$C0)),
                FF = cbind(mod1$FF, mod2$FF),
                V = mod1$V + mod2$V,
                GG = bdiag(list(mod1$GG, mod2$GG)),
                W = bdiag(list(mod1$W, mod2$W)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### Print method for "dlm" objects
print.dlm <- function(x, ...) {
    nmRef <- c("FF","V","GG","W","JFF","JV","JGG","JW","X","m0","C0")
    nm <- names(x)
    what <- !sapply(x, is.null)
    ind <- match(nmRef,nm)
    ind <- ind[!is.na(ind)]
    good <- nm[ind][what[ind]]
    if ( !is.na(match("X",good)) && nrow(x$X) > 2 ) {
        x$X <- rbind(formatC(x$X[1:2,,drop=FALSE],digits=4,format="fg"),
                       c("...",rep("",ncol(x$X)-1)))
    }
    print(x[match(good,nm)], quote=FALSE)
    invisible(x)
}

###### Negative loglikelihood
dlmLL <- function(y, mod, debug=FALSE)
{
    ## calculations based on singular value decomposition
    ## Note: V must be nonsingular
    ## The C code relies on the order of the elements in 'mod'
    ## mod = list(m0, C0, FF, V, GG, W)
    storage.mode(y) <- "double"
    if (!debug) {
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz]) - 1
            storage.mode(mod$JFF) <- "integer"
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz]) - 1
            storage.mode(mod$JV) <- "integer"
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz]) - 1
            storage.mode(mod$JGG) <- "integer"
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz]) - 1
            storage.mode(mod$JW) <- "integer"
        }
        if (any(c(tvFF,tvV,tvGG,tvW))) {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W",
                               "JFF", "JV", "JGG", "JW", "X"), names(mod))]
            return(.Call("dlmLL", y, mod, tvFF, tvV, tvGG, tvW, PACKAGE="dlm"))
        } else {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W"), names(mod))]
            return(.Call("dlmLL0", y, mod, PACKAGE="dlm"))
        }
    }
    else {
        y <- as.matrix(y)
        n <- nrow(y)
        ll <- 0 # negative loglikelihood
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz])
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz])
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
        }
        tvFV <- tvFF || tvV
        ## preliminary calculations, if possible (non time-varying case)
        if ( !tvV ) {
            tmp <- La.svd(mod$V,nu=0)
            Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
            Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
            sqrtVinv <- Dv.inv * t(Uv)
            sqrtV <- sqrt(tmp$d) * Uv # t()%*%() = V 
            if ( !tvFF ) 
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
        }
        if ( !tvW ) {
            svdW <- La.svd(mod$W,nu=0)
            sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W 
        }
        tmp <- La.svd(mod$C0,nu=0)
        Ux <- t(tmp$vt); Dx <- sqrt(tmp$d)
        for (i in seq(length=n)) {
            ## set time-varying matrices
            if ( tvFF ) 
                mod$FF[mod$JFF[,-3,drop=FALSE]] <- mod$X[i,mod$JFF[,3]]
            if ( tvV ) {
                mod$V[mod$JV[,-3,drop=FALSE]] <- mod$X[i,mod$JV[,3]]
                tmp <- La.svd(mod$V,nu=0)
                Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
                Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
                sqrtVinv <- Dv.inv * t(Uv)
                sqrtV <- sqrt(tmp$d) * Uv # t()%*%() = V 
            }
            if ( tvGG ) 
                mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
            if ( tvW ) {
                mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
                svdW <- La.svd(mod$W,nu=0)
                sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W 
            }
            if ( tvFV ) 
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
            ## prior
            a <- mod$GG %*% mod$m0
            tmp <- La.svd(rbind( Dx*t(mod$GG%*%Ux), sqrtW ), nu=0)
            Ux.prior <- t(tmp$vt)
            Dx.prior <- tmp$d
            ## one-step forecast
            f <- mod$FF %*% a
            tmp <- La.svd(rbind( Dx.prior*t(mod$FF%*%Ux.prior), sqrtV ), nu=0)
            Uy <- t(tmp$vt)
            Dy <- tmp$d
            ## posterior
            D.inv <- 1/Dx.prior
            D.inv[abs(D.inv)==Inf] <- 0
            tmp <- La.svd(rbind(sqrtVinv%*%mod$FF%*%Ux.prior,
                                diag(x=D.inv,nrow=length(D.inv))), nu=0) 
            Ux <- Ux.prior %*% t(tmp$vt)
            Dx <- 1/tmp$d 
            Dx[abs(Dx)==Inf] <- 0
            e <- as.matrix(y[i,]-f)
            mod$m0 <- a + crossprod(Dx*t(Ux)) %*%
                tF.Vinv %*% e
            ## update scaled negative loglikelihood
            ll <- ll + 2*sum(log(Dy)) + crossprod(crossprod(Uy,e)/Dy)
        }
    }
    return(drop(ll)*0.5) # negative loglikelihood
}


"dlmMLE" <- function(y, parm, build, method = "L-BFGS-B",
                     ..., debug = FALSE)
{
    storage.mode(y) <- "double"
    logLik <- function(parm, ...)
    {
        mod <- build(parm, ...)
        return(dlmLL(y=y, mod=mod, debug=debug))
    }
    out <- optim(parm, logLik, method=method, ...)
    return(out)
}


dlmFilter <- function(y, mod, debug = FALSE, simplify = FALSE)
{
    ## Note: V must be nonsingular
    storage.mode(y) <- "double"
    mod1 <- mod
    yAttr <- attributes(y)
    ytsp <- tsp(y)
    y <- as.matrix(y)
    timeNames <- dimnames(y)[[1]]
    stateNames <- names(mod$m0)
    if (!debug) {
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz]) - 1
            storage.mode(mod$JFF) <- "integer"
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz]) - 1
            storage.mode(mod$JV) <- "integer"
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz]) - 1
            storage.mode(mod$JGG) <- "integer"
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz]) - 1
            storage.mode(mod$JW) <- "integer"
        }
        if (any(c(tvFF,tvV,tvGG,tvW))) {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W",
                               "JFF", "JV", "JGG", "JW", "X"), names(mod))]
            ans <- .Call("dlmFilter", y, mod, tvFF, tvV, tvGG, tvW, PACKAGE="dlm")
        } else {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W"), names(mod))]
            ans <- .Call("dlmFilter0", y, mod, PACKAGE="dlm")
        }
        names(ans) <- c("m", "U.C", "D.C", "a", "U.R", "D.R", "f")
    }
    else {
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz])
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz])
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
        }
        tvFV <- tvFF || tvV
        m <- rbind(mod$m0,matrix(0,nr=nrow(y),nc=length(mod$m0))) # filtered values
        a <- matrix(0,nr=nrow(y),nc=length(mod$m0))
        f <- matrix(0,nr=nrow(y),nc=ncol(y))
        U.C <- vector(1+nrow(y),mode="list")
        D.C <- matrix(0,1+nrow(y),length(mod$m0))
        U.R <- vector(nrow(y),mode="list")
        D.R <- matrix(0,nrow(y),length(mod$m0))
        ## preliminary calculations, if possible (non time-varying case)
        if ( !tvV ) {
            tmp <- La.svd(mod$V,nu=0)
            Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
            Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
            sqrtVinv <- Dv.inv * t(Uv)
            sqrtV <- Dv * Uv # t()%*%() = V 
            if ( !tvFF ) 
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
        }
        if ( !tvW ) {
            svdW <- La.svd(mod$W,nu=0)
            sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W 
        }
        tmp <- La.svd(mod$C0,nu=0)
        U.C[[1]] <- t(tmp$vt)
        D.C[1,] <- sqrt(tmp$d)
        for (i in seq(length=nrow(y))) {
            ## set time-varying matrices
            if ( tvFF ) 
                mod$FF[mod$JFF[,-3,drop=FALSE]] <- mod$X[i,mod$JFF[,3]]
            if ( tvV ) {
                mod$V[mod$JV[,-3,drop=FALSE]] <- mod$X[i,mod$JV[,3]]
                tmp <- La.svd(mod$V,nu=0)
                Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
                Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
                sqrtVinv <- Dv.inv * t(Uv)
                sqrtV <- sqrt(tmp$d) * Uv # t()%*%() = V 
            }
            if ( tvGG ) 
                mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
            if ( tvW ) {
                mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
                svdW <- La.svd(mod$W,nu=0)
                sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W 
            }
            if ( tvFV ) 
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
            ## prior
            a[i,] <- mod$GG %*% m[i,]
            tmp <- La.svd(rbind( D.C[i,]*t(mod$GG%*%U.C[[i]]), sqrtW ), nu=0)
            U.R[[i]] <- t(tmp$vt)
            D.R[i,] <- tmp$d
            ## one-step forecast
            f[i,] <- mod$FF %*% a[i,]
            ## posterior
            D.Rinv <- 1/D.R[i,]
            D.Rinv[abs(D.Rinv)==Inf] <- 0
            tmp <- La.svd(rbind(sqrtVinv%*%mod$FF%*%U.R[[i]],
                                diag(x=D.Rinv,nrow=length(D.Rinv))), nu=0)
            U.C[[i+1]] <- U.R[[i]] %*% t(tmp$vt)
            foo <- 1/tmp$d; foo[abs(foo)==Inf] <- 0
            D.C[i+1,] <- foo
            m[i+1,] <- a[i,] + crossprod(D.C[i+1,]*t(U.C[[i+1]])) %*%
                tF.Vinv %*% as.matrix(y[i,]-f[i,])
        }
        m <- drop(m); a <- drop(a); f <- drop(f)
        attributes(f) <- yAttr
        ans <- list(m=m,U.C=U.C,D.C=D.C,a=a,U.R=U.R,D.R=D.R,f=f)
    }
    
    ans$m <- drop(ans$m); ans$a <- drop(ans$a); ans$f <- drop(ans$f)
    attributes(ans$f) <- yAttr
    if (!is.null(ytsp)) {
        tsp(ans$a) <- ytsp
        tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
        class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 1) c("mts","ts") else "ts" 
    }
    if (!(is.null(timeNames) && is.null(stateNames))) {
        dimnames(ans$a) <- list(timeNames, stateNames)
        dimnames(ans$m) <- list(if(is.null(timeNames)) NULL else c("",timeNames),
                                stateNames)
    }
    if (simplify)
        return(c(mod=list(mod1), ans))
    else
    {
        attributes(y) <- yAttr
        ans <- c(y=list(y), mod=list(mod1), ans)
        class(ans) <- "dlmFiltered"
        return(ans)
    }
}


dlmSmooth <- function(modFilt, debug = FALSE)
{
    big <- 1 / sqrt(.Machine$double.eps)
    mod <- c(modFilt[match(c("m", "U.C", "D.C", "a", "U.R", "D.R"),names(modFilt))],
             modFilt$mod[match(c("GG", "W", "JGG", "JW", "X"), names(modFilt$mod))]) 
    mAttr <- attributes(mod$m)
    mod$m <- as.matrix(mod$m)
    mod$a <- as.matrix(mod$a)
    if (!debug) {
        ## define flags for time-varying components
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz]) - 1
            storage.mode(mod$JGG) <- "integer"
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz]) - 1
            storage.mode(mod$JW) <- "integer"
        }
        if (tvGG || tvW) {
            ans <- .Call("dlmSmooth", mod, tvGG, tvW, big, PACKAGE="dlm")
            names(ans) <- c("s", "U.S", "D.S")
        }
        else {
            ans <- .Call("dlmSmooth0", mod, big, PACKAGE="dlm")
            names(ans) <- c("s", "U.S", "D.S")
        }
        } else {
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
        }
        n <- length(mod$U.R) # number of obs
        p <- NCOL(mod$m) # dimension of state vector
        s <- rbind(matrix(0,n,p), mod$m[n+1,])
        U.S <- vector("list", length=n+1)
        U.S[[n+1]] <- mod$U.C[[n+1]]
        D.S <- rbind(matrix(0,n,p), mod$D.C[n+1,])
        ## preliminary calculations, if possible (time-invariant case)
        if ( !tvW ) {
            tmp <- La.svd(mod$W,nu=0)
            Dw <- sqrt(tmp$d)
            Dw.inv <- pmin(1/Dw, big)
            sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1) 
        }
        if (n > 0)
            for (i in n:1)
            {
                ## set relevant time-varying matrices
                if ( tvGG )
                    mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
                if ( tvW ) {
                    mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
                    tmp <- La.svd(mod$W,nu=0)
                    Dw <- sqrt(tmp$d)
                    Dw.inv <- pmin(1/Dw, big)
                    sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1) 
                }
                Dinv <- 1/mod$D.R[i,]; Dinv[abs(Dinv)==Inf] <- 0
                H <- crossprod(mod$D.C[i,]*t(mod$U.C[[i]])) %*%
                    t(mod$GG) %*% crossprod(Dinv*t(mod$U.R[[i]]))
                Dinv <- 1/mod$D.C[i,]
                Dinv[abs(Dinv)==Inf] <- 0
                tmp <- La.svd(rbind( sqrtWinv%*%mod$GG, Dinv*t(mod$U.C[[i]])), nu=0)
                Dinv <- 1/tmp$d
                Dinv[abs(Dinv)==Inf] <- 0
                tmp <- La.svd(rbind(Dinv*tmp$vt, D.S[i+1,]*t(H%*%U.S[[i+1]])))
                U.S[[i]] <- t(tmp$vt)
                D.S[i,] <- tmp$d
                s[i,] <- mod$m[i,] + H %*% (s[i+1,]-mod$a[i,])
            }
        ans <- list(s=s, U.S=U.S, D.S=D.S)
    }
    attributes(ans$s) <- mAttr
                                        #     if (!is.null(tsp(mod$m))) {
#         tsp(ans$s) <- tsp(mod$m)
#         class(ans$s) <- if (NCOL(ans$s) > 1) c("mts","ts") else "ts"
#     }
#     dimnames(ans$s) <- dimnames(mod$m)
    return(ans)
}


dlmBSample <- function(modFilt)
{
    eps <- .Machine$double.eps^.4
    mod <- c(modFilt[match(c("m", "U.C", "D.C", "a", "U.R", "D.R"),names(modFilt))],
             modFilt$mod[match(c("GG", "W", "JGG", "JW", "X"), names(modFilt$mod))]) 
    n <- length(mod$U.R) # number of obs
    p <- NCOL(mod$m) # dimension of state vector
    mtsp <- tsp(mod$m)
    if (p==1) {
        mod$m <- as.matrix(mod$m)
        mod$a <- as.matrix(mod$a)
    }
    if (is.null(mod$JGG))
        tvGG <- FALSE
    else {
        tvGG <- TRUE
        nz <- mod$JGG != 0
        mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
    }
    if (is.null(mod$JW))
        tvW <- FALSE
    else {
        tvW <- TRUE
        nz <- mod$JW != 0
        mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
    }
    tvGW <- tvGG || tvW
    theta <- matrix(0,n+1,p)
    ## preliminary calculations, if possible (time-invariant case)
    if ( !tvW ) {
        tmp <- La.svd(mod$W,nu=0)
        Dw <- sqrt(tmp$d) 
        Dw <- pmax(Dw, eps)
        Dw.inv <- 1/Dw
        sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1) 
        if ( !tvGG ) tG.Winv <- t(mod$GG) %*% crossprod(sqrtWinv)
    }
    ## generate last theta
    theta[n+1,] <- mod$m[n+1,] + mod$U.C[[n+1]] %*% matrix(mod$D.C[n+1,]*rnorm(p))
    ## generate all the other theta's
    for (i in (n:1))
    {
        if ( tvGG )
            mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
        if ( tvW ) {
            mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
            tmp <- La.svd(mod$W,nu=0)
            Dw <- sqrt(tmp$d) 
            Dw <- pmax(Dw, eps)
            Dw.inv <- 1/Dw
            sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1) 
        }
        if ( tvGW )
            tG.Winv <- t(mod$GG) %*% crossprod(sqrtWinv)
        D.inv <- 1/mod$D.C[i,]; D.inv[abs(D.inv)==Inf] <- 0
        tmp <- La.svd(rbind(sqrtWinv %*% mod$GG %*% mod$U.C[[i]],
                            diag(x=D.inv,nrow=length(D.inv))), nu=0)
        U.H <- mod$U.C[[i]] %*% t(tmp$vt)
        D.H <- 1/tmp$d; D.H[abs(D.H)==Inf] <- 0
        h <- mod$m[i,] + crossprod(D.H*t(U.H)) %*%
            tG.Winv %*% as.matrix(theta[i+1,]-mod$a[i,])
        theta[i,] <- h + U.H %*% matrix(D.H*rnorm(p))
    }
    if (!is.null(mtsp))
    {
        theta <- drop(theta)
        tsp(theta) <- mtsp
        class(theta) <- if (p > 1) c("mts","ts") else "ts" 
    }
    return(theta=theta)
}



rwishart <- function(df, p = nrow(SqrtSigma), SqrtSigma = diag(p))
{
    ## generate a Wishart-distributed matrix - from S-news, due to B. Venables
    ## note: Sigma = crossprod(SqrtSigma), and df must be integer
    if((Ident <- missing(SqrtSigma)) && missing(p))
        stop("either p or SqrtSigma must be specified")
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
    if(p > 1)
    {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    if(Ident)
        crossprod(Z)
    else
        crossprod(Z %*% SqrtSigma)
}


dlmSvd2var <- function(u, d)
{
    if (is.matrix(u)) {
        if ( nrow(u) != length(d) || ncol(u) != length(d) )
            stop("inconsistent dimensions")
        return( crossprod( d * t(u) ) )
    }
    if (is.list(u)) {
        if ( length(u) != NROW(d) )
            stop("length of 'u' must be equal to the number of rows of 'd'")
        return(lapply(seq(along=u), function(i) crossprod( d[i,] * t(u[[i]]) )))
    }
    stop("wrong argument 'u'")
}

dlmForecast <- function(mod, nAhead=1, method=c("plain","svd"), sampleNew=FALSE) {
    method <- match.arg(method)

    ## NEW
    if (is.dlmFiltered(mod)) {
        modFuture <- mod$mod
        lastObsIndex <- NROW(mod$m)
        modFuture$C0 <- dlmSvd2var(mod$U.C[[lastObsIndex]], mod$D.C[lastObsIndex,])
        if (is.ts(mod$m))
            modFuture$m0 <- window(mod$m, start=end(mod$m))
        else
            modFuture$m0 <- mod$m[lastObsIndex,]
        mod <- modFuture
    }
    ## END NEW

    if (method == "svd") {
        stop("Method \"svd\" is not available yet")
        ans <- .Call("dlmForecast", mod, as.integer(nAhead), PACKAGE="dlm")
        names(ans) <- c("a", "U.R", "D.R", "f", "U.Q", "D.Q")
        ans$a <- ans$a[-1,,drop=FALSE]
        ans$U.R <- ans$U.R[-1]
        ans$D.R <- ans$D.R[-1,,drop=FALSE]
        if ( sampleNew ) {
            newStates <- vector("list", sampleNew)
            newObs <- vector("list", sampleNew)
            newS <- matrix(0, nAhead, p)
            newO <- matrix(0, nAhead, m)
            tmp <- La.svd(mod$V,nu=0)
            Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
            tmp <- La.svd(mod$W,nu=0)
            Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
            for (i in 1:sampleNew) {
                newS[1,] <- ans$U.R[[1]] %*% rnorm(p, sd=ans$D.R[nAhead,]) + ans$a[1,]
                newO[1,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + mod$FF %*% newS[1,]
                if ( nAhead > 1 )
                    for (it in 2:nAhead) {
                        newS[it,] <- crossprod(Ut.W, rnorm(p, sd=D.W)) + mod$GG %*% newS[it-1,] 
                        newO[it,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + mod$FF %*% newS[it,]            
                    }
                newStates[[i]] <- newS
                newObs[[i]] <- newO
            }      
            ans$newStates <- newStates
            ans$newObs <- newObs
        }
    } else {
        ytsp <- tsp(mod$m0)
        p <- length(mod$m0)
        m <- nrow(mod$FF)
        a <- rbind(mod$m0, matrix(0,nAhead,p))
        R <- vector("list",nAhead+1)
        R[[1]] <- mod$C0
        f <- matrix(0,nAhead,m)
        Q <- vector("list", nAhead)
        for (it in 1:nAhead) {
            a[it+1,] <- mod$GG %*% a[it,]
            R[[it+1]] <- mod$GG %*% R[[it]] %*% t(mod$GG) + mod$W
            f[it,] <- mod$FF %*% a[it+1,]
            Q[[it]] <- mod$FF %*% R[[it+1]] %*% t(mod$FF) + mod$V
        }
        a <- a[-1,,drop=FALSE]
        R <- R[-1]
        if ( sampleNew ) {
            newStates <- vector("list", sampleNew)
            newObs <- vector("list", sampleNew)
            newS <- matrix(0, nAhead, p)
            newO <- matrix(0, nAhead, m)
            tmp <- La.svd(mod$V,nu=0)
            Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
            tmp <- La.svd(mod$W,nu=0)
            Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
            for (i in 1:sampleNew) {
                tmp <- La.svd(R[[1]],nu=0)
                newS[1,] <- crossprod(tmp$vt, rnorm(p, sd=sqrt(tmp$d))) + a[1,]
                newO[1,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + mod$FF %*% newS[1,]
                if ( nAhead > 1 )
                    for (it in 2:nAhead) {
                        newS[it,] <- crossprod(Ut.W, rnorm(p, sd=D.W)) + mod$GG %*% newS[it-1,] 
                        newO[it,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + mod$FF %*% newS[it,]            
                    }
                newStates[[i]] <- newS
                newObs[[i]] <- newO
            }      
            if (!is.null(ytsp)) {
                a <- ts(a, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
                f <- ts(f, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
                newStates <- lapply(newStates, function(x)
                                    ts(x, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]))
                newObs <- lapply(newObs, function(x)
                                    ts(x, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]))
            }
            ans <- list(a=a, R=R, f=f, Q=Q, newStates=newStates, newObs=newObs)
        } else {
            if (!is.null(ytsp)) {
                a <- ts(a, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
                f <- ts(f, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
            }
            ans <- list(a=a, R=R, f=f, Q=Q)
        }
    }

    return(ans)
}

###### One-step forecast errors
residuals.dlmFiltered <- function(object, ...,
                                  type=c("standardized", "raw"), sd=TRUE) {
    type <- match.arg(type)
    if (is.null(object$mod$JFF)) tvFF <- FALSE else tvFF <- TRUE
    if (is.null(object$mod$JV)) tvV <- FALSE else tvV <- TRUE
    FF <- object$mod$FF
    if (!( tvFF || tvV )) { ## constant model
        f <- object$a %*% t(FF)
        res <- drop(object$y - f) # one-step forecasting errors
        if (sd || (type == "standardized")) {
            V <- object$mod$V
            SD <- drop(t(sqrt(sapply(seq(along=object$U.R),
                                     function(i)
                                     diag(crossprod(object$D.R[i,] *
                                                    t(FF%*%object$U.R[[i]])) + V)))))
        }
    } else 
    if ( !tvFF ) { ## only V time-varying
        f <- object$a %*% t(FF)
        res <- drop(object$y - f) # one-step forecasting errors
        if (sd || (type == "standardized")) {
            nz <- object$mod$JV != 0
            JV <- cbind(row(object$mod$JV)[nz], col(object$mod$JV)[nz],
                        object$mod$JV[nz])
            V <- object$mod$V
            getSD <- function(i) {
                V[JV[,-3,drop=FALSE]] <- object$mod$X[i,JV[,3]]
                diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V)
            }
            SD <- drop(t(sqrt(sapply(seq(along=object$U.R), getSD))))
        }
    } else 
    if ( !tvV ) { ## only FF time-varying
        if (!(sd || (type == "standardized"))) {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            getFore <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                FF %*% object$a[i,]
            }
            f <- drop(t(sapply(seq(along=object$U.R), getFore)))
            res <- drop(object$y - f) # one-step forecasting errors
        } else {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            V <- object$mod$V
            getBoth <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                c(FF %*% object$a[i,],
                  diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V))
            }
            tmp <- t(sapply(seq(along=object$U.R), getBoth))
            m <- ncol(tmp) / 2
            res <- drop(object$y - tmp[,1:m])
            SD <- drop(sqrt(tmp[,-(1:m)]))
        }
    } else { ## both FF and V time-varying 
        if (!(sd || (type == "standardized"))) {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            getFore <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                FF %*% object$a[i,]
            }
            f <- drop(t(sapply(seq(along=object$U.R), getFore)))
            res <- drop(object$y - f) # one-step forecasting errors
        } else {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            nz <- object$mod$JV != 0
            JV <- cbind(row(object$mod$JV)[nz], col(object$mod$JV)[nz],
                        object$mod$JV[nz])
            V <- object$mod$V
            getBoth <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                V[JV[,-3,drop=FALSE]] <- object$mod$X[i,JV[,3]]
                c(FF %*% object$a[i,],
                  diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V))
            }
            tmp <- t(sapply(seq(along=object$U.R), getBoth))
            m <- ncol(tmp) / 2
            res <- drop(object$y - tmp[,1:m])
            SD <- drop(sqrt(tmp[,-(1:m)]))
        }
    }
    
    if ( type == "standardized" )
        res <- res / SD
    if (sd) {
        if (is.ts(res)) attributes(SD) <- attributes(res) # makes a time series of SD
        return(list(res=res, sd=SD))
    } else
    return(res)
}

###### Diagnostic plots
tsdiag.dlmFiltered <- function (object, gof.lag = 10, ...) {
    stdres <- residuals(object, sd=FALSE)
    if ((ns <- NCOL(stdres)) == 1) {
        oldpar <- par(mfrow = c(3, 1))
        on.exit(par(oldpar))
        plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
        abline(h = 0)
        acf(stdres, plot = TRUE, main = "ACF of Residuals", 
            na.action = na.pass)
        nlag <- gof.lag
        pval <- numeric(nlag)
        for (i in 1:nlag) pval[i] <- Box.test(stdres, i, type = "Ljung-Box")$p.value
        plot(1:nlag, pval, xlab = "lag", ylab = "p value",
             ylim = c(0, 1), main = "p values for Ljung-Box statistic")
        abline(h = 0.05, lty = 2, col = "blue")
    } else {
        ask <- dev.interactive()
        oldpar <- par(mfrow = c(3, 1), oma=c(0, 0, 2, 0), "ask")
        on.exit(par(oldpar))
        hasNames <- !is.null(nm <- attr(stdres,"dimnames")[[2]])
        for (j in 1:ns) {
            plot(stdres[,j], type = "h", main = "Standardized Residuals", ylab = "")
            abline(h = 0)
            if (ask) par(ask=FALSE)
            acf(stdres[,j], plot = TRUE, main = "ACF of Residuals", 
                na.action = na.pass)
            nlag <- gof.lag
            pval <- numeric(nlag)
            for (i in 1:nlag) pval[i] <- Box.test(stdres[,j], i, type = "Ljung-Box")$p.value
            plot(1:nlag, pval, xlab = "lag", ylab = "p value",
                 ylim = c(0, 1), main = "p values for Ljung-Box statistic")
            abline(h = 0.05, lty = 2, col = "blue")
            mtext(if (hasNames) nm[j] else paste("Series",j), line=1, outer=TRUE)
            if (ask) par(ask=TRUE)
        }
    }
}

###### Generating a random DLM
dlmRandom <- function(m, p, nobs = 0, JFF, JV, JGG, JW)
{
    ### Assume for now that each of FF, V, GG, W is either fixed or
    ### time-varying in every entry
    FF <- matrix(rnorm(m*p),m,p)
    V <- rwishart(2*m, m)
    GGtrial <- matrix(rnorm(p*p),p,p)
    e <- eigen(GGtrial)
    if ((ab <- max(abs(e$values))) > 1) {
        r <- runif(1)
        GG <- with(e, Re(vectors %*% (r * values / ab * solve(vectors))))
    }
    else
        GG <- GGtrial
    W <- rwishart(2*p, p)
    m0 <- rep(0,p)
    C0 <- diag(nrow = p) * 100
    if (nobs > 0) {
        count <- 0
        if (hasArg(JGG) && JGG == TRUE) {
            tvGG <- TRUE
            JGG <- structure(1:(p*p), dim=dim(GG))
            count <- p * p
        }
        else {
            tvGG <- FALSE
            JGG <- NULL
        }
        if (hasArg(JW) && JW == TRUE) {
            tvW <- TRUE
            JW <- structure(count + 1:(p*p), dim=dim(W))
            count <- count + p * p
        }
        else {
            tvW <- FALSE
            JW <- NULL
            tmp <- La.svd(W,nu=0)
            Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
        }
        if (hasArg(JFF) && JFF == TRUE) {
            tvFF <- TRUE
            JFF <- structure(count + 1:(m*p), dim=dim(FF))
            count <- count + m * p
        }
        else {
            tvFF <- FALSE
            JFF <- NULL
        }
        if (hasArg(JV) && JV == TRUE) {
            tvV <- TRUE
            JV <- structure(count + 1:(m*m), dim=dim(V))
        }
        else {
            tvV <- FALSE
            JV <- NULL
            tmp <- La.svd(V,nu=0)
            Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
        }
        if (any(c(tvFF, tvV, tvGG, tvW))) {
            newStates <- matrix(0, nr = nobs + 1, nc = p)
            newObs <- matrix(0, nr = nobs, nc = m)
            X <- matrix(0, nr = nobs,
                        nc = m * p * tvFF + m * m * tvV +
                        p * p * tvGG + p * p * tvW)
            tmp <- La.svd(C0,nu=0)
            newStates[1,] <- crossprod(tmp$vt, rnorm(p, sd=sqrt(tmp$d)))
            for (it in 1:nobs) {
                count <- 0
                if (tvGG) {
                    GGtrial <- matrix(rnorm(p*p),p,p)
                    e <- eigen(GGtrial)
                    if ((ab <- max(abs(e$values))) > 1) {
                        r <- runif(1)
                        GG <- with(e, Re(vectors %*% (r * values / ab * solve(vectors))))
                    }
                    else
                        GG <- GGtrial
                    X[it, 1:(p*p)] <- GG
                    count <- p * p
                }
                if (tvW) {
                    W <- rwishart(2*p, p)
                    X[it, count + 1:(p*p)] <- W
                    count <- count + p * p
                    tmp <- La.svd(W,nu=0)
                    Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
                }
                if (tvFF) {
                    FF <- matrix(rnorm(m*p),m,p)
                    X[it, count + 1:(m*p)] <- FF
                    count <- count + m * p
                }
                if (tvV) {
                    V <- rwishart(2*m, m)
                    X[it, count + 1:(m*m)] <- V
                    tmp <- La.svd(V,nu=0)
                    Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
                }
                newStates[it + 1,] <- crossprod(Ut.W, rnorm(p, sd=D.W)) + GG %*% newStates[it,] 
                newObs[it,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + FF %*% newStates[it+1,]
            }
            mod <- dlm(list(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W, JFF=JFF, JV=JV,
                            JGG=JGG, JW=JW, X=X))
            ans <- list(mod = mod, theta = newStates[-1,,drop=FALSE], y = newObs)
        }
        else {
            mod <- dlm(list(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W))
            tmp <- dlmForecast(mod, nAhead = nobs, sampleNew = 1)
            ans <- list(mod = mod, theta = tmp$newStates[[1]], y = tmp$newObs[[1]])
        }
    }
    else
        ans <- dlm(list(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W))
    return(ans)
}
