require(deSolve)
require(primer)
###############################################################################
#### Examples from <A primer of ecology with R>
###############################################################################

### 3.1 Discrete density-dependent (Logistic) growth
dlogistic <- function(alpha = 0.01, rd = 1, N0 = 2, t = 15) {
  N <- c(N0, numeric(t))
  for (i in 1:t) N[i + 1] <- 
    N[i] + rd * N[i] * (1 - alpha * N[i])
  return(N)
}

Nts = dlogistic()
t = 15; a = 0.01
plot(0:t, Nts)
abline(h = 1/a, lty = 3)

total.incr = Nts[1:t + 1] - Nts[1:t]
per.capita.incr = total.incr / Nts[1:t]
plot(Nts[1:t], total.incr)
plot(Nts[1:t], per.capita.incr)

N0s = c(0, runif(30) * 1.1 * 1/a)
N <- sapply(N0s, function(n) dlogistic(N0 = n))
matplot(0:t, N, type = 'l', lty = 1, lwd = 0.75, col = 1)
text(t, 1/a, expression(italic('K') == 1/alpha), adj = c(1, 0))

# Influence of rd on the logistic growth
rd.v = seq(1.3, 2.8, by = 0.3)
t = 15
Ns = sapply(rd.v, function(r) dlogistic(rd = r, t = t))
matplot(0:t, Ns, type = "l", col = 1)

# Presentation of limit cycles
tmp = data.frame(rd = as.factor(rd.v), t(Ns))
Ns2 = reshape(tmp, varying = list(2:ncol(tmp)), idvar = 'rd',
              v.names = 'N', direction = 'long')
str(Ns2)
library(lattice)
print(xyplot(N ~ time | rd, data = Ns2, type = "l", layout = c(3, 2, 1),col = 1))

# Attractors (equilibrium state) as a function of rd (Bifurcation)
num.rd = 201; t = 400
rd.s = seq(1, 3, length = num.rd)
tmp <- sapply(rd.s, function(r) dlogistic(rd = r, N0 = 99, t = t))
tmp.s = stack(as.data.frame(tmp))
names(tmp.s) = c('N', 'Old.Column.ID')
tmp.s$rd = rep(rd.s, each = t + 1)
tmp.s$time = rep(0:t, num.rd)
N.bif = subset(tmp.s, time > 0.5 * t)
plot(N ~ rd, data = N.bif, pch = ".", xlab = quote("r"["d"]))

# Sensitivity to initial conditions
N.init <- c(97, 98, 99); t <- 30
Ns <- sapply(N.init, function(n0) dlogistic(rd = 2.7, N0 = n0, t = t))
matplot(0:t, Ns, type = "l", col = 1)

# Boundedness, variance of species abundance 
# The upper and lower bound increase with rd systematically
# This is the persistence in one specie case.

### 3.2 continuous logistic growth

# local stability around the equilibria
# Growth rate vs. N
pop.growth.rate = expression(r * N * (1 - alpha * N))
r <- 1; alpha <- 0.01; N <- 0:120
plot(N, eval(pop.growth.rate), type = "l",
     ylab = "Population Growth Rate (dN/dt)", xlab = "N")
abline(h = 0); legend("topright", "r=1", lty = 1)
N <- c(0, 10, 50, 100, 115)
points(N, eval(pop.growth.rate), cex = 1.5)
text(N, eval(pop.growth.rate), letters[1:5], adj = c(0.5, 2))
arrows(20, 2, 80, 2, length = 0.1, lwd = 3)
arrows(122, -2, 109, -2, length = 0.1, lwd = 3)

# Symbolic differntiation
dF.dN = deriv(pop.growth.rate, 'N')
N <- c(0, 1/alpha)
eval(dF.dN)

##### 4 Populations in Space
##### 4.1 Source-sink Dynamics
L1 <- 2; L2 <- 0.4
A <- matrix(c(1, 0, L1 -1, L2), nrow = 2, byrow = TRUE)
eigen(A)
L1s <- seq(1, 3, by = 0.01)
p1 <- sapply(L1s, function(l1) {
  A[2, 1] <- l1 - 1
  eigen(A)$vectors[1, 1]/sum(eigen(A)$vectors[, 1])
})
plot(L1s, p1, type = "l", ylab = "Source Population",xlab = expression(lambda[1]))

##### 4.6 Core-Satellite Simulations
out.CS.10 <- MetaSim(method = 'hanski', NSims = 10)
matplot(out.CS.10$t, out.CS.10$Ns, type = "l", xlab = "Time",
        ylab = "Occupancy", sub = out.CS.10$method)
out.CS.Lots <- MetaSim(method = 'hanski', NSims = 50, Time = 1000)
hist(out.CS.Lots$Ns[1001, ], breaks = 10, main = NULL, 
     xlab = expression("Occupancy (" * italic("p") * ")"),
     ylab = "Number of Species",
     sub = paste(out.CS.Lots$method, " Model", sep = ""))
system.time(out.L.Lots <- MetaSim(NSims = 50, Time = 500, method = "levins"))
hist(out.L.Lots$Ns[501,],breaks = 10, 
     xlab = expression("Occupancy(" *italic("p") * ")"), 
     ylab = "Number of Species", main = NULL,
     sub = paste(out.L.Lots$method, " Model", sep = ""))



#### 5.1 Two competition species: discrete and continuous time models
# Discrete model
dlvcomp2 <- function(N, alpha, rd = c(1, 1)) {
  N1.t1 <- N[1] + rd[1] * N[1] * (1- alpha[1, 1] * N[1] - alpha[1, 2] * N[2])
  N2.t1 <- N[2] + rd[2] * N[2] * (1- alpha[2, 1] * N[1] - alpha[2, 2] * N[2])
  c(N1.t1, N2.t1)
}
alphs = matrix(c(0.01, 0.005, 0.008, 0.01), ncol = 2, byrow = T)
t = 20
N = matrix(NA, nrow = t + 1, ncol = 2)
N[1, ] = c(10, 10)
for (i in 1:t) N[i + 1, ] = dlvcomp2(N[i, ], alphs)
matplot(0:t, N, type = "l", col = 1, ylim = c(0, 110))
abline(h = 1/alphs[1, 1], lty = 3)
text(0, 1/alphs[1, 1], "K", adj = c(0, 0))
legend("right", c(expression("Sp.1 " * (alpha[21] == 0.008)),
                  expression("Sp.2 " * (alpha[12] == 0.005))), lty = 1:2, bty = "n")

# Continuous model
lvcomp2 <- function(t, n, parms) {
  with(as.list(parms), {
    dn1dt = r1 * n[1] * (1 -a11 * n[1] - a12 * n[2])
    dn2dt = r2 * n[2] * (1 -a21 * n[1] - a22 * n[2])
    list(c(dn1dt, dn2dt))
  })
}
library(deSolve)
parms <- c(r1 = 1, r2 = 0.1, a11 = 0.2, a21 = 0.1, a22 = 0.02, a12 = 0.01)
initialN <- c(2, 1)
out <- ode(y = initialN, times = 1:100, func = lvcomp2, parms = parms)
matplot(out[, 1], out[, -1], type = "l")


#### 8. Multiple Basins of Attraction
####
lvcompg <- function(t, n, parms) {
  r = parms[[1]]
  a = parms[[2]]
  dns.dt = r * n * (1 - a %*% n)
  list(c(dns.dt))
}

r = c(r1 = 0.6, r2 = 1, r3 = 2)
a <- matrix(c(a11 = 0.001, a12 = 0.002, a13 = 0.002, a21 = 0.002,
              a22 = 0.00101, a23 = 0.002, a31 = 0.002, a32 = 0.002,
              a33 = 0.00102), nrow = 3, ncol = 3)
parms = list(r, a)

t = seq(0,40,by = 0.1)
ni = 200
std = 10
N0 = sapply(1:30, function(i) rnorm(3, mean = ni, sd = std))
N0[, 1] = ni

par(mar = c(2, 2, 1, 1))
layout(matrix(1:30, ncol = 5))
for (i in 1:30) {
  lvout <- ode(N0[, i], t, lvcompg, parms)
  matplot(t, lvout[, 2:4], type = "l", lwd = 1.5, col = 1)
  if (all(N0[, i] == 200)) {
    text(3, 500, 'Equal', srt = 90)
  }
  else {
    text(3, 500, paste('Sp,', which.max(N0[, i])), srt = 90)
  }
  lastN = lvout[nrow(lvout), 2:4]
  text(3, max(lastN), paste("Sp.", which.max(lastN)), adj = c(0,1))
}


p <- c(N = 1, as = 0.01, af = 0.01, b = 0.02, qs = 0.075,
       qf=0.005,hs=0,hf=0.2,ls=0.05,lf=0.05,rs=0.5,
       rf=0.5,W=0)
t <- 1:200
Initial <- c(F = 10, S = 10)
S.out1 = ode(Initial, t, scheffer, p)
matplot(t, S.out1[, -1], type='l')
legend('right', c('F', 'S'), lty= 1:2, bty = 'n' )
p['N'] = 4
S.out2 = ode(Initial, t, scheffer, p)
matplot(t, S.out2[, -1], type='l')

N.s = seq(0.5, 4, by = 0.1)
t = 1:1000
S.s <- t(sapply(N.s, function(x) {
  p["N"] <- x
  ode(Initial, t, scheffer, p)[length(t), 2:3]
}))

matplot(N.s, S.s, type = "l")
legend("right", c("F", "S"), lty = 1:2, bty = "n")

Initial.Eutrophic <- c(F = 600, S = 10)
S.s.E <- t(sapply(N.s, function(x) {
  p["N"] <- x
  ode(Initial.Eutrophic, c(1, 1000), scheffer, p)[2, 2:3]
}))
matplot(N.s, S.s.E, type = "l")
legend("right", c("F", "S"), lty = 1:2, bty = "n")



###############################################################################
#### Examples from Tutorial of <differential equations in R> on RUser!2011 
###############################################################################

library(deSolve)
library(scatterplot3d)

model <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <- r * N * (1 - N / K)
    list(dN)
  })
}
y = c(N = 0.1)
parms = c(r = 0.1, K = 10)
times = seq(from = 0, to = 100, by = 1)
out = ode(y, times, model, parms)
plot(out)

rigidode <- function(t, y, parms) {
  dy1 <- -2 * y[2] * y[3]
  dy2 <- 1.25* y[1] * y[3]
  dy3 <- -0.5* y[1] * y[2]
  list(c(dy1, dy2, dy3))
}
yini <- c(y1 = 1, y2 = 0, y3 = 0.9)
times <- seq(from = 0, to = 20, by = 0.01)
out <- ode (times = times, y = yini, func = rigidode, parms = NULL)
par(mar = c(0, 0, 0, 0))
scatterplot3d(out[,-1], xlab = "", ylab = "", zlab = "", label.tick.marks = FALSE)

################################################################################
##
## Example for setting absolute tolerance (atol) to zero
## - see deSolve vignette:
## Soetaert, K, Petzoldt, T, Setzer, W: "Solving Initial Value Differential
## Equations in R". http://cran.r-project.org/package=deSolve
##
################################################################################

library(deSolve)

## A Lotka-Volterra model with producer (P) and consumer (C)
PCmod <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dP <- c*P   - d*C*P      # producer
    dC <- e*P*C - f*C        # consumer
    res <- c(dP, dC)
    list(res)
  })
}

xstart <- c(P = 0.5, C = 1)
times  <- seq(0, 190, 0.1)
parms  <- c(c = 0.1, d = 0.1, e = 0.1, f = 0.1)

out <- ode(xstart, times, PCmod, parms)
plot(out)
## normal solution

## now we set a rather extreme parameter set
parms  <- c(c = 5, d = 0.1, e = 0.1, f = 0.1)

out <- ode(xstart, times, PCmod, parms)
tail(out)
## we get NaN for both state variables!!!

plot(out)
## shows that the solution was unstable

diagnostics(out)
## the solver method did not even detect that there was something wrong

## reduce absolute tolerance to a very small value (zero)
## so that only relative tolerances are allowed
out <- ode(xstart, times, PCmod, parms, atol = 0)

plot(out)
## shows steep peaks but is solved without problems

## Note that it makes absolutely no sense to set both, atol and rtol to zero




#################Lorenz equations###############################################
chaos <- function(t, state, parameters) {
  with(as.list(c(state)), {
    dx <- -8/3 * x + y * z
    dy <- -10 * (y - z)
    dz <- -x * y + 28 * y - z
    list(c(dx, dy, dz))
  })
}
yini <- c(x = 1, y = 1, z = 1)
yini2 <- yini + c(1e-6, 0, 0)
times <- seq(0, 100, 0.01)
out <- ode(y = yini, times = times, func = chaos, parms = 0)
out2 <- ode(y = yini2, times = times, func = chaos, parms = 0)
plot(out, out2, col = c("blue", "orange"), main = c("Xvalue", "Yvalue", "Zvalue"),
     xlim = list(c(20, 30), c(25, 30), NULL), mfrow = c(1, 3))
plot(out[,"x"], out[,"y"], pch = ".", main = "Lorenz butterfly",
     xlab = "x", ylab = "y")
XY <- subset(out, select = c("x", "y"), subset = y < 10 & x < 40)
plot(XY, main = "Lorenz butterfly", xlab = "x", ylab = "y", pch = ".")

derivs <- function(t, y, parms)
  with (as.list(parms), list(r * y * (1-y/K)))
parms <- c(r = 1, K = 10)
yini <- c(y = 2)
yini2 <- c(y = 12)
times <- seq(from = 0, to = 30, by = 0.1)
out <- ode(y = yini, parms = parms, func = derivs, times = times)
out2 <- ode(y = yini2, parms = parms, func = derivs, times = times)
plot(out, out2, lwd = 2)

outlist <- list()
plist <- cbind(r = runif(30, min = 0.1, max = 5),
               K = runif(30, min = 8, max = 15))
for (i in 1:nrow(plist))
  outlist[[i]] <- ode(y = yini, parms = plist[i,], func = derivs, times = times)
plot(out, outlist)

obs2 <- data.frame(time = c(1, 5, 10, 20, 25), y = c(12, 10, 8, 9, 10))
obs1 <- data.frame(time = c(1, 5, 10, 20, 25), y = c(1, 6, 8, 9, 10))
plot(out, out2, col = c("blue", "red"), lwd = 2,
     obs = list(obs1, obs2),
     obspar = list(col = c("blue", "red"), pch = 18, cex = 2))

# 
require(ReacTran)
N = 1000
Grid = setup.grid.1D(N = N, L = 100000)
Esturay <- function(t, y, parms) {
  NH3 = y[1:N]
  O2 = y[(N+1):(2*N)]
  tranNH3 = tran.1D(C = NH3, D = 1e7, v = 1000,
                    C.up = 500, C.down = 10, dx =Grid)$dC
  tranO2 = tran.1D(C = O2, D = 1e7, v = 1000,
                   C.up = 100, C.down = 250, dx = Grid)$dC
  r_nit = 0.1 * O2 / (O2 + 1) * NH3
  dNH3 = tranNH3 - r_nit
  dO2 = tranO2 - 2 * r_nit + 0.1 * (300 - O2)
  list(c(dNH3, dO2), r_nit = r_nit)
}
system.time(
  std <- steady.1D(y = runif(2 * N), parms = NULL, names = c('NH3', 'O2'),
                   func = Esturay, dimens = N, positive =TRUE)
)
plot(std, which = c("NH3", "O2", "r_nit"), lwd = 2,
     mfrow = c(3,1), grid = Grid$x.mid, xlab = "distance, m",
     ylab = c("mmol m-3", "mmol m-3", "mmol m-3 d-1"))

# rootsolve get Jaconbian matrix from ODE
mod <- function (t = 0, y, parms = NULL) {
  dy1 <- y[1] * (1 + 2 * y[2])
  dy2 <- y[2] * (3 * y[1] + 4 + 5 * y[3])
  dy3 <- y[3] * (6 * y[2] + 7 + 8 * y[4])
  dy4 <- y[4] * (9 * y[3] + 10)
  return(as.list(c(dy1, dy2, dy3, dy4)))
}
jacobian.full(y = c(1, 2, 3, 4), func = mod)

mod.matrix <- function(t = 0, y, parms = parms) {
  r = parms[[1]]
  M = parms[[2]]
  dY <- y * (r + M %*% y)
  list(c(dY))
}
r = c(1, 4, 7, 10)
M = matrix(c(0, 2, 0, 0,
             3, 0, 5, 0,
             0, 6, 0, 8,
             0, 0, 9, 0), ncol = 4, byrow = TRUE)
parms = list(r, M)
jacobian.full(y = c(1, 2, 3, 4), func = mod.matrix, parms = parms)

f2<-function(x) {
  X <- matrix(nrow = 5, x)
  X %*% X %*% X - matrix(nrow = 5, data = 1:25, byrow = TRUE)
}
x <- multiroot(f2, start = 1:25 )$root
X <- matrix(nrow = 5, x)
X%*%X%*%X


#################External variables in dynamic models##############################
times <- seq(0, 100, by = 0.1)
signal <- as.data.frame(list(times = times, import = rep(0, length(times))))
signal$import <- ifelse((trunc(signal$times) %% 2 == 0), 0, 1)
input <- approxfun(signal, rule = 2)
input(seq(from = 0.98, to = 1.01, by = 0.005))

SPCmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    import <- input(t)
    dS <- import - b * S * P + g * C
    dP <- c * S * P - d * C * P
    dC <- e * P * C - f * C
    res <- c(dS, dP, dC)
    list(res, signal = import)
  })
}
parms <- c(b = 0.1, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0)
xstart <- c(S = 1, P = 1, C = 1)
out <- ode(y = xstart, times = times, func = SPCmod, parms)


#########################Events in dynamic models##############################
pharmaco <- function(t, blood, p) {
  dblood = -b * blood
  list(dblood)
}
b = 0.6
yini = c(blood = 0)
injectevents = data.frame(var = 'blood', time = 0:20,
                          value = 40, method = 'add')
times <- seq(from = 0, to = 20, by = 1/24)
outDrug <- ode(func = pharmaco, times = times, y = yini, parms = NULL,
               method = "impAdams",
               events = list(data = injectevents))


##########An event triggered by a root: A Bouncing Ball#########################
ball <- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8
  list(c(dy1, dy2))
}
yini <- c(height = 0, velocity = 10)

rootfunc <- function(t, y, parms) return (y[1])
eventfunc <- function(t, y, parms) {
  y[1] <- 0
  y[2] <- -0.9*y[2]
  return(y)
}
times <- seq(from = 0, to = 20, by = 0.01)
out <- ode(times = times, y = yini, func = ball,
           parms = NULL, rootfun = rootfunc,
           events = list(func = eventfunc, root = TRUE))
attributes(out)$troot
for (i in seq(1, 2001, 10)) {
  plot(out, which = "height", type = "l", lwd = 1,
       main = "", xlab = "Time", ylab = "Height"
  )
  points(t(out[i,1:2]), pch = 21, lwd = 1, col = 1, cex = 2,
         bg = rainbow(30, v = 0.6)[20-abs(out[i,3])+1])
  Sys.sleep(0.01)
}

## =============================================================================
## Logistic growth with harvesting
## =============================================================================

require(deSolve)

derivs <- function(t, y, parms) 
  list(r * y * (1-y/K))

r <- 1 
K <- 10
yini <- c(y = 2)
times <- seq(from = 0, to = 20, by = 0.1)

# First run: unharvested
out1 <- ode(y = yini, times = times, func = derivs, parms = NULL)

# Second run: harvest at preset times

harvest <- data.frame(var = "y",
                      time =  seq(2, 20, by = 2),
                      value = 0.5,
                      method = "multiply")

out2 <- ode(y = yini, times = times, func = derivs, parms = NULL,
            events = list(data = harvest))

# Third run: harvest when critical density is reached

rootfunc  <- function(t, y, p) 
  return(y - 0.8*K)

eventfunc <- function(t, y, p) 
  return(y * 0.5)

out3 <- ode(y = yini, times = times, func = derivs, parms = NULL,
            rootfun = rootfunc, events = list(func = eventfunc, root = TRUE))

# Plot different scenarios

plot(out1, out2, out3, lwd = 2, col = "black")
legend ("bottomright", lwd = 2, lty = 1:3,
        legend = c("unharvested", "2-day harvest", "harvest at 80% of K"))

################################################################################
#### Examples from <A practical guide to ecological modelling using R!>
################################################################################
r1 = 3
r2 = 2
K1 = 1.5
K2 = 2
alpha1 = 2
alpha2 = 1
Lotka <- function(t, N, pars) {
  dN1 <- r1 * N[1] * (1- (N[1] + alpha2 * N[2])/K1)
  dN2 <- r2 * N[2] * (1- (N[2] + alpha1 * N[1])/K2)
  list(c(dN1, dN2))
}


################################################################################
#### Examples from <Ecological models and data in R>
################################################################################

library(emdbook)
library(bbmle)
data(MyxoTiter_sum)
myxdat = subset(MyxoTiter_sum, grade == 1)
gm = mean(myxdat$titer)
cv = var(myxdat$titer)/mean(myxdat$titer)
shape0 = gm/cv
scale0 = cv
shapevec = 10:100
scalevec = seq(0.01, 0.3, by = 0.01)
gammaNLL1 = function(shape, scale) {
  -sum(dgamma(myxdat$titer, shape = shape, scale = scale, log = TRUE))
}
surf = matrix(nrow = length(shapevec), ncol = length(scalevec))
for (i in 1:length(shapevec)) {
  for (j in 1:length(scalevec)) {
    surf[i, j] = gammaNLL1(shapevec[i], scalevec[j])
  }
}
contour(shapevec, scalevec, log10(surf))
curve3d(log10(gammaNLL1(x, y)), from = c(10, 0.01),
        to = c(100, 0.3), n = c(91, 30), sys3d = "image")
gridsearch2d(gammaNLL1, v1min = 10, v2min = 0.01, 
             v1max = 100, v2max = 0.3, logz = TRUE)

library(MCMCpack)
gammaNLL2B <- function(p) {
  sum(dgamma(myxdat$titer, shape = p[1], scale = p[2], log = TRUE))
}
m3 <- MCMCmetrop1R(gammaNLL2B, theta.init = c(shape = 20, scale = 0.05),
                   thin = 30, mcmc = 30000, optim.lower = rep(0.004, 2),
                   optim.method = 'L-BFGS-B', tune = 3)


## using <igraph> package to visualize graph
library(igraph)
A = matrix(c(0, 1, 0, 1, 1,
             1, 0, 1, 1, 0,
             0, 1, 0, 0, 1,
             1, 1, 0, 0, 1,
             1, 0, 1, 1, 0), nrow = 5, byrow = T)
G = graph.adjacency(A)
plot.igraph(G)
