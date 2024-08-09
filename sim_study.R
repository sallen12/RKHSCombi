################################################################################
##### simulation study

N <- 10000

a1 <- 1
a2 <- 1
a3 <- 1.1

X0 <- rnorm(N)
X1 <- rnorm(N)
X2 <- rnorm(N)
X3 <- rnorm(N)
eps <- rnorm(N)

Y <- X0 + a1*X1 + a2*X2 + a3*X3 + eps


################################################################################
##### forecast distributions

G <- function(x, a1, a2, a3) pnorm(x, 0, sqrt(2 + a1^2 + a2^2 + a3^2))

F1 <- function(x, X0, X1, a1, a2, a3) pnorm(x, X0 + a1*X1, sqrt(1 + a2^2 + a3^2))
F2 <- function(x, X0, X2, a1, a2, a3) pnorm(x, X0 + a2*X2, sqrt(1 + a1^2 + a3^2))
F3 <- function(x, X0, X3, a1, a2, a3) pnorm(x, X0 + a3*X3, sqrt(1 + a1^2 + a2^2))


Q1 <- function(u, X0, X1, a1, a2, a3) qnorm(u, X0 + a1*X1, sqrt(1 + a2^2 + a3^2))
Q2 <- function(u, X0, X2, a1, a2, a3) qnorm(u, X0 + a2*X2, sqrt(1 + a1^2 + a3^2))
Q3 <- function(u, X0, X3, a1, a2, a3) qnorm(u, X0 + a3*X3, sqrt(1 + a1^2 + a2^2))


F_LP <- function(x, X0, X1, X2, X3, a1, a2, a3, w1 = 1/3, w2 = 1/3, w3 = 1/3) {
  w1*F1(x, X0, X1, a1, a2, a3) + w2*F2(x, X0, X2, a1, a2, a3) + w3*F3(x, X0, X3, a1, a2, a3)
}
Q_LP <- function(u, X0, X1, X2, X3, a1, a2, a3, w1 = 1/3, w2 = 1/3, w3 = 1/3, x = seq(-10, 10, 0.01)) {
  uu <- sapply(x, function(xx) F_LP(xx, X0, X1, X2, X3, a1, a2, a3, w1, w2, w3))
  if (length(u) > 1) {
    xx <- sapply(1:length(u), function(i) x[which.min(abs(u[i] - uu[i, ]))])
  } else {
    xx <- x[which.min(abs(u - uu))]
  }
  
  return(xx)
}


Q_QA <- function(u, X0, X1, X2, X3, a1, a2, a3, w1 = 1/3, w2 = 1/3, w3 = 1/3) {
  w1*Q1(u, X0, X1, a1, a2, a3) + w2*Q2(u, X0, X2, a1, a2, a3) + w3*Q3(u, X0, X3, a1, a2, a3)
}
F_QA <- function(x, X0, X1, X2, X3, a1, a2, a3, w1 = 1/3, w2 = 1/3, w3 = 1/3, u = seq(0.01, 0.99, 0.01)) {
  xx <- sapply(u, function(uu) Q_QA(uu, X0, X1, X2, X3, a1, a2, a3, w1, w2, w3))
  if (length(x) > 1) {
    uu <- sapply(1:length(x), function(i) u[which.min(abs(x[i] - xx[i, ]))])
  } else {
    uu <- u[which.min(abs(x - xx))]
  }
  
  return(uu)
}


################################################################################
##### probabilistic calibration

library(WeightedForecastVerification)

Z1 <- F1(Y, X0, X1, a1, a2, a3)
Z2 <- F2(Y, X0, X2, a1, a2, a3)
Z3 <- F3(Y, X0, X3, a1, a2, a3)
Z_LP <- F_LP(Y, X0, X1, X2, X3, a1, a2, a3)
Z_QA <- F_QA(Y, X0, X1, X2, X3, a1, a2, a3)


pit1 <- pit_hist(Z1, ranks = F, ymax = 0.3, xlab = NULL, ylab = NULL, 
                 title = expression(f[1]), xticks = F, yticks = F)
pit2 <- pit_hist(Z2, ranks = F, ymax = 0.3, xlab = NULL, ylab = NULL,
                 title = expression(f[2]), xticks = F, yticks = F)
pit3 <- pit_hist(Z3, ranks = F, ymax = 0.3, xlab = NULL, ylab = NULL,
                 title = expression(f[3]), xticks = F, yticks = F)
pit4 <- pit_hist(Z_LP, ranks = F, ymax = 0.3, xlab = NULL, ylab = NULL,
                 title = expression(f[LP]), xticks = F, yticks = F)
pit5 <- pit_hist(Z_QA, ranks = F, ymax = 0.3, xlab = NULL, ylab = NULL,
                 title = expression(f[QA]), xticks = F, yticks = F)

pits <- gridExtra::grid.arrange(pit1, pit2, pit3, pit4, pit5, nrow = 1)
ggsave("Figures/sim_pit_1e6.png", plot = pits, height = 2, width = 10)


################################################################################
##### marginal calibration

U <- runif(N)

V1 <- Q1(U, X0, X1, a1, a2, a3)
V2 <- Q2(U, X0, X2, a1, a2, a3)
V3 <- Q3(U, X0, X3, a1, a2, a3)
V_LP <- Q_LP(U, X0, X1, X2, X3, a1, a2, a3)
V_QA <- Q_QA(U, X0, X1, X2, X3, a1, a2, a3)


pit_hist(G(V1, a1, a2, a3), ranks = F, ymax = 0.3, xlab = NULL, title = expression(f[1]))
pit_hist(G(V2, a1, a2, a3), ranks = F, ymax = 0.3, xlab = NULL, title = expression(f[2]))
pit_hist(G(V3, a1, a2, a3), ranks = F, ymax = 0.3, xlab = NULL, title = expression(f[3]))
pit_hist(G(V_LP, a1, a2, a3), ranks = F, ymax = 0.3, xlab = NULL, title = expression(f[LP]))
pit_hist(G(V_QA, a1, a2, a3), ranks = F, ymax = 0.3, xlab = NULL, title = expression(f[QA]))

mc_plot <- function(x1, x2, z, xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {
  dif <- sapply(z, function(zz) mean(x1 <= zz) - mean(x2 <= zz))
  df <- data.frame(z = z, d = dif)
  ggplot(df) + geom_line(aes(x = z, y = d)) +
    geom_hline(aes(yintercept = 0), linetype = "dotted", color = "red") +
    scale_x_continuous(name = xlab, limits = xlims) +
    scale_y_continuous(name = ylab, limits = ylims) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    ggtitle(title)
}

mcp1 <- mc_plot(V1, Y, z = seq(-10, 10, 0.01), ylims = c(-0.1, 0.1), title = expression(f[1]))
mcp2 <- mc_plot(V2, Y, z = seq(-10, 10, 0.01), ylims = c(-0.1, 0.1), title = expression(f[2]))
mcp3 <- mc_plot(V3, Y, z = seq(-10, 10, 0.01), ylims = c(-0.1, 0.1), title = expression(f[3]))
mcp4 <- mc_plot(V_LP, Y, z = seq(-10, 10, 0.01), ylims = c(-0.1, 0.1), title = expression(f[LP]))
mcp5 <- mc_plot(V_QA, Y, z = seq(-10, 10, 0.01), ylims = c(-0.1, 0.1), title = expression(f[QA]))

mcps <- gridExtra::grid.arrange(mcp1, mcp2, mcp3, mcp4, mcp5, nrow = 1)
ggsave("Figures/sim_mcp_1e6.png", plot = mcps, height = 2, width = 10)

var(Y)
var(V1)
var(V2)
var(V3)
var(V_LP)
var(V_QA)


