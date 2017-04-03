##confidence of class - number of replicates

fallrep <- xtabs(Number.1 ~ paste(Station_code, Replicate, sep = ":") + SpeciesForam, data = forams, drop = TRUE)
fallrep.meta <- as.data.frame(matrix(unlist(strsplit(rownames(fallrep), ":")), ncol = 2, byrow = TRUE))
names(fallrep.meta) <- c("station", "replicate")


fallrep.meta$mO2 <- chem$mO2[sapply(fallrep.meta$station, function(g)which(rownames(chem) == g))]
fallrep.shan <- exp(div(fallrep))

x11()
plot(fallrep.meta$mO2, fallrep.shan, type = "n", xlab = expression(paste("Two Year Minimum ", O[2])), ylab = expression(paste("Foram  ", exp(H * minute[bc]))))
sapply(((1:(length(fallrep.shan) / 3)) - 1) * 3 + 1, function(n)
  points(fallrep.meta$mO2[n:(n + 2)], fallrep.shan[n:(n + 2)], type = "o", pch = 20)
)
points(chem$mO2, exp(div(fall)), col = 2, pch = 4)


repsd <- tapply(exp(div(fallrep)), fallrep.meta$station, sd)
repmean <- tapply(exp(div(fallrep)), fallrep.meta$station, mean)
qrdata <- data.frame(mO2 = chem$mO2, repsd = repsd)[!is.na(chem$mO2), ]
plot(chem$mO2, repsd)
plot(chem$mO2, repmean)

library(quantreg)
mod <- rq(repsd ~ mO2, tau = .95, data = qrdata)    #why won't this work
mod <- lm(repsd ~ mO2, data = qrdata)    #
mod <- lm(repsd ~ repmean)    #

x11(4, 4)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0))
plot(repmean, repsd, xlab = "Mean effective number of species", ylab = "Standard Deviation of effective no. of species")
abline(mod)



boundaries <- c(bottom = 1, bad = 5, poor = 10, moderate = 15, good = 20, top = 25)


pclasses <- function(xs, n) {
  sapply(xs, function(x) {
    exsd <- predict(mod, newdata = data.frame(repmean = x))
    better <- pnorm(boundaries, x, exsd / sqrt(n), lower = F)
    p <- list()
    p$bad = 1 - better["bad"]
    p$poor = 1 - better["poor"] - p$bad
    p$moderate = 1 - better["moderate"] - p$poor - p$bad
    p$good = 1 - better["good"] - p$moderate - p$poor - p$bad
    p$high = 1 - p$good - p$moderate - p$poor - p$bad
    p <- unlist(p)
    p
  })
}
pmisclass <- function(pclass, xs) {
  unlist(lapply(1:5, function(n) {
    if (n == 1) {
      x <- xs >= boundaries[n] & xs <= boundaries[n + 1]
    } else{
      x <- xs > boundaries[n] & xs <= boundaries[n + 1]
    }
    mis <- 1 - pclass[n, x]
  }))
}

xs <- seq(1, 25, .1)
p1 <- pclasses(xs, 1)
q1 <- pmisclass(p1, xs)
p3 <- pclasses(xs, 3)
q3 <- pmisclass(p3, xs)
p6 <- pclasses(xs, 6)
q6 <- pmisclass(p6, xs)


x11()
par(mfrow = c(3, 1), mar = c(0, 3, .5, 1), oma = c(3, 0, 1, 0), mgp = c(1.5, .5, 0))
matplot(xs, t(p1), type = "l", xaxt = "n",  xlab = "", ylab = "Confidence of Class, n=1")
abline(v = boundaries)
matplot(xs, t(p3), type = "l", xaxt = "n", xlab = "", ylab = "Confidence of Class, n=3")
abline(v = boundaries)

matplot(xs, cbind(q1, q3, q6), type = "l", xaxt = "n", xlab = "", ylab = "Missclassification rate")
abline(v = boundaries)
axis(1, outer = TRUE)
title(xlab = "Exponential Shannon", outer = TRUE)
