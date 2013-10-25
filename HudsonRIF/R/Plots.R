

diffCoexPIF.plot = function(hudson, probe, pch.max = 0)
{
  # Filters and normalise the the maximum circle size if required
  if (pch.max == 0) {
    pch.max = max(hudson$PIFi)
  }
  # Scales the data points to the maximal size defined
  PIF.DE = hudson$PIFi[hudson$DElist]/max(hudson$PIFi) * pch.max
  # plots the coexpression values in each condition
  plot(x = hudson$rBij[probe, ], y = hudson$rAij[probe, ], xlim = c(-1, 1),
       ylim = c(-1, 1), cex = PIF.DE, xlab = hudson$conds$B,
       ylab = hudson$conds$A, pch = 20, main = probe)
  # adds the identity line (no coexpression change between conditions)
  abline(coef = c(0, 1), col = "green", lwd = 2)
}

diffCoex.bg.plot = function(hudson)
{
  # plots the coexpression values in each condition
  plot(x=hudson$rBij, y=hudson$rAij, xlim=c(-1,1), ylim=c(-1,1),
       xlab = hudson$conds$B, ylab = hudson$conds$A, pch=20, col="grey")
  # adds the identity line (no coexpression change between conditions)
  abline(coef=c(0,1), col="green", lwd=2)
}

