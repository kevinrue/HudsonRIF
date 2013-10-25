

diffCoexPIF.plot = function(hudson, probe, pch.max = 0)
{
  # Filters and normalise the the maximum circle size if required
  if (pch.max == 0) {
    pch.max = max(hudson$PIFi)
  }
  # Scales the data points to the maximal size defined
  PIF.DE = hudson$PIFi[hudson$DElist]/max(hudson$PIFi) * pch.max
  # plot(x=rAij.coefs, y=rBij.coefs, xlim=c(-1,1), ylim=c(-1,1),  pch=20, col="grey") # background coexpressions (heavy)
  plot(x = hudson$rBij[probe, ], y = hudson$rAij[probe, ], xlim = c(-1, 
                                                                    1), ylim = c(-1, 1), cex = PIF.DE, xlab = hudson$conds$B, 
       ylab = hudson$conds$A, pch = 20, main = probe)
  abline(coef = c(0, 1), col = "green", lwd = 2)
}