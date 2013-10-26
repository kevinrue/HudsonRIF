

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

diffExpr.plot = function(hudson, formula, xlim=0, ylim=0)
{
  # Extracts from the formula the names of the conditions to compare
  probeX = as.character(formula[[3]])
  probeY = as.character(formula[[2]])
  # Broaden the scale to at least one log2 fold-change
  if(xlim == 0) {xlim = range(exprs(hudson$eSet)[probeX,which(
    pData(hudson$eSet)[,hudson$classCol] %in% unlist(hudson$conds))])+c(-0.5,0.5)}
  if(ylim == 0) {ylim = range(exprs(hudson$eSet)[probeY,which(
    pData(hudson$eSet)[,hudson$classCol] %in% unlist(hudson$conds))])+c(-0.5,0.5)}
  # Data points for condition B
  plot(formula = formula, data = as.data.frame(t(exprs(
    hudson$eSet[,which(pData(hudson$eSet)[,hudson$classCol] == hudson$conds$B)]))),
       main = paste(hudson$conds$B,": Black\n",hudson$conds$A,": Red"), xlim = xlim,
       ylim = ylim, pch=16)
  # Linear regression for condition B
  abline(lm(formula = formula,  data = as.data.frame(t(
    exprs(hudson$eSet[,which(pData(hudson$eSet)[,hudson$classCol] == hudson$conds$B)])))))
  # Data points for condition A
  points(formula = formula, data = as.data.frame(t(exprs(
    hudson$eSet[,which(pData(hudson$eSet)[,hudson$classCol] == hudson$conds$A)]))),
         col = "red", pch = 17)
  # Linear regression for condition A
  abline(lm(formula = formula,  data = as.data.frame(t(exprs(
    hudson$eSet[,which(pData(hudson$eSet)[,hudson$classCol] == hudson$conds$A)])))),
         col = "red")
}

MA.plot = function(hudson, symmetric=T)
{
  if(symmetric){plot(x=hudson$Ai, y=hudson$dEi, ylim = rep(max(abs(hudson$dEi)), 2)*c(-1,1),
  xlab = "A", ylab = "M")}
  else{plot(x=hudson$Ai, y=hudson$dEi)}
  abline(h=0, col="red")
}