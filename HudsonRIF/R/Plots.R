

diffCoexPIF.plot = function(hudson, probe, pch.max = 0)
{
  # Check the validity of user-defined variables
  paramCheck.diffCoexPIF.plot(hudson=hudson, probe=probe, pch.max=pch.max)
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

diffCoex.bgplot = function(hudson)
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

MA.plot = function(hudson, symmetric=TRUE)
{
  # Plot with symmetric Y-axis range
  if(symmetric){plot(x=hudson$Ai, y=hudson$dEi, ylim = rep(max(abs(hudson$dEi)), 2)*c(-1,1),
  xlab = "A", ylab = "M")}
  # Plot with Y(axis range fitted to data
  else{plot(x=hudson$Ai, y=hudson$dEi)}
  abline(h=0, col="red")
}

paramCheck.diffCoexPIF.plot = function(hudson, probe, pch.max)
{
  # hudson is a working output of the Hudson wrapper
  if(class(hudson) != "Hudson"){
    stop("\"hudson\" (", class(hudson),") is not an Hudson object.", call.=F)
  }
  # Note: $eSet, $classCol, $regulator.list, EiAB, Ai, dEi, dCij not necessary
  if(any(!c("conds", "DElist", "PIFi", "rAij", "rBij", "RIFi") %in% names(hudson))){
    stop("\"hudson\" (", names(hudson),") is not a valid Hudson object.
         One item of $conds, $DElist, $PIFi, $rAij, $rBij or $RIFi is missing.", call.=F)
  }
  # Check that probe is a valid probeID
  if(!is.character(probe)){
    stop("\"probe\" (", probe,") is not a character object.", call.=F)
  }
  if(! probe %in% rownames(hudson$rAij)){
    stop("\"probe\" (", probe,") is not a valid probe ID in rAij rows.", call.=F)
  }
  if(! probe %in% rownames(hudson$rBij)){
    stop("\"probe\" (", probe,") is not a valid probe ID in rBij rows.", call.=F)
  }
  # pch.max should be a positive integer
  if(!is.numeric(pch.max)){
    stop("\"pch.max\" (", pch.max,") is not a numeric object.", call.=F)
  }
  if(pch.max < 0){
    stop("\"pch.max\" (", pch.max,") should be a positive integer.", call.=F)
  }
}

