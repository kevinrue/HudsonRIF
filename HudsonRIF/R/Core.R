

Hudson = function(eSet, contrast, classCol, DElist, abs.PIF=TRUE, regulator.list="")
{
  # Checks the validity of user-defined variables
  cat("Checking input variables...", fill=T)
  paramCheck.Hudson(eSet, contrast, classCol, DElist, abs.PIF, regulator.list)
  # Extracts the individual components of the contrast
  cat("Identifying conditions to compare...", fill=T)
  conds = list(A=as.character(contrast[[2]]), B=as.character(contrast[[3]]))
  # Identifies the list of regulators to consider (default: all features not in DElist)
  cat("Defining list of potential regulators...", fill=T)
  if(any(regulator.list == "")){regulator.list = rownames(exprs(eSet))[!rownames(exprs(eSet)) %in% DElist]}
  # Calculate the average expression of each gene in each conditions
  cat("Computing average gene expression in each condition...", fill=T)
  EiAB = calculateEiAB(eSet=eSet, classCol=classCol, conds$A, conds$B)
  # Calculates the average abundance of each gene across conditions
  cat("Computing average gene abundance across condition...", fill=T)
  Ai = (EiAB[,conds$A] + EiAB[,conds$B]) / 2
  # Calculate the differential expression log2 fold change for each gene in contrast
  cat("Computing differential expression between conditions...", fill=T)
  dEi = EiAB[,conds$A] - EiAB[,conds$B]
  # Calculate the PIF
  cat("Computing Phenotypic Impact Factor...", fill=T)
  PIFi = Ai*dEi
  # Calculate the coexpression of all pairs of genes in each condition
  cat("Computing all pairwise co-expression values in condition A...(This may take a moment)", fill=T)
  rAij = cor(x=t(exprs(eSet[,which(pData(eSet)[,classCol] == conds$A)])), method="spearman")
  rAij = rAij[regulator.list, DElist]
  cat("Computing all pairwise co-expression values in condition B...(This may take a moment)", fill=T)
  rBij = cor(x=t(exprs(eSet[,which(pData(eSet)[,classCol] == conds$B)])), method="spearman")
  rBij = rBij[regulator.list, DElist]
  # Calculate the difference in coexpression between the two conditions
  cat("Computing differential co-expression between conditions...", fill=T)
  dCij = rAij - rBij
  # Calculate the RIF
  cat("Computing Regulatory Impact Factor...", fill=T)
  RIFi = calculateRIFi(PIFi=PIFi, dCij=dCij, DElist=DElist, abs.PIF=abs.PIF)
  output = list(eSet=eSet, DElist=DElist, classCol=classCol, conds=conds, regulator.list=regulator.list,
                EiAB=EiAB, Ai=Ai, dEi=dEi, PIFi=PIFi, rAij=rAij, rBij=rBij, dCij=dCij, RIFi=RIFi)
  class(output) = "Hudson"
  cat("Done.", fill=T)
  return(output)
}

# subset an expression dataset according to its phenodata
calculateEiAB = function(eSet, classCol, A, B)
{
  # Calculates the mean expression of features in samples from class A
  EiAB = calculateEiA(eSet, classCol, A)
  # Calculates the mean expression of features in samples from class B
  EiAB = cbind(EiAB, calculateEiA(eSet, classCol, B))
  # Relevant column names
  colnames(EiAB) = c(A,B)
  # Returns the 2-column matrix
  return(EiAB)
}

calculateEiA = function(eSet, classCol, A)
{
  # Calculates the mean expression of features in samples from a given class
  return(apply(X=exprs(eSet[,which(pData(eSet)[,classCol] == A)]), MARGIN=1, FUN="mean"))
}

calculateRIFi = function(PIFi, dCij, DElist, abs.PIF=T)
{
  # Subset the PIF values to the DE genes
  DE.PIFi = PIFi[DElist]
  # Squares the values
  dCij.coefs = dCij ^ 2
  # The matrix product return the RIF values (absolute PIF differs from original formula)
  if(abs.PIF == T){return(apply(X=((dCij.coefs %*% abs(DE.PIFi)) / length(DElist)), MARGIN=1, FUN=sum))}
  # Original Hudson formula using non-absolute PIF value
  else{return(apply(X=((dCij.coefs %*% DE.PIFi) / length(DElist)), MARGIN=1, FUN=sum))}
}

paramCheck.Hudson = function(eSet, contrast, classCol, DElist, abs.PIF, regulator.list)
{
  ## eSet is an expressionSet
  if(class(eSet) != "ExpressionSet"){
    stop("\"eSet\" (", class(eSet),") is not an ExpressionSet.", call.=F)}
  # Contrast should be a formula
  if(class(contrast) != "formula"){
    stop("\"contrast\" (", class(contrast),") is not a formula.", call.=F)}
  ## contrast respect format A~B
  if(length(contrast[[2]]) > 1 ){
    stop("\"contrast\" (", paste(as.character(contrast[c(2,3)]), collapse=" ~ ") ,") contains more than one term on the left side.", call.=F)}
  if(length(contrast[[3]]) > 1){
    stop("\"contrast\" (", paste(as.character(contrast[c(2,3)]), collapse=" ~ ") ,") contains more than one term on the right side.", call.=F)}
  ## classCol is a valid column name in pData(eSet)
  if(!classCol %in% names(pData(eSet))){
    stop("\"classCol\" (", classCol ,") is not a valid column name in \"pData(eSet)\"", call.=F)}
  # A and B should be valid levels of classCol
  conds = list(A=as.character(contrast[[2]]), B=as.character(contrast[[3]]))
  if(!conds$A %in% levels(as.factor(pData(eSet)[,classCol]))){
    stop("(", conds$A ,") is not a valid class level in \"pData(eSet)[,", classCol,"]\"", call.=F)}
  if(!conds$B %in% levels(as.factor(pData(eSet)[,classCol]))){
    stop("(", conds$B ,") is not a valid class level in \"pData(eSet)[,", classCol,"]\"", call.=F)}
  ## DElist fully included in list of features
  if(sum(DElist %in% rownames(exprs(eSet))) != length(DElist)){
    stop("\"DElist\ contains ", sum(!DElist %in% rownames(exprs(eSet))) ," feature(s) absent from \"rownames(exprs(eSet))\"", call.=F)
  }
  ## abs.PIF is logical
  if(!is.logical(abs.PIF)){
    stop("\"abs.PIF\" (", class(abs.PIF),") is not a logical", call.=F)
  }
  # If regulator.list was user-defined
  if(any(regulator.list != "")){
    ## regulator.list fully included in list of features
    if(sum(!regulator.list %in% rownames(exprs(eSet))) > 0){
      stop("\"regulator.list\ contains ", sum(!regulator.list %in% rownames(exprs(eSet))) ," feature(s) absent from \"rownames(exprs(eSet))\"", call.=F)
    }
    ## no overlap between DElist and regulator.list
    if(sum(DElist %in% regulator.list) > 0){
      stop("\"DElist\ and regulator.list overlap by (", sum(DElist %in% regulator.list),") feature(s) ", call.=F)
    }
  }
}