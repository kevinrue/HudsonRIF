
print.Hudson = function(x, ...)
{
  # eSet information
  cat("$eSet\n")
  print(x$eSet); cat("\n")
  # DElist information
  cat("$DElist\n")
  cat(paste(cat(x$DElist[1:5]), "\n...\n (Total:", length(x$DElist),"values)\n\n"))
  # conds information
  cat("$conds\n")
  print(unlist(x$conds))
  cat(" (", length(x$conds)," conditions)\n\n", sep="")
  # EiAB information
  cat("$EiAB\n")
  print(head(x$EiAB));
  cat("...\n (",paste(dim(x$EiAB), collapse="x"), " matrix)\n\n", sep="")
  # Ai information
  cat("$Ai\n")
  print(x$Ai[1:5])
  cat("...\n (",paste(length(x$Ai)), " numeric)\n\n", sep="")
  # dEi information
  cat("$dEi\n")
  print(x$dEi[1:5])
  cat("...\n (",paste(length(x$dEi)), " numeric)\n\n", sep="")
  # PIFi information
  cat("$PIFi\n")
  print(x$PIFi[1:5])
  cat("...\n (",paste(length(x$PIFi)), " numeric)\n\n", sep="")
  # rAij information
  cat("$rAij\n")
  print(x$rAij[1:5,1:5])
  cat("...\n (",paste(dim(x$rAij), collapse="x"), " matrix)\n\n", sep="")
  # rBij information
  cat("$rBij\n")
  print(x$rBij[1:5,1:5])
  cat("...\n (",paste(dim(x$rBij), collapse="x"), " matrix)\n\n", sep="")
  # dCij information
  cat("$dCij\n")
  print(x$dCij[1:5,1:5])
  cat("...\n (",paste(dim(x$dCij), collapse="x"), " matrix)\n\n", sep="")
  # RIFi information
  cat("$RIFi\n")
  print(x$RIFi[1:5])
  cat("...\n (",paste(length(x$RIFi)), " numeric)\n\n", sep="")
}