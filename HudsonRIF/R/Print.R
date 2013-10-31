
print.Hudson = function(x, ...)
{
  # eSet information
  cat("$eSet\n")
  print(x$eSet); cat("\n")
  # DElist information
  cat("$DElist\n")
  print(head(x$DElist))
  if(length(x$DElist) > 6){cat("...\n")}
  cat(" (Total:", length(x$DElist),"values)\n\n")
  # conds information
  cat("$conds\n")
  print(unlist(x$conds))
  cat(" (", length(x$conds)," conditions)\n\n", sep="")
  # EiAB information
  cat("$EiAB\n")
  print(head(x$EiAB))
  if(nrow(x$EiAB) > 6){cat("...\n")}
  cat(" (",paste(dim(x$EiAB), collapse="x"), " matrix)\n\n", sep="")
  # Ai information
  cat("$Ai\n")
  print(head(x$Ai))
  if(length(x$Ai) > 6){cat("...\n")}
  cat(" (",paste(length(x$Ai)), " numeric)\n\n", sep="")
  # dEi information
  cat("$dEi\n")
  print(head(x$dEi))
  if(length(x$dEi) > 6){cat("...\n")}
  cat(" (",paste(length(x$dEi)), " numeric)\n\n", sep="")
  # PIFi information
  cat("$PIFi\n")
  print(head(x$PIFi))
  if(length(x$PIFi) > 6){cat("...\n")}
  cat(" (",paste(length(x$PIFi)), " numeric)\n\n", sep="")
  # rAij information
  cat("$rAij\n")
  print(x$rAij[1:min(5,nrow(x$rAij)),1:min(5,ncol(x$rAij))])
  if(nrow(x$rAij) > 6 | ncol(x$rAij) > 6){cat("...\n")}
  cat(" (",paste(dim(x$rAij), collapse="x"), " matrix)\n\n", sep="")
  # rBij information
  cat("$rBij\n")
  print(x$rBij[1:min(5,nrow(x$rBij)),1:min(5,ncol(x$rBij))])
  if(nrow(x$rBij) > 6 | ncol(x$rBij) > 6){cat("...\n")}
  cat(" (",paste(dim(x$rBij), collapse="x"), " matrix)\n\n", sep="")
  # dCij information
  cat("$dCij\n")
  if(nrow(x$dCij) > 6 | ncol(x$dCij) > 6){cat("...\n")}
  print(x$dCij[1:min(5,nrow(x$dCij)),1:min(5,ncol(x$dCij))])
  cat(" (",paste(dim(x$dCij), collapse="x"), " matrix)\n\n", sep="")
  # RIFi information
  cat("$RIFi\n")
  print(head(x$RIFi))
  if(length(x$Ai) > 6){cat("...\n")}
  cat(" (",paste(length(x$RIFi)), " numeric)\n\n", sep="")
}