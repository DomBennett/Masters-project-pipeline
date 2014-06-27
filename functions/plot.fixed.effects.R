plot.fixed.effects <- function(model, main, sub = I(formula(model))){
  stopifnot(class(model)=="mer")
  
  split.effect.names<-function(effect.names, factors){
    #Splits effect names into LHS (factor) and RHS (level) for later labelling
    to.return <- data.frame(full.name=effect.names, factor=rep(NA, length(effect.names)), 
                            level=rep(NA, length(effect.names)))
    for (i in 1:length(factors)){
      matches <- grep(factors[i], effect.names)  
      to.return$factor[matches] <- factors[i]
      to.return$level[matches] <- substring(effect.names[matches], nchar(factors[i]) + 1)
    }
    return(to.return)
  }  
  
  factor_names <- attr(attr(model@frame,"terms"),"term.labels")
  effect.names <- split.effect.names(names(fixef(model)),factor_names)
  if(is.null(model@call$family)){
    y <- fixef(model)[-1]
    ysmall <- fixef(model) - 1.96 * summary(model)@coefs[,2]
    ybig <- fixef(model) + 1.96 * summary(model)@coefs[,2]
    ylab.1 <- ""
    baseline = 0
  }else{	
    if (model@call$family == "poisson"){
      y <- exp(fixef(model))[-1]
      ysmall <- exp(fixef(model) - 1.96 * summary(model)@coefs[,2])
      ybig <- exp(fixef(model) + 1.96 * summary(model)@coefs[,2])
      ylab.1 <- "back-transformed "
      baseline <- 1
    }else{
      print("Not implemented yet")
      stop("Plot not implemented yet for this family")
    }
  }
  max.y <- ceiling(max(y))
  min.y <- floor(min(y))
  baseline.level <- levels(model@frame[,2])[1]
  ylab <- paste(ylab.1, "fits (", baseline.level," = ",baseline,")", sep="")
  sub.parts <- as.character(formula(model))
  sub <- paste(sub.parts[2], sub.parts[1], sub.parts[3])
  
  cex <- 1/sqrt(par("mfrow")[2])
  #cex<-1.5
  plot(y, ylim=c(min.y, max.y), xlim=c(0.75, length(y)+0.25),ylab=ylab,xlab="Fixed effect", pch=19, main=main, sub=sub, cex.lab=cex,
       cex.axis = cex, cex.sub=cex)
  abline(h=baseline,col="red")
  for (i in 2:length(fixef(model))){
    lines(c(i-1,i-1),c(ysmall[i],ybig[i]))
    text(i-1, min.y, effect.names$factor[i], srt=90, adj=0,col="gray", cex=cex)
    text(i-1, max.y, effect.names$level[i], srt=90, adj=1,col="gray", cex=cex)
  }
}
