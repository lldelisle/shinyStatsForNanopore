generateFullSizeResData <- function(dataFrameWithFiles){
  # It will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  myCols <- c("length", "nMapped", "Coverage", "number")
  fullSizeResData <- NULL
  simplifiedNames <- simplifiedNamesByEnd(dataFrameWithFiles$name)
  for (i in 1:nrow(dataFrameWithFiles)){
    sizeResData1 <- read.delim(dataFrameWithFiles$datapath[i])
    sizeResData1$number <- sizeResData1$nMapped + sizeResData1$nUnmapped
    sizeResData1 <- sizeResData1[, myCols]
    sizeResData1$sample <- simplifiedNames[i]
    fullSizeResData <- rbind(fullSizeResData, sizeResData1)
  }
  return(fullSizeResData)
}

getMergedSizeRes <- function(intervalls, fullSizeResData){
  fakeSizeRes <- data.frame(length = intervalls[-c(1)])
  fakeSizeRes$number <- 0
  fakeSizeRes$Coverage <- 0
  fakeSizeRes$sample <- fullSizeResData$sample[1]
  mergedSizeRes <- rbind(fullSizeResData[, c("length", "number",
                                             "Coverage", "sample")],
                         fakeSizeRes)
  mergedSizeRes$interval <- cut(mergedSizeRes$length, intervalls)
  mergedSizeRes$base <- mergedSizeRes$number * mergedSizeRes$length
  return(mergedSizeRes)
}

getDfIntervalsTotal <- function(mergedSizeRes){
  dfIntervalsPerCov <- aggregate(list(nbOfBases = mergedSizeRes$Coverage),
                                 by = list(interval = mergedSizeRes$interval,
                                           sample = mergedSizeRes$sample),
                                 FUN = sum)
  dfIntervalsPerCov$type <- "Coverage"
  dfIntervalsPerBase <- aggregate(list(nbOfBases = mergedSizeRes$base),
                                  by = list(interval = mergedSizeRes$interval,
                                            sample = mergedSizeRes$sample),
                                  FUN = sum)
  dfIntervalsPerBase$type <- "Sequenced"
  dfIntervalsTotal <- rbind(dfIntervalsPerBase, dfIntervalsPerCov)
  return(dfIntervalsTotal)
}

getIndividualMatSizePerBase <- function(dfIntervalsTotal){
  myTest <- cast(dfIntervalsTotal, interval~type~sample, value = "nbOfBases")
  matSizePerBase <- matrix(myTest, nrow = dim(myTest)[1],
                           ncol = prod(dim(myTest)[2:3]))
  matSizePerBase[is.na(matSizePerBase)] <- 0
  rownames(matSizePerBase) <- rownames(myTest)
  colnames(matSizePerBase) <- paste0(rep(dimnames(myTest)$sample,
                                         each = dim(myTest)[2]),
                                     "_",
                                     rep(colnames(myTest), dim(myTest)[3]))
  if (all(grepl("_single_", dimnames(myTest)$sample))){
    mySamplesOrder <- dimnames(myTest)$sample[order(
      sapply(dimnames(myTest)$sample,
             function(s){
               v <- strsplit(s, "_")[[1]]
               return(as.numeric(v[which(v == "single") + 1]))
             }))]
    matSizePerBase <- matSizePerBase[, paste0(rep(mySamplesOrder,
                                                  each = dim(myTest)[2]),
                                              "_",
                                              rep(colnames(myTest),
                                                  dim(myTest)[3]))]
  }
  matSizePerBase <- matSizePerBase[rev(levels(dfIntervalsTotal$interval)), ]
  return(matSizePerBase)
}



plotReadVsBase <- function(intervalls, fullSizeResData, firstLineTitle){
  fullSizeResData <- fullSizeResData[order(fullSizeResData$length), ]
  v1 <- unlist(apply(fullSizeResData, 1, function(v){
    rep(v["length"],
        v["number"])
  }))
  csvr1 <- cumsum(as.numeric(rev(v1)))
  v2 <- unlist(apply(fullSizeResData[, grep("sample",
                                            colnames(fullSizeResData),
                                            invert = T)],
                     1, function(v){
                       c(rep(v["Coverage"] / v["nMapped"], v["nMapped"]),
                         rep(0, v["number"] - v["nMapped"]))
                     }))
  csvr2 <- cumsum(as.numeric(rev(v2)))
  vList <- list(v1, v2)
  names(vList) <- c("Sequenced", "Coverage")
  csvrList <- list(csvr1, csvr2)
  myPValue <- c(0.5, 0.75, 0.9)
  xmax <- length(vList[[1]])
  ymax <- tail(csvrList[[1]], 1)
  layout(mat = matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE),  height = c(1, 8),
         widths = c(8, 2))
  par(mar = c(5, 4, 4, 2) + 0.1)
  plot(csvrList[[1]], type = "l", xlab = "Number of reads",
       ylab = "Number of bases",
       main = paste0(firstLineTitle, "Cumulative bases"))
  for (i in 1:length(vList)){
    if (i != 1){
      lines(csvrList[[i]], col = i)
    }
    for (j in 1:length(myPValue)){
      myValue <- tail(csvrList[[i]], 1) * myPValue[j]
      myLty <- j + 1
      myX <- which.min(abs(csvrList[[i]] - myValue))
      segments(x0 = 0, x1 = myX, y0 = myValue, lty = myLty, col = i)
      segments(myX, myValue, y1 = 0, lty = myLty, col = i)
      myLabel <- paste0("top ", myX, "\n =", round(myX / xmax * 100), "%")
      myAdj <- c(1, 0)
      if (i != 1){
        myLabel <- paste0(myLabel, "\n", rev(vList[[1]])[myX], " bases")
        myAdj <- c(0, 1)
      }
      text(x = myX, y = myValue,
           labels = myLabel,
           adj = myAdj, col = i)
    }
  }
  legend(x = xmax / 2,
         y = ymax / 2,
         legend = names(vList),
         lty = 1,
         col = 1:length(vList))
  mergedSizeRes <- getMergedSizeRes(intervalls, fullSizeResData)
  dfIntervalsPerRead <- aggregate(mergedSizeRes$number,
                                  by = list(interval = mergedSizeRes$interval),
                                  FUN = sum)
  dfIntervalsPerRead <- dfIntervalsPerRead[order(dfIntervalsPerRead$interval,
                                                 decreasing = T), ]
  matSizePerRead <- matrix(dfIntervalsPerRead$x,
                           nrow = length(levels(dfIntervalsPerRead$interval)),
                           ncol = 1)
  rownames(matSizePerRead) <- rev(levels(dfIntervalsPerRead$interval))
  colorUsed <- rainbow(100)[seq(1,
                               100,
                               length.out = nrow(matSizePerRead) +
                                 1)[1:nrow(matSizePerRead)]]
  par(mar = c(1, 4, 0, 2) + 0.1)
  b <- barplot(matSizePerRead, horiz = T,
               col = colorUsed, axes = T,
               xlim = c(-xmax * 0.02, xmax * 1.02))
  dfIntervalsTotal <- getDfIntervalsTotal(mergedSizeRes)
  myTest2 <- cast(dfIntervalsTotal, interval ~ type, value = "nbOfBases",
                  fun.aggregate = sum)
  rownames(myTest2) <- myTest2$interval
  myTest2 <- myTest2[rev(levels(dfIntervalsPerRead$interval)), names(vList)]
  myMat <- as.matrix(myTest2)
  rownames(myMat) <- rownames(myTest2)
  colnames(myMat) <- colnames(myTest2)
  par(mar = c(5, 0, 4, 0) + 0.1)
  b <- barplot(myMat, horiz = F, col = colorUsed,
               axes = T, ylim = c(-ymax * 0.05, ymax * 1.05),
               las = 2)
  cs <- cumsum(myMat[, 1])
  xs <- c(0, cs[-length(cs)]) + myMat[, 1] / 2
  for (i in 1:length(cs)){
    if (myMat[i, 1] / sum(myMat[, 1]) > 0.05){
      text(y = xs[i],
           x = b[1],
           labels = rownames(myMat)[i],
           adj = c(0, 0.5))
    }
  }
}
