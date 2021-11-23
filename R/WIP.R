# # Load/install dependencies
# CRAN_packs <- c("magrittr", "parallel", "optparse", "plyr", "ggplot2", "grid", "gtools", "reshape2",
#                 "EnvStats", "RColorBrewer")
# BIOC_packs <- c("GenomicAlignments", "GenomicRanges", "rtracklayer")
# sapply(CRAN_packs, function(x) installLoad_CRAN(x))
# sapply(BIOC_packs, function(x) installLoad_BioC(x))


# # Functions for plotting.
# # Score RD: Makes an average of the relative differences of neighboring positions
# scoreRD <- function(CHR, d = 6){
#   scoreRD <- rep(0, length(CHR)); names(scoreRD) <- names(CHR)
#   ds <- NULL; us <- NULL
#   for(P in (d+1):(length(CHR)-d-1)){
#     evalPos <- CHR[P]
#     leftFlank <- CHR[(P-d):(P-1)]
#     rightFlank <- CHR[(P+1):(P+d)]
#     if(sum(leftFlank == 0) > 3 | sum(rightFlank == 0) > 3){next()}
#     us <- mean(((leftFlank - evalPos + 1) / (evalPos + 1))) # Up-stream Relative difference
#     ds <- mean(((rightFlank - evalPos + 1) / (evalPos + 1))) # Downstream Relative difference
#     scoreRD[P] <- mean(c(us, ds))
#   }
#   # scoreRD <- (scoreRD - min(scoreRD)) /  (max(scoreRD) - min(scoreRD))
#   return(scoreRD)
# }
#
# #Gives scoreA given a numeric vector containing read ends
# scoreA <- function(CHR, d = 6, minEnv = 15){
#     scoreA <- rep(NA, length(CHR)); names(scoreA) <- names(CHR)
#     for (P in (d+1):(length(CHR)-d)){
#         evalPos <- CHR[P]
#         leftFlank <- CHR[(P-d):(P-1)]
#         rightFlank <- CHR[(P+1):(P+d)]
#         if( min(c(rightFlank, leftFlank)) < minEnv ){next()}
#         m1 <- mean(leftFlank)
#         m2 <- mean(rightFlank)
#         sd1 <- sd(leftFlank)
#         sd2 <- sd(rightFlank)
#         scoreA[P] <- 1 - (((2*evalPos) + 1) /
#                               ((0.5*abs(m1-sd1)) + evalPos + (0.5*abs(m2-sd2)) + 1))
#     }
#     scoreA[scoreA < 0] <- 0
#     return(scoreA)
# }
#
#
# # ScoreMeth: Quantifies methylation as a percentage.Based on a weighted ratio to
# # the "d" positions flanking both sides (Birkedal et al., 2015)
# scoreMeth <- function(CHR, d = 6){
#   scoreMeth <- rep(0, length(CHR)); names(scoreMeth) <- names(CHR)
#   rfWeigth <- seq(1, 0.5, -0.5/(d-1))[1:d]; lfWeigth <- rev(seq(1,0.5,-0.5/(d-1))[1:d])
#   srfw <- sum(rfWeigth); slfw <- sum(lfWeigth)
#   for(P in (d+1):(length(CHR)-d-1)){
#     evalPos <- CHR[P]
#     leftFlank <- CHR[(P-d):(P-1)]
#     rightFlank <- CHR[(P+1):(P+d)]
#     if(sum(leftFlank == 0) > 3 | sum(rightFlank == 0) > 3){next()}
#     if(evalPos >= 0 & (sum(rightFlank) > 0 | sum(leftFlank) > 0) ){
#       scoreMeth[P] <- 1-(evalPos/(0.5*((sum(lfWeigth*leftFlank)/slfw) + sum(rfWeigth*rightFlank)/srfw)))
#     }
#   }
#   scoreMeth[is.na(scoreMeth)] <- 0
#   return(scoreMeth)
# }
#
# #ScoreZ - based on Z.score
# scoreZ <- function(CHR, d = 6){
#     scoreZ <- rep(NA, length(CHR)); names(scoreZ) <- names(CHR)
#     for (P in (d+1):(length(CHR)-d-1)){
#         evalPos= CHR[P]; leftFlank = CHR[(P-d):(P-1)]; rightFlank = CHR[(P+1):(P+d)]
#         meanNeigh <- mean(c(leftFlank, rightFlank))
#         sdNeigh <- sd(c(leftFlank, rightFlank))
#         scoreZ[P] <- (evalPos - meanNeigh) / sdNeigh
#     }
#     scoreZ <- scoreZ * -1
#     return(scoreZ)
# }
#
# # Stops score
# scoreStop <- function(stopReads, readtReads, minCov = 50){
#   res <- stopReads/readtReads
#   autoFail <- readtReads < minCov
#   res[autoFail] <- 0
#   res <- res[-1]; res <- c(res, 0); names(res) <- names(stopReads) # Move one up-stream
#   res[is.na(res)] <- 0
#   res[is.infinite(res)] <- 0
#   return(res)
# }
#
# # read coverage file and generate CHR object
# readCov <- function(file, chromosome){
#   coverage = read.delim(file, header = F, stringsAsFactors = F)
#   CHR = coverage[coverage[,1] == chromosome, 3]
#   names(CHR) <- coverage[coverage[,1] == chromosome, 2]
#   return(CHR)
# }
#
# # read ALLREADENDS zip file and generate CHR object
# zReadCov <- function(file, chromosome, rColumn){
#   coverage = read.delim(gzfile(file), header = T, stringsAsFactors = F)
#   CHR = coverage[coverage[,1] == chromosome, rColumn]
#   names(CHR) <- coverage[coverage[,1] == chromosome, 2]
#   return(CHR)
# }
#
# # Outputs table of confusion, precision, recall, AUC and MCC of predictor variable with binary classification
# modMetrics <- function(RESULTS, successes){
#     evalRES <- matrix(NA, nrow = 1, ncol = 9)
#     res <- RESULTS; zres <- scale(res)
#     MCC <- matthewCorrCoeff(scores = res, successes =  successes)
#     selThr <- res[which.max(MCC)]; zThr <- zres[which.max(MCC)]
#     calls <- as.numeric(res >= selThr)
#     known <- rep(0, length(calls)); known[successes] <- 1
#     TP <-  sum(calls & known)
#     FP <-  sum(calls) - TP
#     TN <- sum(!calls & !known)
#     FN <- sum(!calls) - TN
#     Prec <- TP/(TP+FP)
#     Reca <- TP/(TP+FN)
#     M <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
#     AUC <- signif(auc(roc(predictions = res, labels = as.factor(known))), 4)
#     eval <- c(TP, TN, FP, FN, Prec, Reca, M, zThr, AUC)
#     evalRES[1,] <- eval
#     colnames(evalRES) <- c("TP", "TN", "FP", "FN", "Precision", "Recall", "MCC", "Z-threshold", "AUC")
#     rownames(evalRES) <- names(RESULTS)
#     return(evalRES)
# }
#
# # Matthew Correlation Coefficient
# matthewCorrCoeff <- function(scores, successes){
#   MCC <- NULL
#     known <- rep(0, length(scores)); known[successes] <- 1
#     NAs <- which(is.na(scores))
#     scores[is.na(scores)] <- 0
#     for (i in 1:length(scores)){
#         thr <- scores[i]
#         calls <- as.numeric(scores >= thr)
#         TP <-  sum(calls & known)
#         FP <-  sum(calls) - TP
#         TN <- sum(!calls & !known)
#         FN <- sum(!calls) - TN
#         MCC[i] <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
#     }
#     MCC[is.nan(MCC)] <- 0
#     MCC[is.infinite(MCC)] <- 0
#     MCC[NAs] <- NA
#   return(MCC)
# }
#
# # Paralelized version of Matthew Correlation coefficient calculation
# parMCC <- function(scores, successes, cluster){
#     known <- rep(0, length(scores)); known[successes] <- 1
#     NAs <- which(is.na(scores))
#     scores[is.na(scores)] <- 0
#     MCC <- parSapply(cluster, 1:length(scores), function(i){
#         thr <- scores[i]
#         calls <- as.numeric(scores >= thr)
#         TP <-  sum(calls & known)
#         FP <-  sum(calls) - TP
#         TN <- sum(!calls & !known)
#         FN <- sum(!calls) - TN
#         return((TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
#     })
#     MCC[is.nan(MCC)] <- 0
#     MCC[is.infinite(MCC)] <- 0
#     MCC[NAs] <- NA
#     return(MCC)
# }
#
# # Binomial proportional confidence interval. For Alpha = 0.05 z = 1.96
# binError <- function(p, z, n){
#     bErr <- z * sqrt((p*(1-p))/n)
#     return(bErr)
# }
#
# # Scatter, density function, and ROC curve plots for assesing scores
# modEvalPlots <- function(RESULTS, successes, MCCthr = F, main = "", multiplot = T){
#     if(multiplot == T) {par(mfrow=c(1,3))}
#     res <- RESULTS
#     MCC <- matthewCorrCoeff(res, successes)
#     selThr <- res[which.max(MCC)]
#     colorvector <- rep("black", length(res)); colorvector[successes] <- "red"
#     par(pch = 16); barplot(res, col = colorvector, border = NA)
#     abline(h = res[which.max(MCC)],col= "blue")
#     plot(ecdf(res[-successes]), main = main, col = "black")
#     lines(ecdf(res[successes]), col = "red")
#     abline(v = res[which.max(MCC)], col = "blue")
#     legend("bottomright", legend = paste("max.MCC =", signif(max(MCC), 4)))
#     if(MCCthr){res <- as.numeric(res >= selThr)}
#     known <- rep(0, length(res)); known[successes] <- 1
#     areaUC <- signif(auc(roc(predictions = res, labels = as.factor(known))), 4)
#     plot(roc(predictions = res, labels = as.factor(known)), col = "red")
#     legend("bottomright", legend = paste("AUC =", areaUC))
#     if(multiplot == T) {par(mfrow=c(1,1))}
# }
# modelEvalPlot <- modEvalPlots #deprecated function name
#
# #Calculate AUC for two vectors (scores, successes)
# modAUC <- function(scores, successes, main = ""){
#     known <- as.numeric(successes)
#     areaUC <- signif(auc(roc(predictions = scores, labels = as.factor(known))), 4)
#     plot(roc(predictions = scores, labels = as.factor(known)), col = "red", main= main)
#     legend("bottomright", legend = paste("AUC =", areaUC))
# }
#
# # Fit the fragmentation and signal score with logistic regression
# integrate_scores <- function(RES_A, RES_B, success){
#   RESULTS <- list()
#   for(i in 1:length(RES_A)){
#     res1 <- RES_A[[i]] #Frag Score
#     res2 <- RES_B[[i]] # Stop Score
#     meth <- rep(0, 1, length(res1)); meth[success] <- 1
#     df <- data.frame(METH = meth, FRAG_SCORE = res1, STOP_SCORE = res2)
#     model <- glm(METH ~., family=binomial(link='logit'), data = df)
#     RESULTS[[i]] <- predict(model, newdata= df ,type='response')
#   }
#   nick <- sapply(rRNA_rEnds, function(x) gsub(pattern = "_Aligned.out.sorted.allreadEnds.gz", replacement = "", x))
#   names(RESULTS) <- nick
#   return(RESULTS)
# }
#
# # Fit the fragmentation and signal score with logistic regression
# integrate_3scores <- function(RES_A, RES_B, RES_C, success){
#   RESULTS <- list()
#   for(i in 1:length(RES_A)){
#     res1 <- RES_A[[i]] # Frag Score
#     res2 <- RES_B[[i]] # Stop Score
#     res3 <- RES_C[[i]] # Z-score
#     meth <- rep(0, 1, length(res1)); meth[success] <- 1
#     df <- data.frame(METH = meth, FRAG_SCORE = res1, STOP_SCORE = res2, Zscore = res3)
#     model <- glm(METH ~., family=binomial(link='logit'), data = df)
#     RESULTS[[i]] <- predict(model, newdata= df ,type='response')
#   }
#   nick <- sapply(rRNA_rEnds, function(x) gsub(pattern = "_Aligned.out.sorted.allreadEnds.gz", replacement = "", x))
#   names(RESULTS) <- nick
#   return(RESULTS)
# }
#
# # Get Coefficients of integrated model
# getCoef <- function(scoreFr, scoreSt, success){
#   meth <- rep(0, 1, length(scoreFr)); meth[success] <- 1
#   df <- data.frame(METH = meth, FRAG_SCORE = scoreFr, STOP_SCORE = scoreSt)
#   model <- glm(METH ~., family=binomial(link='logit'), data = df)
#   return(model$coefficients)
# }
#
# # Fit one score to the logistic regression
# fitOneScore <- function(res, success){
#   meth <- rep(0, 1, length(res)); meth[success] <- 1
#   df <- data.frame(METH = meth, SCORE = res)
#   model <- glm(METH ~., family=binomial(link='logit'), data = df)
#   RESULTS <- predict(model, newdata= df ,type='response')
#   return(RESULTS)
# }
#
# # Write table with TAB separated files standards
# writeMyTable <- function(X, filename){
#   write.table(X, file = filename, sep= "\t", col.names = NA, quote = F)
# }
#
# # Extracting the last n characters from a string in R
# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }
#
# # Logistic Regression
# logisticReg <- function(interc, var1, cf1, var2, cf2){
#   res <- 1/(1 + exp(-(var1*cf1 + var2*cf2 + interc))) # Logistic Regression Formula
#   return(res)
# }
#
# # List files in work dir that match pattern pat
# listFilePatt <- function(pattern, path = "."){
#   files <- list.files(path)[grep(pattern = pattern, x = list.files(path))]
#   return(files)
# }
#
# # Object size
# oSize <- function(x){
#   print(object.size(x), units = "auto")
# }
#
# # Define Chr, Start & End of a GENE based on UTR annotations
# defCoor <- function(gene){
#   g <- gene
#   chr <- geneAnnot[g,1]
#   coor <- geneAnnot[g,2:3]
#   strand <- geneAnnot[g, 6]
#   U3 <- UTR3Annot[UTR3Annot[, "gene"] == g, 2:3]
#   U5 <- UTR5Annot[UTR5Annot[, "gene"] == g, 2:3]
#   start <- NULL; end <- NULL
#   if (strand == "+"){
#     if(nrow(U5) > 0){
#       start <- U5[1]
#     } else {
#       start <- coor[1]
#     }
#     if(nrow(U3) > 0){
#       end <- U3[2]
#     } else {
#       end <- coor[2]
#     }
#   } else if (strand == "-"){
#     if(nrow(U5) > 0){
#       start <- U5[2]
#     } else {
#       start <- coor[2]
#     }
#     if(nrow(U3) > 0){
#       end <- U3[1]
#     } else {
#       end <- coor[1]
#     }
#   }
#   start <- as.integer(start) ; end <- as.integer(end)
#   return(c(chr, start, end, strand))
# }
#
# # Calculate adjusted coefficients
# adjCoeff <- function(F18s, F25s, S18s, S25s, R18s, R25s, medGoal, known){
#   adjF <- mean(median(F18s)/medGoal, median(F25s)/medGoal)
#   score.frag <- c(scoreRD(round(F18s/adjF)), scoreRD(round(F25s/adjF)))
#   score.stop <- c(scoreStop(round(S18s/adjF), round(R18s/adjF)), scoreStop(round(S25s/adjF), round(R25s/adjF)))
#   cf <- getCoef(score.frag, score.stop, known)
#   res <- logisticReg(cf[1], score.frag, cf[2], score.stop, cf[3])
#   maxMCC <- max(matthewCorrCoeff(res, known_RRNA_2Ome))
#   modThresh <- res[which.max(matthewCorrCoeff(res, known))]
#   modVals <- c(cf, maxMCC, modThresh); names(modVals)[4:5] <- c("maxMCC", "Thresh")
#   return(modVals)
# }
#
# # Strand Specific coverage
# strSpCov <- function(covStrand, covOppStrand){
#   cov <- as.numeric(covStrand[,3]) - as.numeric(covOppStrand[,3])
#   cov[cov < 0] <- 0
#   cov <- data.frame(chr= covStrand[,1], coor = covStrand[,2], cov)
#   return(cov)
# }
#
# # Smooth y given x,y coordinates
# smoothLine <- function(x, y){
#   sl <- smooth.spline(x, y, spar = 0.7)
#   y <- sl$y; names(y) <- sl$x
#   return(y)
# }
#
# # Turns a character matrix into a dataframe more "smartly"
# asDataFrame <- function(x){
#   tmpCol <- colnames(x)
#   out <- data.frame(lapply(split(x, col(x)), type.convert, as.is = TRUE),
#              stringsAsFactors = FALSE)
#   colnames(out) <- tmpCol
#   return(out)
# }
#
# # Read rRNA data from rEnd files
# rRNAdata <- function(rEndsFile){
#   F18s <- zReadCov(rEndsFile, "RDN18-1", "r3pos")
#   F25s <- zReadCov(rEndsFile, "RDN25-1", "r3pos")
#   S18s <- zReadCov(rEndsFile, "RDN18-1", "r5pos")
#   S25s <- zReadCov(rEndsFile, "RDN25-1", "r5pos")
#   R18s <- zReadCov(rEndsFile, "RDN18-1", "Fpos")
#   R25s <- zReadCov(rEndsFile, "RDN25-1", "Fpos")
#   return(list(F18s = F18s, F25s = F25s, S18s = S18s, S25s = S25s, R18s = R18s, R25s = R25s))
# }
#
# adjCoeff <- function(h_F18s, h_F25s, l_S18s, l_S25s, l_R18s, l_R25s, h_S18s, h_S25s,
#                      h_R18s, h_R25s, medGoal, known){
#     adjF <- mean(median(h_F18s)/medGoal, median(h_F25s)/medGoal)
#     score.frag <- c(scoreRD(round(h_F18s/adjF)), scoreRD(round(h_F25s/adjF)))
#     score.stop <- c(scoreStop(round(l_S18s/adjF), round(l_R18s/adjF)), scoreStop(round(l_S25s/adjF), round(l_R25s/adjF))) -
#         c(scoreStop(round(h_S18s/adjF), round(h_R18s/adjF)), scoreStop(round(h_S25s/adjF), round(h_R25s/adjF)))
#     score.stop[score.stop < 0] <- 0
#     cf <- getCoef(score.frag, score.stop, known)
#     res <- logisticReg(cf[1], score.frag, cf[2], score.stop, cf[3])
#     maxMCC <- max(matthewCorrCoeff(res, known_RRNA_2Ome))
#     modThresh <- res[which.max(matthewCorrCoeff(res, known))]
#     modVals <- c(cf, maxMCC, modThresh); names(modVals)[4:5] <- c("maxMCC", "Thresh")
#     return(modVals)
# }
#
# # Results Dataframe
# resDF_Thresh <- function(lDATA, hDATA, stopThresh = 0, fragThresh = 0, fCHThresh = 0, minMedF = 0, evalGenes){
#     stopRES <- list()
#     for(eGene in evalGenes){
#         frGeneL <- lDATA$DATA[[eGene]]$frags
#         stGeneL <- lDATA$DATA[[eGene]]$stops
#         rtGeneL <- lDATA$DATA[[eGene]]$readt
#         frGeneH <- hDATA$DATA[[eGene]]$frags
#         stGeneH <- hDATA$DATA[[eGene]]$stops
#         rtGeneH <- hDATA$DATA[[eGene]]$readt
#         frGeneT <- frGeneL + frGeneH
#         medBool <- medNeighBool(frGeneT, minMedF)
#         frScore <- scoreA(frGeneT)
#         frScore2 <- scoreRD(frGeneT)
#         zScore <- scoreZ(frGeneT)
#         stScoreL <- scoreStop(stGeneL+1, rtGeneL+1)
#         stScoreH <- scoreStop(stGeneH+1, rtGeneH+1)
#         fCHstop <- stScoreL/stScoreH
#         # Thresholds
#         callsGene <- which(stScoreL >= stopThresh & frScore >= fragThresh & fCHstop >= fCHThresh & medBool)
#         callsGene <- callsGene[callsGene > 6 & (callsGene <= length(rtGeneL) - 7)]
#         chrGene <- geneAnnot[eGene,1]
#         strandGene <- geneAnnot[eGene,6]
#         if(strandGene == "+"){
#             coorsGene <- geneAnnot[eGene,2]:geneAnnot[eGene,3]
#         } else if(strandGene == "-"){
#             coorsGene <- rev(geneAnnot[eGene,2]:geneAnnot[eGene,3])
#         }
#         bedGene <- list() # Defining yet another dummy variable...
#         if(length(callsGene) > 0){
#             for(l in 1:length(callsGene)){
#                 selRes <- callsGene[l]
#                 neighF <- frGeneT[(selRes-6):(selRes+6)]
#                 nS_l <- stGeneL[selRes+1] + 1 # +1 Position Upstream +1 PseudoCount
#                 nR_l <- rtGeneL[selRes+1] + 1
#                 nS_h <- stGeneH[selRes+1] + 1
#                 nR_h <- rtGeneH[selRes+1] + 1
#                 fisher.pVal <- fisher.test(rbind(c(nR_l,nS_l),c(nR_h,nS_h)), alternative ="l")$p.value
#                 FCH.stop <- (nS_l/nR_l)/(nS_h/nR_h)
#                 bedGene[[l]] <- c(chrGene,
#                                   coorsGene[selRes],
#                                   coorsGene[selRes],
#                                   eGene,
#                                   paste(chrGene, coorsGene[selRes], eGene, sep = "."),
#                                   strandGene,
#                                   frScore[selRes],
#                                   frScore2[selRes],
#                                   zScore[selRes],
#                                   stScoreL[selRes],
#                                   stScoreH[selRes],
#                                   FCH.stop,
#                                   fisher.pVal,
#                                   nS_l,
#                                   nR_l,
#                                   nS_h,
#                                   nR_h,
#                                   paste(neighF, collapse = ","))
#             }
#             stopRES[[eGene]] <- t(sapply(bedGene, c))
#         }
#     }
#     stopRES <- stopRES[!sapply(stopRES, is.null)]
#     if(length(stopRES) == 0){
#         stop("No results found with specified thresholds!")
#     }
#     RESMAT <- NULL
#     for(m in 1:length(stopRES)){
#         RESMAT <- rbind(RESMAT, stopRES[[m]])
#     }
#     RESMAT <- asDataFrame(RESMAT)
#     colnames(RESMAT) <- c("chr", "start", "end", "gene", "id", "strand",
#                           "scoreA", "scoreRD", "scoreZ", "stopScoreL", "stopScoreH",
#                           "FCHstop", "p.Value", "stopsL", "readtL", "stopsH", "readtH",
#                           "fragEnvir")
#     return(RESMAT)
# }
#
# # Parallel version of Results Data frame function
# parResDF_Thresh <- function(lDATA, hDATA, stopThresh = 0, fragThresh = 0, fCHThresh = 0, minMedF = 0, evalGenes, cluster, noTresh = F){
#     stopRES <- parLapply(cluster, evalGenes, function(eGene){
#         frGeneL <- lDATA$DATA[[eGene]]$frags
#         stGeneL <- lDATA$DATA[[eGene]]$stops
#         rtGeneL <- lDATA$DATA[[eGene]]$readt
#         frGeneH <- hDATA$DATA[[eGene]]$frags
#         stGeneH <- hDATA$DATA[[eGene]]$stops
#         rtGeneH <- hDATA$DATA[[eGene]]$readt
#         frGeneT <- frGeneL + frGeneH
#         medBool <- medNeighBool(frGeneT, minMedF)
#         frScore <- scoreA(frGeneT)
#         frScore2 <- scoreRD(frGeneT)
#         zScore <- scoreZ(frGeneT)
#         stScoreL <- scoreStop(stGeneL+1, rtGeneL+1)
#         stScoreH <- scoreStop(stGeneH+1, rtGeneH+1)
#         fCHstop <- stScoreL/stScoreH
#         # Thresholds
#         if(noTresh == F){
#             callsGene <- which(stScoreL >= stopThresh & frScore >= fragThresh & fCHstop >= fCHThresh & medBool)
#         }else {callsGene <- 1:length(frGeneL)}
#         callsGene <- callsGene[callsGene > 6 & (callsGene <= length(rtGeneL) - 7)]
#         chrGene <- geneAnnot[eGene,1]
#         strandGene <- geneAnnot[eGene,6]
#         if(strandGene == "+"){
#             coorsGene <- geneAnnot[eGene,2]:geneAnnot[eGene,3]
#         } else if(strandGene == "-"){
#             coorsGene <- rev(geneAnnot[eGene,2]:geneAnnot[eGene,3])
#         }
#         bedGene <- list() # Defining yet another dummy variable...
#         if(length(callsGene) > 0){
#             for(l in 1:length(callsGene)){
#                 selRes <- callsGene[l]
#                 neighF <- frGeneT[(selRes-6):(selRes+6)]
#                 nS_l <- stGeneL[selRes + 1] + 1 # +1 Position Upstream +1 PseudoCount
#                 nR_l <- rtGeneL[selRes + 1] + 1
#                 nS_h <- stGeneH[selRes + 1] + 1
#                 nR_h <- rtGeneH[selRes + 1] + 1
#                 fisher.pVal <- fisher.test(rbind(c(nR_l - nS_l,nS_l),c(nR_h - nS_h,nS_h)), alternative ="l")$p.value
#                 if(nR_l > 100 & nR_h > 100){
#                   FCH.stop <- (nS_l/nR_l)/(nS_h/nR_h)
#                   }else{FCH.stop <- NA}
#                 bedGene[[l]] <- c(chrGene,
#                                   coorsGene[selRes],
#                                   coorsGene[selRes],
#                                   eGene,
#                                   paste(chrGene, coorsGene[selRes], eGene, sep = "."),
#                                   strandGene,
#                                   frScore[selRes],
#                                   frScore2[selRes],
#                                   zScore[selRes],
#                                   stScoreL[selRes],
#                                   stScoreH[selRes],
#                                   FCH.stop,
#                                   fisher.pVal,
#                                   nS_l,
#                                   nR_l,
#                                   nS_h,
#                                   nR_h,
#                                   paste(neighF, collapse = ","))
#             }
#             return(asDataFrame(t(sapply(bedGene, c))))
#         }
#     })
#     stopRES <- stopRES[!sapply(stopRES, is.null)]
#     if(length(stopRES) == 0){
#         stop("No results found with specified thresholds!")
#     }
#     RESMAT <- rbindlist(stopRES)
#     colnames(RESMAT) <- c("chr", "start", "end", "gene", "id", "strand",
#                           "scoreA", "scoreRD", "scoreZ", "stopScoreL", "stopScoreH",
#                           "FCHstop", "p.Value", "stopsL", "readtL", "stopsH", "readtH",
#                           "fragEnvir")
#     return(RESMAT)
# }
#
#
# # Read rRNA data from rEnd files
# rRNAdata <- function(rEndsFile){
#     F18s <- zReadCov(file.path(INPUTDIR, rEndsFile), "RDN18-1", "r3pos")
#     F25s <- zReadCov(file.path(INPUTDIR, rEndsFile), "RDN25-1", "r3pos")
#     S18s <- zReadCov(file.path(INPUTDIR, rEndsFile), "RDN18-1", "r5pos")
#     S25s <- zReadCov(file.path(INPUTDIR, rEndsFile), "RDN25-1", "r5pos")
#     R18s <- zReadCov(file.path(INPUTDIR, rEndsFile), "RDN18-1", "Fpos")
#     R25s <- zReadCov(file.path(INPUTDIR, rEndsFile), "RDN25-1", "Fpos")
#     return(list(F18s = F18s, F25s = F25s, S18s = S18s, S25s = S25s, R18s = R18s, R25s = R25s))
# }
#
# # Median of neighboring positions above threshold
# medNeighBool <- function(CHR, minMed, d = 6){
#     medN <- rep(F, length(CHR)); names(scoreA) <- names(CHR)
#     for (P in (d+1):(length(CHR)-d-1)){
#         medN[P] <- median(c(CHR[(P-d):(P-1)], CHR[(P+1):(P+d)])) >= minMed
#     }
#     return(medN)
# }
#
# # Vector indicating methylation status for known rRNA sites
# knownStatus <- function(res){
#     real <- NULL
#     for (iter in 1:nrow(res)){
#         real[iter] <- paste(res[iter,1:2], collapse = " ") %in% known2Om_coors
#     }
#     falPos <- NULL
#     for (iter in 1:nrow(res)){
#         falPos[iter] <- paste(res[iter,1:2], collapse = " ") %in% FP_2Om_coors
#     }
#     realChr <- as.character(real)
#     realChr[real] <- "Known"
#     realChr[falPos] <- "FalsePositives"
#     realChr[realChr == "FALSE"] <- "Unknown"
#     return(as.factor(realChr))
# }
#
# # Sample read Ends with c(P (fail),Q(success)) probabilities of sampling each read
# sampleREnds <- function(rEnds, prob){
#     sampRE <- NULL
#     for(j in 1:length(rEnds)){
#         sampRE[j] <- sum(sample(0:1, rEnds[j],T, prob))
#     }
#     return(sampRE)
# }
#
# # Sample geneData to desired median depth
# sampDesMed <- function(geneData, desiredMedian){
#     suc <- 1 / (median(geneData$frags)/desiredMedian)
#     fai <- 1 - suc
#     sFrags <- sampleREnds(geneData$frags, prob = c(fai,suc))
#     sStops <- sampleREnds(geneData$stops, prob = c(fai,suc))
#     sReadt <- as.integer(round(geneData$readt*suc))
#     return(list(readt = sReadt, stops = sStops, frags= sFrags))
# }
#
# # Rank of sums of ranks in dataframe X columns (...yeah... bit complicated...)
# rankRanks <- function(x) {
#     rRanks <- apply(x, 2, rank, ties.method = "average") %>%
#         rowSums() %>% rank(ties.method = "average")
#     return(nrow(x) - rRanks + 1)
# }
#
# # Filter Results based on MMC and known rRNA sites
# filtResMMC <- function(res, keep = F, doPlot = T){
#     rRNAranks <- (max(res$rank) + 1 ) - res$rank[res$status == "Known" | res$status == "FalsePositives"]
#     which2Ome <- which(res$status[res$status == "Known" | res$status == "FalsePositives"] == "Known")
#     mcc <- matthewCorrCoeff(rRNAranks, which2Ome)
#     lastTrue <- (max(res$rank) +1 ) - min(rRNAranks[which2Ome])
#     thr <- rRNAranks[which.max(mcc)]
#     thr <- (max(res$rank) +1 ) - thr
#     if(doPlot == T){
#         plot(res$rank,
#              col = res$status,
#              main= "Ranking -- Unknown + True/False rRNA 2'-O-meth. sites",
#              ylab = "Rank in results")
#         abline(c(thr,0), col = "gray50")
#     }
#     if(keep == F){
#         # fRes <- res[res$rank <= thr & res$status == "Unknown",]
#         fRes <- res[res$rank <= thr,]
#     } else if (keep == T){
#         fRes <- res
#         fRes$passThr <- "Fail"
#         fRes$passThr[res$rank <= thr] <- "Pass"
#         fRes$passThr <- as.factor(fRes$passThr)
#         # fRes <- fRes[fRes$rank <= lastTrue & fRes$status == "Unknown",]
#         fRes <- fRes[fRes$rank <= lastTrue,]
#     }
#     fRes$rank <- rank(fRes$rank, ties.method = "random")
#     fRes <- fRes[order(fRes$rank),]
#     return(fRes)
# }
#
# # Generate Random Sites
# randResults <- function(desDepth, randCycles = 1){
#     samRES <- NULL
#     for (j in 1:randCycles){
#         for(i in 1:length(desDepth)){
#             sDATA_h <- list(a=sampDesMed(hDATA$DATA$`RDN18-1`, desDepth[i]),
#                             b=sampDesMed(hDATA$DATA$`RDN25-1`, desDepth[i]))
#             sDATA_l <- list(a=sampDesMed(lDATA$DATA$`RDN18-1`, desDepth[i]),
#                             b=sampDesMed(lDATA$DATA$`RDN25-1`, desDepth[i]))
#             names(sDATA_l) <- c("RDN18-1", "RDN25-1"); names(sDATA_h) <- c("RDN18-1", "RDN25-1")
#             sDATA_l <- list(DATA = sDATA_l)
#             sDATA_h <- list(DATA = sDATA_h)
#             sRes <- resDF_Thresh(lDATA = sDATA_l,
#                                  hDATA = sDATA_h,
#                                  stopThresh = 0.01,
#                                  fragThresh = 0.3,
#                                  fCHThresh = 1,
#                                  minMedF = 10,
#                                  evalGenes = tail(g, 2))
#             sRes$medFR <- sapply(1:nrow(sRes), function(x) median(as.numeric(unlist(strsplit(sRes$fragEnvir[x], ",")))[-7]))
#             sRes$status <- knownStatus(sRes)
#             sRes$log10pVal <- -log10(sRes$p.Value)
#             samRES <- rbind(samRES, sRes)
#         }
#     }
#     return(samRES)
# }
#
# # Calculate Threshold Parameters
# calcThr <- function(res){
#     rRNAranks <- (max(res$rank) + 1 ) - res$rank[res$status == "Known" | res$status == "FalsePositives"]
#     which2Ome <- which(res$status[res$status == "Known" | res$status == "FalsePositives"] == "Known")
#     mcc <- matthewCorrCoeff(rRNAranks, which2Ome)
#     lastTrue <- (max(res$rank) +1 ) - min(rRNAranks[which2Ome])
#     thr <- rRNAranks[which.max(mcc)]
#     thr <- (max(res$rank) +1 ) - thr
#     modPar <- c(max(mcc), thr)
#     names(modPar) <- c("maxMCC", "rankThr")
#     return(modPar)
# }
#
# # Crate BED file with range of positions for sequence retrieval with BEDTOOLS
# bedRange <- function(bedDF, d){
#     data.frame(bedDF[,1], as.integer(bedDF[,3]) - d -1,
#                as.integer(bedDF[,3]) + d, bedDF[,5],
#                bedDF[,4], bedDF[,6], stringsAsFactors = F)
# }
#
# # Function to perform t tests on groups G1 and G2 by score with an alternative hypothesis
# tTestScores <- function(resList, expG1, expG2, score, tTestAlternative){
#     expG1 <- as.integer(expG1)
#     expG2 <- as.integer(expG2)
#     g1 <- sapply(expG1, function(x) resList[[x]][score])
#     g1 <- matrix(unlist(g1), ncol = length(expG1), byrow = F)
#     g2 <- sapply(expG2, function(x) resList[[x]][score])
#     g2 <- matrix(unlist(g2), ncol = length(expG2), byrow = F)
#     tt_score <- NULL
#     for(i in 1:nrow(g1)){
#         sc1 <- g1[i,]; sc2 <- g2[i,]; tRes <- NA
#         if(sum(!is.na(sc1)) >= 2 & sum(!is.na(sc2)) >= 2 ){
#             tRes <- try(t.test(sc1, sc2, alternative = tTestAlternative)$p.value, silent = T)
#         }
#         tt_score[i] <- tRes
#     }
#     return(tt_score)
# }
#
# # Processing readEnds from rEndsFileName and gene annotation
# rEndsProcess <- function(rEndsFileName, geneAnnot){
#     rownames(geneAnnot) <- geneAnnot[,4]
#     rEndsFile <- read.delim(gzfile(rEndsFileName), header = T, stringsAsFactors = F)
#     covFpos <- rEndsFile[,c("X0", "X1", "Fpos")]
#     covFneg <- rEndsFile[,c("X0", "X1", "Fneg")]
#     r5pos <- rEndsFile[,c("X0", "X1", "r5pos")]
#     r5neg <- rEndsFile[,c("X0", "X1", "r5neg")]
#     r3pos <- rEndsFile[,c("X0", "X1", "r3pos")]
#     r3neg <- rEndsFile[,c("X0", "X1", "r3neg")]
#     GENES <- geneAnnot[,4]
#     chrs <- unique(r5pos[,1])
#     coors <- sapply(chrs, function(x) which(r5pos[,1] == x) %>% min()) - 1
#     DATA <- lapply(GENES, function(gene){
#         geneCoors <- c(geneAnnot[gene, 2], geneAnnot[gene, 3]) + coors[geneAnnot[gene, 1]]
#         if (geneAnnot[gene, 6] == "+"){
#             sta <- min(geneCoors)
#             end <- max(geneCoors)
#             readt <- covFpos[sta:end, 3] #Read-through
#             stops <- r5pos[sta:end, 3] #Reverse-Stops
#             frags <- r3pos[sta:end, 3] #Fragmentation
#         } else if (geneAnnot[gene, 6] == "-"){
#             sta <- max(geneCoors)
#             end <- min(geneCoors)
#             readt <- covFneg[sta:end, 3]
#             stops <- r5neg[sta:end, 3]
#             frags <- r3neg[sta:end, 3]
#         }
#         OUT <- list(readt, stops, frags)
#         names(OUT) <- c("readt", "stops", "frags")
#         return(OUT)
#     }
#     )
#     names(DATA) <- GENES
#     medRT <- sapply(DATA, function(x) {
#         res <- median(x[[1]])
#         return(res)
#     }
#     )
#     medST <- sapply(DATA, function(x) {
#         res <- median(x[[2]])
#         return(res)
#     }
#     )
#     medFR <- sapply(DATA, function(x) {
#         res <- median(x[[3]])
#         return(res)
#     }
#     )
#     ALLDATA <- list(DATA = DATA, medRT = medRT, medST = medST, medFR = medFR)
#     return(ALLDATA)
# }
#
# # Read a BED-6 file and transform base-0 coordinates to base-1
# readBED <- function(fileName, idAsRowNames = T){
#     tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
#     if(ncol(tmp) == 6){
#         tmp[,2] <- tmp[,2] + 1
#         colnames(tmp) <- bed6_header
#     }else if(ncol(tmp) == 12){
#         tmp[,2] <- tmp[,2] + 1
#         tmp[,7] <- tmp[,2] + 1
#         colnames(tmp) <- bed12_header
#     }else{stop("Number of columns in file is neither 6 or 12")}
#     if(idAsRowNames){rownames(tmp) <- tmp$name}
#     return(tmp)
# }
#
# # Function write bed file (Input = base1, output = base-0)
# writeBED <- function(bed, filename){
#     if(ncol(bed) == 6){
#         bed[,2] <- bed[,2] - 1
#     }else if(ncol(bed) == 12){
#         bed[,2] <- bed[,2] - 1
#         bed[,7] <- bed[,7] - 1
#     }else{
#         stop("BED is neither BED-6 or BED-12 column length")
#     }
#     write.table(bed, filename, sep= "\t", col.names = F, row.names = F, quote = F)
# }
#
#
# # Crate BED file with range of positions for sequence retrieval with BEDTOOLS
# bedRange <- function(bedDF, d){
#     cbind(as.character(bedDF[,1]), as.integer(bedDF[,3]) - d -1,
#           as.integer(bedDF[,3]) + d,
#           as.character(bedDF[,5]),
#           as.character(bedDF[,4]),
#           as.character(bedDF[,6]))
# }
#
# # Some logical functions
# great <- function(x, Than, orEqual = F, value = F){
#     if(value){
#         if(orEqual){x[x >= Than]}else if(!orEqual){x[x > Than]}
#     }else if(!value){
#         if(orEqual){x >= Than}else if(!orEqual){x > Than}
#     }
# }
#
# less <- function(x, Than, orEqual = F, value = F){
#     if(value){
#         if(orEqual){x[x <= Than]}else if(!orEqual){x[x < Than]}
#     }else if(!value){
#         if(orEqual){x <= Than}else if(!orEqual){x < Than}
#     }
# }
#
#
# equal <- function(x,y, value = F){
#     if(value){x[x == y]}else if(!value){x == y}
# }
#
# # Function to go up along a vector until finding lower values, with a threshold of fails going up
# stepWiseUp <- function(x, failThr, start = 2){
#     best <- start
#     fails <- 0
#     for(i in start:length(x)){
#         if(fails > failThr){break()}
#         if(x[i] > x[i-1]){
#             best <- i
#             fails <- 0
#         }else{fails <- fails + 1}
#     }
#     return(best)
# }
#
# # Genomic to transcriptomic coordinates based on gene annotation (Bed Format - base1)
# genCoor_2_transCorr <- function(bedLikeRow, geneAnnot){
#     tmpGA <- geneAnnot[geneAnnot$chr == as.character(bedLikeRow[1,1]) &
#                            geneAnnot$strand == bedLikeRow[1,6] &
#                            geneAnnot$start <= as.numeric(bedLikeRow[1,3]) &
#                            geneAnnot$end >= as.numeric(bedLikeRow[1,3]),]
#     lGenes <- lapply(tmpGA$name, function(x) exonBlockGen(x, geneAnnot))
#     names(lGenes) <- tmpGA$name
#     tmpGA$relCoor <- sapply(lGenes, function(x) match(as.numeric(bedLikeRow[1,3]), x))
#     return(na.omit(tmpGA[,c("name", "relCoor")]))
# }
#
# genCoor_2_transCorr_V2 <- function(bedLikeRow, geneAnnot){
#     tmpGA <- geneAnnot[geneAnnot$chr == as.character(bedLikeRow[1]) &
#                            geneAnnot$strand == bedLikeRow[6] &
#                            geneAnnot$start <= as.numeric(bedLikeRow[3]) &
#                            geneAnnot$end >= as.numeric(bedLikeRow[3]),]
#     lGenes <- lapply(tmpGA$name, function(x) exonBlockGen(x, geneAnnot))
#     names(lGenes) <- tmpGA$name
#     tmpGA$relCoor <- sapply(lGenes, function(x) match(as.numeric(bedLikeRow[3]), x))
#     return(na.omit(tmpGA[,c("name", "relCoor")]))
# }
#
# # Transcriptomic to genomic coordinates - parallelized
# transCoor_2_genCoor <- function(cl, relCoors, geneAnnot){
#     res <- parSapply(cl, 1:nrow(relCoors), function(i){
#         id <- paste(relCoors[i, ], collapse = ".")
#         gene <- relCoors[i,1]
#         relPos <- relCoors[i,2] %>% as.numeric()
#         chr <- geneAnnot[gene,]$chr
#         strand <- geneAnnot[gene,]$strand
#         if(strand == "+"){
#             genCoor <- (geneAnnot[gene,]$start:geneAnnot[gene,]$end)[relPos]
#         }else if(strand == "-"){
#             genCoor <- (geneAnnot[gene,]$end:geneAnnot[gene,]$start)[relPos]
#         }
#         return(c(chr, genCoor, genCoor, id, 0, strand))
#     }) %>% t()
#     return(res)
# }
#
# # Function to change the range of numeric values from 0 to 1 - Feature Scaling
# range0_1 <- function(x){
#     z <- (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
#     return(z)
# }
#
# # Standardized Mean Difference (2 groups - unbiased)
# smd <- function(g1, g2){
#     meanDiff <- mean(g1) - mean(g2)
#     SDwithin <- sqrt(((length(g1)-1)*sd(g1)^2 + (length(g2)-1)*sd(g2)^2)/(length(g1)+length(g2)-2))
#     smd <- meanDiff/SDwithin
#     return(smd)
# }
#
# # Function for plotting a 2 axis line plot with baseR
# # Based on https://stackoverflow.com/questions/6142944/how-can-i-plot-with-2-different-y-axes#6143251
# secAxisPlot <- function(xAxis, yAxis, yAxis2, xlab = "", ylab = "", y2lab = "", main = "", secCol= 2){
#     ## add extra space to right margin of plot within frame
#     par(mar=c(5, 5, 4, 6) + 0.1)
#     ## Plot first set of data and draw its axis
#     plot(xAxis, yAxis, pch = 16, axes = FALSE, ylim = range(yAxis), xlab= "", ylab = "",
#          type = "b", col = "black", main = main)
#     axis(2, ylim = range(yAxis), col = "black", las=1)  ## las=1 makes horizontal labels
#     mtext(ylab, side = 2, line = 4)
#     box()
#     ## Allow a second plot on the same graph
#     par(new = TRUE)
#     ## Plot the second plot and put axis scale on right
#     plot(xAxis, yAxis2, pch=15, xlab="", ylab = "", ylim = range(yAxis2),
#          axes= FALSE, type="b", col= secCol)
#     ## a little farther out (line=4) to make room for labels
#     mtext(y2lab, side = 4, col = secCol, line=4)
#     axis(4, ylim = range(yAxis2), col = secCol, col.axis= secCol, las=1)
#     ## Draw the xAxis axis
#     axis(1, pretty(range(xAxis), 10))
#     mtext(xlab, side = 1, col = "black", line = 2.5)
# }
#
# #GC content
# GC_content <- function(sequence){
#     tmp <- strsplit(sequence, "") %>% unlist()
#     (tmp == "G" | tmp == "C") %>% sum() / length(tmp)
# }
#
# #Mode function by Nick https://stackoverflow.com/users/5222/nick
# Mode <- function(x) {
#     ux <- unique(x)
#     ux[which.max(tabulate(match(x, ux)))]
# }
#
# # Nucleotide Binary Matrix
# nucBinMat <- function(tmp){
#     binMat <- matrix(0, nrow = length(tmp), ncol = 5)
#     colnames(binMat) <- c("A", "C", "G", "T", "N")
#     for(i in 1:length(tmp)){
#         binMat[i, tmp[i]] <- 1
#     }
#     return(binMat[,1:4])
# }
#
# # t.Test and FCH contrast between two groups (matrix: columns = sample, rows = 'genes')
# tTest_lFCH <- function(g1Dat, g2Dat, na.rm = T){
#     # NAs filter (min 2 samples by group)
#     keep <- rowSums(!is.na(g1Dat)) >= 2 & rowSums(!is.na(g2Dat)) >=2 ; table(keep)
#     g1Dat <- g1Dat[keep,]
#     g2Dat <- g2Dat[keep,]
#     g1_mu <- rowMeans(g1Dat, na.rm = na.rm)
#     g2_mu <- rowMeans(g2Dat, na.rm = na.rm)
#     pV <- sapply(1:nrow(g1Dat), function(i) {
#         tTest <- try(t.test(g1Dat[i,], g2Dat[i,]), silent = T)
#         if(class(tTest) == "try-error"){
#             return(NA)
#         }else return(tTest$p.value)
#     })
#     logFCH <- log2(g1_mu/g2_mu)
#     meanDiff <- g1_mu - g2_mu
#     tmp_RES <- data.frame(g1Mean = g1_mu, g2Mean= g2_mu, logFCH = logFCH, pVal = pV, meanDiff = meanDiff)
#     return(tmp_RES)
# }
#
# #Merge two data frames by rownames and renames rows
# mergebyRowNames <- function(x, y, keepAll = F){
#     tmp <- merge(x, y, by = "row.names", all = keepAll)
#     rownames(tmp) <- tmp[,1]
#     tmp <- tmp[,-1]
#     return(tmp)
# }
#
# # Merge by rownames using dplyr
# fullJoinByRowNames <- function(x, y){
#     x$rowNames <- rownames(x)
#     y$rowNames <- rownames(y)
#     tmpO <- dplyr::full_join(x, y, by = "rowNames")
#     rownames(tmpO) <- tmpO$rowNames
#     return(tmpO[, !colnames(tmpO) == "rowNames"])
# }
#
# leftJoinByRowNames <- function(x, y, keepAll = F){
#     x$rowNames <- rownames(x)
#     y$rowNames <- rownames(y)
#     tmpO <- dplyr::left_join(x, y, by = "rowNames")
#     rownames(tmpO) <- tmpO$rowNames
#     return(tmpO[,!colnames(tmpO) == "rowNames"])
# }
#
#
# # Read a BEDPE file and transform base-0 coordinates to base-1
# readBEDPE <- function(fileName, idAsRowNames = T){
#     tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
#     colnames(tmp) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2")
#     tmp$start1 <- tmp$start1 + 1
#     tmp$start2 <- tmp$start2 + 1
#     if(idAsRowNames){rownames(tmp) <- tmp$name}
#     return(tmp)
# }
#
# # Read BedPE file for getting read counts (strand specific)
# readBEDPE3 <- function(fileName){
#     tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
#     colnames(tmp) <- c("chr", "start_r1", "end_r1", "start_r2", "end_r2", "strand")
#     tmp$start_r1 <- tmp$start_r1 + 1
#     tmp$start_r2 <- tmp$start_r2 + 1
#     return(tmp)
# }
#
# # Read a custom BEDPE file and transform base-0 coordinates to base-1
# #readBEDPE2 <- function(fileName){
# #    tmp <- read.delim(fileName, header = F, stringsAsFactors = F)
# #    colnames(tmp) <- c("chr1", "start1", "end2", "times")
# #    tmp$start1 <- tmp$start1 + 1
# #    return(tmp)
# #}
#
# # Generate exon coordinates block
# exonBlockGen <- function(iGene, geneAnnot){
#     iStart <- geneAnnot[iGene,]$start
#     iEnd <- geneAnnot[iGene,]$end
#     iStrand <- geneAnnot[iGene,]$strand
#     tmpA <- geneAnnot[iGene,]$blockStarts %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
#     tmpB <- geneAnnot[iGene,]$blockSizes %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
#     iBlocks <- sapply(1:length(tmpA), function(i){
#         c(tmpA[i],(tmpA[i] + tmpB[i] -1))
#     }) %>% t
#     iBlocks <- iStart + iBlocks
#     if(iEnd %in% as.vector(iBlocks)){
#         if(iStrand == "+"){
#             sapply(1:nrow(iBlocks), function(k){
#                 iBlocks[k,1]:iBlocks[k,2]
#             }) %>% unlist %>% as.numeric
#         }else if(iStrand == "-"){
#             sapply(1:nrow(iBlocks), function(k){
#                 iBlocks[k,1]:iBlocks[k,2]
#             }) %>% unlist %>% rev %>% as.numeric
#         }
#     }else{stop(paste("Malformed exon structure at gene", iGene))}
# }
#
# # Function paired end BEDfile to transcriptome
# pairEndReadsToGenes <- function(parCluster, bedPEGzFile, gAnnot){
#     inputFile <- readBEDPE3(gzfile(bedPEGzFile))
#     commChrom <- intersect(unique(gAnnot$chr), unique(inputFile$chr))
#     genPReads <- parLapply(parCluster, commChrom, function(iChr){
#         tmpIF <- subset(inputFile, chr == iChr)
#         tmpGA <- subset(gAnnot, chr == iChr)
#         tmpO <- lapply(tmpGA$name, function(iGene){
#             exBlock <- exonBlockGen(iGene, geneAnnot)
#             if(geneAnnot[iGene, "strand"] == "+"){
#                 tmpDF <- subset(tmpIF, strand == "+" & start_r1 %in% exBlock & end_r2 %in% exBlock)
#                 tmpDF$start <- match(tmpDF$start_r1, exBlock)
#                 tmpDF$end <- match(tmpDF$end_r2, exBlock)
#                 tmpDF <- tmpDF[,c("chr", "start", "end")]
#             }else if(geneAnnot[iGene, "strand"] == "-"){
#                 tmpDF <- subset(tmpIF, strand == "-" & start_r2 %in% exBlock & end_r1 %in% exBlock)
#                 tmpDF$start <- match(tmpDF$end_r1, exBlock)
#                 tmpDF$end <- match(tmpDF$start_r2, exBlock)
#                 tmpDF <- tmpDF[,c("chr", "start", "end")]
#             }
#             tmpDF <- ddply(tmpDF,.(chr, start, end), nrow)
#             if(ncol(tmpDF) == 4){
#                 colnames(tmpDF)[4] <- "times"
#                 tmpDF$len <- tmpDF$end - tmpDF$start
#             }
#             rownames(tmpDF) <- NULL
#             return(tmpDF)
#         })
#         names(tmpO) <- tmpGA$name
#         return(tmpO)
#     })
#     # Flatten list
#     tmp <- do.call(c, genPReads)
#     genPReads <- tmp
#     return(genPReads)
# }
#
# # RowMeans by column groups
# rowMeansColG <- function(DF, colGroups, na.rm = T){
#     cG <- unique(colGroups)
#     out <- sapply(cG, function(x){
#         rowMeans(DF[,colGroups == x], na.rm = na.rm)
#     }) %>% as.data.frame()
#     colnames(out) <- cG
#     rownames(out) <- rownames(DF)
#     return(out)
# }
#
# # Extract countData from DATA object
# getCountData <- function(cluster, dat, GENES, geneAnnot, minInserts, maxInsLen){
#     rawData <- parLapply(cluster, GENES, function(iGene) {
#         iStrand <- geneAnnot[iGene, "strand"]
#         tmpD <- dat[[iGene]]
#         tmpD <- subset(tmpD, tmpD$len <= maxInsLen)
#         if(nrow(tmpD) == 0){
#             tmpRE <- list(rE5 = NA, rE3 = NA, cov = NA)
#         }else if(sum(tmpD$times) < minInserts){
#             tmpRE <- list(rE5 = NA, rE3 = NA, cov = NA)
#         }else{
#             tmpRE_5 <- tmpRE_3 <- tmpRE_c <- rep(0L, length(exonBlockGen(iGene, geneAnnot)))
#             for(k in 1:nrow(tmpD)){
#                 tmpRE_5[tmpD[k,]$start] <- tmpRE_5[tmpD[k,]$start] + tmpD[k,]$times
#                 tmpRE_3[tmpD[k,]$end] <- tmpRE_3[tmpD[k,]$end] + tmpD[k,]$times
#                 tmpRE_c[tmpD[k,]$start:tmpD[k,]$end] <- tmpRE_c[tmpD[k,]$start:tmpD[k,]$end] + tmpD[k,]$times
#             }
#             tmpRE <- list(rE5 = tmpRE_5, rE3 = tmpRE_3, cov = tmpRE_c)
#         }
#         return(tmpRE)
#     })
#     return(rawData)
# }
#
#
# # Load/install CRAN packages
# installLoad_CRAN <- function(package){
#     if (!require(package, character.only = T)) {
#         install.packages(package, dependencies = TRUE,
#                          repos = "http://cran.us.r-project.org")
#         library(package, character.only = T, quietly = T)
#     }
# }
#
# # Load/install Bioconductor packages
# installLoad_BioC <- function(package){
#     if (!require(package, character.only = T)) {
#         if(version$minor %>% as.numeric() >= 3.5){
#             if (!requireNamespace("BiocManager", quietly = TRUE)){
#                 installLoad_CRAN("BiocManager")
#             }
#             BiocManager::install(package)
#             library(package, character.only = T, quietly = T)
#         }
#         if(version$minor %>% as.numeric() < 3.5){
#             source("https://bioconductor.org/biocLite.R")
#             biocLite(package, ask = F)
#             library(package, character.only = T, quietly = T)
#         }
#     }
# }
#
# # Number of inserts in sample
# insertsNum <- function(pairEndReads){
#     sapply(pairEndReads, function(x){
#         sum(x$times)
#     })
# }
#
# # Extract rawData from DATA object
# getRawData <- function(dat, GENES, geneAnnot, motifOme, maxL){
#     rawData <- lapply(GENES, function(iGene) {
#         iStrand <- geneAnnot[iGene, "strand"]
#         tmpD <- dat[[iGene]]
#         iCoors <-  motifOme[[iGene]]
#         tmpD$len <- abs(tmpD$start - tmpD$end)
#         tmpD <- tmpD[tmpD[,"len"] <= maxL,]
#         if(nrow(tmpD) == 0){
#             tmpRE <- cbind(NA, NA, NA, NA, NA) %>% as.matrix()
#         }else if(sum(tmpD$times) < minInserts){
#             tmpRE <- cbind(NA, NA, NA, NA, NA) %>% as.matrix()
#         }else{
#             tmpRE_5 <- tmpRE_3 <- tmpRE_c <- rep(0L, geneAnnot[iGene,]$seqs %>% nchar)
#             if(iStrand == "-"){
#                 for(k in 1:nrow(tmpD)){
#                     tmpRE_5[tmpD[k,]$end] <- tmpRE_5[tmpD[k,]$end] + tmpD[k,]$times
#                     tmpRE_3[tmpD[k,]$start] <- tmpRE_3[tmpD[k,]$start] + tmpD[k,]$times
#                     tmpRE_c[tmpD[k,]$end:tmpD[k,]$start] <- tmpRE_c[tmpD[k,]$end:tmpD[k,]$start] + tmpD[k,]$times
#                 }
#             }else if(iStrand == "+"){
#                 for(k in 1:nrow(tmpD)){
#                     tmpRE_5[tmpD[k,]$start] <- tmpRE_5[tmpD[k,]$start] + tmpD[k,]$times
#                     tmpRE_3[tmpD[k,]$end] <- tmpRE_3[tmpD[k,]$end] + tmpD[k,]$times
#                     tmpRE_c[tmpD[k,]$start:tmpD[k,]$end] <- tmpRE_c[tmpD[k,]$start:tmpD[k,]$end] + tmpD[k,]$times
#                 }
#             }else{stop("No valid strand specified in gene annotation")
#             }
#             rowN <- paste(iGene, iCoors, sep = "_")
#             tmpRE <- cbind(rowN, tmpRE_5[iCoors], tmpRE_3[iCoors-1], tmpRE_c[iCoors], tmpRE_c[iCoors -1]) %>% as.matrix()
#             colnames(tmpRE) <- c("coor", "rE5","rE3", "rT5", "rT3")
#         }
#         return(tmpRE)
#     })
#     rawData <- do.call("rbind", rawData) %>% na.omit() %>% asDataFrame()
#     return(rawData)
# }
#
# # Get all the relative coordinates for all genes using a BED-like matrix row
# getRelPos <- function(BEDrow){
#     tmp <- BEDrow
#     iChr <- tmp$chr
#     iStart <- tmp$start
#     iStrand <- tmp$strand
#     tmpA <- subset(geneAnnot, chr == iChr & strand == iStrand & start <= iStart & end >= iStart)
#     tmpCoor <- sapply(tmpA$name, function(x){
#         tmp <- which(exonBlocks[[x]] == iStart)
#         if(!as.logical(length(tmp))){
#             tmp <- NA
#         }
#         return(tmp)
#     }) %>% unlist
#     paste0(tmpA$name,"_", tmpCoor)[!is.na(tmpCoor)]
# }
#
# # From relative to genomic coordinate
# relToGenCoors <- function(relCoor, genAnnot = geneAnnot){
#     splt <- unlist(strsplit(relCoor, split = "_"))
#     gene <- splt[1] %>% as.character()
#     genCoor <- exonBlockGen(gene, genAnnot)[as.numeric(splt[2])]
#     return(c(genAnnot[gene,]$chr, genCoor, genCoor, relCoor, 0, genAnnot[gene,]$strand))
# }
#
# # Extract window of sequences from gene sequence vector
# extracSeq <- function(relCoor, window = 10, fullSeqs = txOmeSEQS){
#     tmp <- unlist(strsplit(relCoor, split = "_"))
#     iGene <- tmp[1]
#     iCoor <- as.numeric(tmp[2])
#     substr(fullSeqs[iGene], start = iCoor - window, stop =  iCoor + window)
# }
#
# # raw cleavage efficiencies in motif (motifOme object required)
# rawClvEff_motif <- function(countDataFile, RTthr = 15, motifOme = motifOme){
#     load(countDataFile)
#     tmpSel <- sapply(countData, function(x) length(x$cov)) > 1
#     countData <- countData[tmpSel]
#     genesI <- intersect(names(countData), names(motifOme))
#     rawClvEff <- lapply(genesI, function(geneI){
#         tmp <- data.frame(coorNames = paste0(geneI,"_", motifOme[[geneI]]),
#                           rE5 = (countData[[geneI]]$rE5)[motifOme[[geneI]]],
#                           rT5 = (countData[[geneI]]$cov)[motifOme[[geneI]]],
#                           rE3 = (countData[[geneI]]$rE3)[motifOme[[geneI]] -1],
#                           rT3 = (countData[[geneI]]$cov)[motifOme[[geneI]] -1], stringsAsFactors = F)
#         tmp$clvEff_5 <- tmp$rE5 / tmp$rT5
#         tmp$clvEff_3 <- tmp$rE3 / tmp$rT3
#         tmp$clvEff_5[tmp$rT5 < RTthr] <- NA
#         tmp$clvEff_3[tmp$rT3 < RTthr] <- NA
#         out <- tmp[(is.na(tmp[,c("clvEff_5", "clvEff_3")]) %>% rowSums()) < 2,]
#         return(out)
#     })
#     rawClvEff <- do.call(rbind, rawClvEff)
#     rawClvEff <- merge(allSites, rawClvEff)
#     return(rawClvEff)
# }
#
# # Optimal distance analysis
# optiMotifDist <- function(rawCE, distLimit){
#     TMP_both <- sapply(1:distLimit, function(thr){
#         TMPdata <- cbind(rawCE$clvEff_5[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= thr],
#                          rawCE$clvEff_3[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= thr]) %>% na.omit
#         c(cor(TMPdata[,1], TMPdata[,2]), nrow(TMPdata))
#     })
#
#     minBoD <- stepWiseUp(TMP_both[1,], failThr = 1, start = 5)
#     TMP_up <- sapply(1:distLimit, function(thr){
#         TMP_up <- cbind(rawCE$clvEff_5[rawCE$doS_motifDist >= minBoD & rawCE$upS_motifDist >= thr],
#                         rawCE$clvEff_3[rawCE$doS_motifDist >= minBoD & rawCE$upS_motifDist >= thr]) %>% na.omit
#         c(cor(TMP_up[,1], TMP_up[,2]), nrow(TMP_up))
#     })
#     TMP_do <- sapply(1:distLimit, function(thr){
#         TMP_do <- cbind(rawCE$clvEff_5[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= minBoD],
#                         rawCE$clvEff_3[rawCE$doS_motifDist >= thr & rawCE$upS_motifDist >= minBoD]) %>% na.omit
#         c(cor(TMP_do[,1], TMP_do[,2]), nrow(TMP_do))
#     })
#     minUpD <- stepWiseUp(TMP_up[1,], failThr = 1, start = 10)
#     minDoD <- stepWiseUp(TMP_do[1,], failThr = 1, start = 10)
#     tmp <- c(minBoD, minUpD, minDoD); names(tmp) <- c("minBoD", "minUpD", "minDoD")
#     return(list(boD = TMP_both, upD = TMP_up, doD = TMP_do, minD = tmp))
# }
#
# # Plotting optimal distance analyisis
# plot_OptimMotifDist <- function(opMoDi, smpName){
#     palette(brewer.pal(9, "Set1")) #RColorBrewer
#     par(mfrow=c(1,3))
#     distLimit <- ncol(opMoDi$boD)
#     secAxisPlot(1:distLimit, opMoDi$boD[2,], opMoDi$boD[1,],
#                 xlab = "Distance to ACA",
#                 ylab = "n ACA sites eval",
#                 y2lab = "Pearson Correlation",
#                 main = "P. Corr ~ ACA dist.-BothS")
#     abline(v = opMoDi$minD["minBoD"], col = "red")
#     abline(h = opMoDi$boD[1, opMoDi$minD["minBoD"]], col = "red")
#
#     secAxisPlot(1:distLimit, opMoDi$upD[2,], opMoDi$upD[1,],
#                 xlab = "Distance to ACA",
#                 ylab = "n ACA sites eval",
#                 y2lab = "Pearson Correlation",
#                 main = "P. Corr ~ ACA dist.-UpS")
#     title(sub = smpName)
#     abline(v = opMoDi$minD["minUpD"], col = "red")
#     abline(h = opMoDi$upD[1, opMoDi$minD["minUpD"]], col = "red")
#
#     secAxisPlot(1:distLimit, opMoDi$doD[2,], opMoDi$doD[1,],
#                 xlab = "Distance to ACA",
#                 ylab = "n ACA sites eval",
#                 y2lab = "Pearson Correlation",
#                 main = "P. Corr ~ ACA dist.-DownS")
#     abline(v = opMoDi$minD["minDoD"], col = "red")
#     abline(h = opMoDi$doD[1, opMoDi$minD["minDoD"]], col = "red")
#     par(mfrow=c(1,1))
# }
#
# # Calculating average cleavage efficiency using rawClvEff_motif output and
# # downStream and upStream thresholds
# avgCleavEff <- function(rawClvEff, doS_thr, upS_thr){
#     tmp <- rawClvEff
#     tmp$clvEff_5[tmp$doS_motifDist < doS_thr] <- NA
#     tmp$clvEff_3[tmp$upS_motifDist < upS_thr] <- NA
#     tmp$avgClvEff <- rowMeans(tmp[,c("clvEff_5", "clvEff_3")], na.rm = T)
#     tmp$avgClvEff[is.na(tmp$avgClvEff)] <- NA
#     out <- rawClvEff
#     out$avgClvEff <- tmp$avgClvEff
#     return(out)
# }
#
# # Cleavage motif frequency (based on first 'N' nucleotides)
# clvNucFreq <- function(cluster, countFile, geneSeqs, nNucs){
#     load(countFile)
#     triNucs <- permutations(n = 4,r = nNucs, v = c("A", "T", "C", "G"),
#                             repeats.allowed = T) %>% apply(MARGIN = 1, function(x){paste(x, collapse = "")})
#     iGenes <- intersect(names(countData)[sapply(countData, function(x) sum(x$rE3)>0)], names(geneSeqs))
#     nucsFreq_5 <- parSapply(cluster, iGenes, function(iG){
#         test <- sapply(which(countData[[iG]]$rE5 > 0), function(i){
#             c(substring(text = geneSeqs[iG], first = i, last = i+nNucs-1), countData[[iG]]$rE5[i])
#         }) %>% t %>% asDataFrame()
#         if(is.null(names(test))){
#             test <- data.frame(a= NA, b = NA)
#         }
#         names(test) <- c("nucs", "times")
#         sapply(triNucs, function(x){
#             sum(subset(test, nucs == x)$times)
#         })
#     }) %>% rowSums()
#     nucsFreq_3 <- parSapply(cluster, iGenes, function(iG){
#         test <- sapply(which(countData[[iG]]$rE3 > 0), function(i){
#             c(substring(text = geneSeqs[iG], first = i+1, last = i+nNucs), countData[[iG]]$rE3[i])
#         }) %>% t %>% asDataFrame()
#         if(is.null(names(test))){
#             test <- data.frame(a= NA, b = NA)
#         }
#         names(test) <- c("nucs", "times")
#         sapply(triNucs, function(x){
#             sum(subset(test, nucs == x)$times)
#         })
#     }) %>% rowSums()
#     list(nucsFreq_5= nucsFreq_5, nucsFreq_3 = nucsFreq_3)
# }
#
# # Which positions are top N
# whichTopN <- function(x, N){
#     names(x) <- 1:length(x)
#     as.numeric(names(tail(sort(x), N)))
# }
# # Which positions are bottom N
# whichBottomN <- function(x, N){
#     names(x) <- 1:length(x)
#     as.numeric(names(head(sort(x), N)))
# }
# # Omit infinite values
# inf.omit <- function(x){
#     x[!is.infinite(x)]
# }
#
# # Keep rows and column names when using normalize.quantiles()
# normalize.quantiles.keepNames <- function(x){
#     x[is.na(x)] <- NA
#     require(preprocessCore)
#     rN <- rownames(x)
#     cN <- colnames(x)
#     tmpO <- data.frame(normalize.quantiles(as.matrix(x)))
#     rownames(tmpO) <- rN; colnames(tmpO) <- cN
#     return(tmpO)
# }
#
# # Clear temporary variables
# clearTmp <- function(){rm(list = grep("tmp", ls(), value = T)); gc(verbose = T)}
#
# # Get model data
# getModeldata <- function(geneAnnot, txOmeSEQS, relCoors){
#     cl <- makeCluster(nCores, type ="FORK")
#     tmpCoors <- parSapply(cl, relCoors, function(x) relToGenCoors(x, geneAnnot)) %>% t %>% asDataFrame()
#     stopCluster(cl)
#     seqData <- data.frame(stringsAsFactors = F, tmpCoors)
#     colnames(seqData) <- c("chr", "start", "end", "id", "score", "strand")
#     class(seqData$start) <- "numeric"
#     class(seqData$end) <- "numeric"
#     seqData <- seqData[!duplicated(seqData),]
#     rownames(seqData) <- seqData$id
#     all_relPos <- sapply(rownames(seqData), function(x) strsplit(x, split = "_")) %>% unlist %>% matrix(ncol = 2, byrow = T)
#     wind <- 30
#     cl <- makeCluster(nCores, type ="FORK")
#     seqData$sequence <- parSapply(cl, 1:nrow(all_relPos), function(i){
#         iGene <- all_relPos[i, 1]
#         iPos <- all_relPos[i, 2] %>% as.numeric()
#         return(txOmeSEQS[[iGene]] %>% substring(first = iPos - wind, last = iPos + wind))
#     })
#     stopCluster(cl)
#     nick <- tempfile()
#     writeLines(seqData$sequence, paste0(nick, ".seqs"))
#     # Free energy estimation (ViennaRNA)
#     comm2 <- paste0("/apps/RH7U2/gnu/ViennaRNA/2.3.1-nogsl/bin/RNAfold -i ", nick,
#                     ".seqs --noconv --noPS -T 30 > ", nick, ".fold.out.fold")
#     system(comm2)
#     seqFold <- readLines(paste0(nick, ".fold.out.fold"))
#     tmp <- seqFold[1:length(seqFold) %% 2 == 0]
#     sep <- grep(tmp[1] %>% strsplit(split = "") %>% unlist, pattern = " ")[1]
#     cl <- makeCluster(nCores, type ="FORK")
#     seqData$fEnergy <- parSapply(cl, 1:length(tmp), function(i) {
#         substring(tmp[i], first = sep+1) %>%
#             gsub(pattern = "\\(", replacement = "") %>%
#             gsub(pattern = "\\)", replacement = "") %>%
#             as.numeric()
#     })
#     stopCluster(cl)
#     dummy <- file.remove(paste0(nick, ".fold.out.fold"))
#
#     #GC content
#     cl <- makeCluster(nCores, type ="FORK")
#     seqData$GCcont <- parSapply(cl, seqData$sequence, GC_content)
#     stopCluster(cl)
#
#     # Relative position in gene
#     cl <- makeCluster(nCores, type ="FORK")
#     seqData$relPosRat <- parSapply(cl, 1:nrow(all_relPos), function(i){
#         gene <- all_relPos[i, 1]
#         relPos <- all_relPos[i, 2] %>% as.numeric()
#         fRelPos <- relPos / length(exonBlockGen(gene, geneAnnot))
#         return(fRelPos)
#     })
#     stopCluster(cl)
#
#     # Nucleotide Matrix
#     wind <- 6
#     cl <- makeCluster(nCores, type ="FORK")
#     tmpSeqs <- parSapply(cl, 1:nrow(all_relPos), function(i){
#         iGene <- all_relPos[i, 1]
#         iPos <- all_relPos[i, 2] %>% as.numeric()
#         return(txOmeSEQS[[iGene]] %>% substring(first = iPos - wind, last = iPos + wind))
#     })
#     stopCluster(cl)
#     selR <- sapply(tmpSeqs, nchar) %>% equal((wind*2)+1)
#     tmpSeqs <- tmpSeqs[selR]
#     seqData <- seqData[selR,]
#     cl <- makeCluster(nCores, type ="FORK")
#     nucMat <- parSapply(cl, tmpSeqs, function(x) {
#         x %>% strsplit("") %>% unlist}) %>% t
#     stopCluster(cl)
#     rownames(nucMat) <- rownames(seqData)
#     seqData$seqs <- tmpSeqs
#     # Nucleotides binary matrix
#     selNucPos <- c(3:wind,(wind+4):((wind*2)+1))
#     selNucMat <- nucMat[, selNucPos]
#     cl <- makeCluster(nCores, type = "FORK")
#     nucBinMatrix <- parLapply(cl, 1:ncol(selNucMat), function(x) nucBinMat(selNucMat[,x])) %>% cbind.data.frame()
#     stopCluster(cl)
#     colnames(nucBinMatrix) <- sapply(paste("POS", selNucPos - wind - 1, sep = "_"),
#                                      function(x) {paste(x, c("A", "C", "G", "T"), sep = "_")}) %>% as.vector()
#     rownames(nucBinMatrix) <- rownames(seqData)
#
#     # Putting everythin together
#     modelDATA <- data.frame(seqData[, c("chr", "start", "end", "id", "score",
#                                         "strand", "fEnergy", "GCcont","relPosRat", "seqs")],
#                             nucBinMatrix)
#     rownames(modelDATA) <- rownames(nucBinMatrix)
#     return(modelDATA)
# }
#
# # Merge lists elements by names
# mergeListsByNames <- function(list){
#     keys <- unique(unlist(lapply(list, names)))
#     setNames(do.call(mapply, c(FUN=c, lapply(list, `[`, keys))), keys)
# }
#
# # Extract random window of consecutive elements in vector
# extractRandWindow <- function(x, p){
#     firstIndex = sample(seq(length(x) - p + 1), 1)
#     x[firstIndex:(firstIndex + p -1)]
# }
#
#
# # get Coverage of range
# getCov_InRange <- function(covObject, gRange){
#     covObject[[as.character(seqnames(gRange))]][(start(gRange)):(end(gRange))]
# }
#
# # Remove alignments that include Soft-clipping, Hard-clipping, deletions and insertions (S|D|I  )
# cleanMapping_PairEnd <- function(pairedEnd_Alignments){
#     tmpCig_f <- cigar(pairedEnd_Alignments@first)
#     tmpCig_l <- cigar(pairedEnd_Alignments@last)
#     pairedEnd_Alignments <- pairedEnd_Alignments[-union(grep(pattern = "D|I|S|H", tmpCig_f),
#                                                         grep(pattern = "D|I|S|H", tmpCig_l))]
# }
#
# # Getting sequence from genome
# getGenSeqInGenome <- function(gene, geneAnnot, genome){
#     tmpA <- geneAnnot[gene,]$blockStarts %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
#     tmpB <- geneAnnot[gene,]$blockSizes  %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
#     iBlocks <- (sapply(1:length(tmpA), function(i){
#         c(tmpA[i],(tmpA[i] + tmpB[i] -1))
#     }) %>% t) + geneAnnot[gene,]$start
#     if(geneAnnot[gene,]$end %in% as.vector(iBlocks)){
#         if(geneAnnot[gene,]$strand == "+"){
#             sapply(1:nrow(iBlocks), function(k){
#                 str_sub(genome[[geneAnnot[gene,]$chr]], iBlocks[k, 1], iBlocks[k,2])
#             }) %>% paste(collapse = "") %>% DNAString()
#         }else if(geneAnnot[gene,]$strand == "-"){
#             sapply(1:nrow(iBlocks), function(k){
#                 str_sub(genome[[geneAnnot[gene,]$chr]], iBlocks[k, 1], iBlocks[k,2])
#             }) %>% paste(collapse = "") %>% DNAString() %>% reverseComplement()
#         }
#     }else{stop(paste("Malformed exon structure at gene", gene))}
# }
#
# # From transcriptomic range to genomic coordinates in BED format
# txRangeToBED <- function(Tgene, TrangeStart, TrangeEnd, geneAnnot){
#     genCoor <- exonBlockGen(Tgene, geneAnnot)[TrangeStart:TrangeEnd]
#     if(geneAnnot[Tgene,]$strand == "-"){
#         genCoor <- rev(genCoor)
#     }
#     gStart <- genCoor[1]
#     gEnd <- genCoor[length(genCoor)]
#     blockCount <- sum(diff(genCoor) !=1) + 1
#     if(sum(diff(genCoor) < 1) > 0){stop("Next exon position cannot be negative")}
#     blockSizes <- paste(diff(which(c(2, diff(genCoor), 2) !=1)), collapse = ",")
#     blockStarts <- paste(c(0, genCoor[which(diff(genCoor) != 1) +1] - genCoor[1]), collapse = ",")
#     c(geneAnnot[Tgene,]$chr, gStart, gEnd, Tgene, 0,
#              geneAnnot[Tgene,]$strand, gStart, gEnd, "0", blockCount,
#              blockSizes, blockStarts)
# }
#
#
# # Shotcuts #####################################################################
# bed6_header <- c("chr", "start", "end", "name", "score", "strand")
# bed12_header <- c("chr", "start", "end", "name", "score", "strand",
#                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
#
#
# # GenomicAlignments  functions #############################################################
# # Merging GAlignmentPairs
# mergeGAPairs <- function(GAPairsList){
#   tmp_first <- lapply(1:length(GAPairsList), function(i){
#     GAPairsList[[i]]@first
#   }) %>% do.call(what = c)
#   tmp_last  <- lapply(1:length(GAPairsList), function(i){
#     GAPairsList[[i]]@last
#   }) %>% do.call(what = c)
#   tmp_names <- lapply(1:length(GAPairsList), function(i){
#     names(GAPairsList[[i]])
#   }) %>% do.call(what = c)
#   GAlignmentPairs(first = tmp_first, last = tmp_last, names = tmp_names)
# }
#
# # Extract Nucleotides from GAPairs object at single nucleotide position
# extract_nuc <- function(GAPairs, genomicPos){
#     ov_f <- findOverlaps(genomicPos, GAPairs@first, ignore.strand = T)@to
#     ov_l <- findOverlaps(genomicPos, GAPairs@last, ignore.strand = T)@to
#     ov_l <- setdiff(ov_l, ov_f)
#     if(length(ov_f) > 0){
#         nucs_f <- str_sub(string = sequenceLayer(mcols(GAPairs@first[ov_f])$seq,
#                                                  mcols(GAPairs@first[ov_f])$cigar),
#                           start = start(genomicPos) - start(GAPairs@first[ov_f]) + 1,
#                           end = start(genomicPos) - start(GAPairs@first[ov_f]) + width(genomicPos))
#     }else{nucs_f <- NULL}
#     if(length(ov_l) > 0){
#         nucs_l <- str_sub(string = sequenceLayer(mcols(GAPairs@last[ov_l])$seq,
#                                                  mcols(GAPairs@last[ov_l])$cigar),
#                           start = start(genomicPos) - start(GAPairs@last[ov_l]) + 1,
#                           end = start(genomicPos) - start(GAPairs@last[ov_l]) + width(genomicPos))
#     }else{nucs_l <- NULL}
#     if(as.logical(strand(genomicPos) == "-")){
#         as.character(complement(DNAStringSet(c(nucs_f, nucs_l))))
#     }else{c(nucs_f, nucs_l)}
# }
#
# # Read counts by range
# readCountsByRange <- function(BAM_file, targetRanges){
#     flag0 <- scanBamFlag(isDuplicate = FALSE, # Base parameters
#                          isNotPassingQualityControls = FALSE,
#                          isPaired = T) #Remove optical duplicates and low quality reads
#     chroms <- seqnames(targetRanges) %>% unique %>% as.character()
#     lapply(chroms, function(iChr){
#         param0 <- ScanBamParam(flag=flag0,
#                                which = targetRanges)
#         rds <- readGAlignmentPairs(BAM_file, use.names = TRUE, param = param0)
#         rds_pos <- rds[strand(rds@first) == "+"]
#         rds_neg <- rds[strand(rds@first) == "-"]
#         tmpRds_pos <- GRanges(seqnames = seqnames(rds_pos),
#                               ranges = IRanges(start = start(rds_pos@first), end = end(rds_pos@last)),
#                               strand = strand(rds_pos@first))
#         tmpRds_neg <- GRanges(seqnames = seqnames(rds_neg),
#                               ranges = IRanges(start = start(rds_neg@last), end = end(rds_neg@first)),
#                               strand = strand(rds_neg@first))
#         rdsExt <- c(tmpRds_pos, tmpRds_neg)
#         tmpCoors <- targetRanges[seqnames(targetRanges) == iChr]
#         tmpOver <- findOverlaps(rdsExt, tmpCoors)
#         readSplit <- split(tmpOver@from, names(tmpCoors)[tmpOver@to])
#         readSplit <- lapply(readSplit, unique) # Remove duplicated reads
#         sapply(readSplit, length) # Counts by
#     }) %>% unlist
# }
#
# # Counting reads by gene in gene annotation as ranges
# readCountsByGene <- function(BAM_file, geneAnnotRanges){
#     flag0 <- scanBamFlag(isDuplicate = FALSE, # Base parameter
#                          isNotPassingQualityControls=FALSE,
#                          isPaired = T) #Remove optical duplicates and low quality reads
#     chroms <- seqnames(geneAnnotRanges) %>% unique %>% as.character()
#     lapply(chroms, function(iChr){
#         param0 <- ScanBamParam(flag=flag0,
#                                which = GRanges(seqnames= iChr,
#                                                ranges=IRanges(start=1,
#                                                               end= length(genome[[iChr]])),
#                                                strand = "*"))
#         rds <- readGAlignmentPairs(BAM_file, use.names = TRUE, param = param0)
#         strand(rds@last) <- strand(invertStrand(rds@last))
#         tmpGA <- geneAnnotRanges[seqnames(geneAnnotRanges) == iChr]
#         fOver <- findOverlaps(rds@first, tmpGA, type = "within")
#         lOver <- findOverlaps(rds@last, tmpGA, type = "within")
#         readSplit_f <- split(fOver@from, tmpGA$name[fOver@to])
#         readSplit_l <- split(lOver@from, tmpGA$name[lOver@to])
#         readSplit <- mergeListsByNames(list(readSplit_f, readSplit_l))
#         readSplit <- lapply(readSplit, unique) # Remove duplicated reads
#         sapply(readSplit, length) # Gene Counts
#     }) %>% unlist
# }
#
# # Counts list to dataframe function
# countsList_ToDataFrame <- function(countsList){
#     tmpNames <- sapply(countsList, names) %>% unlist %>% unique %>% sort
#     tmpDF <- lapply(1:length(countsList), function(i){
#         tmp <- data.frame("x" = countsList[[i]][tmpNames], row.names = tmpNames)
#         tmp[is.na(tmp)] <- 0
#         tmp
#     }) %>% do.call(what = cbind)
#     colnames(tmpDF) <- names(countsList)
#     tmpDF
# }
#
# # TPM normalization
# TPMnorm <- function(tmpDF){
#   cS <- colSums(tmpDF)
#   tmp <- sapply(1:ncol(tmpDF), function(i){
#     ((tmpDF[,i]/cS[i]) * 1000000)
#   })
#   tmp
# }
