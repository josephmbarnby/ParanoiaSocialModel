# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        #####     NETWORK MODELLING   ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# to run this script you must first have run script "hisi_script_OSF"                        
                        
#network stuff
library(psychonetrics)
library(qgraph) # v1.6.5
library(bootnet)

#model selection approach
names(test.full.fits.df2)
comb.modselect <- estimateNetwork(test.full.fits.df2[c(11, 14:19)], 
                                  default = "ggmModSelect",
                                  stepwise = T
)
plot(comb.modselect, 
     layout = 'circle',
     details = T)
layout(1)
centralityPlot(list(EstimateNetwork = comb.modselect,
                    GGMmodel = comb.net,
                    Glasso = comb.glasso),
               include = "all")

#condition on the moderator

combsims.modelling$LevelLoop <- ordered(cut(combsims$GPTSTotal, 5, labels = F))

par(mfrow = c(3,2))
predictors <- list()
graphs  <- list()
names(combsims.modelling)
#for (i in 1:5){

library("mgm")

data <- test.full.fits.df2[c(11, 14:19)]
names(data)
par(mfrow=c(2,3))
par(mar=c(4,4,4,4))
hist(data[2:4])

predictor<- mgm::mgm(data,
                             type = c(rep('g', 7)),
                             level = c(rep(1, 7)),
                             scale = T,
                             lambdaSel = "EBIC",
                             ruleReg = "AND",
                            moderators = c(1))


tb <- table(scale(data$GPTSTotal))
names(tb) <- round(as.numeric(names(tb)), 2)
par(mar=c(4,4,4,4))
barplot(tb, axes=FALSE, xlab="", ylim=c(0, 3000))
axis(2, las=2, c(0, 1000, 2000, 3000))

low <- condition(predictor, values = list("1"= -0.85))
med <- condition(predictor, values = list("1"= 0))
hig <- condition(predictor, values = list("1"= 1))
vhi <- condition(predictor, values = list("1"= 2))
cln <- condition(predictor, values = list("1"= 3))
ucl <- condition(predictor, values = list("1"= 4))

l_cond <- list(low, med, hig, vhi, cln, ucl)

max_val <- max(max(l_cond[[1]]$pairwise$wadj),
               max(l_cond[[2]]$pairwise$wadj),
               max(l_cond[[3]]$pairwise$wadj),
               max(l_cond[[4]]$pairwise$wadj),
               max(l_cond[[5]]$pairwise$wadj),
               max(l_cond[[6]]$pairwise$wadj))

library(qgraph)
library(bootnet)
par(mfrow=c(2,3))
graphs <- list()

for(i in 1:6) graphs[[i]] <- qgraph(l_cond[[i]]$pairwise$wadj, layout="circle", 
                     edge.color=l_cond[[i]]$pairwise$edgecolor, 
                     labels = colnames(data), 
                     label.cex = 1.5,
                     maximum = max_val, 
                     posCol = "blue",
                     theme = "colorblind",
                     edge.labels = TRUE, edge.label.cex=4, 
                     vsize=20, esize=18, border.width = 3,
                     title = paste0("GPTS Score = ", (c("-0.85", "-0", "1", "2", "3", "4"))[i]),
                     title.cex = 1.5)


set.seed(1) 
res_obj <- resample(object = predictor, 
                    data = data, 
                    nB = 50)
plotRes(res_obj, 
        cex.label = 2, cex.bg = 5, cex.mean = 1.3, layout.gap.pw.mod = 0.05,
        axis.ticks.mod = c(-.1, -.05, 0, .05, .1), 
        labels = colnames(data), 
        layout.width.labels = .40, lwd.qtl = 4)

#sanity check

names(data)
data.a <- data
data.a$zPara <- scale(data.a$GPTSTotal)

ppcor::pcor(data, method = "spearman")$estimate

sanEst <- list()

sanEst[[1]] <- ppcor::pcor(data.a[which(data.a$zPara > -0.85 & data.a$zPara < 0),][2:7], method = "spearman")$estimate
sanEst[[2]] <- ppcor::pcor(data.a[which(data.a$zPara > 0 & data.a$zPara < 1),][2:7], method = "spearman")$estimate
sanEst[[3]] <- ppcor::pcor(data.a[which(data.a$zPara > 1 & data.a$zPara < 2),][2:7], method = "spearman")$estimate
sanEst[[4]] <- ppcor::pcor(data.a[which(data.a$zPara > 2 & data.a$zPara < 3),][2:7], method = "spearman")$estimate

ppcor::pcor(data)

par(mfrow=c(2,3))

for(i in 1:4){
  sanEst[[i]] <- as.data.frame(sanEst[[i]])
  sanEst[[i]]$Level <- i
  print(sanEst[[i]])}
do.call(rbind, sanEst)

for(i in 1:4){
  corrplot(as.matrix(sanEst[[i]][1:6]), type="upper", insig = "blank", method = "number",title = sanEst[[i]][7])}

#accuracy
comb.accuracy <- bootnet::bootnet(comb.modselect, 
                         nBoots = 1000, 
                         default = "glasso", 
                         statistics = "edge", 
                         nCores = 8)
plot(comb.accuracy, labels = T, order = "sample")
plot(comb.accuracy, labels = T, order = "sample", plot = "interval")

#stability

comb.stability <- bootnet(comb.modselect, 
                              type = "case", 
                              statistics = c("edge", "strength", "closeness", "betweenness", "expectedInfluence"),
                              nCores = 8)

plot(comb.stability, statistics = c("strength", "closeness", "betweenness", "expectedInfluence"))
corStability(comb.stability, cor = 0.7, verbose = T)

