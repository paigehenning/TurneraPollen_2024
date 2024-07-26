https://rpubs.com/natmurad/WGCNA

library(doParallel)
registerDoParallel(cores=15)
library(WGCNA) 
library(flashClust) 
library(tidyverse) 
library(magrittr)
enableWGCNAThreads(nThreads = 15) 
options(stringsAsFactors = FALSE)

 
expressiondata = read.csv("Combined.tab", header=TRUE, sep="\t")
rownames(expressiondata)= expressiondata$ID
expressiondata = expressiondata[-c(1)]
flipped <- data.frame(t(expressiondata))
flipped=fix(flipped)
expression = flipped[, colSums(flipped != 0) > 0]
rm(flipped)


sampleTree = hclust(dist(expression), method = "average")
nGenes = ncol(expression)
nSamples = nrow(expression)
traitData = read.csv("trait.txt",sep="\t",header=TRUE)
Samples = rownames(expression)
traitRows = match(Samples, traitData$Plant)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()
sampleTree2 = hclust(dist(expression), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")



powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expression, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
	 
	abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=8
adjacency = adjacency(expression, power = softPower, type = "signed") 
TOM = TOMsimilarity(adjacency) 
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = 20
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods) 
MEList = moduleEigengenes(expression, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
	 

MEDissThres = 0.2

merge = mergeCloseModules(expression, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


nGenes = ncol(expression)
nSamples = nrow(expression)

MEs0 = moduleEigengenes(expression, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



##########
#trait
###########

Peso10dias = as.data.frame(datTraits$SC)
names(Peso10dias) = "Mating type"
###################
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(expression, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(expression, Peso10dias, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Peso10dias), sep="")
names(GSPvalue) = paste("p.GS.", names(Peso10dias), sep="")

###########################################

genes = colnames(expression)
inModule = is.finite(match(moduleColors, moduleColors))
modGenes = genes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)
modTOMSignificantes = which(modTOM>0.4)


geneInfo0 = data.frame(ESTs = genes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

modOrder = order(-abs(cor(MEs, Peso10dias, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

moretraits = exportNetworkToCytoscape(modTOM,
                               edgeFile = "AllmodulesEdgeFile.txt",
                               nodeFile = "AllmodulesNodeFile.txt",
                               weighted = TRUE,
                               threshold = 0.4,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
							   
							   
write.table(moduleTraitCor,file="moduletraitcor.txt",sep = "\t")
write.table(moduleTraitPvalue,file="moduletraitpvalue.txt",sep = "\t")

modules = c("magenta","turquoise")

genes = colnames(expression)
inModule = is.finite(match(moduleColors, modules))
modGenes = genes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)
modTOMSignificantes = which(modTOM>0.4)

geneInfo0 = data.frame(ESTs = genes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

modOrder = order(-abs(cor(MEs, Peso10dias, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

SI = exportNetworkToCytoscape(modTOM,
                               edgeFile = "mutant_posthresh0.6EdgeFile.txt",
                               nodeFile = "mutant_posthresh0.6NodeFile.txt",
                               weighted = TRUE,
                               threshold = 0.6,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
							   
							   