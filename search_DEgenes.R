#####
# Gaelle LELANDAIS <gaelle.lelandais@u-psud.fr>
#
# Ce script R est destine aux etudiants de Master. Il a pour objectif de les aider dans la
# comprehension des methodes de selection des genes differentiellement exprimes.
####

####
# Contributeurs de ce script
#
# Melina GALLOPIN; Thomas DENECKER; StÈphane LE CROM; 
# Daniel GAUTHERET
# Les Ètudiants qui ont repÈrÈ (et reprÈreront encore) des
# petites erreurs !
####

####
# Librairies R necessaires : LIMMA, MASS, DESeq2, edgeR 
#### 

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#pdf("graphics.pdf")

# Lecture des donn√©es
countData = read.table("count_dataFile.txt", header = T, row.names = 1)
# --> Ces donnees sont en relation avec l'article :
#
# "Determination of a Comprehensive Alternative Splicing Regulatory Network and 
# Combinatorial Regulation by Key Factors during the Epithelial-to-Mesenchymal Transition"
#
# Yang et al. (2016, Molecular and Cellular Biology)
#
# --> Deux points de temps ont ete etudies par RNAseq : "Day0" et "Day7", 3 replicats sont
# disponibles pour chaque point de temps.

#-----------------------------------------------------------------------------------------
# Pretraitement des donnees entre les experiences (filtrage des faibles et normalisation)
#-----------------------------------------------------------------------------------------

## Suppression des genes sans valeurs (0 pour toutes les experiences)
length(which(rowSums(countData)==0))  ## 2957
countData <- countData[-which(rowSums(countData) == 0),]

# Distribution des valeurs de comptages dans les diff√©rents √©chantillons
boxplot(log(countData + 1), ylab = "#reads (log scale)",
        main = "Boxplot of read counts", col = c(rep("grey", 3), rep("blue", 3)),
        names = c(paste("Day0_", 1:3, sep = ""), paste("Day7_", 1:3, sep = "")))
# --> Les distributions des donnees de comptage sont tres proches entre les differentes
# conditions.

#-----------------------------------------------------------------------------------------
# Etudes des logFC (log Fold Changes)
#-----------------------------------------------------------------------------------------

# Pour eviter les valeurs 0 (problematiques pour certains calculs, 
# les donnees sont legerement modifiees)
countData2 <- countData + 1

# Les replicats du temps T0 sont utilises comme reference pour l'analyse differentielle,
# calculs des logFC pour les 3 replicats du temps Day7
logFCvalues = log2(countData2[,4:6] / countData2[,1:3])
# Valeurs de logFC > 0 --> Genes induit
# Valeurs de logFC < 0 --> Genes reprimes.

# Valeurs moyennes des r√©plicats
logFCmean = apply(logFCvalues, 1, mean)

# Distribution des valeurs
hist(logFCmean, nclass = 100, main = "Distribution of logFC values",
                xlab = "logFC", ylab = "# of genes")
# --> La distribution des valeurs est centree sur 0. Cela signifie que la grande majorite
# des genes a une expression qui n'a pas bouge.
# --> Equilibre entre les genes induits (logFC > 0) et les genes reprimes (logFC < 0).

summary(logFCmean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-6.2960 -0.5017 -0.1355 -0.1038  0.2949  8.6700 

# --> Un log2FC = 6 signifie un Ratio = 2^6 = 64, c'est un differentielle tres important
# entre la condition Day7 et la condition Day0.

# Seuils logFC > 2 et logFC < -2
abline(v = 2, col = "red", lty = "dashed")
abline(v = -2, col = "green", lty = "dashed")

sum(logFCmean > 2)
# --> 467 genes sont selectionnes
sum(logFCmean < -2)
# --> 551 genes sont selectionnes

# Prise en compte de la "reproductibilite" des donnees
# --> Calcul de la deviation standard (ecart type)
logFCsd = apply(logFCvalues, 1, sd)

# Representation des logFC en fonction de SD
plot(logFCsd, logFCmean, pch = 20, 
     xlab = "Standard deviation", ylab = "logFC",
     main = "DE analysis based on logFC AND SD")

# Rappel des seuils logFC > 2 et logFC < -2
abline(h = 2, col = "red", lty = "dashed")
abline(h = -2, col = "green", lty = "dashed")

# Calcul des parametres t (statistique de Student)
tValues = logFCmean/(logFCsd*sqrt(ncol(logFCvalues)))
summary(tValues)
# --> Attention certaines valeurs de T ne sont pas calculees car la SD est nulle ....

# Distribution des valeurs de t
hist(tValues, nclass = 50, main = "Distribution of t-values", 
     xlab = "t-value (Student)", ylab = "# of genes")

# Les valeurs de T en fonction des logFC 
plot(logFCmean, tValues, pch = 20, 
     main = "Drawback associated to the classical t-value",
     xlab = "logFC", ylab = "t-value (Student)")
abline(v = 2, col = "red", lty = "dashed")
abline(v = -2, col = "green", lty = "dashed")
# gene ANGPTL2 a une valeyr de T = 280.93 !!!
countData["ANGPTL2",]
#Counts_SAQuant_Day_0_1_1 Counts_SAQuant_Day_0_2_1 Counts_SAQuant_Day_0_3_1 Counts_SAQuant_Day_7_1_1
#ANGPTL2                      111                      111                      112                      490
#Counts_SAQuant_Day_7_2_1 Counts_SAQuant_Day_7_3_1
#ANGPTL2                      490                      497

# ---> Donnees tres reproductibles !! 
# Mais est-ce suffisant pour considerer ce gene comme le "plus" interessant ?

# CONCLUSION :::: La statistique de STUDENT n'est pas un bon critere de selection des genes
# differentiellement exprimes, la deviation standard (ecrat type) a un poids trop important
# par rapport au logFC.

#-----------------------------------------------------------------------------------------
# Etudes des donnees de comptage (reads)
# Utilisation d'un modele de poisson
#-----------------------------------------------------------------------------------------

library(MASS)
countData3 = as.matrix(countData2)

# Vecteur pour garder les pvalues
pValPoisson = NULL

# L'analyse se fait gene par gene
for(i in 1:nrow(countData3)){

  print(i)
  # H0 : une seule valeur pour LAMBDA
  res = fitdistr(countData3[i, 1:6], "Poisson")
  # H1 : dex valeurs diff√©rentes pour LAMBDA
  res1 = fitdistr(countData3[i, 1:3], "Poisson")
  res2 = fitdistr(countData3[i, 4:6], "Poisson")
  
  # Ratio des vraisemblances, pour evaluer la pertinence de H0
  Dval = -2*(res$loglik - (res1$loglik + res2$loglik))
  
  # Calcul de la p-value
  Pval = 1 - pchisq(Dval, df = 1)
  
  # Cette valeur est gardee dans le vecteur des pvalues
  pValPoisson = c(pValPoisson, Pval)
    
# fin du for()
}

# ajout des noms de genes
names(pValPoisson) = row.names(countData3)

# Selection des genes p-value < 0.01
nrow(countData3[pValPoisson < 0.01,])
# --> 13656 genes sont selectionne. C'est beaucoup !....

#-----------------------------------------------------------------------------------------
# Remise en question du modele ....
#-----------------------------------------------------------------------------------------

## Graphique moyenne variance 
m <- apply(countData, 1, mean)
v <- apply(countData, 1, var)
smoothScatter(log(m+1), log(v+1), xlab="moyenne", 
              ylab="variance", 
              main = "Moyenne vs. variance")
lines(abline(0,1))
lines(lowess(m,v), col=2)

# Le modele de Poisson "predit" moyenne = variance.
# Cela n'est pas le cas ici !
# Au final, ce modele est valable quand des replicats techniques sont etudies. 
# Ici cela n'est pas le cas.

# CONCLUSION :::: Le MODELE DE POISSON n'est pas pertinent losrque des replicats biologiques
# sont etudiees. Il est necessaire d'aborder des modelisations plus complexes (cf. ci dessous).

#-----------------------------------------------------------------------------------------
# Utilisation de la methode DEseq2  
#-----------------------------------------------------------------------------------------

library(DESeq2)

## Formatage des donnees 
conds <- factor(c("Day0","Day0","Day0","Day7","Day7","Day7"))
colData <- data.frame(condition=conds)
ddsObjet <- DESeqDataSetFromMatrix(countData = countData, colData = colData,formula(~ condition))

## Taille de librairie
barplot(colSums(countData), col  = c(rep("grey", 3), rep("blue", 3)),
        main = "Library sizes")

ddsObjet <- estimateSizeFactors( ddsObjet )
sizeFactors(ddsObjet)
#1         2         3         4         5         6 
#1.1836498 1.0345727 0.9822464 0.9160292 0.9393742 0.9850124
# --> La normalisation sera ici tres "faible". Les donnees ne seront pas beaucoup modifiees. 

# Besoin d'aide ?
browseVignettes("DESeq2")

ddsEstim <- DESeq(ddsObjet)
resDESeq <- results(ddsEstim)

# Ecriture des resultats
write.table(resDESeq, row.names = T, quote = F, sep = "\t",
            file = "DESeq2_results.txt")

## Histogramme des p-values
hist(resDESeq$padj, breaks = 25, col = "grey", xlab = "Ajusted p-values", main = "")

# Quelques graphiques a discuter ...
plotDispEsts(ddsEstim)
plotMA(ddsEstim, alpha = 0.001)

#-----------------------------------------------------------------------------------------
# Analyse differentielle avec limma + voom 
#-----------------------------------------------------------------------------------------

library(limma)
browseVignettes("limma")

## On recupere le code disponible sur le userguide
library(edgeR)
y <- DGEList(counts=countData, genes=rownames(countData))
y <- calcNormFactors(y)
designLimma <- model.matrix(~conds)
v <- voom(countData,designLimma,plot=TRUE)
fitlimmaVoom <- lmFit(v,designLimma)
fitlimmaVoom <- eBayes(fitlimmaVoom)

resLimmaVoom = topTable(fitlimmaVoom, number = nrow(countData))[row.names(countData),]

# Ecriture des resultats
write.table(resLimmaVoom, row.names = T, quote = F, sep = "\t",
            file = "LIMMA-Voom_results.txt")

#-----------------------------------------------------------------------------------------
# Comparaison des differentes methodes
#-----------------------------------------------------------------------------------------

# Creation d'une table avec tous les resultats
allDiffRes = cbind(resDESeq[row.names(countData), c("baseMean", "log2FoldChange", "padj")],
                   resLimmaVoom[row.names(countData), c("logFC", "adj.P.Val")])

colnames(allDiffRes) = c("DESeq2_baseMean", "DESeq2_logFC", "DESeq2_padj",
                         "LIMMAvoom_logFC", "LIMMAvoom_padj")

# ecriture des resultats
write.table(allDiffRes, file = "allResults.txt", sep = "\t", quote = F)

focusFewGenes <- c("ESRP1","ESRP2","RBM47","QKI","NOVA1","PTBP2","MBNL1","RBFOX1","ZEB1","TWIST1")
allDiffRes[focusFewGenes,]
countData[focusFewGenes,]

#dev.off()

#####
# TRAVAIL PERSONNEL :
#
# Question 1 : Les methodes donnent-elles des resultats coherents entre elles ?
# --> Comparer les nombres genes selectionnes logFC > 2 ou logFC < -2 et padj < 0.01
# --> Les 50 "meilleurs" genes selectionnes par chaque methode sont-ils les memes ?
#
# Question 2 : Etude des genes du chromosome 18 (fichier "geneList_chrom18.txt")
# --> Avez vous identifie des genes interessants sur le chromosome 18 ?
# --> Si oui, visualisation sur IGV (10 meilleurs genes). Comment les reads sont-ils repartis sur 
# ces genes ?
####
