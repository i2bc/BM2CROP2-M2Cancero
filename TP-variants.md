Exercices Variants
================
Daniel Gautheret - I2BC & Bioinformatics Core, Gustave Roussy
17/10/2018

### I: ANALYSE D'UNE TABLE DE VARIANTS PRODUITE PAR VARSCAN

Lire la table de tous les variants d'un echantillon tumoral produit par Varscan2 (Fichier "ovc-Tpre\_varscan2\_annot.tab") et sauvegarder le résultats dans une variable "var". (Fonction read.table, avec option header=T; sep="").

Quelle est Structure de la table? (fonction str)

Quel est le nombre de variants totaux? (fonction nrow)

Quel est le nombre de mutations par chromosome? (fonction table)

Trier les chromosomes du moins au plus muté. (fonctions table et order)

Tracer la distribution des freq alleliques. (fonction hist)

Tracer les AF le long du genome. (fonction plot) (Note: ce n'est pas exactement une couverture, il faudrait additionner ref et alt.)

Filtrer les couvertures avec profondeur &gt;10. Combien reste-t-il de positions? Tracez à nouveau les AF.

les chromosomes sont-ils affichés dans le bon ordre ? Non car les noms "chrXX" ne sont pas tries par ordre alphabétique ni numérique. Il faudrait passer par une table des chromosomes. Voir script baf.r : identifie les chr et trace les limite des chr.

### II: ANALYSE SIMPLIFIEE DE LA SIGNATURE MUTATIONNELLE DOMINANTE

Dans cette section, nous allons établir la signature mutationnelle de l'echantillon et la comparer aux signatures d'Alexandrov et al. (2013)

Quelles sont les 6 signatures mutationelles possibles (décrites en se referant a la pyrimidine de la paire WC mutee)?

Le script ci-dessous a pour effet de collecter les signatures (X&gt;Y) et

``` r
#Signatures possibles:  C>T (G>A), C>A (G>T), C>G (G>C), T>C (A>G), T>A (A>T), T>G (A>C)
# cree un tableau pour convertir A>T en T>A etc. (12 à 6)

sig12=c("A>T","A>G","A>C","T>A","T>G","T>C","G>A","G>T","G>C","C>A","C>T","C>G")
sig6=c("T>A","T>C","T>G","T>A","T>G","T>C","C>T","C>A","C>G","C>A","C>T","C>G")
convertsig=sig6
names(convertsig)=sig12

# cree le vecteur de signature dans le meme ordre que Alexandrov
allsig=rep(0,6)
names(allsig)=c("C>A","C>G","C>T","T>A","T>C","T>G")

refvar=data.frame(var$reference,var$variant)
n=nrow(refvar)

for (i in (1:n)) {
    sig=paste(unlist(refvar[i,]),collapse=">") 
    if (grepl ("^.>.$", sig)) {  # verifie qu'il s'agit d'une mutation "1>1"
       sig6=convertsig[sig]
       allsig[sig6]=allsig[sig6]+1
    }
}
barplot(allsig[1:6]) 

# Une alternative + facile à comprendre mais moins élégante:
#for (i in (1:n)) {
#  sig=paste(unlist(var2[i,]),collapse=">") 
#  if (sig=="T>A"){sig6="T>A"}
#    else if (sig=="T>A"){sig6="T>A"}
#      else if (sig=="A>T"){sig6="T>A"}
#        else if (sig=="A>G"){sig6="T>C"}
#     etc.  
```

### III. ANALYSE DES VARIANTS SOMATIQUES

Lire dans la variable "som" la table des variants Somatiques produits par Varscan2 (Tumor vs Normal) dans le fichier "ovc-somatic-variants.tsv". (fonction read.table option header=T; sep="")

Afficher la structure de la table. (fonction str)

Combien y a-t-il de variants somatiques prédits chez ce patient?

Lister tous les genes concernés, avec chromosome et position

Lister tous les variants du chr 1. (Fonction: subset ou en appliquant un test sur la colonne)

Quelles sont les valeurs possibles pour SIFT? (fonction table)

Quelles sont les valeurs possibles pour Polyphen2?

Quels sont les genes avec mutations deleteres selon SIFT? Sur quels chromosomes?

Visualiser avec la commande "table" la concordance entre scores SIFT et polyphen.

Ecrire dans un fichier la liste de tous les genes et de leur changement de residu selon le format requis par Cbioportal MutationMapper. Affichez quelques lignes de ce fichier dans Cbioportal MutationMapper.

### IV. VISUALISATION SOUS IGV DE MUTATIONS SOMATIQUES

Dans le tableau ovc-somatic-variants.tsv, extraire les mutations situées sur le chr17. Afficher les informations gène, chromosome, position, somatic status.

Charger sous IGV les fichiers:

-   ovc-N-chr17.bam
-   ovc-T-chr17.bam

Visualiser les mutations délétères sous IGV. Cas intéressants:

-   FGFR3 somatique, TP53 somatique, WNT16 somatique, PDGFD,
-   LOH BRCA1: chr17:41,222,955-41,222,995
