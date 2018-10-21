#' 
#' 
#' ## Premiers pas avec R
#' 
#' 
## ----results="hold"------------------------------------------------------
2+1

3-2

4*8

16/5

7.15 * sqrt(9)

#' 
#' On peut sauvegarder un résultat en affectant une valeur à une variable, à l'aide du symbole "&lt;-". On peut ainsi réutiliser ces variables pour effectuer des calculs plus complexes.
#' 
## ---- results="hold"-----------------------------------------------------
x <- 2+1

y <- 7.15 * sqrt(4)

x + y

#' 
#' 
#' On peut combiner des objects dans une même liste en utilisant la commande c (c comme "combine").
#' 
## ------------------------------------------------------------------------
monObjet <- c(1,2,3,4,10,59)

#' 
#' On peut ainsi avoir accès au contenu de l'objet en tapant le nom de l'object:
## ------------------------------------------------------------------------
monObjet

#' 
#' 
#' On peut effectuer des opérations sur l'objet:
#' 
## ----results="hold"------------------------------------------------------
monObjet + 19

monObjet * 10

#' 
#' D'autres options sont également disponibles:
#' 
## ----results="hold"------------------------------------------------------
sum(monObjet)
mean(monObjet)
max(monObjet)
min(monObjet)
sd(monObjet)
range(monObjet)
length(monObjet)

#' 
#' 
#' 
#' 
#' 
#' ## Types d'objet sous R 
#' 
#' 
## ----results="hold"------------------------------------------------------
monObjet
monObjet[1]        # Accède au premier élement de monObjet
monObjet[5]        # Accède au 5ème élement de monObjet

#' 
#' Une liste peut aussi contenir autre chose que des nombres, par exemple des chaînes de caractères. 
## ------------------------------------------------------------------------
names <- c("Jack", "Alba", "Irene")
names[2]
estUneFille <- c(FALSE, TRUE, TRUE)
estUneFille[2]

#' 
#' 
#' ![](Rclasses.png)
#' 
#' 
#' 
#' 
## ----results="hold"------------------------------------------------------
monObjet

monObjet > 2

#' 
#' 
#' 
## ----results="hold"------------------------------------------------------
sum(monObjet > 2)

monObjet[monObjet > 2]

#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ## Les graphiques 
#' 
## ------------------------------------------------------------------------
plot(monObjet)

#' 
#' 
## ------------------------------------------------------------------------
barplot(monObjet)

pie(monObjet)

#' 
#' 
## ------------------------------------------------------------------------
flowerdata <- iris      

#' 
#' 
## ------------------------------------------------------------------------
class(flowerdata)

#' 
#' Un "data.frame" est une table. 
#' 
#' 
## ------------------------------------------------------------------------
dim(flowerdata)         

#' 
#' 
## ------------------------------------------------------------------------
flowerdata[11,2]        

flowerdata[11,]        

flowerdata[,2]           

#' 
#' 
#' 
## ------------------------------------------------------------------------
plot(flowerdata[,2])

hist(flowerdata[,2])

boxplot(flowerdata, col=rainbow(5))

#' 
#' 
#' ```{}
#' pdf("flower_images.pdf")
#' 
#' boxplot(flowerdata, col=rainbow(5))
#' 
#' dev.off()
#' ```
#' 
#' 
#' 
#' ## L'essentiel à retenir
#' 
#' Pour comprendre le fonctionnement de R, il faut comprendre la phrase suivante: 
#' 
#' >  _To understand computations in R, two slogans are helpful:
#' >       Everything that exists is an object.
#' >       Everything that happens is a function call.
#' >                     (John Chamber)_
#' 
#' Quelques remarques:
#' 
#'   * La commande `ls()` permet d'obtenir la liste des objets existants dans votre environnement.
#'   * Il existe différentes classes d'objet : matrice, vecteur, chaînes de charactères, fonctions... 
#'   * La fonction `class()` permet de connaître la classe d'un objet.
#'   * La fonction `str()` permet de connaître la structure de l'objet.
#'   * De nombreuses fonctions sont disponibles sous R: la fonction cosinus `cos()`, `plot()` pour réaliser des graphiques...
#'   * L'utilisation de certaines fonctions avancées nécéssite l'installation de _packages_.
#' 
#' 
