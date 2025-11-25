# Master -> Detección de anomalías
# Juan Carlos Cubero. Universidad de Granada

# Instale los paquetes tal y como se especifica a continuación:

install.packages("ggplot2")
install.packages("tidyverse")

# ggbiplot:
install.packages("devtools")  # Para poder instalar paquetes de Github
library(devtools)
install_github("vqv/ggbiplot")  # Si pregunta, seleccione la opción 3.

# El directorio en el que se instale no debe tener caracteres como acentos, ?, etc.
# 
# Si no funciona lo anterior, probar alguna de las siguientes:
# devtools::install_github("richardjtelford/ggbiplot", ref = "experimental")
# install_github("ggbiplot", "vqv")
# install_github("vqv/ggbiplot", ref = "experimental")

install.packages("fitdistrplus")  # Ajuste de una distribución
install.packages("outliers")      # Grubbs
install.packages("CerioliOutlierDetection")  #MCD Hardin Rocke estimación robusta de la matriz de covarianzas
install.packages("mvoutlier")     # corr.plot , MCD ChiC 
install.packages("mvnormtest")    # mvn. Test Normalidad multivariante
install.packages("MVN")           # Test de normalidad multivariante     
install.packages("DescTools")     # LOF  (También en Rlof, dbscan, DMwR, etc.)
install.packages("cluster")       # pam
install.packages("ggpubr")        # ggqqplot
install.packages("heplots")       # cqplot qq plot para Mahalanobis

