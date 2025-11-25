
# Máster -> Detección de anomalías
# Juan Carlos Cubero. Universidad de Granada

###########################################################################
# Funciones utilizadas a lo largo del curso
###########################################################################

# rm(list=ls()) 


###########################################################################
# Realiza un plot de todos los registros
# Permite cambiar el color con el que se visualiza un conjunto de registros. 
# Dicho conjunto se especifica en el parámetro claves.a.mostrar 
# y se muestran con el segundo color especificado. 
# El primer color es el usado para los registros que NO están en claves.a.mostrar

plot_2_colores = function (datos, 
                           claves.a.mostrar, 
                           titulo = "",
                           colores = c("black", "red")){
  
  num.datos = nrow(as.matrix(datos))
  seleccionados =  rep(FALSE, num.datos)
  seleccionados[claves.a.mostrar] = TRUE
  colores.a.mostrar = rep(colores[1], num.datos)
  colores.a.mostrar [seleccionados] = colores[2]
  
  plot(datos, col=colores.a.mostrar, main = titulo)
}


claves_outliers_IQR_columna = function(columna,  coef = 1.5){
  son.outliers.IQR = son_outliers_IQR(as.data.frame(columna), 1, coef)
  return (which(son.outliers.IQR  == TRUE))
}


###########################################################################
# Función análoga a son_outliers_IQR, salvo que devuelve un vector
# de claves en vez de un vector de bools

claves_outliers_IQR = function(datos, ind.columna, coef = 1.5){
  columna.datos = datos[,ind.columna]
  return (claves_outliers_IQR_columna(columna.datos, coef))
}



###########################################################################
# Calcula los outliers IQR con respecto a una columna 
# Devuelve un vector de bools indicando si el registro i-ésimo 
# de datos es o no un outlier IQR con respecto a la columna ind.columna
# coef es 1.5 para los outliers normales y hay que pasarle 3 para los outliers extremos

son_outliers_IQR = function (datos, ind.columna, coef = 1.5){
  columna.datos = datos[,ind.columna]
  cuartil.primero = quantile(columna.datos)[2]  
  #quantile[1] es el mínimo y quantile[5] el máximo.
  cuartil.tercero = quantile(columna.datos)[4] 
  iqr = cuartil.tercero - cuartil.primero
  extremo.superior.outlier = (iqr * coef) + cuartil.tercero
  extremo.inferior.outlier = cuartil.primero - (iqr * coef)
  son.outliers.IQR  = columna.datos > extremo.superior.outlier |
    columna.datos < extremo.inferior.outlier
  return (son.outliers.IQR)
}


###########################################################################
# Calcula los outliers IQR con respecto a ALGUNA columna
# Devuelve un vector de claves indicando si el registro i-ésimo 
# de datos es o no un outlier IQR con respecto a ALGUNA columna
# coef es 1.5 para los outliers normales y  3 para los outliers extremos

claves_outliers_IQR_en_alguna_columna = function(datos, coef = 1.5){
  df.clave.columnas = data.frame()
  claves.outliers =  sapply(1:ncol(datos), 
                               function(x) claves_outliers_IQR(datos, x, coef)
  )
  claves.outliers.en.alguna.columna = unlist(claves.outliers)
  return (claves.outliers.en.alguna.columna)
}




#######################################################################
# Devuelve los nombres de aquellas filas especificadas en el parámetro claves
# filas es un vector de bools 

nombres_filas = function (datos, claves) {
  num.claves = length(claves)
  nombres.filas = row.names(as.data.frame(datos))[claves]
  
  return (nombres.filas)
}




#######################################################################
# función base para diag_caja_outliers_IQR y diag_caja

diag_caja_grafico_base = function(datos, indice.columna){
  # Importante: Para que aes busque los parámetros en el ámbito local, 
  # debe incluirse  environment = environment()
  nombre.columna = colnames(datos)[indice.columna]
  ggboxplot = ggplot(data = as.data.frame(datos), 
                     aes(x=factor(""), 
                         y = datos[,indice.columna]) , 
                     environment = environment()) + 
              xlab(nombre.columna) + ylab("") 
  return (ggboxplot)
}

#######################################################################
# Muestra un diagrama de caja
# Calcula los outliers IQR y los muestra como puntos en rojo en un BoxPlot

diag_caja_outliers_IQR = function (datos, ind.columna, coef.IQR = 1.5){
  # Si quisiéramos líneas horizontales en los límites de las cajas
  # habría que añadir 
  # + stat_boxplot(geom = 'errorbar')   
  
   outliers.IQR = son_outliers_IQR(datos, ind.columna, coef = coef.IQR)
   ggboxplot =  diag_caja_grafico_base(datos, ind.columna) + 
                stat_boxplot(coef = coef.IQR) +
                geom_boxplot(coef = coef.IQR, outlier.colour = "red") 
                # Importante: geom_boxplot debe ir después de stat_boxplot
   
   return (ggboxplot)
}



#######################################################################
# Muestra un diagrama de caja
# También muestra las etiquetas de los registros indicados en 
# el parámetro claves.a.mostrar 

diag_caja = function (datos, ind.columna, claves.a.mostrar = c()){
  num.filas = nrow(datos)
  num.claves = length(claves.a.mostrar)
  nombres.filas = vector (mode = "character", length = num.filas)
  nombres.filas = rep("", num.filas)
  nombres.claves = nombres_filas(datos, claves.a.mostrar)
  
  for (i in 1:num.claves)
    nombres.filas[claves.a.mostrar[i]]  = nombres.claves[i]
  
  ggboxplot = diag_caja_grafico_base(datos, ind.columna) + 
    geom_boxplot(outlier.shape = NA) + # Para que no imprima los outliers IQR calculados dentro del mismo geom_boxplot
    geom_text(aes(label = nombres.filas)) 
  
  return (ggboxplot)
}





#######################################################################
# Muestra de forma conjunta todos los diagramas de caja de las variables de datos
# Para ello, normaliza previamente los datos.
# También muestra las etiquetas de los registros indicados en claves.a.mostrar

diag_caja_juntos = function (datos, titulo = "", claves.a.mostrar = c()){  
  # Importante: Para que aes busque los parámetros en el ámbito local, 
  # debe incluirse  environment = environment()
  
  # Para hacerlo con ggplot, lamentablemente hay que construir antes una tabla 
  # que contenga en cada fila el valor que a cada tupla le da cada variable 
  # -> usamos pivot_longer de tidyverse 
  
  
  datos = as.data.frame(scale(datos))
  nombres.de.filas = nombres_filas (datos, claves.a.mostrar)
  
  # Añadimos una columna con el nombre de las filas y pivotamos

  Row = c("Row")  # Nombre que le damos a la columna con los nombres de las filas
  
  datos.pivotados <- 
    datos %>% 
    tibble::rownames_to_column(var = Row) %>% 
    pivot_longer(-Row, names_to = "Variables", values_to = "ZScore")
  
  columna.factor = as.factor(as.data.frame(datos.pivotados)[,Row])
  
  # Eliminamos los que no están en claves.a.mostrar
  levels(columna.factor)[!levels(columna.factor) %in% nombres.de.filas] = ""  
  
  ggplot(data = datos.pivotados, 
         aes(x=Variables, y=ZScore), 
         environment = environment()) + 
    ggtitle(titulo) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_text(aes(label = columna.factor), size = 3) 
}



#######################################################################
# Muestra un biplot del conjunto de datos
# Se muestran los nombres de los registros indicados en claves.a.mostrar
# El color usado para dichos registros es el segundo del parámetro colores
# El título para el grupo de dichos registros es el especificado en titulo.grupo.a.mostrar
# El parámetro titulo especifica el título principal del gráfico

biplot_2_colores = function (datos, 
                             claves.a.mostrar = c(), 
                             titulo = "",
                             titulo.grupo.a.mostrar = "Outliers",
                             colores = c("black","red")){
  nombres = rownames(datos)
  claves.datos = c(1:nrow(datos))
  son.a.mostrar = claves.datos %in% claves.a.mostrar
  nombres[!son.a.mostrar] = ''
  
  PCA.model = princomp(scale(datos))
  outlier.shapes = c(".","x") 
  biplot = ggbiplot(PCA.model,
                    obs.scale = 1,
                    var.scale = 1 ,
                    varname.size = 5,
                    groups =  son.a.mostrar,
                    alpha = 1/2) #alpha = 1/10
  biplot = biplot + labs(color = titulo.grupo.a.mostrar)
  biplot = biplot + scale_color_manual(values = colores)
  biplot = biplot + geom_text(label = nombres,
                              stat = "identity",
                              size = 3,
                              hjust=0,
                              vjust=0)
  biplot = biplot + ggtitle(titulo)
}


#######################################################################
# Muestra un biplot de un conjunto de datos diferenciados por color
# El color lo determina la asignación de cada dato a un cluster 
# Las asignaciones de datos a cluster se indican en asignaciones.clustering
# También se muestran los outliers cuyas claves vienen indicadas en claves.outliers
 
biplot_outliers_clustering = function(datos, 
                                      titulo = "Outliers por el método de Clustering", 
                                      titulo.color = "Asignaciones Clustering",
                                      titulo.outlier = "Outliers",
                                      asignaciones.clustering,
                                      claves.outliers){
  son.outliers = rep(FALSE, nrow(datos))
  son.outliers[claves.outliers] = TRUE
  
  bip = biplot_colores_formas(datos, 
                              titulo, titulo.color, titulo.outlier,
                              asignaciones.clustering,
                              son.outliers,
                              claves.outliers)
  bip 
}

#######################################################################
# Muestra un biplot del conjunto de datos
# Los datos se muestran diferenciados por color y por forma
# Las asignaciones de cada dato a su color y forma vienen dadas por los vectores
# asignaciones.colores y asignaciones.formas 
# También se muestran las etiquetas de los registros indicados
# en el parámetro opcional claves.a.mostrar 

biplot_colores_formas = function (datos, 
                                  titulo, titulo.color = '', titulo.forma = '', 
                                  asignaciones.colores, asignaciones.formas,
                                  claves.a.mostrar = c()){
  PCA.model = princomp(scale(datos))
  
  son.a.mostrar = rep(FALSE, nrow(datos))
  son.a.mostrar[claves.a.mostrar] = TRUE
  nombres.a.mostrar = rownames(datos)
  nombres.a.mostrar[!son.a.mostrar] = ''

  asignaciones.colores = factor(asignaciones.colores)
  asignaciones.formas  = factor(asignaciones.formas)

  
  bip = ggbiplot(PCA.model, obs.scale = 1, var.scale=1 , varname.size = 3, alpha = 0) +              
    geom_point(aes(shape = asignaciones.formas, colour = asignaciones.colores))  +
    labs(shape = titulo.forma) +
    labs(colour = titulo.color) +
    ggtitle(titulo) +
    geom_text(label = nombres.a.mostrar, stat = "identity", size = 3, hjust=0, vjust=0)      
  
  bip
}

#######################################################################
# Calcula las distancias de cada dato al centroide de su cluster
# Las asignaciones de cada dato a su cluster se indican en asignaciones.clustering
# Cada centroide es una fila del data frame datos.centroides.normalizados

distancias_a_centroides = function (datos.normalizados, 
                                    asignaciones.clustering, 
                                    datos.centroides.normalizados){
  
  sqrt(rowSums(   (datos.normalizados 
                   - 
                   datos.centroides.normalizados[asignaciones.clustering,])^2  ))
}


#######################################################################
# Revierte la función de normalización (z-score)

desnormaliza = function(datos, filas.normalizadas){
  medias        = colMeans(datos)
  desviaciones  = apply(datos, 2, sd , na.rm = TRUE)
  
  filas.desnormalizadas  = sweep(filas.normalizadas, 2, desviaciones, "*")
  filas.desnormalizadas  = sweep(filas.desnormalizadas, 2, medias, "+")
  
  filas.desnormalizadas 
}




top_clustering_outliers = function(datos.norm, 
                                   asignaciones.clustering, 
                                   datos.centroides.norm, 
                                   num.outliers){
  
  dist_centroides = distancias_a_centroides (datos.norm, 
                                             asignaciones.clustering, 
                                             datos.centroides.norm)
  
  claves = order(dist_centroides, decreasing=T)[1:num.outliers]
  
  list(distancias = dist_centroides[claves]  , claves = claves)
}





