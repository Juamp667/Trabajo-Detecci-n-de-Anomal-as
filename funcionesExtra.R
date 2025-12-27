library(htmlwidgets)
library(htmltools)
library(IRdisplay)
library(tidyverse)
library(patchwork)
library(GGally)
library(nortest)
library(ggplot2)
library(fitdistrplus)  # Ajuste de una distribución -> denscomp 
library(ggpubr)    # ggqqplot -para los gráficos QQ-
library(ggbiplot)  # biplot
library(outliers)  # Grubbs
library(MVN)       # mvn: Test de normalidad multivariante  
library(CerioliOutlierDetection)  #MCD Hardin Rocke
library(mvoutlier) # corr.plot 
library(DescTools) # lof
library(cluster)   # PAM
library(heplots)
library(plotly)
library(EnvStats)


#' Cambiar las dimensiones de gráficas de salida
#'
#' @param w Ancho (width) de salida
#' @param h Alto (height) de salida
#'
#' @examples
#' plot_dim(10,5)
#'
plot_dim <- \(w,h) options(repr.plot.width=w, repr.plot.height=h)


#' Retorna una lista de gráficas tipo histograma+boxplot con 
#' medias+-std sobre todas las columnas de ds
#'
#' @param ds Dataset cuyas columnas desea graficar
#' @param bins_=50 Número de cajas de los histogramas
#'
#' @examples
#' feat_plots(ds)
#' 
feat_plots <- function(ds, bins_=50) {
  lapply(
    colnames(ds)[1:length(ds)],
    \(col_name){

      # Calcular estadísticas
      x_mean <- mean(ds[[col_name]], na.rm = TRUE)
      x_sd   <- sd(ds[[col_name]], na.rm = TRUE)

      stats_df <- data.frame(
        x = c(x_mean, x_mean - x_sd, x_mean + x_sd),
        type = c("Media", "-1 SD", "+1 SD")
      )


      tema <- theme_minimal(base_size = 16) +
        theme(
          axis.title.x = element_text( size = 18),
          axis.text.x  = element_text(size = 16),
          axis.text.y  = element_text(size = 16),
          legend.text = element_text(size = 18)
        )
      p1 <- ds %>%
        ggplot(aes(x = .data[[col_name]])) +
        geom_histogram(fill = "cadetblue2", color = "black", bins = bins_) +
        
        geom_vline(
          data = stats_df,
          aes(xintercept = x, color = type, linetype = type),
          linewidth = 1
        ) +
        
        scale_color_manual(
          name = "Estadísticos",
          values = c(
            "Media" = "#1F4E5F",
            "-1 SD" = "#5FA4B8",
            "+1 SD" = "#5FA4B8"
          )
        ) +
        scale_linetype_manual(
          name = "Estadísticos",
          values = c(
            "Media" = "solid",
            "-1 SD" = "dashed",
            "+1 SD" = "dashed"
          )
        ) +
        labs( title = col_name, y="Cuentas") +
        tema

      p2 <- ds %>%
        ggplot(aes(x = .data[[col_name]])) +
        geom_boxplot(fill = "cadetblue2", color = "black") +
        scale_y_continuous(breaks = c()) +
        # scale_x_continuous(breaks = c()) +
        # labs(x = "")
        tema

      (p1 / p2) + plot_layout(heights = c(3, 1))
    }
  )
}


#' Modifificación feat_plots con densidades en vez de histogramas
#' 
feat_plots2 <- function(ds, cluster_var_name, bins_=50) {
  lapply(
    colnames(ds),
    \(col_name){
      if (col_name!=cluster_var_name){
        p1 <- ds %>%
              arrange(.data[[cluster_var_name]]) %>%
              ggplot(aes(x = .data[[col_name]], group = factor(.data[[cluster_var_name]]), fill = factor(.data[[cluster_var_name]]), alpha = factor(.data[[cluster_var_name]]))) +
              geom_density(color = "black") +
              scale_alpha_manual(values = c("0" = 1, "1" = 0.5)) +
              labs(x = "", title = col_name) +
              tema1

        p2 <- ds %>%
          ggplot(aes(x = .data[[col_name]], group= .data[[cluster_var_name]], fill = .data[[cluster_var_name]])) +
          geom_boxplot(color = "black") +
          scale_y_continuous(breaks = c()) +
          scale_x_continuous(breaks = c()) +
          labs(x = "")

        (p1 / p2) + plot_layout(heights = c(3, 1))
        }
    }
  )
}


tema <- theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text( size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text = element_text(size = 18)
  )

tema1 <- theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text( size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text = element_text(size = 18)
  )


# Funcion para obtener info IQR sobre un vector
quantiles_IQR_analysis <- function(x){
  quants <- quantile(x, c(0.25,0.75))
  q25 <- as.numeric(quants[1])
  q75 <- as.numeric(quants[2])
  iqr <- q75-q25

  sup_o <- q75 + 1.5*iqr # Límite superior para calificar como outlier
  sup_oe <- q75 + 3.0*iqr # Límite superior para calificar como outlier extremo

  inf_o <- q25 - 1.5*iqr # Límite inferior para calificar como outlier
  inf_oe <- q25 - 3.0*iqr # Límite inferior para calificar como outlier extremo

  son_o <- (x<inf_o)|(x>sup_o)  #Son ourtliers
  son_oe <- (x<inf_oe)|(x>sup_oe)  #Son ourtliers extremos

  #Retorno una lista con todos los datos
  list(
    iqr=iqr,
    cuartil.primero=q25,
    cuartil.tercero=q75,
    extremo.superior.IQR=sup_o,
    extremo.inferior.IQR=inf_o,
    extremo.superior.IQR.extremo=sup_oe,
    extremo.inferior.IQR.extremo=inf_oe,
    son.outliers.IQR=son_o,
    son.outliers.IQR.extremos=son_oe
  )
}



#' Estudiar normalidad de los datos de un dataframe mediante
#' qqplots, test y comparación distrib-histograma.
#' 
#' @param ds data.frame con columnas a analizar
#' @param bins_ número de bins de los histogramas
#'
#' @return lista de gráficos, uno por columna
#' 
#' @examples
#' feat_plots_normality(ds)
#'
feat_plots_normality <- function(ds, bins_ = 50) {
  lapply(
    colnames(ds)[1:length(ds)],
    \(col_name) {
      
      # Vector de datos (sin NAs)
      x <- ds[[col_name]]
      x_no_na <- x[!is.na(x)]
      
      # Calcular estadísticas
      x_mean <- mean(x_no_na)
      x_sd   <- sd(x_no_na)
      
      stats_df <- data.frame(
        x    = c(x_mean, x_mean - x_sd, x_mean + x_sd),
        type = c("Media", "-1 SD", "+1 SD")
      )
      
      # Test de Lilliefors (normalidad)
      lillie_res <- lillie.test(x_no_na)
      pval       <- lillie_res$p.value
      p_label    <- paste0(
        "Lilliefors p-value = ",
        formatC(pval, format = "e", digits = 2)
      )
      
      # Histograma + normal teórica + líneas de media y ±1SD
      p1 <- data.frame(x = x_no_na) %>%
        ggplot(aes(x = x)) +
        geom_histogram(
          aes(y = ..density..),
          fill  = "cadetblue2",
          color = "black",
          bins  = bins_
        ) +
        # Curva normal teórica ajustada (misma media y sd)
        stat_function(
          fun  = dnorm,
          args = list(mean = x_mean, sd = x_sd),
          color = "red",
          linewidth = 1.2
        ) +
        # Líneas verticales de media y ±1 SD
        geom_vline(
          data  = stats_df,
          aes(xintercept = x, color = type, linetype = type),
          linewidth = 1
        ) +
        scale_color_manual(
          name = "Estadísticos",
          values = c(
            "Media" = "#1F4E5F",
            "-1 SD" = "#5FA4B8",
            "+1 SD" = "#5FA4B8"
          )
        ) +
        scale_linetype_manual(
          name = "Estadísticos",
          values = c(
            "Media" = "solid",
            "-1 SD" = "dashed",
            "+1 SD" = "dashed"
          )
        ) +
        labs(
          x        = "",
          y        = "Densidad",
          title    = col_name,
          subtitle = p_label
        ) +
        theme_minimal(base_size = 16) +
        theme(
          axis.title.x = element_text(face = "bold", size = 18),
          axis.text.x  = element_text(size = 16)
        )
      
      # Q–Q plot
      p2 <- ggplot(data.frame(x = x_no_na), aes(sample = x)) +
        stat_qq(color = "cadetblue4") +
        stat_qq_line(color = "red") +
        labs(
          x = "Cuantiles teóricos",
          y = "Cuantiles muestrales"
        ) +
        theme_minimal(base_size = 14)
      
      # Combinación: histograma arriba, QQ abajo
      (p1 / p2) 
    }
  )
}


grubbs_recursive_fwer <- function(x, alpha) {
    x_original <- x
    removed_idx <- c()
    iteration <- 0
    
    repeat {
        iteration <- iteration + 1
        
        # Grubbs test
        gt <- grubbs.test(x, two.sided = TRUE)
        pval <- gt$p.value
        
        # FWER: 1 - (1-alpha)^iteration
        fwer <- 1 - (1 - alpha)^iteration
        
        # Parar si deja de ser significativo o se supera FWER
        if (pval >= alpha || fwer >= 0.05) break
        
        # Identificar el outlier extremo
        outlier <- as.numeric(strsplit(gt$alternative, " ")[[1]][3])
        idx <- which(x == outlier)[1]
        
        # Guardar índice global
        global_idx <- which(x_original == outlier)[1]
        removed_idx <- c(removed_idx, global_idx)
        
        # Eliminar
        x <- x[-idx]
        
        # Parar si quedan < 3 valores
        if (length(x) < 3) break
    }
    
    return(list(clean = x,
                removed = removed_idx,
                n_outliers = length(removed_idx)))
}


pca_3d <- function(data, op1=0.5, op2=1,
                   labels = NULL,
                   clusters = NULL,
                   scale_data = TRUE,
                   width = 1000, height = 800,
                   label_clusters = NULL,
                   marker_size = 4,
                   cluster_titles = NULL,
                   cluster_opacity = NULL,   # <- NUEVO: opacidad por cluster
                   title=""
                   ) {

  data_num <- data[, sapply(data, is.numeric), drop = FALSE]
  pca <- prcomp(data_num, scale. = scale_data)

  pcs <- as.data.frame(pca$x[, 1:3])
  colnames(pcs) <- c("PC1", "PC2", "PC3")

  if (!is.null(labels)) {
    if (length(labels) != nrow(pcs)) stop("labels debe tener la misma longitud que las filas del dataset.")
    pcs$label <- labels
  } else pcs$label <- ""

  var_exp <- round(100 * summary(pca)$importance[2, 1:3], 2)

  loadings <- as.data.frame(pca$rotation[, 1:3])
  loadings$varname <- rownames(loadings)

  scale_factor <- max(abs(as.matrix(pcs[,1:3]))) * 0.8
  loadings[,1:3] <- loadings[,1:3] * scale_factor

  # --- helpers ---
  cluster_name <- function(lev, levs) {
    if (is.null(cluster_titles)) return(paste0("Cluster ", lev))

    if (!is.null(names(cluster_titles)) && lev %in% names(cluster_titles))
      return(as.character(cluster_titles[[lev]]))

    pos <- match(lev, levs)
    if (!is.na(pos) && pos <= length(cluster_titles))
      return(as.character(cluster_titles[[pos]]))

    paste0("Cluster ", lev)
  }

  cluster_alpha <- function(lev, levs, show_text) {
    if (!is.null(cluster_opacity)) {
      # caso vector nombrado
      if (!is.null(names(cluster_opacity)) && lev %in% names(cluster_opacity))
        return(as.numeric(cluster_opacity[[lev]]))

      # caso vector en orden según levels()
      pos <- match(lev, levs)
      if (!is.na(pos) && pos <= length(cluster_opacity))
        return(as.numeric(cluster_opacity[[pos]]))
    }
    # fallback: tu lógica original
    if (show_text) op2 else op1
  }

  # --- FIGURA ---
  if (!is.null(clusters)) {

    if (length(clusters) != nrow(pcs)) stop("clusters debe tener misma longitud que las filas del dataset.")
    pcs$cluster <- factor(clusters)
    levs <- levels(pcs$cluster)

    if (is.null(label_clusters)) label_clusters <- character(0)
    label_clusters <- as.character(label_clusters)

    # init
    first_lev <- levs[1]
    df0 <- pcs[pcs$cluster == first_lev, , drop = FALSE]

    show_text0 <- first_lev %in% label_clusters
    alpha0 <- cluster_alpha(first_lev, levs, show_text0)

    fig <- plot_ly(
      data = df0,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = if (show_text0) "markers+text" else "markers",
      text = if (show_text0) df0$label else "",
      textposition = "top center",
      hoverinfo = if (show_text0) "text" else "none",
      marker = list(size = marker_size, opacity = alpha0),
      name = cluster_name(first_lev, levs),
      width = width, height = height
    )

    # resto
    if (length(levs) > 1) {
      for (lev in levs[-1]) {
        dfk <- pcs[pcs$cluster == lev, , drop = FALSE]
        show_textk <- lev %in% label_clusters
        alphak <- cluster_alpha(lev, levs, show_textk)

        fig <- fig %>% add_trace(
          data = dfk,
          x = ~PC1, y = ~PC2, z = ~PC3,
          type = "scatter3d",
          mode = if (show_textk) "markers+text" else "markers",
          text = if (show_textk) dfk$label else "",
          textposition = "top center",
          hoverinfo = if (show_textk) "text" else "none",
          marker = list(size = marker_size, opacity = alphak),
          name = cluster_name(lev, levs)
        )
      }
    }

  } else {

    fig <- plot_ly(
      data = pcs,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = "markers+text",
      text = ~label,
      textposition = "top center",
      hoverinfo = "text",
      marker = list(size = 3, opacity = 0.5),
      width = width, height = height
    )
  }

  # loadings
  for (i in 1:nrow(loadings)) {
    fig <- fig %>% add_trace(
      x = c(0, loadings$PC1[i]),
      y = c(0, loadings$PC2[i]),
      z = c(0, loadings$PC3[i]),
      type = "scatter3d",
      mode = "lines+text",
      text = loadings$varname[i],
      textposition = "top center",
      line = list(width = 6),
      showlegend = FALSE,
      inherit = FALSE
    )
  }

  fig %>% layout(
    title = ifelse(title!="", title, "PCA en 3 dimensiones"),
    scene = list(
      xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
      yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
      zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
    )
  )
}



encontrar_k_optimo <- function(data, k_min = 2, k_max = 10, nstart = 10, graficar = TRUE) {
  
  # Guardamos los valores de silueta por cada k
  silhouette_scores <- numeric(k_max)
  
  for (k in k_min:k_max) {
    km <- kmeans(data, centers = k, nstart = nstart)
    sil <- silhouette(km$cluster, dist(data))
    silhouette_scores[k] <- mean(sil[, 3])
  }
  
  # Determinar el k óptimo
  k_optimo <- which.max(silhouette_scores)
  
  # Entrenar el modelo final
  modelo_final <- kmeans(data, centers = k_optimo, nstart = nstart)
  
  # Graficar si se desea
  if (graficar) {
    plot(k_min:k_max,
         silhouette_scores[k_min:k_max],
         type = "b",
         xlab = "Número de clusters (k)",
         ylab = "Silhouette promedio",
         main = "Selección óptima de k con Silhouette")
    abline(v = k_optimo, col = "red", lty = 2)
  }
  
  return(list(
    k_optimo = k_optimo,
    silhouette = silhouette_scores,
    modelo = modelo_final
  ))
}


encontrar_k_optimo_kmedoids <- function(data, 
                                        k_min = 2, 
                                        k_max = 10, 
                                        graficar = TRUE) {
  
  # Calcular distancias una sola vez
  d <- dist(data)

  # Vector para guardar siluetas
  silhouette_scores <- numeric(k_max)

  for (k in k_min:k_max) {
    pam_model <- pam(d, k = k, diss = TRUE)
    sil <- silhouette(pam_model$clustering, d)
    silhouette_scores[k] <- mean(sil[, 3])
  }

  # Elegir el k óptimo
  k_optimo <- which.max(silhouette_scores)

  # Modelo final con k óptimo
  modelo_final <- pam(d, k = k_optimo, diss = TRUE)

  # Gráfica opcional
  if (graficar) {
    plot(k_min:k_max,
         silhouette_scores[k_min:k_max],
         type = "b",
         xlab = "Número de clusters (k)",
         ylab = "Silhouette promedio",
         main = "Selección óptima de k con Silhouette (k-medoids)")
    abline(v = k_optimo, col = "red", lty = 2)
  }
  
  return(list(
    k_optimo = k_optimo,
    silhouette = silhouette_scores,
    modelo = modelo_final
  ))
}


rosner_cotas <- function(rosner_obj, x, only_detected = TRUE) {
  if (is.null(rosner_obj$all.stats)) {
    stop("El objeto no tiene $all.stats. ¿Seguro que viene de EnvStats::rosnerTest?")
  }

  st <- rosner_obj$all.stats
  nms <- names(st)

  # 1) localizar Lambda (crítico) en all.stats
  lambda_col <- {
    cand <- grep("lambda|crit|critical|cv|cutoff|threshold", nms, ignore.case = TRUE, value = TRUE)
    if (length(cand) == 0) stop("No encuentro la columna Lambda/critical en all.stats.")
    cand[1]
  }

  # 2) localizar el índice de observación (Obs.Num) o, si no, el valor candidato (Value)
  obsnum_col <- {
    cand <- grep("obs\\.num|obsnum|obs_num|obs", nms, ignore.case = TRUE, value = TRUE)
    if (length(cand) == 0) NA_character_ else cand[1]
  }
  value_col <- {
    cand <- grep("^value$|value\\b|x\\.value|candidate", nms, ignore.case = TRUE, value = TRUE)
    if (length(cand) == 0) NA_character_ else cand[1]
  }

  k <- nrow(st)
  m <- rosner_obj$n.outliers
  k_use <- if (isTRUE(only_detected)) max(0L, min(k, m)) else k
  if (k_use == 0L) return(st[0, , drop = FALSE])

  st_use <- st[1:k_use, , drop = FALSE]

  # índice del punto eliminado en cada paso (preferimos Obs.Num; si no, lo buscamos por Value)
  removed_idx <- integer(k_use)

  if (!is.na(obsnum_col)) {
    removed_idx <- as.integer(st_use[[obsnum_col]])
  } else if (!is.na(value_col)) {
    # OJO: si hay valores repetidos, esto puede elegir el primero que coincida.
    vals <- st_use[[value_col]]
    for (i in seq_len(k_use)) {
      idx <- which(x == vals[i])[1]
      if (is.na(idx)) stop("No pude mapear un Value de all.stats a un índice en x (¿NA o no coincide?).")
      removed_idx[i] <- idx
    }
  } else {
    stop("No encuentro Obs.Num ni Value en all.stats para saber qué observación se elimina en cada paso.")
  }

  # Recalcular mean/sd iterativas y cotas
  x_work <- x
  idx_work <- seq_along(x)

  out <- data.frame(
    step = seq_len(k_use),
    n = integer(k_use),
    mean = numeric(k_use),
    sd = numeric(k_use),
    lambda = as.numeric(st_use[[lambda_col]]),
    lower = numeric(k_use),
    upper = numeric(k_use),
    removed_index = integer(k_use),
    removed_value = numeric(k_use)
  )

  for (i in seq_len(k_use)) {
    mu <- mean(x_work, na.rm = TRUE)
    s  <- sd(x_work, na.rm = TRUE)

    out$n[i] <- length(x_work)
    out$mean[i] <- mu
    out$sd[i] <- s
    out$lower[i] <- mu - out$lambda[i] * s
    out$upper[i] <- mu + out$lambda[i] * s

    # localizar en el "work set" el índice original a eliminar
    orig_idx <- removed_idx[i]
    pos <- which(idx_work == orig_idx)
    if (length(pos) == 0) {
      stop("El índice a eliminar ya no está en el conjunto de trabajo (¿inconsistencia en all.stats?).")
    }

    out$removed_index[i] <- orig_idx
    out$removed_value[i] <- x_work[pos[1]]

    # eliminar para el siguiente paso
    x_work <- x_work[-pos[1]]
    idx_work <- idx_work[-pos[1]]
  }

  out
}


# interactive_plot <- function(plot, plot_idx, width = 1000, height = 800){

#     plot_file = paste("interactivePlots/p",plot_idx,".html", sep="") 

#     saveWidget(
#         widget = plot,  
#         file = plot_file, 
#         selfcontained =   TRUE)

#     html <- paste(readLines(plot_file, warn = FALSE), collapse = "\n")

#   IRdisplay::display_html(
#     as.character(
#       htmltools::tags$iframe(
#         srcdoc = html,
#         style = paste0(
#           "width:", width, "px; height:", height, "px; border:0;",
#           " display:block; margin:0 auto;"
#         )
#       )
#     )
#   )

# }

interactive_plot <- function(plot, plot_idx, width = 1000, height = 800){

  plot_file <- paste0("interactivePlots/p", plot_idx, ".html")

  plot <- plotly::config(plot, displayModeBar = FALSE)
  plot <- plotly::toWebGL(plot)
  plot <- plotly::style(plot, hoverinfo = "none")

  htmlwidgets::saveWidget(
    widget = plot,
    file = plot_file,
    selfcontained = TRUE
  )


  html <- paste(readLines(plot_file, warn = FALSE), collapse = "\n")

  IRdisplay::display_html(
    as.character(
      htmltools::tags$iframe(
        srcdoc = html,
        style = paste0(
          "width:", width, "px; height:", height, "px; border:0;",
          " display:block; margin:0 auto;"
        )
      )
    )
  )

}