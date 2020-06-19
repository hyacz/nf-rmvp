#---------------------------------------------------------------------
# Script Name:   Report.R
# Description:   Some visualization functions and auxiliary data 
#                cleaning functions
# Author:        Haohao Zhang
# Created:       2018-12-20 at 20:53
# Last Modified: 2020-06-04 at 22:01
#--------------------------------------------------------------------

filter.points <- function(x, y, w, h, dpi=300, scale=1) {
    x <- ceiling((x - min(x)) / (max(x) - min(x)) * w * dpi / scale)
    y <- ceiling((y - min(y)) / (max(y) - min(y)) * h * dpi / scale)
    index <- !duplicated(cbind(x, y))
    return(index)
}

classify.points <- function(x, threshold) {
    style.points <- matrix(NA, length(x), length(threshold))  # identify which style points belong to.
    th <- c(Inf, threshold)
    for (t in 2:length(th)) {
        style.points[, t - 1] <- (x >= th[t] & x < th[t - 1])
    }
    return(style.points)
}

#' Title
#'
#' @param Pmap CHOR, POS, P-Value
#' @param taxa.name 
#' @param pch 
#' @param cex 
#' @param cex.axis 
#' @param cex.lab 
#' @param col 
#' @param title the title of the image, if NULL, will be set to "Distribution 
#'    of \code{taxa.name}"
#' @param xlab The title of the x axis.
#' @param ylab The title of the y axis.
#' @param density 
#' @param density.bin.size 
#' @param density.bin.max 
#' @param density.col 
#' @param signal 
#' @param signal.pch 
#' @param signal.cex 
#' @param signal.col 
#' @param threshold 
#' @param threshold.col 
#' @param threshold.lty 
#' @param threshold.lwd 
#' @param vline 
#' @param vline.col 
#' @param vline.lty 
#' @param vline.lwd 
#' @param file.memo A string that is the prefix of the output image file.
#' @param file.type A string or NULL is used to determine the type of output 
#'    file. Can be "jpg", "pdf", "tiff". If it is NULL, it will use 
#'    \code{\link[grDevices:dev]{dev.new()}} to create a new graphics device 
#'    in the current environment, which may be RStudioGD or the default 
#'    device of the system.
#' @param file.width The width of the output file, in inches.
#' @param file.height The height of the output file, in inches.
#' @param file.dpi The resolution of the image, specifying how many pixels 
#'    per inch.
#'
#' @return
#' @export
#'
#' @examples
drawManhattanPlot <-
function(Pmap,
         taxa.name = "rMVP",
         pch = 19,
         cex = 0.5,
         cex.axis = 1,
         cex.lab = 1.2,
         col = c("dodgerblue4", "deepskyblue"),
         title = NULL,
         xlab = NULL,
         ylab = NULL,
         density = TRUE,
         density.bin.size = 1e6,
         density.bin.max = NULL,
         density.col = c("darkgreen", "yellow", "red"),
         signal = TRUE,
         signal.pch = 19,
         signal.cex = 0.8,
         signal.col = "red",
         threshold = NULL,
         threshold.col = "grey45",
         threshold.lty = 2,
         threshold.lwd = 1,
         vline = NULL,
         vline.col = "red",
         vline.lty = 2,
         vline.lwd = 1,
         file.memo = NULL,
         file.type = "jpg",
         file.width = 14,
         file.height = 5,
         file.dpi = 300
) {
  # === Argument parsing & Definition =============
  cat(paste0("Rectangular_Manhattan Plotting ", taxa.name, "...\n"))

  den.fold <- 20   # Ratio of Manhattan plot height to SNP density plot height

  # Text objects in figure
  if (is.null(title)) { title <- paste("Manhattan plot of", taxa.name) }
  if (is.null(xlab)) { xlab <- "Chromosome" }
  if (is.null(ylab)) { ylab <- expression(-log[10](italic(p))) }

  # Threshold
  if (0 %in% threshold) { threshold <- threshold[threshold != 0]}
  threshold.n <- length(threshold)
  if (threshold.n > length(threshold.col)) { threshold.col <- rep_len(threshold.col, threshold.n)}
  if (threshold.n > length(threshold.lty)) { threshold.lty <- rep_len(threshold.lty, threshold.n)}
  if (threshold.n > length(threshold.lwd)) { threshold.lwd <- rep_len(threshold.lwd, threshold.n)}

  # Draw Vertical line
  vline.n <- length(vline)
  if (vline.n > length(vline.col)) { vline.col <- rep_len(vline.col, vline.n)}
  if (vline.n > length(vline.lty)) { vline.lty <- rep_len(vline.lty, vline.n)}
  if (vline.n > length(vline.lwd)) { vline.lwd <- rep_len(vline.lwd, vline.n)}

  # Setup output device
  if (!is.null(file.type)) {
    name <- paste("Rectangular-Manhattan", taxa.name, sep = ".")
    if (!is.null(file.memo)) { name <- paste(file.memo, name, sep = ".")}
    switch(
      file.type,
      jpg = jpeg(paste0(name, ".jpg"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi, quality = 100),
      pdf = pdf(paste0(name, ".pdf"), width = file.width, height = file.height),
      tiff = tiff(paste0(name, ".tiff"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi)
    )
  } else {
    if (is.null(dev.list())) { dev.new(width = file.width, height = file.height) }
  }
  par(mar = c(5, 6, 4, 3), xaxs = "i", yaxs = "r", xpd = TRUE)

  # === Data preprocessing ========================
  P.values <- suppressWarnings(as.numeric(Pmap[, 3]))
  is.valid <- !is.na(P.values) & P.values > 0 & P.values <= 1
  Pmap     <- Pmap[is.valid, ]

  # Get digitally encoded chromosome
  chr       <- as.character(Pmap[, 1])
  is.num    <- suppressWarnings(!is.na(as.numeric(chr)))
  chr.max   <- max(as.numeric(chr[is.num]))
  chr.max   <- ifelse(is.infinite(chr.max), 0, chr.max)
  chr.other <- sort(unique(chr[!is.num]))
  chr.n     <- length(unique(chr))

  for (i in 1:length(chr.other)) {
    chr <- replace(chr, chr == chr.other[i], chr.max + i)
  }
  Pmap[, 1] <- as.numeric(chr)
  
  # Order Pmap
  Pmap <- apply(Pmap, 2, as.numeric)
  Pmap <- Pmap[order(Pmap[, 1], Pmap[, 2]), ]
  
  # Convert chromosome position to absolute coordinates
  pos.list  <- tapply(Pmap[, 2], Pmap[, 1], list)
  chr.len   <- sapply(pos.list, max)
  band      <- ceiling((sum(chr.len) / 100))
  Pmap[, 2] <- unlist(lapply(1:length(pos.list), function(x) {
    return(pos.list[[x]] + sum(chr.len[seq_len(x - 1)]) + x * band)
  }))
  
  # -log10
  Pmap[, 3] <- -log10(Pmap[, 3])
  if (!is.null(threshold)) { threshold <- -log10(sort(threshold)) }
  
  # 
  pos.list      <- tapply(Pmap[, 2], Pmap[, 1], list)
  chr.label     <- c(sort(as.numeric(unique(chr[is.num]))), chr.other)
  chr.label.pos <- sapply(pos.list, function(x) {(max(x) + min(x)) / 2})
  
  # vline
  if (!is.null(vline)) { vline.x <- sapply(vline, function(x) { Pmap[x, 2] })}
  
  # Filter points
  is.visable <- filter.points(Pmap[, 2], Pmap[, 3], file.width, file.height, file.dpi, scale = 1)
  Pmap <- Pmap[is.visable, ]
  print(paste(sum(!is.visable), "points has been filter."))
  
  # Get point style
  point.style <- NULL
  style.pch   <- NULL
  style.cex   <- NULL
  style.col   <- NULL
  is.signal   <- rep(FALSE, nrow(Pmap))
  
  # get signal points style
  if (signal && !is.null(signal.col) && !is.null(threshold)) {
      point.style <- classify.points(Pmap[, 3], threshold)
      is.signal   <- apply(point.style, 1, any)
      style.pch   <- rep_len(signal.pch, length(threshold))
      style.cex   <- rep_len(signal.cex, length(threshold))
      style.col   <- rep_len(signal.col, length(threshold)) 
  }
  
  # get normal points style, different chromosomes can have different styles.
  point.style <- cbind(
    point.style,
    sapply(unique(Pmap[, 1]), function(x) {
      i <- rep(FALSE, nrow(Pmap))
      i[(Pmap[, 1] == x) & !is.signal ] <- TRUE
      return(i)
    })
  )
  
  style.pch <- c(style.pch, rep_len(pch, chr.n))
  style.cex <- c(style.cex, rep_len(cex, chr.n))
  style.col <- c(style.col, rep_len(col, chr.n))
  
  # === Drawing ===================================
  # Calculate canvas size
  xlim.max <- max(Pmap[, 2])
  ylim.max <- ceiling(max(Pmap[, 3]))
  if (!is.null(threshold)) { ylim.max <- ceiling(max(ylim.max, threshold)) }
  
  xlim <- c(0, xlim.max)
  ylim <- c(0, ylim.max)
  if (density) {
    xlim <- c(0, 1.01 * xlim.max)
    ylim <- c(-ylim.max / den.fold, ylim.max)
  }
  
  # Draw canvas and axis
  plot(
    NULL,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    cex.axis = cex.axis,
    cex.lab = cex.lab,
    font = 2,
    axes = FALSE,
    main = title
  )
  
  y.at     <- seq(0, ylim.max, ceiling(ylim.max / 10))
  y.legend <- tail(y.at, 1)
  
  axis(1, at = c(0, chr.label.pos), cex.axis = cex.axis, font = 2, labels = c("Chr", chr.label))
  axis(2, at = y.at, cex.axis = cex.axis, font = 2, labels = y.at)
  
  # Draw SNP density
  if (density) {
    # Get density data
    density.list <- Densityplot(
      map = Pmap[, 1:3],
      col = density.col,
      plot = FALSE,
      bin = density.bin.size,
      legend.max = density.bin.max
    )
    # draw chromosome grey background
    for (chr.pos in pos.list) {
      x  <- c(min(chr.pos), min(chr.pos), max(chr.pos), max(chr.pos))
      y0 <- -0.5 * ylim.max / den.fold
      y1 <- -1.5 * ylim.max / den.fold
      y  <- c(y0, y1, y1, y0)
      polygon(x, y, col = "grey", border = "grey")
    }
    segments(
      x0 = Pmap[, 2],
      y0 = y0,
      x1 = Pmap[, 2],
      y1 = y1,
      col = density.list$den.col,
      lwd = 0.1
    )
    legend(
      x = max(Pmap[, 2]) + band,
      y = y.legend,
      title = "",
      legend = density.list$legend.y,
      pch = 15,
      pt.cex = 2.5,
      col = density.list$legend.col,
      cex = 0.8,
      bty = "n",
      y.intersp = 1,
      x.intersp = 1,
      yjust = 1,
      xjust = 0,
      xpd = TRUE
    )
  }
  
  # Draw threshold
  if (!is.null(threshold)) {
    for (thr in 1:length(threshold)) {
      abline(
        h = threshold[thr],
        col = threshold.col[thr],
        lty = threshold.lty[thr],
        lwd = threshold.lwd[thr],
        xpd = FALSE
      )
    }
  }
  
    # Draw points
  for (s in 1:ncol(point.style)) {
    p <- point.style[, s]
    points(
      Pmap[p, 2],
      Pmap[p, 3],
      pch = style.pch[s],
      cex = style.cex[s],
      col = style.col[s]
    )
  }   
    
    # Draw Vertical line
  if (!is.null(vline)) {
    for (v in 1:length(vline)) {
      abline(
        v = vline.x[v],
        col = vline.col[v],
        lty = vline.lty[v],
        lwd = vline.lwd[v],
        xpd = FALSE
      )
    }
  }
    
  if (!is.null(file.type)) { invisible(dev.off()) }
}


#' QQ Plot
#'
#' @param P.val 
#' @param taxa.name 
#' @param pch Either an integer specifying a symbol or a single character to be 
#'    used as the default in plotting points. See \code{\link[graphics]{points}} 
#'    for possible values and their interpretation. Note that only integers and
#'    single-character strings can be set as a graphics parameter (and not NA 
#'    nor NULL).
#' @param cex A numerical value giving the amount by which plotting text and
#'    symbols should be magnified relative to the default. This starts as 1 
#'    when a device is opened, and is reset when the layout is changed, e.g. 
#'    by setting mfrow. see \code{\link[graphics]{par}}. 
#' @param cex.axis The magnification to be used for axis annotation relative
#'    to the current setting of cex.
#' @param cex.lab The magnification to be used for x and y labels relative to
#'    the current setting of cex.
#' @param col A specification for the default plotting color. 
#' @param title the title of the image, if NULL, will be set to "QQ Plot 
#'    of \code{taxa.name}"
#' @param xlab The title of the x axis.
#' @param ylab The title of the y axis.
#' @param box A Boolean value that controls whether to draw a box around
#'    QQplot.
#' @param signal A Boolean value that controls whether to draw points that
#'    exceed the threshold using different styles. Note: this option needs to
#'    be used with \code{threshold}. see bellow.
#' @param signal.pch 
#' @param signal.cex 
#' @param signal.col 
#' @param conf.int A Boolean value that controls whether to draw a confidence
#'    interval.
#' @param conf.int.col The color of confidence interval. The default is "grey".
#' @param threshold 
#' @param threshold.show A Boolean value that controls whether to draw
#'    threshold lines.
#' @param threshold.col 
#' @param threshold.lty 
#' @param threshold.lwd 
#' @param reference 
#' @param reference.col 
#' @param reference.lty 
#' @param reference.lwd 
#' @param file.memo the prefix of the output image file.
#' @param file.type A string or NULL is used to determine the type of output 
#'    file. Can be "jpg", "pdf", "tiff". If it is NULL, it will use 
#'    \code{\link[grDevices:dev]{dev.new()}} to create a new graphics device 
#'    in the current environment, which may be RStudioGD or the default 
#'    device of the system.
#' @param file.width The width of the output file, in inches.
#' @param file.height The height of the output file, in inches.
#' @param file.dpi The resolution of the image, specifying how many pixels 
#'    per inch.
#'
#' @export
#'
#' @examples
#' 
drawQQPlot <-
function(P.val,
         taxa.name,
         pch = 19,
         cex = 0.5,
         cex.axis = 1,
         cex.lab = 1.2,
         col = c("dodgerblue4"),
         title = NULL,
         xlab = NULL,
         ylab = NULL,
         box = FALSE,
         signal = TRUE,
         signal.pch = 19,
         signal.cex = 0.8,
         signal.col = "red",
         conf.int = TRUE,
         conf.int.col = "grey",
         threshold = NULL,
         threshold.show = FALSE,
         threshold.col = "grey45",
         threshold.lty = 2,
         threshold.lwd = 1,
         reference = TRUE,
         reference.col = "red",
         reference.lty = 1,
         reference.lwd = 2,
         file.memo = NULL,
         file.type = "jpg",
         file.width = 5,
         file.height = 5,
         file.dpi = 300
) {
  # === Argument parsing & Definition =============
  cat(paste0("Q_Q Plotting ", taxa.name, "...\n"))
  
  # Text objects in figure
  if (is.null(title)) { title <- paste("QQplot of", taxa.name) }
  if (is.null(xlab)) { xlab <- expression(Expected~~-log[10](italic(p))) }
  if (is.null(ylab)) { ylab <- expression(Observed~~-log[10](italic(p))) }
  
  # Threshold
  if (0 %in% threshold) { threshold <- threshold[threshold != 0]}
  threshold.n <- length(threshold)
  if (threshold.n > length(threshold.col)) { threshold.col <- rep_len(threshold.col, threshold.n)}
  if (threshold.n > length(threshold.lty)) { threshold.lty <- rep_len(threshold.lty, threshold.n)}
  if (threshold.n > length(threshold.lwd)) { threshold.lwd <- rep_len(threshold.lwd, threshold.n)}
  
  # setup output device
  if (!is.null(file.type)) {
    name <- paste("QQplot", taxa.name, sep = ".")
    if (!is.null(file.memo)) { name <- paste(file.memo ,  name, sep = ".") }
    switch(
      file.type,
      jpg = jpeg(paste0(name, ".jpg"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi, quality = 100),
      pdf = pdf(paste0(name, ".pdf"), width = file.width, height = file.height),
      tiff = tiff(paste0(name, ".tiff"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi)
    )
  } else {
    if (is.null(dev.list())) { dev.new(width = file.width, height = file.height) }
  }
  par(mar = c(5, 5, 4, 2), xpd = TRUE)
  
  # === Data preprocessing ========================
  # filter & order P-values
  P.val     <- suppressWarnings(as.numeric(P.val))
  is.valid  <- !is.na(P.val) & P.val > 0 & P.val <= 1
  P.val     <- P.val[is.valid]
  P.val     <- -log10(P.val[order(P.val)])
  P.val.n   <- length(P.val)
  quantiles <- -log10((1:P.val.n) / (P.val.n + 1))
  
  if (!is.null(threshold)) { threshold <- -log10(sort(threshold)) }
  
  # filter
  is.visable <- filter.points(quantiles, P.val, file.width, file.height, file.dpi)
  
  # calculate the confidence interval of QQ-plot
  if (conf.int) {
    xi  <- ceiling((10 ^ -quantiles[is.visable]) * P.val.n)
    c95 <- -log10(qbeta(0.95, xi, P.val.n - xi + 1))
    c05 <- -log10(qbeta(0.05, xi, P.val.n - xi + 1))
  }
  
  # Get point style
  point.style <- NULL
  style.pch   <- NULL
  style.cex   <- NULL
  style.col   <- NULL
  is.signal   <- rep(FALSE, length(P.val))
  
  # get signal points style
  if (signal && !is.null(signal.col) && !is.null(threshold)) {
    point.style <- classify.points(P.val, threshold)
    is.signal   <- apply(point.style, 1, any)
    style.pch   <- rep_len(signal.pch, length(threshold))
    style.cex   <- rep_len(signal.cex, length(threshold))
    style.col   <- rep_len(signal.col, length(threshold)) 
  }
  
  # get normal points style
  point.style <- cbind(point.style, !is.signal)
  style.pch <- c(style.pch, pch)
  style.cex <- c(style.cex, cex)
  style.col <- c(style.col, col)
  
  # === Drawing ===================================
  # Calculate canvas size
  xlim.max  <- ceiling(max(quantiles))
  ylim.max <- ceiling(max(P.val))
  if (conf.int) {
    ylim.max <- ceiling(max(ylim.max, c95, c05))
  }
  if (!is.null(threshold) && threshold.show) {
    ylim.max <- ceiling(max(ylim.max, threshold))
  }
  
  # Draw canvas and axis
  plot(
    NULL,
    xlim = c(0, xlim.max),
    ylim = c(0, ylim.max),
    axes = FALSE,
    cex.axis = cex.axis,
    cex.lab = cex.lab,
    xlab = xlab,
    ylab = ylab,
    main = title
  )
  x.at <- seq(0, xlim.max, ceiling(xlim.max / 10))
  y.at <- seq(0, ylim.max, ceiling(ylim.max / 10))
  axis(1, at = x.at, labels = x.at, cex.axis = cex.axis)
  axis(2, at = y.at, labels = y.at, cex.axis = cex.axis)
  
  # Draw the confidence interval of QQ-plot
  if (conf.int) {
    x = quantiles[is.visable]
    polygon(
      x = c(rev(x), x),
      y = c(rev(c05), c95),
      col = conf.int.col,
      border = conf.int.col
    )
  }
  
  # Draw threshold
  if (!is.null(threshold) && threshold.show) {
    for (thr in 1:length(threshold)) {
      abline(
        h = threshold[thr],
        col = threshold.col[thr],
        lty = threshold.lty[thr],
        lwd = threshold.lwd[thr],
        xpd = FALSE
      )
    }
  }
  
  # Draw reference line (y = x)
  if (reference) {
    lines(
      x = c(-0.1, max(quantiles)),
      y = c(-0.1, max(quantiles)),
      col = reference.col,
      lty = reference.lty,
      lwd = reference.lwd
    )
  }
  
  # Draw points
  for (s in 1:ncol(point.style)) {
    p <- point.style[, s]
    points(
      quantiles[p & is.visable],
      P.val[p & is.visable],
      pch = style.pch[s],
      cex = style.cex[s],
      col = style.col[s]
    )
  }
  
  # Draw box
  if (box) {
    box()
  }
  
  # Close device
  if (!is.null(file.type)) {
    invisible(dev.off())
  }
}


Densityplot <- function(map, col = c("darkgreen", "yellow", "red"), main = "SNP Density", bin = 1e6,
                       band = 3, width = 5, legend.len = 10, legend.max = NULL, legend.pt.cex = 3,
                       legend.cex = 1, legend.x.intersp = 1, legend.y.intersp = 1, plot = TRUE) {
    
    ## Step 1: preprocess map
    # filter map
    map <- as.matrix(map)
    map <- map[!is.na(map[, 2]), ]
    map <- map[!is.na(map[, 3]), ]
    map <- map[map[, 2] != 0, ]

    # get the number of chromosomes
    max.chr <- max(as.numeric(map[, 2]), na.rm = TRUE)
    if (is.infinite(max.chr)) { max.chr <- 0 }
    
    # deal with x,y
    map.xy.index <- suppressWarnings(which(!as.numeric(map[, 2]) %in% c(0:max.chr)))
    if (length(map.xy.index) != 0) {
        chr.xy <- unique(map[map.xy.index, 2])
        for (i in 1:length(chr.xy)) {
            map[map[, 2] == chr.xy[i], 2] <- max.chr + i
        }
    }
    
    # sort map
    map <- map[order(as.numeric(map[, 2]), as.numeric(map[, 3])), ]
    chr <- as.numeric(map[, 2])
    pos <- as.numeric(map[, 3])
    chr.num <- unique(chr)
    chorm.maxlen <- max(pos)
    
    # Step 2: count SNP
    pos.x      <- list()
    col.index  <- list()
    maxbin.num <- NULL
    for (i in 1:length(chr.num)) {
        pos.x[[i]] <- pos[which(chr == chr.num[i])]
        cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
        if (cut.len <= 1) {
            col.index[[i]] = 1
        } else {
            cut.r          <- cut(pos.x[[i]], cut.len, labels = FALSE)
            eachbin.num    <- table(cut.r)
            maxbin.num     <- c(maxbin.num, max(eachbin.num))
            col.index[[i]] <- rep(eachbin.num, eachbin.num)
        }
    }
    Maxbin.num <- max(maxbin.num)
    maxbin.num <- Maxbin.num
    if (!is.null(legend.max)) {
        maxbin.num <- legend.max
    }
    col = colorRampPalette(col)(maxbin.num)
    col.seg = NULL
    for (i in 1:length(chr.num)) {
        if (!is.null(legend.max)) {
            if (legend.max < Maxbin.num) {
                col.index[[i]][col.index[[i]] > legend.max] <- legend.max
            }
        }
        col.seg <- c(col.seg, col[round(col.index[[i]] * length(col) / maxbin.num)])
    }
    if (length(map.xy.index) != 0) {
        for (i in 1:length(chr.xy)) {
            chr.num[chr.num == max.chr + i] <- chr.xy[i]
        }
    }
    chr.num <- rev(chr.num)

    # image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
    # round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
    # t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
    
    ## Step 3: Deal with legend label and color.
    len      <- round(seq(0, maxbin.num, length = legend.len))[2]
    legend.y <- seq(0, maxbin.num, len)

    if (!is.null(legend.max)) {
        if (legend.max < Maxbin.num) {
            if (!maxbin.num %in% legend.y) {
                legend.y <- c(legend.y, paste0(">=", maxbin.num))
            } else {
                legend.y[length(legend.y)] <- paste0(">=", maxbin.num)
            }
            legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
        } else {
            if (!maxbin.num %in% legend.y) { 
                legend.y <- c(legend.y, maxbin.num) 
            }
            legend.y.col <- c(legend.y[-1])
        }
    } else {
        if (!maxbin.num %in% legend.y) {
            legend.y     <- c(legend.y, paste0(">", max(legend.y)))
            legend.y.col <- c(legend.y[c(-1, -length(legend.y))], maxbin.num)
        } else {
            legend.y.col <- c(legend.y[-1])
        }
    }
    
    legend.y.col <- as.numeric(legend.y.col)
    legend.col   <- c("grey", col[round(legend.y.col * length(col) / maxbin.num)])
    
    ## Step 4: Plot or Return
    if (plot) {
        # draw a canvas
        plot(
            NULL,
            xlim = c(0, chorm.maxlen + chorm.maxlen / 10),
            ylim = c(0, length(chr.num) * band + band),
            main = main,
            axes = FALSE,
            xlab = "",
            ylab = "",
            xaxs = "i",
            yaxs = "i"
        )
        
        # draw each Chr
        for (i in 1:length(chr.num)) {
            # draw bg
            polygon(
                x = c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])),
                y = c(
                    -width / 5 - band * (i - length(chr.num) - 1),
                     width / 5 - band * (i - length(chr.num) - 1),
                     width / 5 - band * (i - length(chr.num) - 1),
                    -width / 5 - band * (i - length(chr.num) - 1)
                ),
                col = "grey",
                border = "grey"
            )
            # draw fg
            segments(
                x0 = pos.x[[i]],
                y0 = -width / 5 - band * (i - length(chr.num) - 1),
                x1 = pos.x[[i]],
                y1 = width / 5 - band * (i - length(chr.num) - 1),
                col = col[round(col.index[[i]] * length(col) / maxbin.num)],
                lwd = 1
            )
        }
        # draw Chr label
        mtext(
            at = seq(band, length(chr.num) * band, band),
            text = paste0("Chr", chr.num),
            side = 2,
            las = 2,
            font = 1,
            cex = 0.6,
            line = 0.2
        )
        
        # draw physical distance ruler
        axis(
            side = 3,
            at = seq(0, chorm.maxlen, length = 10),
            labels = c(NA, paste0(round((
                seq(0, chorm.maxlen, length = 10)
            )[-1] / 1e6, 0), "Mb")),
            font = 1,
            cex.axis = 0.8,
            tck = 0.01,
            lwd = 2,
            padj = 1.2
        )
        
        # draw legend
        legend(
            x = (chorm.maxlen + chorm.maxlen / 100),
            y = (-width / 2.5 - band * (length(chr.num) - length(chr.num) - 1)),
            legend = legend.y,
            col = legend.col,
            pch = 15,
            bty = "n",
            cex = legend.cex,
            pt.cex = legend.pt.cex,
            xjust = 0,
            yjust = 0,
            x.intersp = legend.x.intersp,
            y.intersp = legend.y.intersp,
            title = "",
            xpd = TRUE
        )
    } else {
        return(list(den.col = col.seg, legend.col = legend.col, legend.y = legend.y))
    }
}

#' SNP Density
#'
#' @param Pmap P value Map
#' @param col The color vector
#' @param dpi Number. Dots per inch for .jpg and .tiff files
#' @param bin.size the window size for counting SNP number
#' @param bin.max maximum SNP number, for winows, which has more SNPs than bin.max, will be painted in same color
#' @param memo Character. A text marker on output files
#' @param outpath Only when file.output = TRUE, determines the path of the output file
#' @param file.type format of output figure
#' @param file.output Whether to output the file
#' @param verbose whether to print detail.
#' 
#' @export
#' @return 
#' Output file:
#' <memo>.SNP_Density.<type>
#'
#' @examples
#' data(pig60K, package = "rMVP")
#' 
#' MVP.Report.Density(pig60K, file.output=FALSE)
#' 
drawDensityPlot <- function(Pmap, col = c("darkgreen", "yellow", "red"), dpi = 300, outpath = getwd(), memo = 'MVP',
                            bin.size = 1e6, bin.max = NULL, file.type = "jpg", file.output = TRUE, verbose = TRUE) {
    cat("SNP_Density Plotting", "\n", verbose = verbose)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    
    w <- 9
    h <- 7
    if (file.output) {
        name <- file.path(outpath, paste0(memo, ".SNP_Density"))
        switch(file.type,
               jpg = jpeg(paste0(name, ".jpg"), width = w * dpi, height = h * dpi, res = dpi, quality = 100),
               pdf = pdf(paste0(name, ".pdf"), width = w, height = h),
               tiff = tiff(paste0(name, ".tiff"), width = w * dpi, height = h * dpi, res = dpi)
        )
        par(xpd = TRUE)
    } else {
        if (is.null(dev.list())) { dev.new(width = w, height = h) }
        par(xpd = TRUE)
    }
    
        Densityplot(map = Pmap[,c(1:3)], col = col, bin = bin.size, legend.max = bin.max, main = paste0("The number of SNPs within ", bin.size / 1e6, "Mb window size"))
    
    if (file.output) { dev.off() }
}

#' Phenotype Histogram
#'
#' @description This function is to draw a single phenotype into a histogram.
#' 
#' @param phe Phenotype data vector, any value that cannot be converted to a 
#'    number will be ignored as NA.
#' @param taxa.name The identifier of the phenotype will be used to generate a
#'    portion of the image file name. If the title parameter is NULL, it will 
#'    also be part of the title.
#' @param break.n A single number giving the number of cells for the histogram.
#'    The default value is 15.
#' @param col The color vector of the histogram. If the number of colors is 
#'    less than break.n, the color will be reused. If the number of colors is 
#'    greater than break.n, only the previous break.n colors will be used.
#' @param title the title of the image, if NULL, will be set to "Distribution 
#'    of \code{taxa.name}"
#' @param xlab The title of the x axis.
#' @param ylab The title of the y axis.
#' @param statistics.show A Boolean value that determines whether to display a 
#'    series of statistics for the phenotype, including mean, standard 
#'    deviation, and P value for Shapiro-Wilk Normality tests or 
#'    Kolmogorov-Smirnov tests. Shapiro-Wilk Normality tests are used when the
#'    sample size is less than 5000, and Kolmogorov-Smirnov tests are used when
#'    the sample size is greater than 5000.
#' @param file.memo A string that is the prefix of the output image file.
#' @param file.type A string or NULL is used to determine the type of output 
#'    file. Can be "jpg", "pdf", "tiff". If it is NULL, it will use 
#'    \code{\link[grDevices:dev]{dev.new()}} to create a new graphics device 
#'    in the current environment, which may be RStudioGD or the default 
#'    device of the system.
#' @param file.width The width of the output file, in inches.
#' @param file.height The height of the output file, in inches.
#' @param file.dpi The resolution of the image, specifying how many pixels 
#'    per inch.
#' 
#' @export
#' 
#' @examples
#' # read the example phenotype file
#' phenoPath <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")
#' phe <- read.table(file = phenoPath, header = TRUE, na.strings = "-999")
#' 
#' # draw a single phenotype histogram
#' drawPhenotypeHistogram(phe[, 2], colnames(phe)[2])
drawPhenotypeHistogram <-
function(phe,
         taxa.name,
         break.n = 15,
         col = c("dodgerblue4", "olivedrab4", "violetred", "darkgoldenrod1", "purple4"),
         title = NULL,
         xlab = NULL,
         ylab = NULL,
         statistics.show = TRUE,
         file.memo = NULL,
         file.type = "jpg",
         file.width = 5,
         file.height = 5,
         file.dpi = 300
) {
  # === Argument parsing & Definition =============
  cat(paste0("Phenotype distribution plotting...\n"))
  
  if (is.null(taxa.name)) { taxa.name <- colnames(phe) }
  
  # Text objects in figure
  if (is.null(title)) { title <- paste("Distribution of", taxa.name) }
  if (is.null(xlab)) { xlab <- "" }
  if (is.null(ylab)) { ylab <- "Density" }
  
  if (!is.null(file.type)) {
    name <- paste("Phe_Distribution", taxa.name, sep = ".")
    if (!is.null(file.memo)) { name <- paste(file.memo ,  name, sep = ".") }
    switch(
      file.type,
      jpg = jpeg(paste0(name, ".jpg"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi, quality = 100),
      pdf = pdf(paste0(name, ".pdf"), width = file.width, height = file.height),
      tiff = tiff(paste0(name, ".tiff"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi)
    )
  } else {
    if (is.null(dev.list())) { dev.new(width = file.width, height = file.height) }
  }
  par(mar = c(5, 5, 4, 2), xpd = TRUE)
  
  # === Data preprocessing ========================
  phe    <- suppressWarnings(as.numeric(phe))
  phe    <- phe[!is.na(phe)]
  breaks <- seq(min(phe), max(phe), length = break.n)
  col    <- rep_len(col, break.n)
  
  # hist
  xx <- hist(phe, plot = FALSE, breaks = breaks)
  ylim.max <- max(xx$density, density(phe)$y) * 1.25  # Leave space for the text
  
  # === Drawing ===================================
  hist(
    x = phe,
    breaks = breaks,
    xlab = xlab,
    ylab = ylab,
    ylim = c(0, ylim.max),
    freq = FALSE,
    col = col,
    font = 2,
    font.lab = 2,
    main = title
  )
  
  lines(density(phe), lwd = 2, xpd = FALSE)
  
  if (statistics.show) {
    # test normal distribution
    if (length(phe) <= 5000) {
      norm.p <- round(shapiro.test(phe)$p, 4)
      test.method <- "Shapiro-Wilk"
    } else {
      norm.p <- round(ks.test(phe,"pnorm")$p, 4)
      test.method <- "Kolmogorov-Smirnov"
    }
    
    # draw text
    text(xx$breaks[1], y = ylim.max * 0.975, labels = paste0("Mean: ", round(mean(phe), 2)), font = 2, adj = 0)
    text(xx$breaks[1], y = ylim.max * 0.92, labels = paste0("Sd: ", round(sd(phe), 2)), font = 2, adj = 0)
    text(xx$breaks[1], y = ylim.max * 0.865, labels = paste0(test.method, ": ", norm.p), font = 2, adj = 0)
  }
  # close device
  if (!is.null(file.type)) { invisible(dev.off()) }
}


#' @rdname drawPhenotypeHistogram
#' 
#' @param phe a phenotype data frame with colnames
#' @param ... A list of name-value pairs. Will be passed to 
#'    drawPhenotypeHistogram.
#' 
#' @export
#' 
#' @examples 
#' # Convenient visualization function to draw each phenotype 
#' # histogram separately
#' drawMultiplePhenotypeHistograms(phe)
drawMultiplePhenotypeHistograms <- function(phe, ...) {
  for (i in seq_len(ncol(phe))) {
    drawPhenotypeHistogram(phe[, i], colnames(phe)[i], ...)
  }
}


drawPCAplot <- 
function(PCA,
         class = NULL,
         cluster.n = 3,
         pch = 19,
         cex = 0.5,
         cex.axis = 1,
         cex.lab = 1.2,
         col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1", "red"),
         title = NULL,
         xlab = NULL,
         ylab = NULL,
         box = FALSE,
         legend.pos = "topright",
         file.memo = NULL,
         file.type = "jpg",
         file.width = 5,
         file.height = 5,
         file.dpi = 300
) {
  # === Argument parsing & Definition =============
  cat(paste0("PCA plotting...\n"))
  
  # Text objects in figure
  if (is.null(title)) { title <- paste("PCA Plot") }
  if (is.null(xlab)) { xlab <- "PC1" }
  if (is.null(ylab)) { ylab <- "PC2" }
  
  # setup output device
  if (!is.null(file.type)) {
    name <- "PCA"
    if (!is.null(file.memo)) { name <- paste(file.memo ,  name, sep = ".") }
    switch(
      file.type,
      jpg = jpeg(paste0(name, ".jpg"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi, quality = 100),
      pdf = pdf(paste0(name, ".pdf"), width = file.width, height = file.height),
      tiff = tiff(paste0(name, ".tiff"), width = file.width * file.dpi, height = file.height * file.dpi, res = file.dpi)
    )
  } else {
    if (is.null(dev.list())) { dev.new(width = file.width, height = file.height) }
  }
  par(mar = c(5, 5, 4, 2), xpd = TRUE)
  
  if (!is.null(class)) {
    if (length(class) != nrow(PCA)) {
      stop("the length of 'class' differs from the row of 'PCA'.")
    } 
    cluster.n <- length(unique(class))
  }
  
  col <- rep_len(col, cluster.n)
  pch <- rep_len(pch, cluster.n)

  # === Data preprocessing ========================
  if (is.null(class)) {
    kc <- kmeans(PCA[, c(1, 2)], cluster.n)
  } else {
    kc <- list(cluster = as.numeric(as.factor(class)))
  }
  
  # === Drawing ===================================
  plot(
    PCA[, 1],
    PCA[, 2],
    pch = pch[kc$cluster],
    col = col[kc$cluster],
    cex = cex,
    cex.axis = cex.axis,
    cex.lab = cex.lab,
    font = 2,
    font.lab = 2,
    main = title,
    xlab = xlab,
    ylab = xlab,
    axes = FALSE
  )
  axis(1, cex.axis = cex.axis)
  axis(2, cex.axis = cex.axis)
  
  if (!is.null(class)) {
    legend(legend.pos, levels(as.factor(class)), col = col, pch = pch, box.col = NA)
  }
  
  # Draw box
  if (box) { box() }
  
  # Close device
  if (!is.null(file.type)) { invisible(dev.off()) }
}

drawLDdecayPlot <- function() {
  # pass
}