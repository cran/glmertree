node_glmertree <- function (mobobj, which = NULL, id = TRUE, pop = TRUE, pointcol = "black", 
                            pointcex = 0.5, boxcol = "black", boxwidth = 0.5, boxfill = "lightgray", 
                            bg = "white", fitmean = TRUE, linecol = "red", 
                            cdplot = FALSE, fivenum = TRUE, breaks = NULL, ylines = NULL, 
                            xlab = FALSE, ylab = FALSE, margins = rep(1.5, 4), mainlab = NULL,
                            ranef = "fixed", fixef = "fixed", plot.obs = TRUE, 
                            observed = TRUE, ...) {
  
  mf <- model.frame.modelparty(mobobj)  
  y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, rhs = 0L)
  if (isTRUE(ylab)) ylab <- names(y)
  if (identical(ylab, FALSE)) ylab <- ""
  if (is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 0, 2)
  y <- y[[1L]]
  X <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, rhs = 1L)
  fitted <- mobobj$fitted[["(fitted)"]] # this identifies node membership
  if (ncol(X) == 0) {
    rval <- switch(class(y)[1L], 
                   Surv = node_surv(mobobj, id = id, mainlab = mainlab, ...), 
                   factor = node_barplot(mobobj, id = id, 
                                         mainlab = mainlab, ...), 
                   ordered = node_barplot(mobobj, id = id, 
                                          mainlab = mainlab, ...), 
                   node_boxplot(mobobj, ...))
    return(rval)
  }
  if (is.factor(y)) y <- factor(y, levels = rev(levels(y)))
  if (is.null(which)) which <- 1L:NCOL(X)
  X <- X[, which, drop = FALSE]
  k <- NCOL(X)
  xlab <- if (!identical(xlab, FALSE)) {
    if (isTRUE(xlab)) colnames(X) else rep(xlab, length.out = k)
  } else rep("", k)
  if (is.factor(y)) {
    if (!requireNamespace("vcd")) 
      stop(sprintf("Package %s is required for spine/CD plots", sQuote("vcd")))
    if (cdplot) {
      num_fun <- function(x, y, yfit, i, name, ...) {
        vcd::cd_plot(x, y, xlab = xlab[i], ylab = ylab, 
                     name = name, newpage = FALSE, margins = margins, 
                     pop = FALSE, ...)
        if (fitmean) {
          grid.lines(x, yfit, default.units = "native", 
                     gp = gpar(col = linecol))
          if (pop) popViewport() else upViewport()
        } else {
          if (pop) popViewport() else upViewport()
        }
      }
    } else {
      xscale <- if (is.null(breaks)) {
        if (fivenum) {
          lapply(X, function(z) {
            if (is.factor(z)) 1 else fivenum(z)
          })
        } else {
          lapply(X, function(z) {
            if (is.factor(z)) 1 else hist(z, plot = FALSE)$breaks
          })
        }
      } else {
        if (is.list(breaks)) breaks else list(breaks)
      }
      num_fun <- function(x, y, yfit, i, name, ...) {
        vcd::spine(x, y, xlab = xlab[i], ylab = ylab, 
                   name = name, newpage = FALSE, margins = margins, 
                   pop = FALSE, breaks = xscale[[i]], ...)
        if (fitmean) {
          xaux <- cut(x, breaks = xscale[[i]], include.lowest = TRUE)
          yfit <- unlist(tapply(yfit, xaux, mean))
          xaux <- prop.table(table(xaux))
          xaux <- cumsum(xaux) - xaux/2
          grid.lines(xaux, yfit, default.units = "native", 
                     gp = gpar(col = linecol))
          grid.points(xaux, yfit, default.units = "native", 
                      gp = gpar(col = linecol, cex = pointcex), 
                      pch = 19)
          if (pop) popViewport() else upViewport()
        } else {
          if (pop) popViewport() else upViewport()
        }
      }
    }
    cat_fun <- function(x, y, yfit, i, name, ...) {
      vcd::spine(x, y, xlab = xlab[i], ylab = ylab, name = name, 
                 newpage = FALSE, margins = margins, pop = FALSE, ...)
      if (fitmean) {
        yfit <- unlist(tapply(yfit, x, mean))
        xaux <- prop.table(table(x))
        xaux <- cumsum(xaux + 0.02) - xaux/2 - 0.02
        grid.lines(xaux, yfit, default.units = "native", 
                   gp = gpar(col = linecol))
        grid.points(xaux, yfit, default.units = "native", 
                    gp = gpar(col = linecol, cex = pointcex), pch = 19)
        if (pop) popViewport() else upViewport()
      } else {
        if (pop) popViewport() else upViewport()
      }
    }
  } else {
    xscale <- sapply(X, function(z) {
      if (is.factor(z)) c(1, length(levels(z))) else range(z)
    })
    yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
    num_fun <- function(x, y, yfit, i, name, ...) {
      xscale[, i] <- xscale[, i] + c(-0.1, 0.1) * diff(xscale[, i])
      pushViewport(plotViewport(margins = margins, name = name, 
                                yscale = yscale, xscale = xscale[, i]))
      if (observed) {
        grid.points(x, y, gp = gpar(col = pointcol, cex = pointcex))
      }
      if (fitmean) {
        grid.lines(x, yfit, default.units = "native", gp = gpar(col = linecol))
      }
      grid.xaxis(at = c(ceiling(xscale[1L, i] * 10), 
                        floor(xscale[2L, i] * 10))/10)
      grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
      grid.rect(gp = gpar(fill = "transparent"))
      if (ylab != "")
        grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, "lines"), rot = 90)
      if (xlab[i] != "") 
        grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, "lines"))
      if (pop) popViewport() else upViewport()
    }
    cat_fun <- function(x, y, yfit, i, name, ...) {
      xlev <- levels(x)
      pushViewport(plotViewport(margins = margins, name = name, yscale = yscale, 
                                xscale = c(0.3, xscale[2L, i] + 0.7)))
      for (j in seq_along(xlev)) {
        by <- boxplot(y[x == xlev[j]], plot = FALSE)
        if (observed) {
          xl <- j - boxwidth/4
          xr <- j + boxwidth/4
          grid.lines(unit(c(xl, xr), "native"), 
                     unit(by$stats[1L], "native"), gp = gpar(col = boxcol))
          grid.lines(unit(j, "native"), 
                     unit(by$stats[1L:2L], "native"), gp = gpar(col = boxcol, lty = 2))
          grid.rect(unit(j, "native"), 
                    unit(by$stats[2L], "native"), width = unit(boxwidth, "native"), 
                    height = unit(diff(by$stats[2:3]), "native"), 
                    just = c("center", "bottom"), gp = gpar(col = boxcol, 
                                                            fill = boxfill))
          grid.rect(unit(j, "native"), unit(by$stats[3L], "native"), 
                    width = unit(boxwidth, "native"), 
                    height = unit(diff(by$stats[3L:4L]), "native"), 
                    just = c("center", "bottom"), gp = gpar(col = boxcol, 
                                                            fill = boxfill))
          grid.lines(unit(j, "native"), unit(by$stats[4L:5L], "native"), 
                     gp = gpar(col = boxcol, lty = 2))
          grid.lines(unit(c(xl, xr), "native"), unit(by$stats[5L], "native"), 
                     gp = gpar(col = boxcol))
          n <- length(by$out)
          if (n > 0L) {
            grid.points(unit(rep.int(j, n), "native"), unit(by$out, "native"), 
                        size = unit(0.5, "char"), gp = gpar(col = boxcol))
          }
        }
      }
      if (fitmean) {
        yfit <- unlist(tapply(yfit, x, mean))
        grid.lines(seq_along(xlev), yfit, default.units = "native", 
                   gp = gpar(col = linecol))
        grid.points(seq_along(xlev), yfit, default.units = "native", 
                    gp = gpar(col = linecol, cex = pointcex), pch = 19)
      }
      grid.rect(gp = gpar(fill = "transparent"))
      grid.xaxis(at = 1L:length(xlev), label = xlev)
      grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
      if (ylab != "") 
        grid.text(ylab, y = unit(0.5, "npc"), x = unit(-3, "lines"), rot = 90)
      if (xlab[i] != "") 
        grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, "lines"))
      if (pop) popViewport() else upViewport()
    }
  }
  rval <- function(node) {
    nid <- id_node(node)
    ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)
    y <- y[ix]
    top_vp <- viewport(
      layout = grid.layout(nrow = k, ncol = 2, 
                           widths = unit(c(ylines, 1), c("lines", "null")), 
                           heights = unit(k, "null")), width = unit(1, "npc"), 
      height = unit(1, "npc") - unit(2, "lines"), 
      name = paste("node_mob", nid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if (id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, info_node(node)$nobs)
    }
    grid.text(mainlab, y = unit(1, "npc") - unit(0.75, "lines"))
    popViewport()
    for (i in 1L:k) {
      xi <- X[ix, i]
      o <- order(xi)
      yi <- y[o]
      xi <- xi[o]
      if (fixef == "constant" && k > 1L) {
        yfit <- node$info$fitted[, names(X)[i]][o]
      } else {
        yfit <- node$info$fitted[o]
      }
      plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
      pushViewport(plot_vpi)
      if (is.factor(xi)) {
        cat_fun(xi, yi, yfit, i, paste("node_mob", nid, "-", i, sep = ""), ...)
      } else {
        num_fun(xi, yi, yfit, i, paste("node_mob", nid, "-", i, sep = ""), ...)
      }
      if (pop) popViewport() else upViewport()
    }
    if (pop) popViewport() else upViewport()
  }
  return(rval)
}
class(node_glmertree) <- "grapcon_generator"


## TODO: Do we need this one? Or can we just use model.frame.glmertree above?
model.frame.modelparty <- function (formula, ...) {
  mf <- formula$data
  if (nrow(mf) > 0L) return(mf)
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- formula$info$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), 
                       names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf[names(nargs)] <- nargs
  if (is.null(env <- environment(formula$info$terms))) env <- parent.frame()
  mf$formula <- Formula::Formula(as.formula(mf$formula))
  eval(mf, env)
}


simple_terminal_func <- function(node, minlength = 15) {
  c("Fixed effects:", 
    paste(abbreviate(names(node$coefficients), minlength = minlength), " ",
          formatC(node$coefficient, digits = 3, format = "f")))
}






node_terminal_glmertree <- function(obj, digits = 3, abbreviate = FALSE,
                                    fill = c("lightgray", "white"), id = TRUE,
                                    just = c("center", "top"),
                                    top = 0.85, align = c("center", "left", "right"),
                                    gp = NULL, FUN = NULL, height = NULL, width = NULL)
{
  nam <- names(obj)
  extract_label <- function(node) formatinfo_node(node, FUN = FUN,
                                                  default = c("terminal", "node"))
  maxstr <- function(node) {
    lab <- extract_label(node)
    klab <- if (is.terminal(node)) {
        ""
      } else {
        unlist(lapply(kids_node(node), maxstr))
      }
    lab <- c(lab, klab)
    lab <- try(unlist(lapply(lab, function(x) strsplit(x, "\n"))), silent = TRUE)
    if (inherits(lab, "try-error")) {
      paste(rep("a", 9L), collapse = "")
    } else {
      return(lab[which.max(nchar(lab))])
    }
  }
  nstr <- if (is.null(width)) {
    maxstr(node_party(obj))
  } else {
    paste(rep("a", width), collapse = "")
  }
  just <- match.arg(just[1L], c("center", "centre", "top"))
  if (just == "centre") {
    just <- "center"
  }
  align <- match.arg(align[1L], c("center", "centre", "left", "right"))
  if (align == "centre") {
    align <- "center"
  }
  rval <- function(node) {
    fill <- rep(fill, length.out = 2)
    lab <- extract_label(node)
    if (!is.null(gp)) {
      outer_vp <- viewport(gp = gp)
      pushViewport(outer_vp)
    }
    if (is.null(height)) {
      height <- length(lab) + 1L
    }
    node_vp <- viewport(x = unit(0.5, "npc"),
                        y = unit(if (just == "top") top else 0.5, "npc"),
                        just = c("center", just),
                        width = unit(1, "strwidth", nstr) * 1.1,
                        height = unit(height, "lines"),
                        name = paste("node_terminal", id_node(node), sep = ""),
                        gp = if (is.null(gp)) gpar() else gp)
    pushViewport(node_vp)
    grid.rect(gp = gpar(fill = fill[1]))
    for (i in seq_along(lab)) {
      grid.text(x = switch(
      align, center = unit(0.5, "npc"), left = unit(1, "strwidth", "a"),
      right = unit(1, "npc") - unit(1, "strwidth", "a")),
      y = unit(length(lab) - i + 1, "lines"), lab[i], just = align)
    }
    if (id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"),
                           y = unit(1.25, "npc"),
                           width = max(unit(1, "lines"),
                                       unit(1.3, "strwidth", nam[id_node(node)])),
                           height = max(unit(1, "lines"),
                                        unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      nid <- id_node(node)
      mainlab <- function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      mainlab <- mainlab(nid, info_node(node)$nobs)
      grid.text(mainlab, y = unit(1, "npc") - unit(0.75, "lines"))
      popViewport()
    }
    if (is.null(gp)) upViewport() else upViewport(2)
  }
  return(rval)
}
class(node_terminal_glmertree) <- "grapcon_generator"







node_growthplot <- function(mobobj, cluster = NULL,
                            id = TRUE, pop = TRUE, bg = "white",
                            curvecol = gray(0.2, alpha = 1/6), curvelwd = 2,
                            linecol = 2, linelwd = 3,
                            ylines = NULL, xlab = FALSE, ylab = FALSE,
                            mainlab = NULL, observed = TRUE, fitmean = TRUE,
                            xscale = NULL, xaxis.at = NULL, xaxis.labs = NULL,
                            gp = gpar(...), ...)
{
  ## obtain dependent variable
  mf <- model.frame(mobobj)
  y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, rhs = 0L)
  if(isTRUE(ylab)) ylab <- names(y)
  if(identical(ylab, FALSE)) ylab <- ""
  if(is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 3, 5)
  y <- y[[1L]]
  
  ## obtain explanatory variables
  X <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, rhs = 1L)
  X <- X[,NCOL(X),drop=FALSE]
  xlab <- if(!identical(xlab, FALSE)) { if(isTRUE(xlab)) colnames(X) else xlab[1L] } else ""
  
  ## axis scaling
  if (is.null(xscale)) xscale <- range(X) + c(-0.1, 0.1) * diff(range(X))
  yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
  
  ## other information
  fitted <- mobobj$fitted[["(fitted)"]]
  cf <- coef(mobobj)
  
  rval <- function(node) {
    
    ## node index
    nid <- id_node(node)
    ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)
    
    ## subsets
    y <- y[ix]
    x <- X[ix, 1L]
    cf <- cf[as.character(nid), ]
    
    ## set up top viewport
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_growthplot", nid, sep = ""),
                       gp = gp)
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    
    ## main title
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    
    if (is.null(mainlab)) { 
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, info_node(node)$nobs)
    }
    grid.text(mainlab)
    popViewport()
    
    plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                     xscale = xscale, yscale = yscale,
                     name = paste0("node_growthplot", nid, "plot"),
                     clip = FALSE)
    pushViewport(plot)
    
    if (is.null(xaxis.at)) {
      grid.xaxis(at = c(ceiling(xscale[1L] * 10), floor(xscale[2L] * 10))/10)
    } else {
      if (is.null(xaxis.labs)) xaxis.labs <- xaxis.at
      grid.xaxis(at = xaxis.at, label = xaxis.labs)
    }
    grid.yaxis(at = c(ceiling(yscale[1L] * 10), floor(yscale[2L] * 10))/10)
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    
    ## plot individual growth curves
    if (observed) {
      if(!is.null(cluster)) {
        cluster <- factor(cluster[ix])
        for(j in levels(cluster)) {
          jx <- which(cluster == j)
          grid.lines(unit(x[jx], "native"), unit(y[jx], "native"), gp = gpar(col = curvecol, lwd = curvelwd))
        }
      }
    }
    
    ## plot model-implied average curve
    if (fitmean) {
      grid.abline(unit(cf[1], "native"), unit(cf[2], "native"), gp = gpar(col = linecol, lwd = linelwd))
    }
    
    if(ylab != "") grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, "lines"), rot = 90)
    if(xlab != "") grid.text(xlab, x = unit(0.5, "npc"), y = unit(-2, "lines"))                
    if(pop) popViewport(2L) else upViewport(2L)
    
  }
  
  return(rval)
}
class(node_growthplot) <- "grapcon_generator"