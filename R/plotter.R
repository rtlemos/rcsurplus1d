#' rcplotter: reference class that plots results of surplus model fits
#'
#' @field buffer matrix. 
#' @field palettes list. 
#' @field default_palettes list. 
#'
#' @import rcvirtual
#' @import grid
#' @importFrom methods setRefClass
#'
# #' @export rcplotter
# #' @exportClass rcplotter
rcplotter <- setRefClass(
    Class = 'rcplotter',
    contains = 'rcvirtual.plotter',
    fields = list(buffer = 'matrix', palettes = 'list',
                  default_palettes = 'list'),
    methods = list(
        initialize     = function(){
            "Initializes the printer object"
            
            #default palettes are colour-blind friendly 
            dfp <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
            .self$default_palettes <- list( 
                line = dfp, 
                fill = dfp,
                CI = c(0.8,0.7),
                neg_zero_pos = c('blue','white','red'),
                zero_pos = c('white', 'black'))
            .self$palettes <- list(
                line = .self$default_palettes$line,
                fill = .self$default_palettes$fill,
                CI = .self$default_palettes$CI,
                neg_zero_pos = .self$default_palettes$neg_zero_pos,
                zero_pos = .self$default_palettes$zero_pos
            )
        },
        
        set_buffer_size = function(nr, nc){
            "Sets up plotter's buffer size (number of rows and columns)"
            
            .self$buffer <- array(list(),c(nr,nc))
        },
        
        set_palette = function(argument, value){
            "Sets the printer's colour palettes"
            
            if (argument == 'all') {
                nm <- names(.self$palettes)
                if (value[1] == 'default') {
                    lapply(seq(along = nm), function(i) {
                        .self$set_palette(nm[i], .self$default_palettes[[i]])
                    })
                } else {
                    lapply(seq(along = nm), function(i) {
                        .self$set_palette(nm[i], value)
                    })
                }
            } else {
                pos <- which(names(.self$palettes) == argument)
                if (length(pos) == 0) {
                    stop('Palette not found: ', argument)
                } else if (value[1] == 'default') {
                    .self$palettes[[pos]] <- .self$default_palettes[[pos]]
                } else {
                    .self$palettes[[pos]] <- value
                }
            }
        },
        
        set_in_buffer = function(myplot, xpos, ypos){
            'Places a plot in the buffer'
            
            .self$buffer[xpos,ypos] <- list(myplot)
        },
        
        get_palette = function(type){
            if (type == 'all') type <- names(.self$palettes)
            out <- mapply(1:length(type), FUN = function(i) {
                pos <- which(names(.self$palettes) == type[i])
                if (length(pos) == 0) {
                    return(type[i])
                } else {
                    print(paste0(type[i],' = c("', 
                         paste0(.self$palettes[[pos]], collapse = '", "') ,'")' 
                    ), quote = FALSE)
                    return('found')
                }
            })
            if (any(out != 'found')) {
                stop('Palette(s) not found: ', out[out != 'found'])
            }
        },
        
        get_buffer_plot = function(){
            'Prints the plots in the buffer'
            
            vplayout <- function(x, y) {
                viewport(layout.pos.row = x, layout.pos.col = y)
            }
            nr <- nrow(.self$buffer)
            nc <- ncol(.self$buffer)
            grid::grid.newpage()
            grid::pushViewport(grid::viewport(layout = grid.layout(nr, nc)))
            for (ii in 1:nr) for (jj in 1:nc) {
                myplot <- .self$buffer[ii,jj][[1]]
                if (!is.null(myplot)) {
                    print(myplot, vp = vplayout(ii,jj))
                }
            }
        },
        
        get_one_density_plot = function(out, mytitle, xpos=1, ypos=1){
            'Plots the posterior density of one hyperparameter'
            
            df <- data.frame(x = out)
            myplot <- ggplot(data = df, aes(x = x)) + geom_density() +
                xlab(mytitle) + ylab('density') +
                scale_colour_manual(values = .self$palettes$line, guide = FALSE)
            .self$set_in_buffer(myplot,xpos,ypos)
        },
        
        get_density_plot = function(out, mytitle, dolegend = TRUE,
                                    xpos = 1, ypos = 1){
            'Plots the posterior density of >= 1 hyperparameters'
            nl <- nlevels(out$model)
            fc <- .self$palettes$fill[1:nl]
            auxplot <- ggplot(data = out, aes(x = x, fill = model)) + 
                geom_density(alpha = 0.2) +
                xlab(mytitle) + ylab('density')
            if (dolegend) {
                myplot <- auxplot + 
                    scale_fill_manual(values = .self$palettes$fill) +
                    theme(legend.position = "top")
            } else {
                myplot <- auxplot + 
                    scale_fill_manual(values = .self$palettes$fill,
                                      guide = FALSE)
            }
            .self$set_in_buffer(myplot,xpos,ypos)
        },
        
        get_scatterplot = function(out, mytitle, xpos = 1, ypos = 1){
            'Plots a scatterplot of two hyperparameters'
            nl <- nlevels(out$model)
            fc <- .self$palettes$fill[1:nl]
            myplot <- ggplot(data = out, aes(x = x, y = y, colour = model)) + 
                geom_point(cex = 2, alpha = 0.2) +
                xlab(mytitle[1]) + ylab(mytitle[2]) +
                scale_colour_manual(values = .self$palettes$fill, guide = FALSE)
            .self$set_in_buffer(myplot, xpos, ypos)
        },
        
        get_ts_fit_plot = function(out, mytitle, mylabs, xpos = 1, ypos = 1){
            'Plots a time series plot with observed and fit (median plus 95 CI)'
            myplot <- ggplot(data = out,aes(x = year, y = obs, z = model)) + 
                geom_ribbon(data = out,
                            aes(x = year, ymin = low95, ymax = high95), 
                            fill = grey(.self$palettes$CI[1]),
                            alpha = .self$palettes$CI[2]) +
                geom_line(data = out,
                          aes(x = year, y = median, color = model)) + 
                geom_point() +
                scale_color_manual(values = .self$palettes$line) +
                xlab(mylabs[1]) + ylab(mylabs[2]) + labs(title = mytitle)
            .self$set_in_buffer(myplot, xpos, ypos)
        },
        
        get_data_plot = function(dt, do_catch = TRUE, do_effort = FALSE,
                                 do_cpue = FALSE, one_row = TRUE){
            nplots <- sum(c(do_catch, do_effort, do_cpue))
            if (nplots == 0) return()
            nr     <- if (one_row) 1 else nplots
            nc     <- if (one_row) nplots else 1
            bufr   <- array(list(), c(nr, nc))
            ridx   <- 1
            cidx   <- 1
            mtitle <- if (one_row) c('Time series of catch',
                                     'Time series of effort',
                                     'Time series of CPUE'
            ) else c('Time series of catch',
                     'Time series of effort', 
                     'Time series of CPUE')
            yl     <- if (one_row) c('catch', 'effort', 'CPUE [ton/h]'
            ) else c('catch', 'effort', 'CPUE')
            
            if (do_catch) {
                catch_plot <- ggplot( dt, aes(x = year, y = catch) ) +
                    geom_line() +
                    geom_point() +
                    ylab(yl[1]) +
                    ggtitle(mtitle[1])
                bufr[ridx,cidx] <- list(catch_plot)
                if (one_row) cidx <- cidx + 1 else ridx <- ridx + 1
            }
            
            if (do_effort) {
                effort_plot <- ggplot(dt, aes(x = year, y = effort)) +
                    geom_line() +
                    geom_point() +
                    ylab(yl[2]) +
                    ggtitle(mtitle[2])
                bufr[ridx,cidx] <- list(effort_plot)
                if (one_row) cidx <- cidx + 1 else ridx <- ridx + 1
            }
            
            if (do_cpue) {
                cpue_plot <- ggplot(dt, aes(x = year, y = catch / effort)) +
                    geom_line() +
                    geom_point() +
                    ylab(yl[3]) +
                    ggtitle(mtitle[3])
                bufr[ridx,cidx] <- list(cpue_plot)
                if (one_row) cidx <- cidx + 1 else ridx <- ridx + 1
            }
            
            vplayout <- function(x, y) {
                viewport(layout.pos.row = x, layout.pos.col = y)
            }
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nr, nc)))
            for (ii in 1:nr) for (jj in 1:nc) {
                myplot <- bufr[ii,jj][[1]]
                if (!is.null(myplot)) {
                    print(myplot, vp = vplayout(ii,jj))
                }
            }
        }
    )
)