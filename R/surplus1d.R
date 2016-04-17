#' rcsurplus
#'
#' @field fit_counter numeric. NUmber of times the model has been fitted.
#' @field models list. List of models handled by the daemon.
#' @field coda_models list. List of CODA output.
#' @field model_names character. List with names of fitted models.
#' @field par_titles data.frame. Parameter names, as they appear in plots.
#' @field IO .IAT. Reference class object that handles OpenBUGS output.
#'
#' @importFrom methods setRefClass
#' @import rcvirtual
#' @import ggplot2
#' @import R2OpenBUGS
#' @import shiny
#' @importFrom shinyBS bsAlert
#' @importFrom shinyBS bsButton
#' @importFrom shinyBS createAlert
#' 
#' @export rcsurplus1d
#' @exportClass rcsurplus1d
rcsurplus1d <- setRefClass(
    Class = "rcsurplus1d",
    fields = list(
        fit_counter = "numeric",
        models = "list",
        coda_models = "list",
        model_names = "character",
        par_titles = "data.frame",
        IO = "rcvirtual.plotter"
    ),
    methods = list(
        initialize = function(){
            .self$fit_counter <- 0
            .self$models <- vector("list", length = 2)
            .self$IO <- rcplotter()
            .self$model_names <- c("Pella-Tomlinson", "Schaefer",
                                   "Fox", "Alternative")
            .self$par_titles <- data.frame(
                pella_tomlinson = c("K", "r", "phi", "log(q)",
                   "log(sigma)","MSY","B(MSY)","F(MSY)"),
                schaefer = c("K", "r", "phi", "log(q)",
                   "log(sigma)","MSY","B(MSY)","F(MSY)"),
                fox = c("K", "r", "phi", "log(q)",
                    "log(sigma)","MSY","B(MSY)","F(MSY)"),
                alternative = c("exp(rho)", "r", "phi", "chi",
                   "log(sigma)","MSY","B(MSY)","F(MSY)"),
                stringsAsFactors = FALSE)
        },
        
        get_data = function(input){
            if (is.null(input$mydata)) {
                mydata <- hake
            } else {
                mydata <- read.csv(input$mydata$datapath)
            }
            return(mydata)
        },
        
        fit_models = function(input){
            "Function to fit surplus production models"
            
            if (any(input$spm == "Pella-Tomlinson")) {
                .self$get_pella_tomlinson(input)
            } else {
                .self$models[[1]] <- NA
            }
            if (any(input$spm == "Schaefer")) {
                .self$get_schaefer(input)
            } else {
                .self$models[[2]] <- NA
            }
            if (any(input$spm == "Fox")) {
                .self$get_fox(input)
            } else {
                .self$models[[3]] <- NA
            }
            if (any(input$spm == "Alternative")) {
                .self$get_alternative(input)
            } else {
                .self$models[[4]] <- NA
            }
            .self$fit_counter <- .self$fit_counter + 1
        },
        
        get_pella_tomlinson = function(input){
            mydata <- .self$get_data(input)
            model_file <- file.path(tempdir(), "pella_tomlinson_model.txt") 
            nobs     <- nrow(mydata)
            C        <- mydata$catch
            Ef       <- mydata$effort
            Y        <- log(C/Ef)
            lsprior  <- input$sprior
            rprior   <- input$rprior
            PHIprior <- input$PHIprior
            Kprior   <- input$Kprior
            lqprior  <- input$qprior
            mymodel  <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                r ~ dunif(rprior[1],rprior[2])
                PHI ~ dunif(PHIprior[1],PHIprior[2])
                lq ~ dunif(lqprior[1],lqprior[2])
                K ~ dunif(Kprior[1],Kprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                q   <- exp(lq)
                lK  <- log(K)
                #######################################
                # Process equation (P[t] = B[t]/K) ####
                #######################################
                meanP[1] <- 1
                for (t in 1:(nobs - 1)) {
                    meanP[t + 1] <- max(P[t] + (r / PHI) * P[t] * 
                                            (1 - pow(P[t], PHI)) - 
                                            C[t] / K, 0.001)
                }
                for (t in 1:nobs) {
                    lmP[t] <- log(meanP[t]) 
                    P[t] ~ dlnorm(lmP[t], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- lq + log(P[t]) + lK
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                MSY  <- r * K * exp(-1 * (1 / PHI + 1) * log(1 + PHI)) 
                BMSY <- K * exp( -1/PHI *log(1 + PHI)) 
                FMSY <- r / (1 + PHI) 
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef","lsprior","rprior",
                           "PHIprior", "Kprior","lqprior")
            params <- c("K", "r", "PHI", "lq", "lsigma", "MSY",
                        "BMSY", "FMSY", "mean")
            inits  <- function() { list(  r = mean(input$rprior),
                                          PHI = mean(input$PHIprior),
                                          lq = mean(input$qprior), 
                                          K = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1 - input$MCMCb)))
            .self$models[[1]] <- bugs(data = data, inits = inits,
                                      parameters.to.save = params,  
                                      model.file = model_file,
                                      n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE,
                                      over.relax = input$overrelax )
            if (input$do_coda) {
                outS <- bugs(data = data, inits = inits, 
                             parameters.to.save = c("K", "r", "PHI",
                                                    "lq", "lsigma"),  
                             model.file = model_file, n.chains = input$MCMCc,
                             n.iter = lchain, n.thin = input$MCMCt, 
                             n.burnin = floor(lchain * input$MCMCb),
                             DIC = FALSE, codaPkg = TRUE,
                             over.relax = input$overrelax )
                .self$coda_models[[1]] <- read.bugs(outS)
            }
        },
        
        get_schaefer = function(input){
            mydata <- .self$get_data(input)
            model_file <- file.path(tempdir(), "schaefer_model.txt") 
            nobs    <- nrow(mydata)
            C       <- mydata$catch
            Ef      <- mydata$effort
            Y       <- log(C/Ef)
            lsprior <- input$sprior
            rprior  <- input$rprior
            Kprior  <- input$Kprior
            lqprior <- input$qprior
            mymodel <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                r ~ dunif(rprior[1],rprior[2])
                lq ~ dunif(lqprior[1],lqprior[2])
                K ~ dunif(Kprior[1],Kprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                q   <- exp(lq)
                lK  <- log(K)
                #######################################
                # Process equation (P[t] = B[t]/K) ####
                #######################################
                meanP[1] <- 1
                for (t in 1:(nobs - 1)) {
                    meanP[t + 1] <- max((P[t] + r * P[t] * (1 - P[t]) - 
                                             C[t] / K), 0.001) #eps=0.001
                }
                for (t in 1:nobs) {
                    lmP[t] <- log(meanP[t]) 
                    P[t] ~ dlnorm(lmP[t], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- lq + log(P[t]) + lK
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                MSY <- r * K / 4
                BMSY <- K / 2
                FMSY <- r/2
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef", "lsprior",
                           "rprior", "Kprior", "lqprior")
            params <- c("K", "r", "lq", "lsigma", "MSY", "BMSY", "FMSY", "mean")
            inits  <- function() {list(r = mean(input$rprior),  
                                       lq = mean(input$qprior), 
                                       K = mean(input$Kprior),
                                       lsigma = mean(input$sprior))}
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1 - input$MCMCb)))
            .self$models[[2]] <- bugs(data = data, inits = inits,
                                      parameters.to.save = params,
                                      model.file = model_file,
                                      n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE,
                                      over.relax = input$overrelax )
            if (input$do_coda) {
                outS <- bugs(data = data, inits = inits, 
                             parameters.to.save = c("K", "r", "lq", "lsigma"),  
                             model.file = model_file, n.chains = input$MCMCc,
                             n.iter = lchain, n.thin = input$MCMCt, 
                             n.burnin = floor(lchain * input$MCMCb),
                             DIC = FALSE, codaPkg = TRUE,
                             over.relax = input$overrelax )
                .self$coda_models[[2]] <- read.bugs(outS)
            }
        },
        
        get_fox = function(input){
            mydata <- .self$get_data(input)
            model_file <- file.path(tempdir(), "fox_model.txt") 
            nobs    <- nrow(mydata)
            C       <- mydata$catch
            Ef      <- mydata$effort
            Y       <- log(C/Ef)
            lsprior <- input$sprior
            rprior  <- input$rprior
            Kprior  <- input$Kprior
            lqprior <- input$qprior
            mymodel <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                r ~ dunif(rprior[1],rprior[2])
                lq ~ dunif(lqprior[1],lqprior[2])
                K ~ dunif(Kprior[1],Kprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                q   <- exp(lq)
                lK  <- log(K)
                #######################################
                # Process equation (P[t] = B[t]/K) ####
                #######################################
                meanP[1] <- 1
                for (t in 1:(nobs - 1)) {
                    meanP[t + 1] <- max((P[t] + r * P[t] * (1 - log(P[t] * K) /
                                         log(K)) - C[t] / K ), 0.001)
                }
                for (t in 1:nobs) {
                    lmP[t] <- log(meanP[t]) 
                    P[t] ~ dlnorm(lmP[t], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- lq + log(P[t]) + lK
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                MSY <- r*K / (exp(1)*log(K))
                BMSY <- K / exp(1)
                FMSY <- r/log(K)
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef",
                           "lsprior", "rprior", "Kprior", "lqprior")
            params <- c("K", "r", "lq", "lsigma", "MSY", "BMSY", "FMSY", "mean")
            inits  <- function() { list(  r = mean(input$rprior),  
                                          lq = mean(input$qprior), 
                                          K = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1 - input$MCMCb)) )
            .self$models[[3]] <- bugs(data = data, inits = inits,
                                      parameters.to.save = params,  
                                      model.file = model_file,
                                      n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE,
                                      over.relax = input$overrelax )
            if (input$do_coda) {
                outS <- bugs(data = data, inits = inits, 
                             parameters.to.save = c("K", "r", "lq", "lsigma"),  
                             model.file = model_file, n.chains = input$MCMCc,
                             n.iter = lchain, n.thin = input$MCMCt, 
                             n.burnin = floor(lchain * input$MCMCb),
                             DIC = FALSE, codaPkg = TRUE,
                             over.relax = input$overrelax )
                .self$coda_models[[3]] <- read.bugs(outS)
            }
        },
        
        get_alternative = function(input){
            mydata <- .self$get_data(input)
            model_file <- file.path(tempdir(), "alternative_model.txt") 
            nobs       <- nrow(mydata)
            C          <- mydata$catch
            Ef         <- mydata$effort
            Y          <- log(C/Ef)
            lsprior    <- input$sprior
            PHIprior   <- input$PHIprior
            eRHOprior  <- input$Kprior
            CHIprior   <- input$qprior
            mymodel    <- function(){
                #######################################
                # Priors ##############################
                #######################################
                lsigma ~ dunif(lsprior[1],lsprior[2])
                PHI ~ dunif(PHIprior[1],PHIprior[2])
                CHI ~ dunif(CHIprior[1],CHIprior[2])
                eRHO ~ dunif(eRHOprior[1],eRHOprior[2])
                #######################################
                # Transformations #####################
                #######################################
                tau <- 1/exp(lsigma) #precision
                eCHI  <- exp(CHI)
                RHO   <- log(eRHO)
                #######################################
                # Process equation: A[t] = log(B[t]/K)
                #######################################
                A[1] ~ dnorm(0,tau)
                for (t in 1:(nobs - 1)) {
                    meanA[t + 1] <- (1 - PHI) * A[t] - eCHI * Ef[t]
                    A[t + 1] ~ dnorm(meanA[t + 1], tau)
                }
                #######################################
                # Observation equation ################
                #######################################
                for (t in 1:nobs ) {
                    mean[t] <- CHI + A[t] + RHO
                    Y[t] ~ dnorm(  mean[t], tau )
                }
                #######################################
                # Derived quantities ##################
                #######################################
                BMSY <- exp( (1.0 / PHI) * log(1.0 - PHI) + RHO )
                MSY <-  BMSY * PHI / (1.0 - PHI)
                FMSY <- -log(1.0 - PHI)
            }
            write.model(mymodel, model_file)
            data   <- list("nobs", "Y", "C", "Ef","lsprior",
                           "PHIprior", "eRHOprior", "CHIprior")
            params <- c("eRHO", "PHI", "CHI", "lsigma", "MSY",
                        "BMSY", "FMSY", "mean")
            inits  <- function() { list(  PHI = mean(input$PHIprior),  
                                          CHI = mean(input$qprior), 
                                          eRHO = mean(input$Kprior),
                                          lsigma = mean(input$sprior) ) }
            lchain <- ceiling(input$MCMCn / (input$MCMCc * (1 - input$MCMCb)))
            .self$models[[4]] <- bugs(data = data, inits = inits,
                                      parameters.to.save = params,  
                                      model.file = model_file,
                                      n.chains = input$MCMCc,
                                      n.iter = lchain, n.thin = input$MCMCt, 
                                      n.burnin = floor(lchain * input$MCMCb),
                                      DIC = TRUE, codaPkg = FALSE,
                                      over.relax = input$overrelax )
            if (input$do_coda) {
                outA <- bugs(data = data, inits = inits, 
                            parameters.to.save = c("eRHO", "PHI",
                                                   "CHI", "lsigma"),  
                            model.file = model_file, n.chains = input$MCMCc,
                            n.iter = lchain, n.thin = input$MCMCt, 
                            n.burnin = floor(lchain * input$MCMCb),
                            DIC = FALSE, codaPkg = TRUE,
                            over.relax = input$overrelax )
                .self$coda_models[[4]] <- read.bugs(outA)
            }
        },
        plot_hyperparameters = function(input){
            nmod <- length(input$spm_hyper)
            all_titles <- c("Carrying capacity", "Growth rate",
                            "Elasticity", "Log-catchability",
                            "Log-variance", "MSY", "BMSY", "FMSY")
            pars <- c(input$plot_KeRHO, input$plot_r, input$plot_PHI,
                      input$plot_qCHI, input$plot_lSIGMA,
                      input$plot_msy,input$plot_bmsy,input$plot_fmsy)
            ploc <- list(pella_tomlinson = c(1:8), schaefer = c(1, 2, NA, 3:7),
                         fox = c(1, 2, NA, 3:7), alternative = c(1, NA, 2:7))
            plot_titles <- list(pella_tomlinson = all_titles, 
                                schaefer = c(all_titles[1:2], all_titles[4:8]), 
                                fox = c(all_titles[1:2], all_titles[4:8]), 
                                alternative = c(all_titles[1], all_titles[3:8]))
            if (nmod > 0) for (kk in 1:nmod) {
                k <- which(input$spm_hyper[kk] == .self$model_names)
                for (i in 1:length(pars)) {
                    pars[i] <- (pars[i] & !is.na(ploc[[k]][i]) )
                }
            }
            npars <- sum(pars)
            if (npars == 0 | nmod == 0) return()
            .self$IO$set_buffer_size(nr = npars, nc = npars)
            i <- 1
            for (ii in 1:length(pars)) if (pars[ii]) {
                j <- 1 #column index
                for (jj in 1:ii) if (pars[jj]) {
                    x <- y <- mytx <- myty <- nm <- NULL
                    for (kk in 1:nmod) {
                        k <- which(input$spm_hyper[kk] == .self$model_names)
                        ib <- ploc[[k]][ii]
                        jb <- ploc[[k]][jj]
                        if (length(.self$models[[k]]) > 1) {
                            aux <- .self$models[[k]]$sims.matrix[,ib]
                            x <- c(x,aux)
                            if (!is.na(ib)) mytx <- plot_titles[[k]][ib]
                            nmaux <- rep(input$spm_hyper[kk], length(aux))
                            nm <- c(nm,nmaux)
                            if (jj < ii) {
                                aux <-  .self$models[[k]]$sims.matrix[,jb]
                                y <- c(y,aux)
                                if (!is.na(jb)) myty <- plot_titles[[k]][jb]
                            }
                        }
                    }
                    if (jj < ii) {
                        out <- data.frame(x = y, y = x, model = as.factor(nm))
                        mytitle <- c(myty,mytx)
                        .self$IO$get_scatterplot(out, mytitle,
                                                 xpos = i, ypos = j)
                    } else {
                        out <- data.frame(x = x, model = as.factor(nm))
                        dolegend <- (i == 1 & j == 1 & 
                                     length(input$spm_hyper) > 1)
                        .self$IO$get_density_plot(
                            out, mytitle = mytx, dolegend = dolegend,
                            xpos = i, ypos = i)
                    }
                    j <- j + 1
                }
                i <- i + 1
            }
            .self$IO$get_buffer_plot()
        },
        
        plot_data = function(input){
            mydata <- .self$get_data(input)
            do_catch  <- input$do_catch
            do_effort <- input$do_effort
            do_cpue   <- input$do_cpue
            one_row   <- input$one_row
            .self$IO$get_data_plot(mydata, do_catch, do_effort, do_cpue,
                                   one_row)
        },
        
        plot_fitted_cpue = function(input){
            mydata <- .self$get_data(input)
            nobs <- nrow(mydata)
            fn <- function(data_vec){
                return(exp(quantile(data_vec,c(0.025, 0.5, 0.975))))
            }
            rn <- c(1:nobs)
            .self$IO$set_buffer_size(nr = 1, nc = 1)
            year <- low95 <- median <- high95 <- mname <- obs <- NULL
            obs_aux <- (mydata$catch / mydata$effort)[rn]
            nmod <- length(input$spm_fit)
            if (nmod > 0) {
                for (kk in 1:nmod) {
                    k <- which(input$spm_fit[kk] == .self$model_names)
                    if (length(.self$models[[k]]) > 1) {
                        lcpue  <- .self$models[[k]]$sims.list[["mean"]]
                        cpue   <- apply(X = lcpue, MARGIN = 2, FUN = fn)
                        year   <- c(year, rn + mydata$year[1] - 1)
                        mname  <- c(mname, rep(input$spm_fit[kk], length(rn)))
                        low95  <- c(low95,  cpue[1,rn])
                        median <- c(median, cpue[2,rn])
                        high95 <- c(high95, cpue[3,rn])
                        obs    <- c(obs, obs_aux)
                    }
                }
                out <- data.frame(year = year, low95 = low95,
                                  median = median, high95 = high95, 
                                  obs = obs, model = as.factor(mname))
                .self$IO$get_ts_fit_plot(
                    out, mytitle = "Observations vs. median and 95% CI", 
                    mylabs = c("year", "cpue"), xpos = 1, ypos = 1)
                .self$IO$get_buffer_plot()
            }
        },
        
        plot_fitted_cpue_old = function(input){
            mydata <- .self$get_data(input)
            nobs <- nrow(mydata)
            fn <- function(data_vec) {
                return(exp(quantile(data_vec,c(0.025, 0.5, 0.975))))
            }
            if (input$spm_fit == "Schaefer and Alternative" &
                !input$fit_one_row) {
                nr <- 2
            } else {
                nr <- 1
            }
            if (input$spm_fit == "Schaefer and Alternative" &
               input$fit_one_row) {
                nc <- 2
            } else {
                nc <- 1
            }
            k  <- 1
            rn <- c(1:nobs)
            mt <- c("Schaefer", "Alternative")
            .self$IO$set_buffer_size(nr = nr, nc = nc)
            for (i in 1:nr) for (j in 1:nc) {
                if (length(.self$models[[k]]) > 1) {
                    lcpue <- .self$models[[k]]$sims.list[["mean"]]
                    cpue  <- apply(X = lcpue, MARGIN = 2, FUN = fn)
                    out   <- data.frame(
                        year = rn + mydata$year[1] - 1,
                        obs = (mydata$catch / mydata$effort)[rn],
                        low95 = cpue[1, rn], median = cpue[2, rn],
                        high95 = cpue[3, rn])
                    .self$IO$get_ts_fit_plot(
                        out, mytitle = mt[k], mylabs = c("year", "cpue"),
                        xpos = i, ypos = j)
                }
                k <- k + 1
            }
            .self$IO$get_buffer_plot()
        },
        
        gui = function(){
            m <- .self
            server <- function(input, output, session) {
                
                output$plot_data  <- renderPlot({
                    m$plot_data(input)
                })
                
                output$model_settings <- renderPrint({
                    cat(floor(input$MCMCt * input$MCMCn / input$MCMCb ))
                })
                
                output$model_status <- renderPrint({
                    if(any(input$Kprior>0) | any(input$rprior>0) | any(input$qprior>0) | any(input$sprior>0) |
                       input$MCMCn>0 | input$MCMCb>0 | input$MCMCt>0 | input$MCMCc>0 | 
                       length(input$spm) == 0 | input$do_coda){
                        shinyBS::createAlert(
                            session, anchorId = "model_inputId", alertId = "model_alertId", 
                            content = "Model(s) not fitted.", 
                            style = "danger", dismiss = FALSE, append = FALSE )
                    }
                    if(input$fit_model > m$fit_counter){
                        shinyBS::createAlert(
                            session, anchorId = "model_inputId", alertId = "model_alertId", 
                            content = "Fitting model(s), please wait.", 
                            style = "warning", dismiss = FALSE, append = FALSE )
                        ok <- m$fit_models(input)
                        shinyBS::createAlert(session, anchorId = "model_inputId", alertId = "model_alertId", 
                                    content = "Model(s) fitted.", 
                                    style = "success", dismiss = FALSE, append = FALSE )
                    }
                })
                
                output$fit_plot  <- renderPlot(
                    m$plot_fitted_cpue(input)
                )
                
                output$plot_hyperparameters <- renderPlot({
                    m$plot_hyperparameters(input)
                })
                
                output$summary_dev <- renderPrint({
                    sep <- c('*** Schaefer ****************************************\n',
                             '*** Alternative *************************************\n')
                    if(input$diag == 'Summary & Deviance'){
                        if(any(input$spm_diag == 'Pella-Tomlinson')){
                            cat('*** Pella-Tomlinson *********************************\n')
                            print(m$models[[1]])
                            cat('\n')
                        }
                        if(any(input$spm_diag == 'Schaefer')){
                            cat('*** Schaefer ****************************************\n')
                            print(m$models[[2]])
                            cat('\n')
                        }
                        if(any(input$spm_diag == 'Fox')){
                            cat('*** Fox *********************************************\n')
                            print(m$models[[3]])
                            cat('\n')
                        }
                        if(any(input$spm_diag == 'Alternative')){
                            cat('*** Alternative *************************************\n')
                            print(m$models[[4]])
                            cat('\n')
                        }
                    } else if(input$diag == 'CODA - Gelman' & input$do_coda){
                        if(input$spm_diag != 'Alternative') {
                            cat(sep[1])
                            print(gelman.diag(m$coda_models[[1]]))
                            cat('\n')
                        }
                        if(input$spm_diag != 'Schaefer'){
                            cat(sep[2])
                            print(gelman.diag(m$coda_models[[2]]))
                            cat('\n')
                        }
                    } else if(input$diag == 'CODA - Geweke' & input$do_coda){
                        if(input$spm_diag != 'Alternative') {
                            cat(sep[1])
                            print(geweke.diag(m$coda_models[[1]], frac1=0.1, frac2=0.5))
                            cat('\n')
                        }
                        if(input$spm_diag != 'Schaefer'){
                            cat(sep[2])
                            print(geweke.diag(m$coda_models[[2]], frac1=0.1, frac2=0.5))
                            cat('\n')
                        }
                    } else if(input$diag == 'CODA - Raftery & Lewis' & input$do_coda){
                        if(input$spm_diag != 'Alternative') {
                            cat(sep[1])
                            print(raftery.diag(data=m$coda_models[[1]], 
                                               q=0.025, r=0.005, s=0.95, converge.eps=0.001))
                            cat('\n')
                        }
                        if(input$spm_diag != 'Schaefer'){
                            cat(sep[2])
                            print(raftery.diag(data=m$coda_models[[2]], 
                                               q=0.025, r=0.005, s=0.95, converge.eps=0.001))
                            cat('\n')
                        }
                    }
                })
            }
            
            ui <- shiny::navbarPage(
                title = 'RC SURPLUS',
                tabPanel("About",
                         shiny::navbarPage('About',
                                    tabPanel( 'Introduction', fluidRow(helpText(
                                        rcsurplus1d.help$Introduction )) ),
                                    tabPanel( 'Model',
                                              fluidRow(
                                                  withMathJax(''),
                                                  helpText( rcsurplus1d.help$Model1 ),
                                                  helpText( rcsurplus1d.help$Model2 ),
                                                  helpText( rcsurplus1d.help$Model3 ),
                                                  helpText( rcsurplus1d.help$Model4 ),
                                                  helpText( rcsurplus1d.help$Model5 ),
                                                  helpText( rcsurplus1d.help$Model6 )
                                              )),
                                    tabPanel( 'UI', fluidRow(helpText( rcsurplus1d.help$UI )) ),
                                    tabPanel( 'Fitting', fluidRow(helpText( rcsurplus1d.help$Fitting )) )
                         )
                ),
                tabPanel("Input",
                         shiny::navbarPage('Input',
                                    tabPanel('Data',
                                             fluidPage(
                                                 sidebarLayout(
                                                     sidebarPanel(
                                                         fileInput('mydata', 'Data file',
                                                                   accept = c(
                                                                       'text/csv',
                                                                       'text/comma-separated-values',
                                                                       'text/tab-separated-values',
                                                                       'text/plain',
                                                                       '.csv',
                                                                       '.tsv'
                                                                   )
                                                         ),
                                                         checkboxInput("do_catch", "Plot catch", TRUE),
                                                         checkboxInput("do_effort", "Plot effort", FALSE),
                                                         checkboxInput("do_cpue", "Plot cpue", FALSE),
                                                         checkboxInput("one_row", "Plot in one row", TRUE)
                                                     ),
                                                     mainPanel(
                                                         plotOutput("plot_data")
                                                     )
                                                 )
                                             )
                                    )
                         )
                ),
                tabPanel('Model',
                         fluidPage(
                             br(),
                             fluidRow(
                                 column(4,
                                        h4("Model(s) to fit (multi-choice)"),
                                        selectizeInput("spm", "",choices =
                                                           c('Pella-Tomlinson', 'Schaefer', 'Fox', 'Alternative'),
                                                       multiple=TRUE),
                                        h4("Bounds of uniform priors"),
                                        sliderInput("Kprior", "Carrying capacity: \\(K\\) and \\(e^\\rho\\)",
                                                    min = myglobal$priorK[1],
                                                    max = myglobal$priorK[2],
                                                    value = c(myglobal$priorK[1],myglobal$priorK[2]),
                                                    step = 100 ),
                                        sliderInput("rprior", "Growth rate: \\(r\\)",
                                                    min = myglobal$priorr[1],
                                                    max = myglobal$priorr[2],
                                                    value = c(myglobal$priorr[1],myglobal$priorr[2]),
                                                    step = 0.01 ),
                                        sliderInput("PHIprior", "Elasticity: \\(\\phi\\)",
                                                    min = myglobal$priorPHI[1],
                                                    max = myglobal$priorPHI[2],
                                                    value = c(myglobal$priorPHI[1],myglobal$priorPHI[2]),
                                                    step = 0.01 ),
                                        sliderInput("qprior", "Log-catchability: \\(\\log (q)\\) and \\(\\chi\\)",
                                                    min = myglobal$priorq[1],
                                                    max = myglobal$priorq[2],
                                                    value = c(myglobal$priorq[1],myglobal$priorq[2]),
                                                    step = 0.01 ),
                                        sliderInput("sprior", "Model log-variance: \\(\\log (\\sigma)\\)",
                                                    min = myglobal$priors[1],
                                                    max = myglobal$priors[2],
                                                    value = c(myglobal$priors[1],myglobal$priors[2]),
                                                    step = 0.01 )
                                 ),
                                 column(4,
                                        h4("MCMC settings"),
                                        sliderInput("MCMCn", "Number of iterations for output:",
                                                    min = myglobal$mcmc_n[1],
                                                    max = myglobal$mcmc_n[2],
                                                    value = myglobal$mcmc_n[1],
                                                    step = 100 ),
                                        sliderInput("MCMCb", "Burn-in ratio:",
                                                    min = myglobal$mcmc_b[1],
                                                    max = myglobal$mcmc_b[2],
                                                    value = myglobal$mcmc_b[2],
                                                    step = 0.05 ),
                                        sliderInput("MCMCt", "Thinning factor:",
                                                    min = myglobal$mcmc_t[1],
                                                    max = myglobal$mcmc_t[2],
                                                    value = myglobal$mcmc_t[1],
                                                    step = 1 ),
                                        sliderInput("MCMCc", "Number of chains:",
                                                    min = myglobal$mcmc_c[1],
                                                    max = myglobal$mcmc_c[2],
                                                    value = myglobal$mcmc_c[1],
                                                    step = 1 ),
                                        checkboxInput("overrelax", "Over-relax", FALSE),
                                        p('Total no. MCMC iterations:'),
                                        verbatimTextOutput('model_settings')
                                 ),
                                 column(4,
                                        h4("Convergence diagnostics"),
                                        checkboxInput("do_coda", "CODA", FALSE),
                                        h4("Model Status"),
                                        shinyBS::bsAlert("model_inputId"),
                                        shinyBS::bsButton('fit_model', label = "Fit model(s)", style = "info", size="mini"),
                                        verbatimTextOutput('model_status')
                                 )
                             )
                         )
                ),
                tabPanel('Output',
                         shiny::navbarPage('Output',
                                    tabPanel('Hyperparameters',
                                             fluidPage(
                                                 sidebarLayout(
                                                     sidebarPanel(
                                                         h4('Plot:'),
                                                         selectizeInput("spm_hyper", "",choices = c(
                                                             'Pella-Tomlinson',
                                                             'Schaefer',
                                                             'Fox',
                                                             'Alternative'), multiple=TRUE),
                                                         checkboxInput("plot_KeRHO", "Carrying capacity: \\(K\\) and/or \\(e^\\rho\\)", FALSE),
                                                         checkboxInput("plot_r", "Growth rate: \\(r\\)", FALSE),
                                                         checkboxInput("plot_PHI", "Elasticity: \\(\\phi\\)", FALSE),
                                                         checkboxInput("plot_qCHI", "Log-catchability: \\(\\log(q)\\) and/or \\(\\chi\\)", FALSE),
                                                         checkboxInput("plot_lSIGMA", "Model log-variance: \\(\\log(\\sigma)\\)", FALSE),
                                                         checkboxInput("plot_msy", "MSY", FALSE),
                                                         checkboxInput("plot_bmsy", "B(MSY)", FALSE),
                                                         checkboxInput("plot_fmsy", "F(MSY)", FALSE)
                                                     ),
                                                     mainPanel(
                                                         plotOutput("plot_hyperparameters")
                                                     )
                                                 )
                                             )
                                    ),
                                    tabPanel('Fit',
                                             fluidPage(
                                                 sidebarLayout(
                                                     sidebarPanel(
                                                         h4('Observed vs. fitted CPUE:'),
                                                         selectizeInput("spm_fit", "",
                                                                        choices = c('Pella-Tomlinson',
                                                                                    'Schaefer',
                                                                                    'Fox',
                                                                                    'Alternative'),
                                                                        multiple=TRUE),
                                                         checkboxInput("fit_one_row", "Plot in one row", TRUE)
                                                     ),
                                                     mainPanel(
                                                         plotOutput("fit_plot")
                                                     )
                                                 )
                                             )
                                    ),
                                    tabPanel('Diagnostics',
                                             fluidPage(
                                                 sidebarLayout(
                                                     sidebarPanel(
                                                         h4('Model:'),
                                                         selectizeInput("spm_diag", "",choices = c('Pella-Tomlinson',
                                                                                                   'Schaefer',
                                                                                                   'Fox',
                                                                                                   'Alternative'), multiple=TRUE),
                                                         h4('Diagnostic:'),
                                                         selectInput("diag", "",choices = c('Summary & Deviance',
                                                                                            'CODA - Gelman',
                                                                                            'CODA - Geweke',
                                                                                            'CODA - Raftery & Lewis'))
                                                     ),
                                                     mainPanel(
                                                         verbatimTextOutput("summary_dev")
                                                     )
                                                 )
                                             )
                                    )
                         )
                )
            )
            
            shiny::shinyApp(ui = ui, server = server)
        }
    )
)