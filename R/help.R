#' Global specs for rcsurplus
#' 
#' \code{rcsurplus1d.help}
# #' @export rcsurplus1d.help
#' 
rcsurplus1d.help <- list(
  
  Introduction = "Package rcsurplus lets you explore a novel surplus
production model for fisheries stock assessment. The model employs a production
function that differs from the canonical logistic (Schaefer) and Gompertz (Fox)
functions, but can still be related to the Pella-Tomlinson formulation.
This function is embedded in a state-space model, where observed 
catch-per-unit-effort indices and measures of fishing effort are used as input. 
You may define Bayesian prior densities for all model hyperparameters 
(carrying capacity, catchability, growth rate and error variance), 
as well as the state (annual stock biomass). 
The well-studied Namibian hake fishery is provided in a Shiny app as a case
study, via which you can compare the canonical models (Pella-Tomlinson,
Schaefer and Fox) with the new model.",
  
  Model1 = "In the classical approach of Pella and Tomlinson (1969), the
production function \\(h(B_t)\\) is defined as
$$ h(B_t) = \\frac{r}{\\phi} B_t 
\\left[1-\\left(\\frac{B_t}{K}\\right)^{\\phi}\\right],$$
where  \\(B_t\\)  denotes biomass for year \\(t\\), 
\\(r\\) is the intrinsic rate of population growth, \\(K\\) 
is the carrying capacity of the environment and \\(\\phi\\in]0,1]\\) is a shape
parameter.",
  
  Model2 = "The production function is employed to describe the stochastic
evolution of the biomass of exploited stocks,
$$B_{t+1} = \\left\\{B_t + h(B_t) - C_t\\right\\} e^{\\xi_t},$$
where \\(C_t\\) denotes annual catch  and \\(\\xi_t\\sim N[0,\\sigma] \\) 
is Gaussian noise with variance \\(\\sigma\\). Thus, the (process) equation for
the Pella-Tomlinson formulation is
$$B_{t+1} = \\left\\{ B_t \\left[ 1+ \\frac{r}{\\phi} 
\\left(1-\\left(\\frac{B_t}{K}\\right)^{\\phi}\\right)\\right] - 
C_t\\right\\} e^{\\xi_t}.$$",
  
  Model3 = "The Schaefer process equation is obtained by fixing
\\(\\phi = 1\\),
$$B_{t+1} = \\left\\{ B_t \\left[ 1+ \\frac{r}{\\phi}
\\left(1-\\frac{B_t}{K}\\right)\\right] 
- C_t \\right\\} e^{\\xi_t},$$
whereas the Fox process equation results from considering the limit
\\(\\phi \\rightarrow 0\\),
$$B_{t+1} = \\left\\{ B_t \\left[ 1+ r \\left(1-\\frac{\\log B_t}
{\\log K}\\right)\\right] - C_t \\right\\} e^{\\xi_t}.$$",
  Model4 = "Fisheries data generally consist of commercial catch and effort
records, from which Catch Per Unit Effort (CPUE) indices may be constructed
(\\(i_t, t = 1,\\ldots,n\\)). 
The latter are assumed proportional to current biomass, that is,
\\(i_t = q B_t\\),  with \\(q\\) being a constant catchability over time. 
In other words, the value of \\(q\\) defines how many fish are caught
(in the appropriate units), on average, for each unit of fishing effort. 
Many authors build upon this structure by inserting log-normally
distributed errors, to reach the stochastic observation equation
$$ i_t = q_t B_t e^{\\omega _t},$$
where \\(\\omega_t \\sim N[0,\\sigma]\\). 
Note that the variance of \\(\\omega_t\\), \\(\\sigma\\), matches that of
\\(\\xi_t\\). Although such equivalence is not required in state-space models,
it is used in order to compare 
results among models cited in the literature (viz., Parent and Rivot, 2012).",
  Model5 = "We believe there are two problems in the three process equations
(Pella-Tomlinson, Schaefer and Fox) described above. First,
they are hard to work with, as they describe a nonlinear relationship between
biomasses at consecutive time instants. Second, they include observed catch,
which should only appear in the observation equation via the CPUE index
\\(i_t\\). We solve the latter issue by replacing 
\\(C_t\\) with \\((1-e^{-F_t})\\times (B_t+h(B_t))\\), 
where \\(F_t\\) denotes the fishing mortality rate, which can be 
expressed as \\(F_t=qE_t\\). In this equation, \\(q\\) represents catchability
(see above) and \\(E_t\\) stands for reported fishing effort. Unlike total
catch, effort can be controlled by management, so we can more readily accept
it as a fixed parameter.",
  
  Model6 = "As for the tractability of the process equation, we introduce a
novel approach. The production function we propose is
$$ h(B_t) = \\frac{r_t}{\\phi} B_t \\left[ 1 - \\left( \\frac{B_t}{K} \\right)
^ {\\phi}\\right],$$
with $$r_t = \\phi \\left( \\frac{B_t}{K} \\right)^{-\\phi}.$$
This leads to the process equation
$$ B_{t+1} = B_t^{1-\\phi}K^{\\phi}e^{-F_t + \\xi_t},$$
which, unlike the three canonical models, yields linear relationships between
log-biomasses at consecutive time instants.",
  
  UI = "This user interface was designed to help understand the features of
package rcsurplus. In section 'About', you can learn about the motivation for
this package, the equations underlying the four available models, and the
fitting procedure. In section 'Input', you can examine the Namibian hake dataset
(catches, effort, and CPUE), decide which model(s) to fit, define the bounds for
the uniform priors associated with the unknown parameters in these models, and 
set up the Markov Chain Monte Carlo fitting method. In section 'Output' you can
explore the results of the model fit(s) you have just defined. You can inspect
posterior density plots and scatterplots, you can compare model estimates
(posterior means and credibility intervals) with actual observations, and you
    can analyse MCMC convergence diagnostics.",
  
  Fitting = "Exploration of posterior distributions of unknown parameters is
accomplished via Markov Chain Monte Carlo methods. For this you will need to
have the software OpenBUGS installed on your computer. The interface between R
and OpenBUGS is provided by the package R2OpenBUGS. Before you decide whether
or not a model fits a data set, be sure to check that the MCMC chains converged
and passed most diagnostics."   
  
)