######### Stuff for the Engineergin Symposium 2019 ##################
### libraries ####
library(rstanarm)
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
library(mstate)
library(bayesplot)
library(dplyr)
library(dirichletprocess)
library(qwraps2)
roxygen2::roxygenise(clean = TRUE)
source("~/aeim/R/simulation.R") # for simulation

### numerical identifiability of parametric model ####
# ---------------------------------------------------------#

lambdas01_t = 0.325
lambdas02_t = 0.36
lambdas12_t = 0.39
gammas01_t = 1.9
gammas02_t = 2.2
gammas12_t = 2.5
cens = c(4.5, 5.5)

set.seed(9911)
covs <- data.frame(id = 1:30000, trt = stats::rbinom(30000, 1L, 0.5))

sim_wei <- rsimid(
  dist01 = "weibull",
  dist02 = "weibull",
  dist12 = "weibull",
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  gammas01 = gammas01_t,
  gammas02 = gammas02_t,
  gammas12 = gammas12_t,
  x = covs,
  cens = cens
)

tmat_mst <- trans.illdeath(names=c("diagnosis","relapse","death"))
sim_wei_mstate <- msprep(time=c(NA,"df_time","os_time"),
                         status=c(NA,"df_event","os_event"),
                         data = sim_wei,
                         trans=tmat_mst)
sim_wei_mstate <- dplyr::left_join(sim_wei_mstate,
                                   sim_wei[ , c("id", "trt")])

prior_intercept = lapply(1:3, function(x)
  rstanarm::normal() )

prior_aux = lapply(1:3, function(x)
  rstanarm::normal() )

formula = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status) ~ 1) )

basehaz = lapply(1:3, function(x)
  "weibull")

sim_recov_1 <- mstte_stan(formula = formula,
                      data = sim_wei_mstate,
                      transition_labels = c("DP", "DX", "DPDX"),
                      basehaz = basehaz,
                      prior           = lapply(1:3, function(x)
                        rstanarm::normal() ),
                      prior_intercept = prior_intercept,
                      prior_aux       = prior_aux,
                      iter = 2000,
                      chains = 4,
                      control = list(adapt_delta = 0.99)
)

saveRDS(sim_recov_1, "~/rfactory/SympEng2019/sim_recov_1.RDS")


# exp_out <- as.data.frame(summary(sim_recov_1))
# exp_out <- exp_out[1:6, ]
#
# actual <- c(
#   "(Intercept)_1"  = log(lambdas01_t),
#   "weibull-shape_1" = gammas01_t,
#
#   "(Intercept)_2" = log(lambdas02_t),
#   "weibull-shape_2" = gammas02_t,
#
#   "(Intercept)_3" = log(lambdas12_t),
#   "weibull-shape_3" = gammas12_t)
#
# fit_1 <- cbind(actual, exp_out)
#
# write.csv(fit_1, "ni_parametric_weibull.csv")
## ----plot-alpha-vs-test-------------------------------------------------- BAYES PLOT -

# posterior <- as.array(stanfit)
# dim(posterior)
# dimnames(posterior)
#
# color_scheme_set("gray")
# p1 <- mcmc_hex(posterior, pars = c("(Intercept)_1", "weibull-shape_1")) +
#   geom_point(aes(x = log(lambdas01_t), y = gammas01_t ), shape = 25, size = 3, fill = "black") + ggtitle('Transition 0 -> 1')  +
#   labs(x = "scale", y = "shape")
#
# p2 <- mcmc_hex(posterior, pars = c("(Intercept)_2", "weibull-shape_2")) +
#   geom_point(aes(x = log(lambdas02_t), y = gammas02_t ), shape = 25, size = 3, fill = "black") + ggtitle('Transition 0 -> 2')  +
#   labs(x = "scale", y = "shape")
#
# p3 <- mcmc_hex(posterior, pars = c("(Intercept)_3", "weibull-shape_3")) +
#   geom_point(aes(x = log(lambdas12_t), y = gammas12_t ), shape = 25, size = 3, fill = "black") + ggtitle('Transition 1 -> 2')  +
#   labs(x = "scale", y = "shape")
#
#
# library(cowplot)
# title <- ggdraw() + draw_label("Posterior joint distribution of Weibull-shape and scale showing true parameter values", fontface='bold')
#
# p <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
# pp <- cowplot::plot_grid(title, p, ncol = 1, rel_heights=c(0.1, 1))
# pp
#
# ggsave(plot = pp, filename = "nipp_wei.png", height = 4, width = 12, units = "in", dpi = 600)

#### numerical identifiability of proportional hazard #######
# ---------------------------------------------------------#
betas01_t = c(X1 = -0.22, X2 = 0.11)
betas02_t = c(X1 =-0.33, X2 = 0.22)
betas12_t = c(X1 =-0.44, X2 = 0.33)

covs <- data.frame(id = 1:30000,
                   X1 = stats::rbinom(30000, 1L, 0.5),
                   X2 = rnorm(30000, 0, 1))

sim_wei <- rsimid(
  dist01 = "weibull",
  dist02 = "weibull",
  dist12 = "weibull",
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  gammas01 = gammas01_t,
  gammas02 = gammas02_t,
  gammas12 = gammas12_t,
  x = covs,
  cens = cens
)

tmat_mst <- trans.illdeath(names=c("diagnosis","relapse","death"))
sim_wei_mstate <- msprep(time=c(NA,"df_time","os_time"),
                         status=c(NA,"df_event","os_event"),
                         data = sim_wei,
                         trans=tmat_mst)
sim_wei_mstate <- dplyr::left_join(sim_wei_mstate,
                                   sim_wei[ , c("id", "X1", "X2")])


formula = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status) ~ X1 + X2) )

basehaz = lapply(1:3, function(x)
  "weibull")

sim_recov_2 <- mstte_stan(formula = formula,
                      data = sim_wei_mstate,
                      transition_labels = c("DP", "DX", "DPDX"),
                      basehaz = basehaz,
                      prior           = lapply(1:3, function(x)
                        rstanarm::normal() ),
                      prior_intercept = prior_intercept,
                      prior_aux       = prior_aux,
                      iter = 2000,
                      chains = 4,
                      control = list(adapt_delta = 0.99)
)

saveRDS(sim_recov_2, "~/rfactory/SympEng2019/sim_recov_2.RDS")

# ---------- Aplication to real dataset --------------------#

lusc = readr::read_tsv("~/rfactory/mstte-data/lusc/lusc_tcga_clinical_data.tsv")
na_count <-sapply(lusc, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

lusc = lusc %>%
  select(pt_id = `Patient ID`,
         age = `Diagnosis Age`,
         sex = `Sex`,
         history = `Prior Cancer Diagnosis Occurence`,
         smoke = `Patient Smoking History Category`,
         stage = `American Joint Committee on Cancer Tumor Stage Code`,
         metastasis = `American Joint Committee on Cancer Metastasis Stage Code`,
         lymph = `Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`,
         burden = `Fraction Genome Altered`,

         pd_y = `Disease Free (Months)`,
         pd_d = `Disease Free Status`,
         os_y = `Overall Survival (Months)`,
         os_d = `Overall Survival Status`)

summary(lusc)
lusc = lusc[!is.na(lusc$os_y), ]
na_count <-sapply(lusc, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

lusc$stage[lusc$stage == "T1b" | lusc$stage == "T1a"] <- "T1"
lusc$stage[lusc$stage == "T2b" | lusc$stage == "T2a"] <- "T2"
lusc$history[lusc$history == "Yes, History of Prior Malignancy" | lusc$history == "Yes, History of Synchronous/Bilateral Malignancy"] <- "Yes"

lusc$metastasis[lusc$metastasis == "M1b" | lusc$metastasis == "M1a"] <- "M1"
lusc$metastasis[lusc$metastasis == "M2b" | lusc$metastasis == "M2a"] <- "M2"

our_summary1 <-
  list("Diagnosis age" =
         list("min" = ~ min(.data$age, na.rm = T),
              "max" = ~ max(.data$age, na.rm = T),
              "mean (sd)" = ~ qwraps2::mean_sd(.data$age,  na_rm = T)),
       "Sex" =
         list("Male" = ~ qwraps2::n_perc0(.data$sex == "Male"),
              "Female"  = ~ qwraps2::n_perc0(.data$sex == "Female") ),
       "Previous malignancy" =
         list("Yes" = ~ qwraps2::n_perc0(.data$history == "Yes",  na_rm = T),
              "No"  = ~ qwraps2::n_perc0(.data$history == "No",  na_rm = T)),
       "Tumor Stage" =
         list("T1" = ~ qwraps2::n_perc0(.data$stage == "T1",  na_rm = T),
              "T2"  = ~ qwraps2::n_perc0(.data$stage == "T2",  na_rm = T),
              "T3"  = ~ qwraps2::n_perc0(.data$stage == "T3",  na_rm = T),
              "T4"  = ~ qwraps2::n_perc0(.data$stage == "T4",  na_rm = T)),
       "Metastasis Stage" =
         list("MX" = ~ qwraps2::n_perc0(.data$metastasis == "MX",  na_rm = T),
              "M0" = ~ qwraps2::n_perc0(.data$metastasis == "M0",  na_rm = T),
              "M1"  = ~ qwraps2::n_perc0(.data$metastasis == "M1",  na_rm = T)),
       "Lymph Node Stage" =
         list("NX" = ~ qwraps2::n_perc0(.data$lymph == "NX",  na_rm = T),
              "N0"  = ~ qwraps2::n_perc0(.data$lymph == "N0",  na_rm = T),
              "N1" = ~ qwraps2::n_perc0(.data$lymph == "N1",  na_rm = T),
              "N2"  = ~ qwraps2::n_perc0(.data$lymph == "N2",  na_rm = T),
              "N3"  = ~ qwraps2::n_perc0(.data$lymph == "N3",  na_rm = T)),
       "Mutation burden (% altered)" =
         list("min" = ~ min(.data$burden, na.rm = T),
              "median" = ~ median(.data$burden, na.rm = T),
              "max" = ~ max(.data$burden, na.rm = T),
              "mean (sd)" = ~ qwraps2::mean_sd(.data$burden,  na_rm = T))
  )
### Overall
whole <- summary_table(lusc, our_summary1)

library(kableExtra)
kable(whole)


library(mice)

lusc$burden[is.na(lusc$burden)] = mice.impute.mean(y = lusc$burden,
                       ry = !is.na(lusc$burden))

sum(is.na(lusc$burden))
library(caret)
lusc$pd_d[is.na(lusc$pd_d)] <- "DiseaseFree"
set.seed(1234)
train.index <- createDataPartition(lusc$pd_d == "Recurred/Progressed" & lusc$os_d == "DECEASED", p = 3/4, list = FALSE)
train <- lusc[ train.index,]
test  <- lusc[-train.index,]

# Suppose you have $d$ dimensional parameters, the optimal scale is approximate $2.4d^(−1/2)$ times the scale of the target distribution, which implies optimal acceptance rates of 0.44 for $d = 1$ and 0.23 for $d$ goes to \infinity.
# You are shooting for an acceptance rate of 25-50% for the Metropolis algorithm.
# reference: Automatic Step Size Selection in Random Walk Metropolis Algorithms, Todd L. Graves, 2011.

dp <- DirichletProcessBeta(train$burden, 1, mhStep = c(0.01, 0.01), mhDraws = 250)
# sAccept Ratio:  0.482 
dp <- Fit(dp, 24000, updatePrior = TRUE)
dp$numberClusters
plot(dp)

