rstantools::rstan_package_skeleton("~/rfactory/mstte")


lusc = readr::read_tsv("/home/mtr/Downloads/lusc_tcga_clinical_data.tsv")
density(lusc$`Fraction Genome Altered`[!is.na(lusc$`Fraction Genome Altered`)])
hist(lusc$`Fraction Genome Altered`[!is.na(lusc$`Fraction Genome Altered`)])
y <- lusc$`Fraction Genome Altered`[!is.na(lusc$`Fraction Genome Altered`)]
library(dirichletprocess)
library(ggplot2)
dp <- DirichletProcessBeta(y, 1)
dp <- Fit(dp, 2000, updatePrior = TRUE)
plot(dp)




library(car)
Y <- car::logit(y)
Y <- scale(Y)
dp_normal <- DirichletProcessGaussian(Y)
dp_normal <- Fit(dp_normal, 10000)
plot(dp)
dp$numberClusters




posteriorFrame <- PosteriorFrame(dp, ppoints(100), ci_size = 0.05)
trueFrame <- data.frame(x=ppoints(100),
                        y=

ggplot() +
  geom_ribbon(data=posteriorFrame,
               aes(x=x, ymin=X2.5., ymax=X97.5.),
               alpha=0.2,
              colour=NA,
               fill="red") +
 geom_line(data=posteriorFrame, aes(x=x, y=Mean), colour="red") +
 geom_line(data=trueFrame, aes(x=x, y=y))
