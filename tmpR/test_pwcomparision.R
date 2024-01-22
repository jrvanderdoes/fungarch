nSims <- 500

results_oneStep <- results_fromLast <-
  data.frame('OR'=rep(NA,nSims),'PW'=NA, 'OP'=NA,
             'MLE'=NA, 'Mean'=NA,'Prev'=NA)

cat('Fit (',nSims,'): ')
set.seed(123456789)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # tmp <- suppressMessages(
  #   generateFitComparePWARCH(N_fit = 500, N_pred = 250,
  #                            alpha = function(t,s){ 12 * t * (1-t) * s * (1-s) },
  #                            delta = function(t){ rep(0.5,length(t)) })
  # )
  tmp <- suppressMessages(
    generateFitComparePWARCH(N_fit = 500, N_pred = 250,
                             alpha = function(t,s){ (1-cos(2*pi*t)^2) * (1-cos(2*pi*s)^2) },
                             delta = function(t){ rep(0.5,length(t)) })
  )

  results_oneStep[i,] <- tmp
}
# saveRDS(results_oneStep,
#         'C://Users//jerem//Downloads//ShowFigures//pw_onestep.rds')


colMeans(results_oneStep)
#   OR       PW       OP      MLE     Mean     Prev
# 188.1037 190.0013 179.4277 188.4361 196.9432 223.4342
table(apply(results_oneStep,1,function(x){which.min(x)}))
#   OR       PW       OP      MLE     Mean     Prev
#   12       28       428     30      2         0
rowMeans(apply(results_oneStep,1,function(x){rank(x)}))
#   OR    PW    OP   MLE  Mean  Prev
#  2.914 3.272 1.344 2.754 4.742 5.974

# results_oneStep
library(tidyverse)
library(ggplot2)
results_oneStep_plot <- results_oneStep %>%
  pivot_longer(cols = OR:Prev)
ggplot(results_oneStep_plot,
       aes(value,color=name,fill=name,group=name)) +
  geom_histogram(aes(y=after_stat(density)),
                 bins=75,alpha=0.25, position="identity") +  # scale histogram y
  geom_density(alpha=0.05,linewidth=1) +
  theme_bw() +
  ggplot2::scale_color_manual(
    values=c("#000000", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2")) +
  ggplot2::scale_fill_manual(
    values=c("#000000", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2"))
