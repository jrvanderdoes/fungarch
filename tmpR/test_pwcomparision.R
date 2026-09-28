nSims <- 100

results_oneStep <- results_fromLast <-
  data.frame('OR'=rep(NA,nSims),'PW'=NA, 'OP'=NA,
             'MLE'=NA, 'Mean'=NA,'Prev'=NA)

cat('Fit (',nSims,'): ')
set.seed(123456789)
for(i in 1:nSims){
  cat(paste0(i,', '))

  tmp <- suppressMessages(
    generateFitComparePWARCH(N_fit = 500, N_pred = 250,
                             alpha = function(t,s){ 20 * t * (1-t) * s * (1-s) },
                             delta = function(t){ rep(0.5,length(t)) })
  )

  results_oneStep[i,] <- tmp
}
# saveRDS(results_oneStep,
#         'C://Users//jerem//Downloads//ShowFigures//pw_onestep.rds')


colMeans(results_oneStep, na.rm=T)
#   OR       PW       OP      MLE     Mean     Prev
# t(1-t)s(1-s)
# 4.419488 4.499896 5.211625 4.441921 4.761606 5.146680
# (1-cos(2pi t)^2)(1-cos(2pi s)^2)
# 169.2658 171.9417 12359.9025 170.3583 175.1066 205.3207
table(apply(results_oneStep,1,function(x){which.min(x)}))
#   OR       PW       OP      MLE     Mean     Prev
#   19       2        463     16      0         0
#   263      68       0       145     24        0
rowMeans(apply(results_oneStep,1,function(x){rank(x)}))
#   OR    PW    OP   MLE  Mean  Prev
#  2.824 3.352 1.266 2.892 4.850 5.816
#  1.886 2.540 5.998 2.126 3.456 4.994

#######################################################################

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
