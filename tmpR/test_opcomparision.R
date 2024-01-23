nSims <- 500

results_oneStep <- results_fromLast <-
  data.frame('OR'=rep(NA,nSims),
             'PW'=rep(NA,nSims), 'OP'=NA,
             'MLE'=NA, 'Mean'=NA,'Prev'=NA)

cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250)
  )

  results_oneStep[i,] <- tmp
}
# saveRDS(results_oneStep,
#         'C://Users//jerem//Downloads//ShowFigures//op_onestep.rds')


colMeans(results_oneStep)
#   OR       PW       OP      MLE     Mean     Prev
# 168.4386 221.4287 163.6856 255.4684 179.1755 196.4754
table(apply(results_oneStep,1,function(x){which.min(x)}))
#   OR       PW       OP      MLE     Mean     Prev
#   63      17        382     0       36        2
rowMeans(apply(results_oneStep,1,function(x){rank(x)}))
#   OR    PW    OP   MLE  Mean  Prev
# 2.082 4.556 1.360 5.710 2.928 4.364

# results_oneStep
library(tidyverse)
library(ggplot2)
results_oneStep_plot <- results_oneStep %>%
  pivot_longer(cols = PW:Prev)
ggplot(results_oneStep_plot,
       aes(value,color=name,fill=name,group=name)) +
  geom_histogram(aes(y=after_stat(density)),
                 bins=75,alpha=0.25, position="identity") +  # scale histogram y
  geom_density(alpha=0.05,linewidth=1) +
  theme_bw() +
  ggplot2::scale_color_manual(
    values=c(#"#000000",
             "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2")) +
  ggplot2::scale_fill_manual(
    values=c(#"#000000",
             "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2"))

