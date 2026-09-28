nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPGARCH(N_fit = 500, N_pred = 250,
                              alpha = function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) },
                              beta = function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v)},
                              delta = function(s,t){ s*(1-s)*t*(1-t) },
                              resolution=20)
  )


  results_oneStep[i,] <- tmp
}
# saveRDS(results_oneStep, 'C://Users//jerem//Downloads//op_onestep.rds')
results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

##    ROUND 1
round(colMeans(results_oneStep,na.rm = T),3)
#   OR      OP      PW     MLE     Mean    Prev    Zero
apply(results_oneStep,MARGIN = 2,
      function(X){round(mean(X,trim = 0.05),3)})
#   OR      OP      PW     MLE     Mean    Prev    Zero
apply(results_oneStep,MARGIN = 2,
      function(X){round(median(X),3)})
#   OR      OP      PW     MLE     Mean    Prev    Zero
table(apply(results_oneStep,1,function(x){which.min(x)}))
#   OR    OP    PW    MLE   Mean    Prev    Zero
round(rowMeans(apply(results_oneStep,1,function(x){rank(x)})),3)
#   OR      OP   PW     MLE  Mean  Prev   Zero

# results_oneStep
library(tidyverse)
library(ggplot2)
results_oneStep_plot <- results_oneStep %>%
  pivot_longer(cols = OR:Zero)
ggplot(results_oneStep_plot,
       aes(x=name,y=value,color=name,fill=name,group=name)) +
  # geom_histogram(aes(y=after_stat(density)),
  #                bins=75,alpha=0.25, position="identity") +  # scale histogram y
  # geom_density(alpha=0.05,linewidth=1) +
  geom_boxplot(alpha=0.05) +
  stat_summary(geom = "errorbar", #mapping=aes(),
               fun.min = mean, fun = mean,
               fun.max = mean, linewidth = 0.75,
               linetype='dotted') +
  theme_bw() +
  scale_y_log10() +
  ggplot2::scale_color_manual(
    values=c(scales::hue_pal()(7))) +
  ggplot2::scale_fill_manual(
    values=scales::hue_pal()(7))

