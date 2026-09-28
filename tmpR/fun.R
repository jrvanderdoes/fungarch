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
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               f1 <- function(x){
                                 # int (x(1-x))^2 dx = 1/30
                                 x *(1-x) * 30
                               }
                               1/400 *
                                 (f1(s)*f1(t)*f1(u)*f1(v) ) },
                             delta = function(s,t){10*s*(1-s)*t*(1-t) },
                             resolution=20)
  )


}
# saveRDS(results_oneStep,
#         'C://Users//jerem//Downloads//ShowFigures//op_onestep.rds')
results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

##    ROUND 1
# tmp <- suppressMessages(
#   generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
#                            alpha = function(s,t,u,v){
#                              f1 <- function(x){
#                                # int (x(1-x))^2 dx = 1/30
#                                x *(1-x) * 30
#                              }
#                              1/650 *
#                                (f1(s)*f1(t)*f1(u)*f1(v) ) },
#                            delta = function(s,t){s*(1-s)*t*(1-t) },
#                            resolution=20)
# )
## Round 2
# tmp <- suppressMessages(
#   generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
#                            alpha = function(s,t,u,v){
#                              f1 <- function(x){
#                                # int (x(1-x))^2 dx = 1/30
#                                x *(1-x) * 30
#                              }
#                              1/600 *
#                                (f1(s)*f1(t)*f1(u)*f1(v) ) },
#                            delta = function(s,t){10*s*(1-s)*t*(1-t) },
#                            resolution=20)
# )
## Round 3
# tmp <- suppressMessages(
#   generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
#                            alpha = function(s,t,u,v){
#                              f1 <- function(x){
#                                # int (x(1-x))^2 dx = 1/30
#                                x *(1-x) * 30
#                              }
#                              1/750 *
#                                (f1(s)*f1(t)*f1(u)*f1(v) ) },
#                            delta = function(s,t){10*s*(1-s)*t*(1-t) },
#                            resolution=20)
# )

round(colMeans(results_oneStep,na.rm = T),2)
#   OR        PW        OP       MLE     Mean     Prev
table(apply(results_oneStep,1,function(x){which.min(x)}))
#   OR      PW   OP     MLE     Mean    Prev
round(rowMeans(apply(results_oneStep,1,functio

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

