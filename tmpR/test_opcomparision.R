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
                               f2 <- function(x){
                                 # int (I(x>0.5))^2 dx = 8/315
                                 x^2 *(1-x^2) * 315/8
                               }
                               f3 <- function(x){
                                 # int (x(1-x))^2 dx = 9/455
                                 x^3 *(1-x^3) * 455/9
                               }
                               f4 <- function(x){
                                 # int (x(1-x))^2 dx = 32/1989
                                 x^4 *(1-x^4) * 1989/32
                               }
                               # g1 <- function(x){
                               #   # int (x(1-x))^2 dx = ??
                               #   (cos(2*pi*x))^2
                               # }
                               1/700 * 1/4 *
                                 (f1(s)*f1(t)*f1(u)*f1(v) +
                                    f2(s)*f2(t)*f2(u)*f2(v) +
                                    f3(s)*f3(t)*f3(u)*f3(v) +
                                    f4(s)*f4(t)*f4(u)*f4(v)) },
                             delta = function(s,t){s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
# saveRDS(results_oneStep, 'C://Users//jerem//Downloads//op_onestep.rds')
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
#                              f2 <- function(x){
#                                # int (I(x>0.5))^2 dx = 8/315
#                                x^2 *(1-x^2) * 315/8
#                              }
#                              f3 <- function(x){
#                                # int (x(1-x))^2 dx = 9/455
#                                x^3 *(1-x^3) * 455/9
#                              }
#                              f4 <- function(x){
#                                # int (x(1-x))^2 dx = 32/1989
#                                x^4 *(1-x^4) * 1989/32
#                              }
#                              # g1 <- function(x){
#                              #   # int (x(1-x))^2 dx = ??
#                              #   (cos(2*pi*x))^2
#                              # }
#                              1/650 * 1/4 *
#                                (f1(s)*f1(t)*f1(u)*f1(v) +
#                                   f2(s)*f2(t)*f2(u)*f2(v) +
#                                   f3(s)*f3(t)*f3(u)*f3(v) +
#                                   f4(s)*f4(t)*f4(u)*f4(v)) },
#                            delta = function(s,t){s*(1-s)*t*(1-t) },
#                            resolution=20)
# )
round(colMeans(results_oneStep,na.rm = T),2)
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

