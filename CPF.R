#' Function for creating CPF (conditional probability function) graphs
#' with the uncertainty intervals.
#' 
#' The user input is a data.frame containing wind direction, wind speed and
#' pollutant concentration data. The wind is split into sectors based on
#' direction and speed, according to "sectors" value set by the user. 
#' 
#' CPF value is calculated for each wind sector as m/n, where n stands for total
#' number of observations in chosen wind sector, while m stands for the number
#' of wind observations above chosen percentile of pollutant concentrations, for
#' the given wind sector. CPF values are shown as black arches, while CPF
#' threshold is shown as a red circle.
#' 
#' The user can choose between uncertainties calculated from binomial ratio and
#' from bootstrapping, by setting the variable "intervals". If bootstrapping is
#' chosen, the number of bootstrap replicates "R" can also be set (default is
#' 1000).
#' 
#' The function returns a CPF plot with uncertainty intervals, if plot = TRUE
#' (default). In the other case, the function returns a data.frame containing
#' wind directions (wd), number of total observations (n), number of
#' observations above threshold (m), CPF values, lower and upper limits for each
#' wind segment.
#'
#' @param data data.frame including wind direction, wind speed and pollutant
#'   concentration columns; mandatory
#'   
#' @param sectors number of wind direction sectors, default is 16
#' 
#' @param ws name for wind speed column in data; default is "ws"
#' 
#' @param wd name for wind direction column in data; default is "wd"
#' 
#' @param pollutant name for pollutant concentrations in data; mandatory
#' 
#' @param percentile percentile for threshold; defaults is 75
#' 
#' @param conf confidence level for calculating confidence intervals; default is
#'   0.68
#'   
#' @param intervals method for interval calculation, c("br", "bootstrap"), "br"
#'   meaning binomial ratio; default is "br"
#'   
#' @param R number of bootstrap replicates, in case that bootstrap confidence
#'   intervals are chosen; default is 1000
#'   
#' @param remove_speed remove speed below this value; default is 1
#' 
#' @param plot if TRUE, a plot is created; if FALSE, a table is produced
#'   containing wind directions (wd), number of total observations (n), number
#'   of observations above threshold (m), CPF values, lower and upper limits for
#'   each wind direction; default is TRUE
#'   
#' @return if plot==TRUE: CPF plot with confidence intervals; if plot==FALSE:
#'   table containing wind directions (wd), number of total observations (n),
#'   number of observations above threshold (m), CPF values, lower and upper
#'   limits for each wind direction
#'   
#' @export
#' 
#' @examples
CPF <- function(data, sectors=16, ws='ws', wd='wd', 
                pollutant=NULL, percentile=75, 
                conf=0.68, intervals="br", R=1000,
                remove_speed=1, plot=TRUE){
  
  # removing unnecessary columns and missing data (NA):
  data <- data[c(ws, wd, pollutant)]
  data <- na.omit(data)
  
  colnames(data) <- c("ws", "wd", "pollutant")
  
  # removing missing data (-999 or alike) and speeds less than 
  # remove_speed value:
  data <- data[data$ws > remove_speed & data$wd >= 0 & data$pollutant >= 0, ]
  
  # angle for wind sectors:
  angle = 360/sectors
  
  # threshold, calculated from the chosen percentile:
  threshold <- quantile(data$pollutant, probs=percentile/100)
  
  
  # wind directions, starting at 0; their number is set with sectors:
  CPF <- data.frame(wd=seq(0, 360-angle, by=angle))
  
  # all wind directions from data are transformed into wind directions from 
  # above
  data$wd[data$wd <= (angle/2)] <- 0
  data$wd[data$wd > 360-(angle/2)] <- 0
  
  for (i in CPF$wd) {
    data$wd[data$wd <= i+(angle/2) & data$wd > i-(angle/2)] <- i
  }
  
  
  # number of events in each sector:
  for (i in 1:nrow(CPF)) {
    CPF$n[i] = nrow(data[data$wd==CPF$wd[i],])
  }
  
  # data above the threshold:
  data_t <- data[data$pollutant>threshold,]
  
  # number of events above threshold in each sector:
  for (i in 1:nrow(CPF)) {
    CPF$m[i] = nrow(data_t[data_t$wd==CPF$wd[i],])
  }
  
  # CPF for each direction:
  CPF$CPF <- CPF$m/CPF$n
  
  # total number of samples:
  N_samples <- nrow(data)
  
  # total number of events above threshold:
  N_t <- nrow(data_t)
  
  # confidence intervals:
  
  # binomial ratio:
  if(intervals == "br"){
    # phi:
    phi <- data.frame(wd=CPF$wd, x1=CPF$m, n1=N_t, 
                      x2=CPF$n-CPF$m, n2=N_samples-N_t)
    
    # confidence intervals for phi, according to Koopman (1984):
    CI <- as.data.frame(BinomRatioCI(phi$x1, phi$n1, phi$x2, phi$n2, 
                                     method = "koop", conf.level = conf))
    phi$phi <- CI$est
    phi$lower <- CI$lwr.ci
    phi$upper <- CI$upr.ci
    
    # confidence intervals for CPF:
    CPF$lower <- phi$lower / (phi$lower + (N_samples-N_t)/N_t)
    CPF$upper <- phi$upper / (phi$upper + (N_samples-N_t)/N_t)
    CPF$upper <- replace(CPF$upper, is.nan(CPF$upper), 1)
  }
  
  # bootstrap:
  else{
    
    CPF_boot <- function(data, indices){
      
      CPF_boot <- data.frame(wd=CPF$wd)
      # drawing samples from data:
      data2 <- data[indices,]
      
      # number of events in each direction:
      for (i in 1:nrow(CPF_boot)) {
        
        CPF_boot$n[i] = nrow(data2[data2$wd==CPF_boot$wd[i],])
      }
      
      # threshold:
      threshold <- quantile(data2$pollutant, probs=percentile/100)
      
      data_t <- data2[data2$pollutant>threshold,]
      
      # number of events above threshold in each direction:
      for (i in 1:nrow(CPF_boot)) {
        
        CPF_boot$m[i] = nrow(data_t[data_t$wd==CPF_boot$wd[i],])
      }
      
      # CPF for each direction:
      CPF_boot$CPF <- CPF_boot$m/CPF_boot$n
      
      return(CPF_boot$CPF)
    }
    
    # bootstrap results:
    boot_results <- boot(data, CPF_boot, R=R)
    
    # bootstrap confidence intervals:
    
    CI <- data.frame(boot_results$t0)
    
    for (i in 1:length(boot_results$t0)) {
      tryCatch(
        {
          CI$lower[i] <- boot.ci(boot_results, conf=conf, type="perc", index=i)$percent[1,4]
        },
        error = function(e){
          message("no interval calculated")
          print(e)
          return(NA)
        }
        
      )
    }
    
    for (i in 1:length(boot_results$t0)) {
      tryCatch(
        {
          CI$upper[i] <- boot.ci(boot_results, conf=conf, type="perc", index=i)$percent[1,5]
        },
        error = function(e){
          message("no interval calculated")
          print(e)
          return(NA)
        }
        
      )
    }
    
    # correcting confidence intervals outside [0,1] interval:
    CI$lower <- replace(CI$lower, CI$lower<0, 0)
    CI$upper <- replace(CI$upper, CI$upper>1, 1)
    
    CPF$lower <- CI$lower
    CPF$upper <- CI$upper
    
  }
  

  # CPF threshold (0.25 for 75th percentile):
  a=(100-percentile)/100
  aa <- data.frame(x=-angle/2, y=a, xend=360-(angle/2), yend=a)
  
  # ticks on y axis:
  my_breaks <- seq(0, max(CPF$upper), by=0.1)
  
  # label for the confidence interval:
  conf_perc <- conf*100
  
  # plot:
  p <- ggplot(CPF)+
    # lower interval:
    geom_rect(aes(xmin=wd-(angle/2), xmax=wd+(angle/2), 
                  ymin=lower, ymax=CPF, fill=paste(conf_perc, "%", " CI",
                                                   sep="")))+
    # upper interval:
    geom_rect(aes(xmin=wd-(angle/2), xmax=wd+(angle/2), 
                  ymin=CPF, ymax=upper), fill="grey")+
    # a:
    geom_segment(data=aa, aes(x=x, y=y, xend=xend,
                              yend=yend, color=as.character(a)))+
    # CPF:
    geom_segment(aes(x=wd-(angle/2), y=CPF, xend=wd+(angle/2),
                     yend=CPF, color="CPF"))+
    coord_polar(start=-angle*pi/360)+
    # broj događaja u segmentu:
    #annotate("text", x=CPF$wd-angle/2, y=0.6, label=CPF$n)+
    scale_x_continuous(breaks=seq(0, 270, by=90),
                       labels=c("N","W",
                                "S","E"),
                       limits=c(-angle/2,(360-(angle/2))))+
    scale_y_continuous(breaks=my_breaks, limits=c(0,NA))+
    scale_color_manual(name="",breaks=c(as.character(a), "CPF"),
                       values=c("red", "black"))+
    scale_fill_manual(name="",values=c("grey"))+
    labs(x="", y="")+
    theme_minimal()+
    annotate('text', x =315, y = my_breaks, label = my_breaks, size = 9*5/14) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())

  # the function returns the plot or the table, as chosen by the user:
  if(plot==TRUE)  return(p)
  else return(CPF)
}
