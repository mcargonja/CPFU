#' Function for plotting bivariate CPF graph with separate graphs with
#' confidence intervals for each wind direction
#'
#' The user input is a data.frame containing wind direction, wind speed and
#' pollutant concentration data. The wind is split into sectors based on
#' direction and speed, according to sectors and speed_interval values set by
#' the user. 
#' 
#' CPF value is calculated for each wind sector as m/n, where n stands for total
#' number of observations in chosen wind sector, while m stands for the number
#' of wind observations above chosen percentile of pollutant concentrations, for
#' the given wind sector. CPF values are shown in color between blue (CPF = 0)
#' and red (CPF = 1), where grey is for CPF = (100 - percentile)/100 (eg. CPF =
#' 0.25 for percentile = 75).
#' 
#' The user can choose between uncertainties calculated from binomial ratio and
#' from bootstrapping, by setting the variable "intervals". If bootstrapping is
#' chosen, the number of bootstrap replicates "R" can also be set (default is
#' 1000).
#' 
#' The function returns a bivariate CPF plot with separate plots with
#' uncertainty intervals, if plot = TRUE (default). In the other case, the
#' function returns a data.frame containing wind directions (wd), number of
#' total observations (n), number of observations above threshold (m), CPF
#' values, lower and upper limits for each wind segment.
#'
#' @param data data.frame including wind direction, wind speed and pollutant
#'   concentration columns; mandatory
#'   
#' @param sectors Number of wind direction sectors; default is 16
#' 
#' @param ws Name for wind speed column in data; default is "ws"
#' 
#' @param wd Name for wind direction column in data; default is "wd"
#' 
#' @param pollutant Name for pollutant concentrations in data; mandatory
#' 
#' @param percentile Percentile for threshold; default is 75
#' 
#' @param remove_speed Remove speed below this value; default is 1
#' 
#' @param speed_interval Interval for creating wind speed segments; default is 1
#' 
#' @param speed_unit Measuring unit for wind speed; default is "m/s"
#' 
#' @param conf Confidence level; default is 0.68
#' 
#' @param intervals Method for interval calculation, c("br", "bootstrap"), "br"
#'   meaning binomial ratio; default is "br"
#'   
#' @param R Number of bootstrap replicates; default is 1000
#' 
#' @param plot If TRUE, a plot is created; if FALSE, a table is produced
#'   containing wind directions (wd), number of total observations (n), number
#'   of observations above threshold (m), CPF values, lower and upper limits for
#'   each wind segment; default is TRUE
#'   
#' @return If plot = TRUE: CPF plot with confidence intervals; if plot = FALSE:
#'   data.frame containing wind directions (wd), wind speed (ws), number of total
#'   observations (n), number of observations above threshold (m), CPF values,
#'   lower and upper limits for each wind segment
#'   
#' @export
#'
#' @examples
CBPF_intervals <- function(data, sectors=16, ws='ws', wd='wd', 
                           pollutant=NULL, percentile=75, 
                           remove_speed=1, speed_interval=1, 
                           speed_unit="m/s", conf=0.68, 
                           intervals="br", R=1000, plot=TRUE){
  
  # removing unnecessary columns and missing data (NA):
  data <- data[c(ws, wd, pollutant)]
  data <- na.omit(data)
  colnames(data) <- c("ws", "wd", "pollutant")
  
  # removing missing data (-999 or alike) and speeds less than 
  # remove_speed value:
  data <- data[data$ws > remove_speed & data$wd >= 0 & data$pollutant >= 0, ]
  
  # total number of samples:
  N_samples <- nrow(data)
  
  # angle for wind sectors:
  angle = 360/sectors
  
  
  
  # threshold, calculated from the chosen percentile:
  threshold <- quantile(data$pollutant, probs=percentile/100)
  
  
  # wind directions, starting at 0; their number is set with sectors:
  CPF_all_directions <- data.frame(wd=seq(0, 360-angle, by=angle))
  
  
  # all wind directions from data are transformed into wind directions from 
  # above
  data$wd[data$wd <= (angle/2)] <- 0
  data$wd[data$wd > 360-(angle/2)] <- 0
  
  for (i in CPF_all_directions$wd) {
    data$wd[data$wd <= i+(angle/2) & data$wd > i-(angle/2)] <- i
  }
  # wind speed is segmented according to speed_interval:
  data$ws <- ceiling(data$ws/speed_interval)*speed_interval
  
  # finding sections with data:
  CPF <- unique(data[c("wd", "ws")])
  
  # number of events in each sector:
  for (i in 1:nrow(CPF)) {
    CPF$n[i] = nrow(data[data$wd==CPF$wd[i] & data$ws==CPF$ws[i],])
  }
  
  # data above the threshold:
  data_t <- data[data$pollutant>threshold,]
  
  # total number of events above threshold:
  N_t <- nrow(data_t)
  
  # number of events above threshold in each sector:
  for (i in 1:nrow(CPF)) {
    CPF$m[i] = nrow(data_t[data_t$wd==CPF$wd[i] & data_t$ws==CPF$ws[i],])
  }
  
  # CPF:
  CPF$CPF <- CPF$m/CPF$n
  
  # confidence intervals:
  
  # binomial ratio:
  if(intervals == "br"){
    # phi:
    phi <- data.frame(wd=CPF$wd, ws=CPF$ws, x1=CPF$m, n1=N_t, 
                      x2=CPF$n-CPF$m, n2=N_samples-N_t)
    
    # confidence intervals for phi, according to Koopman (1984):
    CI <- as.data.frame(BinomRatioCI(phi$x1, phi$n1, phi$x2, phi$n2, method = "koop"))
    phi$phi <- CI$est
    phi$lower <- CI$lwr.ci
    phi$upper <- CI$upr.ci
    
    # confidence intervals for CBPF:
    CPF$lower <- phi$lower / (phi$lower + (N_samples-N_t)/N_t)
    CPF$upper <- phi$upper / (phi$upper + (N_samples-N_t)/N_t)
    CPF$upper <- replace(CPF$upper, is.nan(CPF$upper), 1)
  }
  
  # bootstrap:
  else {
    
    CPF_boot <- function(data, indices){
      
      CPF_boot <- data.frame(wd=CPF$wd, ws=CPF$ws)
      # drawing samples from data:
      data2 <- data[indices,]
      
      # number of events in each sector:
      for (i in 1:nrow(CPF_boot)) {
        CPF_boot$n[i] = nrow(data2[data2$wd==CPF_boot$wd[i] & data2$ws==CPF_boot$ws[i],])
      }
      
      threshold <- quantile(data2$pollutant, probs=percentile/100)
      data_t <- data2[data2$pollutant>=threshold,]
      
      # number of events above threshold in each sector:
      for (i in 1:nrow(CPF_boot)) {
        CPF_boot$m[i] = nrow(data_t[data_t$wd==CPF_boot$wd[i] & data_t$ws==CPF_boot$ws[i],])
      }
      
      
      # CPF for each sector:
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
          CI$lower[i] <- boot.ci(boot_results, conf=conf, type="perc", 
                                 index=i)$percent[1,4]
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
          CI$upper[i] <- boot.ci(boot_results, conf=conf, type="perc", 
                                 index=i)$percent[1,5]
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
  
  # ticks on y axis:
  my_breaks <- seq(0, max(CPF$ws)+2, by=speed_interval*2)
  
  # facet labels:
  CPF$labels <- factor(paste(CPF$wd, "°", sep=""),
                       levels = paste(sort(unique(CPF$wd)), "°", sep=""))
  
  # label for the confidence interval:
  conf_perc <- conf*100
  
  # average value of CPF (0.25 for 75th percentile):
  a=(100-percentile)/100
  aa <- data.frame(x=min(data$ws), y=a, xend=max(data$ws), yend=a)
  
  # plot:
  # bivariate CPF plot:
  p1 <- ggplot(data=CPF)+
    geom_rect(aes(xmin=wd-(angle/2), xmax=wd+(angle/2), 
                  ymin=ws-speed_interval, ymax=ws, fill=CPF))+
    coord_polar(start=-angle*pi/360)+
    scale_x_continuous(breaks=seq(0,315,by=45),
                       labels=c(0,45,90,135,180,225,270,
                                paste("v/",speed_unit,sep="")),
                       limits=c(-angle/2,(360-(angle/2))))+
    scale_y_continuous(breaks=my_breaks, limits=c(0,NA))+
    scale_fill_gradientn(
      colors=c("blue", "grey", "firebrick1", "red4"),
      limits=c(0,1),
      values=c(0,a,2*a,1))+
    theme_minimal()+
    annotate('text', x =315, y = my_breaks, label = my_breaks, 
             size = 3)+ 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())+
    labs(x="", y="")
  
  # subplots with intervals:
  p2 <- ggplot(data=CPF)+
    # intervals:
    geom_rect(aes(xmin=ws-(speed_interval*0.45), ymin=lower, 
                  xmax=ws+(speed_interval*0.45), ymax=upper, 
                  fill=paste(conf_perc, "%", " CI",
                             sep="")))+
    # a:
    geom_segment(data=aa, aes(x=x-speed_interval, y=y, 
                              xend=xend+speed_interval, yend=yend, 
                              color=as.character(a)))+
    # CPF values:
    geom_segment(aes(x=ws-(speed_interval*0.45), y=CPF, 
                     xend=ws+(speed_interval*0.45), yend=CPF, color="CPF"))+
    facet_wrap(~ labels)+
    scale_color_manual(name="",values=c("black", "red"),
                       breaks=c("CPF", as.character(a)))+
    scale_fill_manual(name="",values=c("grey"), breaks=c(paste(conf_perc, "%", " CI",
                                                                    sep="")))+
    labs(x=paste("v/", speed_unit, sep=""), y="")+
    theme_minimal()+
    theme(legend.position="bottom")
  
  # wind directions with no data are also included in the final table 
  # (with CPF=NA):
  CPF <- left_join(CPF_all_directions, CPF, by=wd)
  
  # ordering the data:
  CPF <- CPF[order(CPF$wd, CPF$ws),]
  rownames(CPF) <- NULL

  # the function returns the plot or the table, as chosen by the user:
  if(plot==TRUE){
    return(ggarrange(p1, p2, ncol=1, heights=c(2, 360/(4*angle))))
  } 
  else return(CPF[,1:7])
}
