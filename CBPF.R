#' Function for plotting the bivariate conditional
#' probability function graph, without the confidence intervals.
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
#' The function returns a bivariate CPF plot, if plot = TRUE (default). In the
#' other case, the function returns a data.frame containing wind directions
#' (wd), number of total observations (n), number of observations above
#' threshold (m), and CPF values for each wind segment.
#' 
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
#' @param percentile Percentile for threshold; defaults is 75
#' 
#' @param remove_speed Remove speed below this value; default is 1
#' 
#' @param speed_interval Interval for creating wind speed sections; default is 1
#' 
#' @param speed_unit Measuring unit for wind speed; default is "m/s"
#' 
#' @param plot If TRUE, a plot is created; if FALSE, a table is produced
#'   containing wind directions (wd), number of total observations (n), number
#'   of observations above threshold (m), CPF values, lower and upper limits for
#'   each wind direction; default is TRUE
#'   
#' @return If plot = TRUE: CPF plot with confidence intervals; if plot = FALSE:
#'   table containing wind directions (wd), wind speed (ws), number of total
#'   observations (n), number of observations above threshold (m), CPF values
#'   
#' @export
#'
#' @examples
CBPF <- function(data, sectors=16, ws='ws', wd='wd', 
                 pollutant=NULL, percentile=75, remove_speed=1,
                 speed_interval=1, speed_unit="m/s", plot=TRUE){
  
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
  
  # finding wind sectors with data:
  CPF <- unique(data[c("wd", "ws")])
  
  # number of events in each sector:
  for (i in 1:nrow(CPF)) {
    CPF$n[i] = nrow(data[data$wd==CPF$wd[i] & data$ws==CPF$ws[i],])
  }
  
  # data above the threshold:
  data_t <- data[data$pollutant>threshold,]
  
  # number of events above threshold in each sector:
  for (i in 1:nrow(CPF)) {
    CPF$m[i] = nrow(data_t[data_t$wd==CPF$wd[i] & data_t$ws==CPF$ws[i],])
  }
  
  # CPF:
  CPF$CPF <- CPF$m/CPF$n
  
  # average value of CPF (0.25 for 75th percentile):
  a=(100-percentile)/100
  
  # ticks on y axis:
  my_breaks <- seq(0, max(CPF$ws)+2, by=speed_interval*2)
  
  # plot:
  p <- ggplot(data=CPF)+
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
  
  # wind directions with no data are also included in the final table 
  # (with CPF=NA):
  CPF <- left_join(CPF_all_directions, CPF, by=wd)
  
  # ordering the data:
  CPF <- CPF[order(CPF$wd, CPF$ws),]
  rownames(CPF) <- NULL

  # the function returns the plot or the table, as chosen by the user:
  if(plot==TRUE) return(p)
  else return(CPF)
}
