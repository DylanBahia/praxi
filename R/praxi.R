#' @importFrom Rcpp evalCpp
#' @import ggplot2
NULL

.praxi.class <- setClass("praxi.class",representation(y="numeric",p="numeric",b="numeric",cost="numeric",res="matrix"))

praxi.class <- function(y,p,b,cost,res)
{
	.praxi.class(y=y,p=p,b=b,cost=cost,res=res)	
}

#' Detects point and collective anomalies in correlated time series data.
#'
#' @name ouralg
#'
#' @description Runs ouralg on a vector of observations.
#'
#' @param y A numeric vector of observations..
#' @param x A positive integer containing the order of the AR model.
#' @param b A positive numeric constant for the penalty. 
#'
#' @return An instance of an S4 class of type \code{\link{praxi}}.
#'
#' @rdname praxi
#'
#' @examples
#' library(praxi)
#' 
#'#simulate an AR(1) time series with two changes in mean
#'
#'gam <- c(0.7); n <- 10000
#'y <- arima.sim(list(ar=gam),n=n)
#'y[300] <- y[300]+10
#'y[5000:5100] <- y[5000:5100]+2
#'
#'#using custom penalty
#'res <- praxi(y=y,p=length(gam),b=5*log(n))
#'#using default penalty
#'res <- praxi(y=y,p=length(gam))
#'anoms <- anomalies(res)
#'pt <- plot(res)
#'print(res)
#'
#'
#' @export

praxi <- function(y,p,b=NULL)
{
  
  if(is.null(b)){
    b <- 4*log(length(y))
  }
  
	result <- ar_alg_call(y,p,b)
	rlist <- praxi.class(y,p,b,result[[1]],result[[2]])
	return(rlist)
}


#' Produces a plot showing anomalies detected by ouralg.
#'
#' @name plot
#'
#' @description A plot method for an S4 object returned by \code{\link{praxi}}.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{praxi}}.
#'
#' @return A ggplot object. 
#'
#' @rdname plot-methods
#'
#' @aliases plot,praxi.class-method
#' 
#' @export

setMethod("plot",signature=list("praxi.class"),function(x)
{
  
  df <- data.frame(t=1:length(x@y),y=x@y,panom=rep(0,length(x@y)))
  
  coll_times <- data.frame(st=numeric(),nd=numeric())
  
  matr <- unname(x@res)
  
  if(nrow(matr)!=0){
    for(i in 1:nrow(matr)){
      curr_row <- matr[i,1:2]
      if(curr_row[1]==curr_row[2]){
        df$panom[curr_row[1]] <- 1
      }else{
        coll_times[nrow(coll_times)+1,] <- c(curr_row[1],curr_row[2])
      }
    }
  }
  
  out <- ggplot2::ggplot(df,ggplot2::aes(t,y))+ggplot2::geom_point(ggplot2::aes(color=factor(panom)),size=1.5)+
    ggplot2::scale_color_manual(values=c("0"="black","1"="red"))+ggplot2::xlab("Index")+ggplot2::ylab("Value")+ggplot2::theme(legend.position = "none")
  
  if(nrow(coll_times)!=0){
    out <- out+ggplot2::geom_rect(data=coll_times,ggplot2::aes(xmin=st,xmax=nd,ymin=-Inf,ymax=Inf),inherit.aes = FALSE,fill="blue",alpha=0.3)
  }
  
  return(out)
})

#' Produces a plot showing anomalies detected by ouralg.
#'
#' @name summary
#'
#' @description A summary method for an S4 object returned by \code{\link{praxi}}.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{praxi}}.
#'
#' @return A numeric matrix containing the respective positions and change in mean of each detected anomaly. 
#'
#' @rdname summary-methods
#'
#' @aliases summary,praxi.class-method
#' 
#' @export

setMethod("summary",signature=list("praxi.class"),function(object)
{
  return(object@res)
})