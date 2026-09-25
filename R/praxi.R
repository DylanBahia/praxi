#' @useDynLib(mypackage, .registration = TRUE)
#' @importFrom(Rcpp, evalCpp)
NULL

.praxi.class <- setClass("praxi.class",representation(y="numeric",p="numeric",b="numeric",cost="numeric",res="matrix"))

praxi.class <- function(y,p,b,cost,res)
{
	.praxi.class(y=y,p=p,b=b,cost=cost,res=res)	
}

#' Find the best segmentation of data for a change-in-slope model
#'
#' @name ouralg
#'
#' @description The CPOP algorithm fits a change-in-slope model to data.
#'
#' @param y A vector of length n containing the data.
#' @param x A positive integer for the order of the AR model.
#' @param b A positive real number for the penalty 
#'
#' @return An instance of an S4 class of type praxi.class.
#'
#' @details \loadmathjax{} The CPOP algorithm fits a change-in-slope model to data. It assumes that we have we have data points \mjeqn{(y_1,x_1),\ldots,(y_n,x_n)}{(y_1,x_1),...,(y_n,x_n)}, ordered
#' so that \mjteqn{x_1 < x_2 < ... < x_n}{x_1 \lt x_2 \lt ... \lt x_n}{x_1 < x_2 < ... < x_n}. For example  \mjeqn{x_i}{x_i} could be a time-stamp of when response \mjeqn{y_i}{y_i} is obtained. We model the response, \mjeqn{y}{y}, as a signal
#' plus noise where the signal is modelled as a continuous piecewise linear function of \mjeqn{x}{x}. That is \mjdeqn{y_i=f(x_i)+\epsilon_i}{y_i=f(x_i)+e_i} where \mjeqn{f(x)}{f(x)} is a continuous
#' piecewise linear function.
#'
#' To estimate the function \mjeqn{f(x)}{f(x)} we specify a set of \mjeqn{N}{N} grid points, \mjeqn{g_{1:N}}{g_{1:N}} with these ordered so that \mjteqn{g_i < g_j}{g_i \lt g_j}{g_i < g_j} if and only if \mjteqn{i < j}{i \lt j}{i < j},
#' and allow the slope of \mjeqn{f(x)}{f(x)} to only change at these grid points. We then estimate the number of changes, the location of the changes, and hence the resulting
#' function \mjeqn{f(x)}{f(x)} by minimising a penalised weighted least squares criteria. This criteria includes a term which measures fit to the data, which is calculated as the sum of the
#'square residuals of the fitted function, scaled by the variance of the noise. The penalty is proportional to the number of changes.
#'
#' That is our estimated function will depend on \mjeqn{K}{K}, the number of changes in slope, their locations, \mjeqn{\tau_1,\ldots,\tau_K}{tau_1,...,tau_K}, and the value of the function
#' \mjeqn{f(x)}{f(x)} at these change points, \mjeqn{\alpha_1,\ldots,\alpha_K}{alpha_1,...,alpha_K}, and its values, \mjeqn{\alpha_0}{alpha_0} at  \mjteqn{\tau_0 < x_1}{tau_0 \lt x_1}{tau_0 < x_1} and
#' \mjeqn{\alpha_{K+1}}{alpha_{K+1}} at some \mjeqn{\tau_{K+1} > x_N}{tau_{K+1}>x_N}. The CPOP algorithm then estimates \mjeqn{K}{K}, \mjeqn{\tau_{1:K}}{tau_{1:K}} and \mjeqn{\alpha_{0:K+1}}{alpha_{0:K+1}} 
#' as the values that solve the following minimisation problem
#'\mjdeqn{\min_{K,\tau_{1:K}\in g_{1:N}, \alpha_{0:K+1} }\left\lbrace\sum_{i=1}^n \frac{1}{\sigma^2_i} \left(y_i -  \alpha_{j(i)}-(\alpha_{j(i)+1}-  \alpha_{j(i)})\frac{x_i-\tau_{j(i)}}{\tau_{j(i)+1}-\tau_{j(i)}}\right)^2+K\beta\right\rbrace}{\min_{K,\tau_{1:K}\in g_{1:N}, \alpha_{0:K+1} }\left\lbrace\sum_{i=1}^n \frac{1}{\sigma^2_i} \left(y_i -  \alpha_{j(i)}-(\alpha_{j(i)+1}-  \alpha_{j(i)})\frac{x_i-\tau_{j(i)}}{\tau_{j(i)+1}-\tau_{j(i)}}\right)^2+K\beta\right\rbrace}
#'
#' where \mjeqn{\sigma^2_1,\ldots,\sigma^2_n}{sigma^2_1,...,sigma^2_n} are the variances of the noise \mjeqn{\epsilon_i}{\epsilon_i} for \mjeqn{i=1,\ldots,n}{i=1,...,n}, and \mjeqn{\beta}{beta} is the
#' penalty for adding a changepoint. The sum in this expression is the weighted residual sum of squares, and the \mjeqn{K\beta}{K*beta} term is the penalty for having \mjeqn{K}{K} changes. 
#'
#' If we know, or have a good estimate of, the residual variances, and the noise is (close to) independent over time then an appropriate choice for the penalty is
#' \mjeqn{\beta=2 \log n}{beta=2log(n)}, and this is the default for CPOP. However in many applications these assumptions will not hold and it is advised to look at segmentations for
#' different value of \mjeqn{\beta}{beta} -- this is possible
#' using CPOP with the CROPS algorithm \code{\link{cpop.crops}}. Larger values of \mjeqn{\beta}{beta} will lead to functions with fewer changes. Also there is a trade-off between the variances of the residuals
#' and \mjeqn{\beta}{beta}: e.g. if we double the variances and half the value of \mjeqn{\beta}{beta} we will obtain the same estimates for the number and location of the changes and the
#' underlying function.
#'
#' @rdname praxi
#'
#' @references \insertRef{doi:10.1080/10618600.2018.1512868}{cpop}
#' @references \insertRef{cpop-jss-article-2024}{cpop}
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
#'#analyse time series using praxi
#'res <- praxi::praxi(y=y,p=length(gam),4*log(n))
#'anoms <- anomalies(praxi)
#'print(anoms)
#'pt <- plot(praxi)
#'print(praxi)
#'
#'
#' @export

praxi <- function(y,p,b=NULL)
{
  
  if(is.null(b)){
    b <- 4*log(length(y))
  }
  
	result <- ar_alg_call(y,p,b)
	#rlist <- list("cpts"=result,"data"=y,"order"=p,"penalty"=b)
	rlist <- praxi.class(y,p,b,result[[1]],result[[2]])
	return(rlist)
}

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


setGeneric("anomalies",function(object,...)standardGeneric("anomalies"))

setMethod("anomalies",signature=list("praxi.class"),function(object)
{
  return(object@res)
})

#' @export
crops <- function(y,p,b_min,b_max){
  
  func <- function(b){
    result <- praxi(y,p,b)
    return(list(result@cost,rep(1,nrow(result@res)),result@res))
  }
  
  result <- praxi::unique(crops::crops(func,b_min,b_max))
  
  return(result)
}

setMethod("summary",signature=list("crops.class"),function(object){
  cat("crops analysis",sep="")
  cat('\n',sep="")
  cat('\n',sep="")
  cat("minimum penalty value = ",min(object@betas)," : maximum penalty value = ",max(object@betas),sep="")
  cat('\n',sep="")
  segs <- crops::segmentations(object)
  if(is.null(segs))
  {
    cat("no anomalies found in the penalty interval [",min(object@betas),",",max(object@betas),"]",sep="")
    cat('\n',sep="")
  }
  else
  {
    cat("number of segmentations calculated : ",nrow(segs),sep="")
    cat('\n',sep="")	    
    cat("least number of anomalies  = ",min(segs$m), " : maximum number of anomalies = ",max(segs$m),sep="")
    cat('\n',sep="")
  }
  invisible()
})

setMethod("unique",signature=list("crops.class"),function(x)
{
  # appease package checks
  . <- NULL
  object<-x
  hash_map <- new.env()
  keys <- object@betas %>% 
    unlist %>% 
    Map(object@method,.) %>% 
    Map(function(.) .[2],.) %>% 
    Map(as.character,.)
  key_value_pairs <- Map(tuple,keys,object@betas %>% unlist)
  hash_map <-    
    key_value_pairs %>%  
    Reduce(function(pair,map) {map[[pair[[1]]]] <- pair[[2]]; return(map);},
           .,
           hash_map,right=TRUE)
  object@betas <-    
    Map(function(key) hash_map[[key]],
        hash_map %>% ls) %>%
    unname %>% 
    as.set
  return(object)
})

setMethod("subset",signature=list("crops.class"), function(x,beta_min=0,beta_max=Inf)
{
  # appease package checks
  . <- NULL
  object <- x
  object@betas %<>% 
    unlist %>% 
    Filter(function(.) . <= beta_max & . >= beta_min,.) %>% 
    as.set
  return(object)            
})

setGeneric("segmentations",function(object) {standardGeneric("segmentations")})
setMethod("segmentations",signature=list("crops.class"),
          function(object)
          {
            # appease package checks
            . <- NULL
            segs <- Map(object@method,unlist(object@betas))
            valid_segs <- Filter(function(x) length(x[[2]]) > 1, segs)
            if(length(valid_segs) == 0)
            {
              return(NULL)
            }
            n <- segs %>% Map(function(.) .[[2]],.) %>% Map(length,.) %>% unlist %>% max
            mat <- segs %>% 
              Map(function(.) paste0("(",.[[3]][1,],",", .[[3]][2,],")"),.) %>% 
              Map(function(.) c(.,rep(NA,n-length(.))),.) %>%
              Reduce(rbind,.,matrix(nrow=0,ncol=n),right=TRUE)
            colnames(mat) <- Map(function(.) paste("anom.",.,sep=""),1:n) %>% unlist      
            return(      
              Map(function(beta,seg) tibble(beta=beta,Qm=seg[[1]],Q=seg[[1]]+beta*length(seg[[2]]),
                                            m=length(seg[[2]])),
                  unlist(object@betas),
                  segs) %>%                                          
                Reduce(add_row,.,tibble(beta=numeric(),Qm=numeric(),Q=numeric(),m=numeric())) %>%
                cbind(.,mat)
            )           
          })

setMethod("plot",signature=list("crops.class"),
          function(x)
          {
            # appease ggplot and tidyverse
            . <- Q <- Qm <- m <- value <- dummy <- NULL
            object <- x
            df <- segmentations(object)
            if(is.null(df))
            {
              return(NULL)
            }
            df <- cbind(df,data.frame("dummy"=1:nrow(df)))
            p <- df %>%
              subset(.,select = -c(beta,Q,Qm,m)) %>%
              melt(., id=c("dummy")) %>% 
              .[complete.cases(.), ] %>%
              ggplot(.,aes(x=value,y=dummy)) %>% 
              add(geom_point()) %>%
              add(labs(x="location",y="penalty")) %>%
              add(geom_hline(aes(yintercept=dummy))) %>%
              add(scale_y_continuous(breaks = seq(1:nrow(df)),labels=signif(df$beta,digits=3),sec.axis = sec_axis( ~.,breaks = seq(1:nrow(df)),labels=signif(df$Qm,digits=4),name="unpenalised cost"))) %>%
              add(theme_bw())
            return(p)       
          })