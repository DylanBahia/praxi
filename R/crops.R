#' @useDynLib praxi, .registration = TRUE
#' @import magrittr
#' @importFrom sets set tuple as.tuple set_union set_symdiff set_is_empty as.set
#' @import tibble
#' @import ggplot2
#' @import reshape
#' @import crops
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
NULL


#' @export
crops <- function(y,p,b_min,b_max){
  
  func <- function(b){
    result <- praxi(y,p,b)
    return(list(result@cost-nrow(result@res)*b,nrow(result@res),result@res))
  }
  
  result <- unique(crops::crops(func,b_min,b_max))
  
  return(result)
}

#' @export
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
  key_value_pairs <- Map(sets::tuple,keys,object@betas %>% unlist)
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

#' @export
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

#' @export
setGeneric("segmentations",function(object) {standardGeneric("segmentations")})
setMethod("segmentations",signature=list("crops.class"),
          function(object)
          {
            # appease package checks
            . <- NULL
            segs <- Map(object@method,unlist(object@betas))
            valid_segs <- Filter(function(x) x[[2]] > 1, segs)
            if(length(valid_segs) == 0)
            {
              return(NULL)
            }
            n <- segs %>% Map(function(.) .[[2]],.) %>% unlist %>% max
            mat <- segs %>% 
              Map(function(.) if(.[[2]]!=0){paste0("(",.[[3]][,1],",", .[[3]][,2],")")} else{NA},.) %>% 
              Map(function(.) c(.,rep(NA,n-length(.))),.) %>%
              Reduce(rbind,.,matrix(nrow=0,ncol=n),right=TRUE)
            colnames(mat) <- Map(function(.) paste("anom.",.,sep=""),1:n) %>% unlist      
            return(      
              Map(function(beta,seg) tibble(beta=beta,Qm=seg[[1]],Q=seg[[1]]+beta*seg[[2]],
                                            m=seg[[2]]),
                  unlist(object@betas),
                  segs) %>%                                          
                Reduce(add_row,.,tibble(beta=numeric(),Qm=numeric(),Q=numeric(),m=numeric())) %>%
                cbind(.,mat)
            )           
          })

#' Visualisation of data, costs, penalty values and anomaly locations.
#'
#' @name plot
#'
#' @description Plot methods for an S4 object returned by \code{\link{crops}}. The plot can also be combined with the original data if required.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{crops}}.
#' @param y A dataframe containing the locations and values of the data points. The data plot is plotted below, and is aligned with, the changepoint locations.
#' 
#' @return A ggplot object. Note - if no changepoints are detected in the penalty interval [beta_min,beta_max], then the value returned is NULL. 
#'
#' @rdname plot-methods
#'
#' @aliases plot,crops.class,data.frame-method
#' 
#' @seealso \code{\link{crops}}.
#'
#' @examples
#' # see the crops example
#'
#' @export  
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
              separate(value, into = c("st", "nd"), sep = ",") %>%
              mutate(
                st = as.numeric(gsub("\\(", "", st)),
                nd = as.numeric(gsub("\\)", "", nd))
              ) %>% 
              ggplot(.) %>% 
              add(geom_hline(aes(yintercept=dummy),linewidth=0.5)) %>%
              add(geom_segment(aes(x=st,xend=nd,y=dummy,yend=dummy),colour="blue",linewidth=1.5)) %>%
              add(geom_point(aes(x=ifelse(st==nd,st,NA),y=ifelse(st==nd,dummy,NA)),colour="red",na.rm=TRUE)) %>%
              add(labs(x="Index",y="Penalty")) %>%
              add(scale_y_continuous(breaks = seq(1:nrow(df)),labels=signif(df$beta,digits=3),sec.axis = sec_axis( ~.,breaks = seq(1:nrow(df)),labels=signif(df$Qm,digits=4),name="Unpenalised cost"))) %>%
              add(theme_bw())
            return(p)       
          })