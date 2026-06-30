.praxi.class <- setClass("praxi.class",representation(y="vector",p="numeric",b="numeric"))

praxi.class <- function(y,p,b)
{
	.praxi.class(y=y,p=p,b=b)	
}

praxi <- function(y,p,b)
{
	result <- ar_alg_call(y,p,b)
	#rlist <- list("cpts"=result,"data"=y,"order"=p,"penalty"=b)
	rlist <- praxi.class(y,p,b)
	return(rlist)
}

