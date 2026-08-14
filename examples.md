#Installing praxi

```{r}
if(system.file(package = "praxi") != ""){
  remove.packages("praxi")
}
library("remotes")
install_github("DylanBahia/praxi",force=TRUE)
```

#Example

## Generate AR time-series

```{r}
set.seed(1)
#time series length
n <- 20000
#AR coefficients
gam <- c(0.5,0.3)

y <- as.numeric(arima.sim(model=list(ar=gam),n=n))
```

## Add anomalies to time-series

```{r}
#Collective anomaly
y[5000:6000] <- y[5000:6000]+2
#Point anomaly
y[11000] <- y[11000]+10
```

## Run algorithm and get results

```{r}
#Model order
p <- length(gam)
#Penalty
b <- 4*log(n)

#Run algorithm
res <- praxi::praxi(y,p=p,b=b)

#Table of anomalies detected (if any)
anomalies(res)
#Plot showing point and collective anomalies detected on dara
plot(res)
```