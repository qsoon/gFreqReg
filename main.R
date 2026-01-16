library(RColorBrewer)
library(colourvalues)
library(grDevices)
# library(SDMTools)
library(network)
library(ggraph)
library(tidygraph)
library(dplyr)
library(gasper)
library(readxl)
library(forecast)
library(ggfortify)
library(Metrics)
library(GNAR)
library(DTWBI)
# library(vars)
library(geosphere)
# library(xlsx)
library(scales)
library(igraph)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)
library(expm)
library(Hmisc)
library(cccd)
library(zoo)
library(expm)
library(MASS)

source("utils.R")


# load trading network in 2023
load("tradingnw.RData")
plot_graph(ITW_G20_2023)


# make trading network in 2019
ITW_G20_2019 <- ITW_G20_2023

trade2019 <- read.csv("BACI_HS92_V202501/BACI_HS92_Y2019_V202501.csv", header=TRUE)
trade2019$q <- as.numeric(trade2019$q)
trade2019 <- na.omit(trade2019)

trade2019$value <- trade2019$v * trade2019$q


identical(sort(unique(trade2019$i)), sort(unique(trade2019$j)))
head(trade2019)

countrycode <- read.csv("BACI_HS92_V202501/country_codes_V202501.csv", header=TRUE)
countrycode_used <- countrycode[countrycode$country_iso3 %in% rownames(ITW_G20_2023$xy),]
countrycode_used <- countrycode_used[countrycode_used$country_code!=280,]

trade2019_used <- trade2019[(trade2019$i %in% countrycode_used$country_code) & (trade2019$j %in% countrycode_used$country_code),]


trade2019_final <- trade2019_used %>%
  # Create new variables i_new and j_new to treat (i, j) and (j, i) as the same pair
  mutate(i_new = pmin(i, j), j_new = pmax(i, j)) %>%
  # Group by these new variables
  group_by(i_new, j_new) %>%
  # Summarize by summing the total
  summarise(
    total = sum(value) / ifelse(sum(i != i_new | j != j_new) > 0, 2, 1),  # Divide by 2 if both directions exist
    .groups = 'drop'
  ) %>%
  # Rename columns back to original names if desired
  rename(i = i_new, j = j_new)


trade2019_final <- trade2019_final %>% arrange(desc(total))
# trade2019_final <- trade2019_final[1:round(0.03*nrow(trade2019_final)),]

trade2019_final$total <- trade2019_final$total/10^9
trade2019_final$total_log <- log(trade2019_final$total)

rm(trade2019, trade2019_used)


N.ITW_G20 <- nrow(ITW_G20_2019$xy)
ITW_G20_2019$A <- matrix(0,0, nrow=N.ITW_G20, ncol=N.ITW_G20)
colnames(ITW_G20_2019$A) <- colnames(ITW_G20_2023$A)
rownames(ITW_G20_2019$A) <- rownames(ITW_G20_2023$A)
e.sp.weight <- NULL
for(k in 1:nrow(trade2019_final)){
  i <- which(countrycode_used$country_code==as.numeric(trade2019_final[k,1]))
  j <- which(countrycode_used$country_code==as.numeric(trade2019_final[k,2]))
  # e.weight[i,j] <- as.numeric(trade2022_final2[k,4])
  # e.weight[j,i] <- as.numeric(trade2022_final2[k,4])
  ITW_G20_2019$A[i,j] <- as.numeric(trade2019_final[k,4])
  ITW_G20_2019$A[j,i] <- as.numeric(trade2019_final[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(min(i,j),max(i,j), as.numeric(trade2019_final[k,4])))
}

ITW_G20_2019$sA <- e.sp.weight[nrow(e.sp.weight):1,]

plot_graph(ITW_G20_2019)


L.ITW_G20_2019 <- laplacian_mat(ITW_G20_2019$A) # laplacian matrix
val1 <- eigensort(L.ITW_G20_2019)
evalues.ITW_G20_2019 <- val1$evalues
evectors.ITW_G20_2019 <- val1$evectors
# largest eigenvalue
lmax.ITW_G20_2019 <- max(evalues.ITW_G20_2019)

N.ITW_G20 <- nrow(L.ITW_G20_2019)


# load economic data in 2019
# load("X.ITW_G20_2023.RData")

#####################
### economic data 
#####################
X.ITW_G20_2019 <- NULL
# p.ITW_G20 <- 11

files <- list.files(path = "./economic/Indicators", pattern = ".csv")
# [1] "balance.csv"            "export.csv"             "foreign_invest.csv"    
# [4] "GDP_growth.csv"         "GDP_per_cap_growth.csv" "GDP_per_cap.csv"       
# [7] "GNI_per_cap.csv"        "gross_capital.csv"      "import.csv"            
# [10] "inflation.csv"          "PLR.csv"   

# explanatory variables
# file_idx <- c(2,8,9,10)
file_idx <- c(8,9)
p.ITW_G20 <- length(file_idx)
for(i in 1:p.ITW_G20){
  tmp <- read.csv(paste("economic/Indicators",files[file_idx[i]], sep="/"), skip=4)
  X.ITW_G20_2019 <- cbind(X.ITW_G20_2019, tmp[tmp$Country.Code %in% rownames(ITW_G20_2023$xy), "X2019"])
}

# for (i in 1:nrow(X.ITW_G20_2019)) {
#   v <- X.ITW_G20_2019[i,]
#   X.ITW_G20_2019[i,] <- (v - mean(v)) / sd(v)
# }

# gross_capital
# → 실물 투자 요인 (growth theory에서 핵심)
# foreign_invest
# → 대외 자본 유입, gross_capital과는 성격이 다름
# inflation
# → 거시 안정성, 통화 환경
# export (또는 import 중 하나만)
# → 개방도 채널


# response variable
Y.ITW_G20_2019 <- NULL

file_idx <- 4
tmp <- read.csv(paste("economic/Indicators",files[file_idx], sep="/"), skip=4)
Y.ITW_G20_2019 <- tmp[tmp$Country.Code %in% rownames(ITW_G20_2019$xy), "X2019"]


gfreqreg.res <- gFreqReg(X=X.ITW_G20_2019, Y=Y.ITW_G20_2019, S=L.ITW_G20_2019, 
                         evalues=evalues.ITW_G20_2019, evectors=evectors.ITW_G20_2019,
                         g=NULL, M=50, sigma=0.3, method="random", seed=1, lambda=0.01)
  


mu.hat <- evectors.ITW_G20_2019 %*% gfreqreg.res$mu.tilde.hat

B.mat <- array(0, dim=c(N.ITW_G20, N.ITW_G20, p.ITW_G20))

for(i in 1:p.ITW_G20){
  B.mat[,,i] <- evectors.ITW_G20_2019 %*% diag(gfreqreg.res$beta.hat[,i]) %*% Conj(t(evectors.ITW_G20_2019))
}

Y.hat.manual <- mu.hat
for(i in 1:p.ITW_G20){
  Y.hat.manual <- Y.hat.manual + B.mat[,,i] %*% X.ITW_G20_2019[,i]
}

Y.hat <- evectors.ITW_G20_2019 %*% gfreqreg.res$Y.tilde.hat

Y.hat.manual - Y.hat

Y.ITW_G20_2019 - Y.hat

gfreqreg.res$beta.hat
gfreqreg.res$rho.hat

par(mar = c(5, 5, 4, 10))
matplot(
  x = 1:nrow(gfreqreg.res$beta.hat),
  y = gfreqreg.res$beta.hat,
  col = rainbow(ncol(gfreqreg.res$beta.hat)),
  type = "b",        # both points and lines
  lty = 1,
  pch = 16,
  xlab = expression(Graph~frequency~index~"\u2113"),
  ylab = expression(hat(beta)(lambda["\u2113"])),
  main = "Graph Frequency-Domain Regression Coefficients"
)

legend(
  "topright",
  inset = c(-0.25, 0),
  legend = expression(X[1], X[2]),
  col = rainbow(ncol(gfreqreg.res$beta.hat)),
  lty = 1,
  pch = 16,
  cex = 1,
  bty = "n",
  xpd = TRUE
)



par(mar = c(5, 5, 4, 10))
matplot(
  x = 1:ncol(gfreqreg.res$beta.hat),
  y = t(gfreqreg.res$beta.hat),
  type = "b",        # both points and lines
  lty = 1,
  pch = 16,
  xlab = "Coefficient index",
  ylab = "Estimated beta",
  main = "Frequency-wise coefficient estimates"
)

legend(
  "topright",
  inset = c(-0.35, 0),
  legend = paste("Frequency", 1:nrow(gfreqreg.res$beta.hat)),
  col = rainbow(nrow(gfreqreg.res$beta.hat)),
  lty = 1,
  pch = 16,
  cex = 0.7,
  bty = "n",
  xpd = TRUE
)




alpha <- 0.05
ind.accept.H0 <- which((gfreqreg.res$Fstat > qf(1-alpha, p.ITW_G20, (50 - p.ITW_G20 - 1)))==FALSE)

ind.accept.H0

range_min <- range(evectors.ITW_G20_2019[,c(2,3,4)])[1]
range_max <- range(evectors.ITW_G20_2019[,c(2,3,4)])[2]

p1 <- plot_graph_custom4(ITW_G20_2019, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2019[,2], value="Value", ratio=0.6,
                         min=range_min, max=range_max, mg=c(4,4,4,4), title=expression(v[2]), main.title.size = 20)
p2 <- plot_graph_custom4(ITW_G20_2019, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2019[,3], value="Value", ratio=0.6,
                         min=range_min, max=range_max, mg=c(4,4,4,4), title=expression(v[3]), main.title.size = 20)

grid.arrange(p1,p2, nrow=1)



ols.fit_economy <- lm(Y.ITW_G20_2019 ~ X.ITW_G20_2019)

ols.fit_economy$coefficients





#####################
### irregular graph 
#####################

## simulation study

x = rep(0:14, 15) #define coordinates of vertices
y = rep(0:14, rep(15,15))

irregular01 <- list()
set.seed(1)
irregular01$xy <- data.frame(x,y) + rnorm(225,0,0.7)# define 'xy' attributes of graph g which represents coordinates
n.irregular01 <- nrow(irregular01$xy)
rownames(irregular01$xy) <- 1:n.irregular01

distmat.irregular01 <- as.matrix(dist(irregular01$xy))
A.irregular01 <- c()
for(i in 1:(nrow(distmat.irregular01)-1)){
  for(j in (i+1):ncol(distmat.irregular01)){
    val <- distmat.irregular01[i,j]
    A.irregular01 <- rbind(A.irregular01, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.irregular01, k=5), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.irregular01[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=n.irregular01, ncol=n.irregular01)

colnames(wmat) <- 1:n.irregular01
rownames(wmat) <- 1:n.irregular01

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}


# sparse weight matrix
# weight matrix
irregular01$A <- wmat

# sparse weight matrix
irregular01$sA <- sp.wmat

irregular01$dist <- distmat.irregular01
irregular01$sdist <- A.irregular01

L.irregular01 <- gasper::laplacian_mat(irregular01$A)
eigenres.irregular01 <- eigen(L.irregular01)

val2 <- eigensort(L.irregular01)
evalues.irregular01 <- val2$evalues
evectors.irregular01 <- val2$evectors


plot_graph(irregular01)
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = evectors.irregular01[,2],
                   min=-0.16, max=0.16, value="value")
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = evectors.irregular01[,3],
                   min=-0.16, max=0.16, value="value")




gfreqreg_all <- NULL
ols_all <- NULL
mse.gfreqreg <- NULL
mse.ols <- NULL
for(iter in 1:100){
  R <- 100
  X.irreg <- array(0, dim=c(n.irregular01,1,R))
  Y.irreg <- array(0, dim=c(n.irregular01,R))
  
  a <- vector(length=R)
  b <- vector(length=R)
  e <- array(0, dim=c(n.irregular01,R))
  
  for(r in 1:R){
    set.seed(100*iter+r)
    a[r] <- rnorm(1, 1, 0.5)
    b[r] <- rnorm(1, 1, 0.5)
    e[,r] <- rnorm(n.irregular01, 0, 0.3)
    
    X.irreg[,1,r] <- a[r]*evectors.irregular01[,2] + b[r]*evectors.irregular01[,3]
    Y.irreg[,r] <- a[r]*evectors.irregular01[,2] - b[r]*evectors.irregular01[,3] + e[,r]
  }
  
  
  gfreqreg.res.irregular01 <- gFreqReg(X=X.irreg, Y=Y.irreg, S=L.irregular01, 
                                       evalues=evalues.irregular01, evectors=evectors.irregular01,
                                       g=NULL, M=50, sigma=0.5, method=NULL, seed=1, lambda=0.01)
  
  mse.gfreqreg <- c(mse.gfreqreg, mean((Y.irreg - gfreqreg.res.irregular01$Y.tilde.hat)^2))
  
  ols.fit <- lm(as.numeric(Y.irreg) ~ as.numeric(X.irreg))
  
  mse.ols <- c(mse.ols, mean((as.numeric(Y.irreg) - ols.fit$fitted.values)^2))
  
  gfreqreg_all <- rbind(gfreqreg_all, as.numeric(gfreqreg.res.irregular01$beta.hat))
  ols_all <- rbind(ols_all, as.numeric(ols.fit$coefficients))
}

colMeans(gfreqreg_all) ; apply(gfreqreg_all, 2, sd)
colMeans(ols_all) ; apply(ols_all, 2, sd)

mean(mse.gfreqreg) ; sd(mse.gfreqreg)
mean(mse.ols) ; sd(mse.ols)


