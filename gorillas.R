##================================================================================================
##
## Spatial occupancy model for Grauer's gorillas with the hSDM R package for hierachical Bayesian 
## species distribution models.
##
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## Andrew Plumptre <aplumptre@wcs.org>
##
## August 2016
##
## Note: Eastern lowland gorilla (Gorilla beringei graueri) is also known as the Grauer's gorilla.
##================================================================================================

## This script uses version v1.4 of the hSDM R package (http://dx.doi.org/10.5281/zenodo.48470). 
## If this version has been replaced by a more recent one on CRAN, you can still install it from GitHub:
## library(devtools)
## install_github(repo="ghislainv/hSDM",ref="v1.4")

## Start the clock
t.start <- Sys.time()

## Load libraries
## "coda" for analysing MCMC
## "hSDM" for hierarchical Bayesian species distribution models
## "sp" for spatial objects
## "raster" for manipulating rasters
## "foreach" for parallel computation
## "doParallel" for parallel computation
## "knitr" to write reports with R Markdown
## "rmarkdown" to process the document, first with knitr and then with pandoc
## "ggplot2" to plot results
## "ggmap" to plot maps including a GoogleEarth background
pkg <- c("coda","hSDM","sp","raster","foreach","doParallel","knitr","rmarkdown","ggplot2","ggmap")
for (i in pkg) {
  if (!require(i,character.only=TRUE)) {
    install.packages(i,dependencies=TRUE)
    require(i,character.only=TRUE)
  }
}

## Create directory to save results
dir.create("results")

## Importing raster stack of explicative variables (5km res)
G.stk <- stack("data/rasters/environ.tif")

## Rename layers
names(G.stk) <- c("bio2","bio12","bio17","dem","treecov","disforlos","minedis","rivdis",
                  "roaddis","rugged","slope","stslopdis","villdis")

## Environmental variables
## =======================
## bio2: mean diurnal temperature range
## bio12: mean annual precipitation
## bio17: precipitation of driest quarter
## dem: elevation above sea level
## rugged: ruggedness of topography
## rivdis: distance to river
## slope: slope, calculated from DEM layer
## stslopdis: distance to steep slopes
## treecov: percentage of tree cover

## Human variables
## ===============
## disforlos: distance to recent deforestation
## roaddis: distance to roads
## villdis: distance to villages
## minedis: distance to mines

## Import observations
Gorilla.obs <- read.csv("data/gorillas/gorillas.csv", header=TRUE)
G.coords <- cbind(Gorilla.obs$EW,Gorilla.obs$NS)
Gorilla.obs.sp <- SpatialPointsDataFrame(coords=G.coords,
                                          data=Gorilla.obs)
projection(Gorilla.obs.sp) <- CRS("+init=epsg:32735")

## Function to plot the observation points on raster maps
fun.obs <- function() {
    plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences==0,],pch=".",add=TRUE)
    plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences>0,],pch=3,add=TRUE)
}
## Plot rasters of covariates
pdf(file="results/covariates.pdf",width=12,height=10)
plot(G.stk,addfun=fun.obs,maxnl=20,nc=4,nr=5,legend.mar=15)
dev.off()

## Plot observations and presence sites
pdf(file="results/observations.pdf")
par(mar=c(3,3,1,1))
plot(G.stk$dem) ## Background=altitude
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences==0,],pch=".",col=grey(0.5),add=TRUE)
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences>0,],pch=3,add=TRUE)
dev.off()

## Extract environmental values and cell number for observations
df <- as.data.frame(extract(G.stk,Gorilla.obs.sp,cellnumber=TRUE))
df$Trials <- Gorilla.obs.sp$Trials
df$Presences <- Gorilla.obs.sp$Presences
df.cc <- df[complete.cases(df),] 

## Normalized continuous covariates
df.norm <- df.cc
Mean <- vector()
Sd <- vector()
for (i in c(2:14)) {
    m <- mean(df.cc[,i],na.rm=TRUE)
    s <- sd(df.cc[,i],na.rm=TRUE)
    Mean <- c(Mean,m)
    Sd <- c(Sd,s)
    df.norm[,i] <- (df.cc[,i]-m)/s
}
## Data-frame with mean and sd for each variable
df.mean.sd <- as.data.frame(rbind(Mean,Sd))
names(df.mean.sd) <- names(df.norm)[c(2:14)]

## Raster stack for predictions (with normalized covariates)
G.stk.pred <- G.stk
for (i in c(2:14)) {
    var.name <- names(df.norm)[i] ## Variable name
    w <- which(names(G.stk.pred)==var.name) ## Position in the stack 
    m <- df.mean.sd[1,var.name] ## Mean
    s <- df.mean.sd[2,var.name] ## Sd
    orig <- values(subset(G.stk.pred,w)) ## Original values
    trans <- (orig-m)/s ## Transformed values
    G.stk.pred[[w]][] <- trans
}

## Plot transformed covariates
pdf(file="results/covariates_transformed.pdf",width=12,height=10)
plot(G.stk.pred,maxnl=20,nc=4,nr=5,legend.mar=10)
dev.off()

## Select only grid cells with no NA
G.df.pred <- as.matrix(G.stk.pred)
w <- complete.cases(G.df.pred) ## Note: w will be used to obtain the cell identifier for predictions in iCAR model
G.df.pred.complete <- as.data.frame(G.df.pred[w,])

## Make a cluster for parallel MCMCs
nchains <- 2
ncores <- nchains ## One core for each MCMC chains
clust <- makeCluster(ncores)
registerDoParallel(clust)

## Starting values and random seed
seed <- 1234
set.seed(seed)
beta.start <- runif(nchains,-1,1)
gamma.start <- runif(nchains,-1,1)
Vrho.start <- runif(nchains,0,10)
seed.mcmc <- round(runif(nchains,0,1e6))

##===============================================
##
## 1. Model with environmental variables
##
##===============================================

## bio2: mean diurnal temperature range
## bio12: mean annual precipitation
## bio17: precipitation of driest quarter
## dem: elevation above sea level
## rugged: ruggedness of topography
## rivdis: distance to river
## slope: slope, calculated from DEM layer
## stslopdis: distance to steep slopes
## treecov: percentage of tree cover

## hSDM model using Zero Inflated Binomial (ZIB)
mod.ZIB.env <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=df.norm$Presences,
                  trials=df.norm$Trials,
                  suitability=~bio2+bio12+bio17+dem+rugged+rivdis+slope+stslopdis+treecov,     
                  observability=~1,
                  data=df.norm,
                  suitability.pred=G.df.pred.complete,
                  burnin=5000,
                  mcmc=5000, thin=5,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
ZIB.env.mcmc <- mcmc.list(lapply(mod.ZIB.env,"[[","mcmc"))

## Outputs summary
ZIB.env.stat <- summary(ZIB.env.mcmc)$statistics
dir.create("results/1.ZIB.env")
sink(file="results/1.ZIB.env/mcmc_summary.txt")
ZIB.env.stat
cat(rep("\n",3))
gelman.diag(ZIB.env.mcmc)
sink()
## Deviance
deviance.ZIB.env <- ZIB.env.stat["Deviance","Mean"]
## Detection probability
gamma.hat <- ZIB.env.stat["gamma.(Intercept)","Mean"]
delta.est <- inv.logit(gamma.hat) ## 0.105
## Plot trace and posterior distributions
pdf("results/1.ZIB.env/mcmc_trace.pdf")
plot(ZIB.env.mcmc)
dev.off()

## Interpretation
## ==============
## The detection probability is 0.105: ~10% of chance of observing a Gorilla if the habitat is suitable.
## Based on *statistical significance* (95% credible intervals) and *biological relevance*,
## we kept only the following environmental variables: dem and tree over. 

##===================================================================
##
## 2. Model with only significant and relevant environmental variables
##
##===================================================================

## hSDM model using Zero Inflated Binomial (ZIB)
mod.ZIB.envsign <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=df.norm$Presences,
                  trials=df.norm$Trials,
                  suitability=~dem+treecov, # envsign covariates
                  observability=~1,
                  data=df.norm,
                  suitability.pred=G.df.pred.complete,
                  burnin=5000,
                  mcmc=5000, thin=5,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
ZIB.envsign.mcmc <- mcmc.list(lapply(mod.ZIB.envsign,"[[","mcmc"))

## Outputs summary
ZIB.envsign.stat <- summary(ZIB.envsign.mcmc)$statistics
dir.create("results/2.ZIB.envsign")
sink(file="results/2.ZIB.envsign/mcmc_summary.txt")
ZIB.envsign.stat
cat(rep("\n",3))
gelman.diag(ZIB.envsign.mcmc)
sink()
## Deviance
deviance.ZIB.envsign <- ZIB.envsign.stat["Deviance","Mean"]
## Plot trace and posterior distributions
pdf("results/2.ZIB.envsign/mcmc_trace.pdf")
plot(ZIB.envsign.mcmc)
dev.off()

## Prediction on the landscape
prob.p <- subset(G.stk,1) ## create a raster for predictions
values(prob.p)[w] <- mod.ZIB.envsign[[1]]$prob.p.pred ## assign predicted values
values(prob.p)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/2.ZIB.envsign/predictions.pdf")
plot(prob.p)
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences==0,],pch=".",col=grey(0.5),add=TRUE)
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p,filename="results/2.ZIB.envsign/pred_mod_ZIB_envsign.tif",overwrite=TRUE)

##================================================================
##
## 3. Model with environmental variables AND human variables
##
##================================================================

## disforlos: distance to recent deforestation
## roaddis: distance to roads
## villdis: distance to villages
## minedis: distance to mines

## hSDM model using Zero Inflated Binomial (ZIB)
mod.ZIB.envhum <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=df.norm$Presences,
                  trials=df.norm$Trials,
                  suitability=~dem+treecov+ # envsign covariates
                    disforlos+roaddis+villdis+minedis, # human covariates
                  observability=~1,
                  data=df.norm,
                  suitability.pred=G.df.pred.complete,
                  burnin=5000,
                  mcmc=5000, thin=5,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
ZIB.envhum.mcmc <- mcmc.list(lapply(mod.ZIB.envhum,"[[","mcmc"))

## Outputs summary
ZIB.envhum.stat <- summary(ZIB.envhum.mcmc)$statistics
dir.create("results/3.ZIB.envhum")
sink(file="results/3.ZIB.envhum/mcmc_summary.txt")
ZIB.envhum.stat
cat(rep("\n",3))
gelman.diag(ZIB.envhum.mcmc)
sink()
## Deviance
deviance.ZIB.envhum <- ZIB.envhum.stat["Deviance","Mean"]
## Plot trace and posterior distributions
pdf("results/3.ZIB.envhum/mcmc_trace.pdf")
plot(ZIB.envhum.mcmc)
dev.off()

## Interpretation
## Based on *statistical significance* (95% credible intervals) and *biological relevance*,
## we kept only the following human variables: disforlos

## The probability of presence increases strongly with proximity to mines
## but this is not a biologically relevant factor (mines dont attract
## Gorilla) so we have to remove this explicative variable from the model

## The probability of presence increases slightly with distance to recent deforestation

##======================================================================
##
## 4. Model with significant environmental variables AND human variables
##
##======================================================================

## hSDM model using Zero Inflated Binomial (ZIB)
mod.ZIB.envhumsign <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=df.norm$Presences,
                  trials=df.norm$Trials,
                  suitability=~dem+treecov+ # envsign covariates
                    disforlos, # humansign covariates
                  observability=~1,
                  data=df.norm,
                  suitability.pred=G.df.pred.complete,
                  burnin=5000,
                  mcmc=5000, thin=5,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
ZIB.envhumsign.mcmc <- mcmc.list(lapply(mod.ZIB.envhumsign,"[[","mcmc"))

## Outputs summary
ZIB.envhumsign.stat <- summary(ZIB.envhumsign.mcmc)$statistics
dir.create("results/4.ZIB.envhumsign")
sink(file="results/4.ZIB.envhumsign/mcmc_summary.txt")
ZIB.envhumsign.stat
cat(rep("\n",3))
gelman.diag(ZIB.envhumsign.mcmc)
sink()
## Deviance
deviance.ZIB.envhumsign <- ZIB.envhumsign.stat["Deviance","Mean"]
## Plot trace and posterior distributions
pdf("results/4.ZIB.envhumsign/mcmc_trace.pdf")
plot(ZIB.envhumsign.mcmc)
dev.off()

## Prediction on the landscape
prob.p <- subset(G.stk,1) ## create a raster for predictions
values(prob.p)[w] <- mod.ZIB.envhumsign[[1]]$prob.p.pred ## assign predicted values
values(prob.p)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/4.ZIB.envhumsign/predictions.pdf")
plot(prob.p)
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences==0,],pch=".",col=grey(0.5),add=TRUE)
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p,filename="results/4.ZIB.envhumsign/pred_mod_ZIB_envhumsign.tif",overwrite=TRUE)

##===============================================
##
## 5. Model with spatial autocorrelation
##
##===============================================

## Landscape and neighbors
ncells <- ncell(G.stk)
neighbors.mat <- adjacent(G.stk, cells=c(1:ncells), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]
cells.pred <- which(w) ## Vector w indicates the cells with environmental information (without NA)

## ZIB.iCAR model
mod.ZIB.iCAR <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB.iCAR(
    ## Observations
    presences=df.norm$Presences,
    trials=df.norm$Trials,
    ## Habitat
    suitability=~dem+treecov+disforlos,
    observability=~1,
    ## Data-set
    data=df.norm,
    ## Spatial structure
    spatial.entity=df.norm$cells,
    n.neighbors=n.neighbors,
    neighbors=adj,
    ## Predictions
    suitability.pred=G.df.pred.complete,
    spatial.entity.pred=cells.pred,
    ## Chains
    burnin=10000, mcmc=10000, thin=10,
    ## Starting values
    beta.start=beta.start[i],
    gamma.start=gamma.start[i],
    Vrho.start=Vrho.start[i],
    ## Priors
    priorVrho="Uniform",
    #priorVrho=10,
    mubeta=0, Vbeta=1.0E6,
    mugamma=0, Vgamma=1.0E6,
    Vrho.max=10,
    ## Various
    seed=seed.mcmc[i], verbose=1,
    save.rho=0, save.p=1) ## Set save.p=1 to save predictive posterior for each spatial cell
  return(mod)
}

## Extract list of MCMCs from output
ZIB.iCAR.mcmc <- mcmc.list(lapply(mod.ZIB.iCAR,"[[","mcmc"))

## Outputs summary
ZIB.iCAR.stat <- summary(ZIB.iCAR.mcmc)$statistics
dir.create("results/5.ZIB.iCAR")
sink(file="results/5.ZIB.iCAR/mcmc_summary.txt")
ZIB.iCAR.stat
cat(rep("\n",3))
gelman.diag(ZIB.iCAR.mcmc)
sink()
## Deviance
deviance.ZIB.iCAR <- ZIB.iCAR.stat["Deviance","Mean"]
## Plot trace and posterior distributions
pdf("results/5.ZIB.iCAR/mcmc_trace.pdf")
plot(ZIB.iCAR.mcmc)
dev.off()

## Spatial random effects
rho <- subset(G.stk,1) ## create a raster
values(rho) <- mod.ZIB.iCAR[[1]]$rho.pred
pdf(file="results/5.ZIB.iCAR/random_effects.pdf")
plot(rho)
dev.off()

## Prediction on the landscape
prob.p <- subset(G.stk,1) ## create a raster for predictions
values(prob.p)[w] <- apply(mod.ZIB.iCAR[[1]]$prob.p.pred,2,mean) ## assign predicted values
values(prob.p)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/5.ZIB.iCAR/predictions.pdf")
plot(prob.p,zlim=c(0,1))
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences==0,],pch=".",add=TRUE)
plot(Gorilla.obs.sp[Gorilla.obs.sp$Presences>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p,filename="results/5.ZIB.iCAR/pred_mod_ZIB_iCAR.tif",overwrite=TRUE)

##===============================================
##
## Model comparison based on deviance
##
##===============================================

## Null model
mod.ZIB.null <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=df.norm$Presences,
                             trials=df.norm$Trials,
                             suitability=~1,
                             observability=~1,
                             data=df.norm,
                             suitability.pred=G.df.pred.complete,
                             burnin=5000,
                             mcmc=5000, thin=5,
                             beta.start=beta.start[i],
                             gamma.start=gamma.start[i],
                             mubeta=0, Vbeta=1.0E6,
                             mugamma=0, Vgamma=1.0E6,
                             seed=seed.mcmc[i], verbose=1,
                             save.p=0)
  return(mod)
}

## Stop cluster
stopCluster(clust)

## Extract list of MCMCs from output
ZIB.null.mcmc <- mcmc.list(lapply(mod.ZIB.null,"[[","mcmc"))

## Deviance
ZIB.null.stat <- summary(ZIB.null.mcmc)$statistics
deviance.null <- ZIB.null.stat["Deviance","Mean"]

## Full or saturated model
## see here: http://www.stat.wisc.edu/courses/st849-bates/lectures/GLMDeviance.pdf
w1 <- which(df.norm$Presences>0)
logL.full <- sum(dbinom(df.norm$Presences[w1],df.norm$Trials[w1],mod.ZIB.null[[1]]$prob.q.latent[w1],log=TRUE))
deviance.full <- -2*logL.full

##= Table of deviance
dev.tab <- data.frame(Model=rep(NA,5),Deviance=rep(0,5),Perc=rep(0,5))
dev.tab$Model <- c("NULL","env","env+hum","env+hum+iCAR","FULL")
dev.tab$Deviance <- c(deviance.null,deviance.ZIB.envsign,deviance.ZIB.envhumsign,
                      deviance.ZIB.iCAR,deviance.full)
dev.tab$Perc <- round(100*(dev.tab$Deviance[1]-dev.tab$Deviance)/(dev.tab$Deviance[1]-dev.tab$Deviance[5]))
##= Export
sink(file="results/deviance.txt")
dev.tab
sink()

##================================================================
##
## TSS (True Skill Statistics) and SDA (Species Distribution Area)
##
##================================================================

## Function to compute threshold dependent indexes
Index.fun <- function(obs,prob.p,thresh) {
    ## Transform probabilities into {0,1} given threshold
    pred <- ifelse(prob.p>=thresh,1,0)
    ## Contingency table (pred/obs)
    n00 <- sum(pred==0 & obs==0,na.rm=TRUE)
    n11 <- sum(pred==1 & obs==1,na.rm=TRUE)
    n01 <- sum(pred==0 & obs==1,na.rm=TRUE)
    n10 <- sum(pred==1 & obs==0,na.rm=TRUE)
    ## Threshold  dependent indexes
    OA <- (n11+n00)/(n11+n10+n00+n01) ## Overall accuracy
    Sensitivity <- n11/(n11+n01)
    Specificity <- n00/(n00+n10)
    TSS <- Sensitivity+Specificity-1
    return(list(OA=OA,TSS=TSS,Sens=Sensitivity,Spe=Specificity))
}

## Suitable sites
df$Suit <- 0
df$Suit[df$Presences>0] <- 1

## Extract predicted probability of presence
df$prob.p <- extract(prob.p,Gorilla.obs.sp)

## TSS as a function of the threshold
OA <- vector()
TSS <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
    Index <- Index.fun(df$Suit,df$prob.p,thresh.seq[i])
    OA[i] <- Index$OA
    TSS[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS
maxTSS <- max(TSS,na.rm=TRUE) ## >> 0.80 => Very good model thanks to CAR
w.th <- which(TSS==maxTSS)
thresh.maxTSS <- mean(thresh.seq[w.th],na.rm=TRUE)
OA.thresh.maxTSS <- Index.fun(df$Suit,df$prob.p,thresh.maxTSS)$OA 
tss.df <- data.frame(masTSS=maxTSS,OA=OA.thresh.maxTSS,prob=thresh.maxTSS)

## Plot evolution of TSS with threshold
pdf(file="results/TSS.pdf")
plot(thresh.seq,TSS,type="l",xlab=c("probability threshold"),ylab="TSS")
abline(v=thresh.maxTSS)
dev.off()

## SDA based on maxTSS
SDA <- prob.p
SDA[SDA>=thresh.maxTSS] <- 1
SDA[SDA<thresh.maxTSS] <- 0
pdf(file="results/SDA.pdf")
plot(SDA, addfun=fun.obs,legend=FALSE)
dev.off()

## Export the result
writeRaster(SDA,filename="results/SDA.tif",overwrite=TRUE)
write.table(round(tss.df,2),file="results/TSS.txt",row.names=FALSE,sep="\t")

## Estimating SDA area (in km2)
n.pix <- sum(values(SDA),na.rm=TRUE)
area.SDA <- n.pix*res(SDA)[1]*res(SDA)[2]/1.0e+6

##==========================================================================
##
## Final SDA map with Google Map background (with ggmap::get_map() function)
##
##==========================================================================

## Reproject in Lat/Long (epsg:4326)
GM.crs <- CRS("+init=epsg:4326")
## SDA
SDA.GM <- projectRaster(from=SDA,crs=GM.crs,method="ngb")
SDA.df <- as.data.frame(SDA.GM,xy=TRUE,na.rm=TRUE)
loc <- extent(SDA.GM)[c(1,3,2,4)]
names(SDA.df) <- c("x","y","pres")
sda.df <- SDA.df[SDA.df$pres==1,]
## Gorilla presence data
Gorilla.GM <- spTransform(Gorilla.obs.sp,GM.crs)
Gorilla.df <- cbind(as.data.frame(Gorilla.GM))
names(Gorilla.df)[c(6,7)] <- c("Long","Lat")
gorilla.df.1 <- Gorilla.df[Gorilla.df$Presences>=1,]
gorilla.df.0 <- Gorilla.df[Gorilla.df$Presences==0,]
## Plot with ggplot2
bg <- ggmap(get_map(location=loc,zoom=7,maptype="terrain",source="google",color="bw"))
g.sda <- geom_raster(mapping=aes(x,y,fill=pres),data=sda.df,alpha=0.5)
g.pres <- geom_point(aes(Long,Lat),data=gorilla.df.1,color="black",size=1,shape=16,alpha=0.5)
g.abs <- geom_point(aes(Long,Lat),data=gorilla.df.0,color="black",size=1,shape=".",alpha=0.5)
gg.plot <- bg + coord_cartesian() + coord_equal() +
  g.sda + scale_fill_continuous(guide = FALSE) + 
  g.pres ## + g.abs ## If we want the absences to be plot
## Save as png image file
ggsave(filename="SDA_ggmap.png",plot=gg.plot,device="png",path="results/",width=10,height=10,units="cm",dpi=300)

# ## Alternative plot with raster::gmap() function
# library(rgdal)
# library(dismo)
# ## Reproject SDA in Google PseudoMercator (epsg:3857)
# GM.crs <- CRS("+init=epsg:3857")
# SDA.GM <- projectRaster(from=SDA,crs=GM.crs,method="ngb")
# Gorilla.GM <- spTransform(Gorilla.obs.sp,GM.crs)
# ## Plot
# png(file="results/SDA_gmap.png",height=1000,width=800)
# gorilla.gmap <- gmap(x=SDA.GM,exp=1,type="terrain",zoom=7,scale=2)
# plot(gorilla.gmap)
# plot(SDA.GM,add=TRUE,legend=FALSE,alpha=0.35)
# plot(Gorilla.GM[Gorilla.GM$Presences>0,],pch=3,add=TRUE)
# dev.off()

##========================================================
##
## Uncertainty for SDA
##
##========================================================

## Matrix with confidence interval fo each pixel
prob.p.quant <- apply(mod.ZIB.iCAR[[1]]$prob.p.pred,2,quantile,c(0.025,0.975))
prob.p.m <- apply(mod.ZIB.iCAR[[1]]$prob.p.pred,2,mean)

## Maps
prob.p.025 <- prob.p.mean <- prob.p.975 <- subset(G.stk,1) ## create rasters for predictions
## Assign predicted values
values(prob.p.025)[w] <- prob.p.quant[1,]
values(prob.p.mean)[w] <- prob.p.m
values(prob.p.975)[w] <- prob.p.quant[2,]
## Set NA where no environmental data
values(prob.p.025)[!w] <- values(prob.p.mean)[!w] <- values(prob.p.975)[!w] <- NA
## Stack
prob.p <- stack(prob.p.025,prob.p.mean,prob.p.975)
names(prob.p) <- c("lower bound","mean","upper bound")
## Plot the predictions
pdf(file="results/proba_uncertainty.pdf",width=10,height=4)
plot(prob.p,addfun=fun.obs,legend=TRUE,zlim=c(0,1),nc=3)
dev.off()

## Extract predicted probability of presence
df$prob.p.025 <- extract(prob.p.025,Gorilla.obs.sp)
df$prob.p.mean <- extract(prob.p.mean,Gorilla.obs.sp)
df$prob.p.975 <- extract(prob.p.975,Gorilla.obs.sp)

##================
## SDA
##================

##========
## For 025

## TSS.025 as a function of the threshold
TSS.025 <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
  Index <- Index.fun(df$Suit,df$prob.p.025,thresh.seq[i])
  TSS.025[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS.025
maxTSS.025 <- max(TSS.025,na.rm=TRUE)
w.th <- which(TSS.025==maxTSS.025)
thresh.maxTSS.025 <- mean(thresh.seq[w.th],na.rm=TRUE)

## SDA.025 based on maxTSS.025
SDA.025 <- prob.p.025
SDA.025[SDA.025>=thresh.maxTSS.025] <- 1
SDA.025[SDA.025<thresh.maxTSS.025] <- 0

## Estimating SDA.025 area
n.pix <- sum(values(SDA.025),na.rm=TRUE)
area.SDA.025 <- n.pix*res(SDA.025)[1]*res(SDA.025)[2]/1.0e+6

##=========
## For mean

## TSS.mean as a function of the threshold
TSS.mean <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
  Index <- Index.fun(df$Suit,df$prob.p.mean,thresh.seq[i])
  TSS.mean[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS.mean
maxTSS.mean <- max(TSS.mean,na.rm=TRUE)
w.th <- which(TSS.mean==maxTSS.mean)
thresh.maxTSS.mean <- mean(thresh.seq[w.th],na.rm=TRUE)

## SDA.mean based on maxTSS.mean
SDA.mean <- prob.p.mean
SDA.mean[SDA.mean>=thresh.maxTSS.mean] <- 1
SDA.mean[SDA.mean<thresh.maxTSS.mean] <- 0

## Estimating SDA.mean area
n.pix <- sum(values(SDA.mean),na.rm=TRUE)
area.SDA.mean <- n.pix*res(SDA.mean)[1]*res(SDA.mean)[2]/1.0e+6

##========
## For 975

## TSS.975 as a function of the threshold
TSS.975 <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
  Index <- Index.fun(df$Suit,df$prob.p.975,thresh.seq[i])
  TSS.975[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS.975
maxTSS.975 <- max(TSS.975,na.rm=TRUE)
w.th <- which(TSS.975==maxTSS.975)
thresh.maxTSS.975 <- mean(thresh.seq[w.th],na.rm=TRUE)

## SDA.975 based on maxTSS.975
SDA.975 <- prob.p.975
SDA.975[SDA.975>=thresh.maxTSS.975] <- 1
SDA.975[SDA.975<thresh.maxTSS.975] <- 0

## Estimating SDA.975 area
n.pix <- sum(values(SDA.975),na.rm=TRUE)
area.SDA.975 <- n.pix*res(SDA.975)[1]*res(SDA.975)[2]/1.0e+6

## Stack
SDA <- stack(SDA.025,SDA.mean,SDA.975)
names(SDA) <- c("lower bound","mean","upper bound")
## Plot the predictions
pdf(file="results/SDA_uncertainty.pdf",width=10,height=4)
plot(SDA, addfun=fun.obs,legend=FALSE,nc=3)
dev.off()

## Export SDA results
SDA.ci <- c(area.SDA.mean,area.SDA.025,area.SDA.975)
SDA.df <- data.frame(SDA=SDA.ci)
row.names(SDA.df) <- c("mean","lower","upper")
write.table(SDA.df,file="results/SDA_uncertainty.txt",row.names=TRUE,sep="\t")

##==================================================================
##
## Relationship between gorilla density and nest site encounter rate
##
##==================================================================

## Load encounter rate and density data
erate.dens <- read.csv("data/gorillas/gorillas-erate-dens.csv", header=TRUE)
names(erate.dens) <- c("Site","erate","Density")
## Regression
lmfit <- lm(Density~erate, data=erate.dens)
model.error <- predict(object=lmfit, interval="confidence")
## Plot
pdf(file="results/erate_density_relationship.pdf")
par(cex=1.4, lwd=2, mar=c(4,4,1,1))
plot(erate.dens$Density~erate.dens$erate,
     xlim=c(0,2.5),ylim=c(0,1.5),
     axes=FALSE,
     xlab="Nest site encounter rate (nb. per km walked)", 
     ylab='Gorilla density (nb. per km2)',pch=19)
lines(erate.dens$erate,predict(lmfit))
lines(erate.dens$erate, model.error[,2], lty=2)
lines(erate.dens$erate, model.error[,3], lty=2)
axis(side=1,at=seq(0,2.5,by=0.5),labels=seq(0,2.5,by=0.5))
axis(side=2,at=seq(0,1.5,by=0.5),labels=seq(0,1.5,by=0.5))
dev.off()

##=============================================================
##
## Weighted mean gorilla density from nest site encounter rates
##
##=============================================================

## Encounter rate data (nest encountered per km walked)
erate.df <- read.csv("data/gorillas/gorillas-erate.csv")
erate.df$erate <- erate.df$NestSiteErate
dens.df <- as.data.frame(predict(lmfit,newdata=erate.df,interval="confidence"))
names(dens.df) <- c("dens.fit","dens.lwr","dens.upr")
dens.df$dens.lwr[dens.df$dens.lwr<0] <- 0.0 ## Correcting for negative densities
erate.dens.df <- cbind(erate.df[-c(7)],dens.df)

## Weighted mean
dens.ci <- apply(dens.df,2,weighted.mean,w=erate.df$Area,na.rm=TRUE)

## Export results
write.table(round(dens.ci,3),file="results/density_uncertainty.txt")
write.table(erate.dens.df,file="results/density_uncertainty_site.txt")

##========================================================
##
## Uncertainty for gorilla population size
##
##========================================================

## Population 95% confidence interval
pop.ci <- SDA.ci*dens.ci

## Export population results
pop.df <- data.frame(pop=pop.ci)
row.names(pop.df) <- c("mean","lower","upper")
write.table(round(pop.df),file="results/pop_uncertainty.txt",row.names=TRUE,sep="\t")

## Decrease in population since 1994
decrease.ci <- 1-(pop.df/16902) ## 16902 gorillas in 1994 from Hall et al. 1998 Oryx Tab. 1)
decrease.df <- data.frame(decrease=decrease.ci)
row.names(decrease.df) <- c("mean","max","min")
names(decrease.df) <- c("decrease")
write.table(round(decrease.df,2),file="results/pop_decrease_uncertainty.txt",row.names=TRUE,sep="\t")

## Time computation
t.stop <- Sys.time() ## Stop the clock
t.diff <- difftime(t.stop,t.start,units="min")
cat(paste0("Computation time (min): ",round(t.diff,2)),file="results/computation_time.txt")

## Save objects
save(list=c("erate.dens.df","dens.ci","SDA.df","pop.df","decrease.df"),file="results/gorillas.rda")

##==================
## Knit short report

## Set knitr chunk default options
opts_chunk$set(echo=FALSE, cache=FALSE,
               results="hide", warning=FALSE,
               message=FALSE, highlight=TRUE,
               fig.show="hide", size="small",
               tidy=FALSE)

## Knit and translate to pdf
dir.create("report")
render("gorillas.Rmd",output_format=c("pdf_document"),output_dir="report") 

##===========================================================================
## End of script
##===========================================================================
