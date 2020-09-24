library(EasyABC)
set.seed(1)

## artificial example to show how to use the 'ABC_sequential' function.
## defining a simple toy model:

toy_model<-function(par){ 
val = c(0,0,0,0)
x = c(1,2,3,4)
rep=10
for (indx in seq(1, rep, by=1)) {
 val[1] =  val[1] + x[1] * par[1] + par[2] + rnorm(1,0,0.1)
 val[2] =  val[2] + x[2] * par[1] + par[2] + rnorm(1,0,0.1)
 val[3] =  val[3] + x[3] * par[1] + par[2] + rnorm(1,0,0.1)
 val[4] =  val[4] + x[4] * par[1] + par[2] + rnorm(1,0,0.1)
}
val[1] = val[1]/rep
val[2] = val[2]/rep
val[3] = val[3]/rep
val[4] = val[4]/rep
val }

PhysiCel_model<-function(par){
	cmd = "./motility2D.exe 1 12.3 0.5 0.2"
	system(cmd)
}

cmd = "./motility2D.exe 1 12.3 0.5 0.2"
system(cmd, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
 
## define prior information
toy_prior=list(c("unif",0,20),c("unif",0,1),c("unif",0,2))

## define the targeted summary statistics
sum_stat_obs= c(0.357654864,0.197460223,0.279777057,0.247467603,0.284797896,0.418199135,0.26002244 ,0.315156145,0.274495159,0.250480047,0.321092513,0.364444009,0.417047732,0.249118033,0.400001263,0.364108632,0.431024186,0.405475457,0.301590111,0.371688764,0.661872856,0.347930995,0.323111452,0.401155024,0.309753994,0.337195674,0.407336223,0.29868571 ,0.359246692,0.454747911,0.436872206,0.258024436,0.36680235 ,0.319961835,0.351699627,0.37039482 ,0.341236939,0.314352688,0.444439582,0.291770135,0.387477023,0.364976286,0.292397876,0.442622385,0.30894623 ,0.380256461,0.333063164,0.406446208,0.356285967,0.326961093,0.359871903,0.304948147,0.397588977,0.293278911,0.291908886,0.275263093,0.290145419,0.383916117,0.356008951,0.43765396)
	 
numSamp = 200

## to perform the Del Moral et al. (2012)'s method:
##
alpha_delmo=0.5
tolerance=0.1
ABC_Delmoral<-ABC_sequential(method="Delmoral", model=toy_model, prior=toy_prior, nb_simul=numSamp, summary_stat_target=sum_stat_obs, alpha=alpha_delmo, tolerance_target=tolerance)

d1=density(ABC_Delmoral$param[1:numSamp,1],weights=ABC_Delmoral$weights)
d2=density(ABC_Delmoral$param[1:numSamp,2],weights=ABC_Delmoral$weights)
d3=density(ABC_Delmoral$param[1:numSamp,3],weights=ABC_Delmoral$weights)
hist(ABC_Delmoral$param[1:numSamp,],breaks=10)
plot(d2,type="l",colblue",lwd=2,xlab="Posterior Distribution",main="",cex.lab=1.5)	 
hist,breaks=10)

# Writing mtcars data
write.table(ABC_Delmoral, file = "OutCalib.txt", sep = "\t",row.names = TRUE, col.names = TRUE)


