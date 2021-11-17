setwd("C:/Users/issa/Desktop/TA/2ata/STA203/projet") #à personnaliser
rm(list=objects())
graphics.off()

library(ggplot2)
library(FactoMineR)
library(glmnet)
library(MASS)
library(pls)
library(ROCR)
library(caret)


##########################################################
#______________Partie 2 Analyse exploratoire______________#
##########################################################


###quest1

## importation des données d'apprentissage
load("cookie.app.RData")
## importation des données de test
load("cookie.val.RData")


#variables explicatives
xtrain=cookie.app[,-1]    #apprentissage
xtest=cookie.val[,-1]     #test

#réponses
ytrain=cookie.app[,1]     #apprentissage
ytest=cookie.val[,1]      #test

#nombre d'observations
n=nrow(xtrain)            

#données centrées réduites
Xtrain=scale(xtrain)/sqrt((n-1)/n)
apply(Xtrain^2,2,mean)

#boxplot des variables explicatives
boxplot(xtrain,border=c(1:ncol(xtrain)))

#traçage des courbes
matplot(t(xtrain),type='l')

#corrélation
library(corrplot)
corrplot(cor(xtrain),diag = FALSE,
         tl.pos = "td", tl.cex = 0.5, method = "color", type = "upper")





###quest2


res=PCA(xtrain,ncp=39,graph = FALSE)

val_prop=res$eig[,1]
sum(val_prop)
barplot(val_prop)
abline(h=1/700,col="red")

which(val_prop>1/700) 
#une prise en compte des 29 axes principaux est suffisante
#pour la reconstruction du nuage

#présentation du nuage sur les 6 premiers axes principaux
par(mfrow=c(3,2))

for (i in c(1:6)){
  
  #plot.new()
  abs=res$var$coord[,i]
  ord=rep(0,700)
  plot(abs,ord,pch=3,col="red",xlim=c(-1,1),xlab=NA, ylab=NA, main=paste("projection du nuage 
  sur l'axe principal N°,",i))
  abline(h=0)
 
}


###quest3


reconstruct=function(res,nr,Xm,Xsd){
  #res: objet retourné de la fonction PCA
  #nr: nombre d'axes à utiliser pour la reconstruction
  #Xm: vecteur des moyennes des variables
  #Xsd: vecteur des corrélations entre les variables
  vp=res$eig[c(1:nr),1]         #les valeurs propres des nr premiers axes principaux
  CP=res$ind$coord[,c(1:nr)]   #les nr premières composantes principales
  l=matrix(rep(1,40),nrow=40)%*%sqrt(1/vp) 
  vectp=l*CP
  coord=as.matrix(res$var$coord[,c(1:nr)],nrow=700)
  X=vectp%*%t(coord)
  X=t(apply(X,1,function(x){x*Xsd+Xm}))
  return(X)
}


Xm=apply(xtrain,2,mean)
Xsd=apply(xtrain,2,sd)

#reconstruction totale du nuage
par_default <- par(mfrow = c(3, 2))
for (nr in c(1,2,3,4,5,39)){
  r=reconstruct(res,nr,Xm,Xsd)
  RMSE=sqrt(sum((xtrain-r)**2)/(700*40)) 
  MAE=sum(abs(xtrain-r))/(700*40)
  matplot(t(r),type='l',main=paste("nr=",nr ,"RMSE=",round(RMSE,4), "MAE",round(MAE,3)),xlab=NA,ylab=NA)
}


#reconstruction de la variable X24
for (nr in c(1,2,3,4,5,39)){
  r=reconstruct(res,nr,Xm,Xsd)[,24] #reconstruction de X24
  RMSE=sqrt(sum((xtrain[,24]-r)**2)/(40)) 
  MAE=sum(abs(xtrain[,24]-r))/(40)
  matplot(r,type='l',main=paste("nr=",nr ,"RMSE=",round(RMSE,4), "MAE",round(MAE,3)),xlab=NA,ylab=NA,col="red")
}
par(par_default)





##########################################################
#______________Partie 3 Régression pénalisée_____________#
##########################################################

### quest 1
par(mfrow=c(1,1))

grid=10^seq(6,-10,length=100)
ridge.fit=glmnet(xtrain,ytrain,alpha=0,lambda=grid)
plot(grid,coef(ridge.fit)[1,],type='l')
plot(coef(ridge.fit)[1,],type='l')


## Calcul de l'estimé de l'intercept par la formule de la section 1.2
theta=coef(ridge.fit)[-1,] ## les valeurs de theta (sans intercept) pour les valeurs de grid
X_bar=as.matrix(apply(xtrain,2,mean),nrow=1)
inters=mean(ytrain)-t(X_bar)%*%theta
lines(inters,type="b",col="brown",pch=4)
#On remarque une superposition des deux courbes


ridge.fit1=glmnet(xtrain,ytrain-mean(ytrain),alpha=0,lambda=grid)
ridge.fit2=glmnet(scale(xtrain,center = T,scale=F),ytrain,alpha=0,lambda=grid)
ridge.fit3=glmnet(scale(xtrain,center = T,scale=F),scale(ytrain,center=T,scale=F),alpha=0,lambda=grid)

plot(coef(ridge.fit)[1,],type='l',ylab="intercept")
lines(coef(ridge.fit1)[1,],type='l',col="red")
lines(coef(ridge.fit2)[1,],type='l',col="blue")
lines(coef(ridge.fit3)[1,],type='l',col="green")
legend("topright", legend=c("modèle initial", "ytrain centré","xtrain centré","xtrain et ytrain centrés"),
      col=c("black","red", "blue","green"), lty=1, cex=0.8)



#estimation de theta quand lambda tend vers 0
lm=0.001
mod=glmnet(Xtrain,(ytrain-mean(ytrain))/sd(ytrain),alpha=0,lambda=lm)
A=function(x,lambda) 
{
  e=eigen(t(x)%*%x)
  INV=e$vectors%*%((1/(e$values+lambda))*t(e$vectors))
  if (!is.element(Inf,INV)) return (INV%*%t(x))
}

A0=A(as.matrix(Xtrain),lm)
theta_lim=A0%*%((ytrain-mean(ytrain))/sd(ytrain))
sum(abs(coef(mod)[-1,]-theta_lim))/700





### quest 2

model=lm.ridge(ytrain~.,data=xtrain,lambda=grid)


#calcul à la main
M=as.matrix(Xtrain)
ridg=function(s){
  solve(t(M)%*%M+s*diag(rep(1,700)))%*%t(M)%*%as.matrix((scale(ytrain)))
}
res=sapply(grid,ridg)

#différence en valeur absolue entre les valeurs
mean(abs(model$coef-coef(ridge.fit)[-1,]))

## quest 3
#on choisit de s'arrêter à grid[60]~0

set.seed(123)
grid=grid[1:60]
B = 4
folds=cvsegments(n,B, type="random")

# matrice des résultats Bx10
cv.errors=matrix(NA,B,60, dimnames=list(NULL, paste0("Dim",1:60)))
cv.bestmod=matrix(NA,B,60,dimnames=list(NULL, paste0("Dim",1:60)))

for(b in 1:B){ 
  subsetb=unlist(folds[[b]])
  #X=scale(xtrain[-subsetb,],center = T,scale = F)  
  #Y=ytrain[-subsetb]-mean(ytrain[-subsetb])
  #model_b=glmnet(X,Y,alpha=0,lambda=grid)
  #model_b=tapply(as.data.frame(x=X,y=Y,alpha=0),grid,glmnet)
  X=xtrain[-subsetb,]
  Y=ytrain[-subsetb]
  for(j in 1:60){
    model=glmnet(X,Y,alpha=0,lambda=grid[j])
    pred=predict(model,s=grid[j],newx=as.matrix(xtrain[subsetb,]))  
    cv.errors[b,j]=mean( (ytrain[subsetb]-pred)^2)
    }
}

# on moyenne colonne par colonne, sur l'ensemble des B lignes
mean.cv.errors = apply(cv.errors,2,mean)
sd.cv.errors = apply(cv.errors,2,sd)
plot(sd.cv.errors)
par(mfrow=c(1,1))
plot(log(grid),mean.cv.errors,ylim=c(4,18))

q=qnorm(0.975)
IC_sup=mean.cv.errors+sd.cv.errors
IC_inf=mean.cv.errors-sd.cv.errors

segments(x0=log(grid),y0=IC_inf,x1=log(grid),y1=IC_sup,col='red')

jbest = which.min(mean.cv.errors )  #l'indice de la plus petite erreur

mean.cv.errors[jbest]  
log(c(mean.cv.errors[jbest]))  




#### avec cv.glmnet() ####
set.seed(123)
cv.out=cv.glmnet(as.matrix(xtrain),ytrain,alpha=0,nfolds = 4,lambda=grid) 
plot(cv.out)

optim_grid=cv.out$lambda.min   #la valeur optimale de k
log(optim_grid)

#réajustement sur les données d'apprentissage
out=glmnet(xtrain,ytrain,alpha=0,lambda=optim_grid)

#erreur de généralisation
ridge.pred=predict(out,s=optim_grid,newx=as.matrix(xtest))
mean((ridge.pred-as.matrix(ytest))^2)


par(mfrow=c(1,1))





##########################################################
#_____Partie 4 Régression logistique pénalisée___________#
##########################################################

### quest1

# étude en données individuelles
# Z_i ~B(1,p(x_i)) indep 
# logit(p(x_i))=x_i beta

ztrain=ifelse(ytrain>18,1,0)
ztest=ifelse(ytest>18,1,0)

mean(ztrain)
mean(ztest)
#les données se trouvent être équilibrés. Il y a presque
#autant de proportion de 0 que de 1 dans ztrain et ztest. 
### quest2

#ridge
mod_ridge=cv.glmnet(as.matrix(xtrain),ztrain,alpha=0,type.measure = "auc")
plot(mod_ridge)
log(mod_ridge$lambda.min)
#lasso
mod_lasso=cv.glmnet(as.matrix(xtrain),ztrain,alpha=1,type.measure = "auc")
plot(mod_lasso)
log(mod_lasso$lambda.min)

#prédictions sur les données test
pred_ridge=predict(mod_ridge,s=mod_ridge$lambda.min,newx=as.matrix(xtest))
z_ridge=pred_ridge>=0.5
err_ridge=mean(z_ridge!=ztest)

pred_lasso=predict(mod_lasso,s=mod_lasso$lambda.min,newx=as.matrix(xtest))
z_lasso=pred_lasso>=0.5
err_lasso=mean(z_lasso!=ztest)

#err_lasso(0.03125)<err_ridge(0.15625)




### quest 3

par(mfrow=c(1,2))

# ridge
pred_test=prediction(pred_ridge,ztest)
pred_train=prediction(predict(mod_ridge,s=mod_ridge$lambda.min,newx=as.matrix(xtrain)),ztrain)
plot(performance(pred_test,"sens","fpr"),xlab="",col=2,main="ROC pour ridge") 
plot(performance(pred_train,"sens","fpr"),xlab="",col=3,add=T) 

#lasso
pred_test=prediction(pred_lasso,ztest)
pred_train=prediction(predict(mod_lasso,s=mod_lasso$lambda.min,newx=as.matrix(xtrain)),ztrain)
plot(performance(pred_test,"sens","fpr"),xlab="",col=2,main="ROC pour lasso") 
plot(performance(pred_train,"sens","fpr"),xlab="",col=3,add=T) 

#On remarque que les courbes ROC du modèle logistique
#pénalisé par lasso est plus performant sur les données
#test. Ses courbes ROC sont plus proches du 
#modèle parfait que celles du modèle pénalisée par ridge.

#On ne peut pas tester l'adéquation ici car on est 
#en données individuelles.








