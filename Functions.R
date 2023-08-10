### Additional used functions in "A novel platform to attract wild Tanimbar corellas (Cacatua goffiniana) for behavioural research"
### provided and written by Roger Mundry


###fe.re.tab - function (Lines: 19 - 269)
###diagnostics.plot - function (Lines: 275 - 307)
###ranef.diagn.plot - function (Lines: 313 - 366)
###glmm.model.stab - function (Lines: 372 - 570)
###m.stab.plot - function (Lines: 576 - 610
###boot.glmm.pred - function (Lines: 616 - 837)
###overdisp.test - function (Lines: 843 - 942)
###factor.int.plot - function (Lines: 948 - 1191)
###glmmTMB_stability - function (Lines: 1199 - 1417)
###glmmTMB_stability - function (Lines: 1422 - 1819)

###########################################################################################################
# fe.re.tab - function
###########################################################################################################
fe.re.tab<-function(fe.model, re, other.vars=NULL, data, treat.covs.as.factors=c(NA, F, T)){
	print("please read the documentation in the beginning of the script")
	#function helping in determinining which random slopes are needed
	#last updated: 2018, June 6
	#latest updates:
		#major revision of how interactions are treated; reveals in the summary
			#for combinations factors whether the number of observations is >1 (seperately for each level of the random effect)
			#for combination of covariates whether the number of unique combinations per level of random effect is >2
			#for combinations of factors and covariates whether the number unique values per covariate is larger than 2 (when treated as covariate)
				#or larger than 1 (when treated as a factor), separately for each commbination of levels of fixed and random effects factors
		#treat.covs.as.factors can be NA in which case covariates aretreates as such and as factors
			#this is the default now
		#output gives information about which fixed effects are covariates and factors aso for interactions
			#and also how covariates were treated
	#input/arguments:
		#data: a dataframe with all relevant variables (including the response)
		#fe.model: character; the model wrt the fixed effect (including the response); e.g., "r~f1*c*f2"
		#re: character; either a vector with the names of the random effects (e.g., c("re1", "re2")) or a random intercepts expression (e.g., "(1|re1)+(1|re2)")
		#other.vars: character, optional; a vector with the names of variables which are to be kept in the data considered and returned
		#treat.covs.as.factors: logical, when set to TRUE covariates will be treated like factors (see value/summary for details)
	#value: list with the following entries:
		#detailed: list with cross-tabulations ffor each combination of (main) fixed and random effect 
		#summary: list tables...
			#telling for each combination of (main) fixed and random effect...
				#the number of levels of the random effect with a given number of unique values of the fixed effect (in case of a covariate)
				#the number of levels of the random effect with a given number of levels of the fixed effect for which at least two cases do exist (in case of a factor)
			#telling for each combination of interaction and random effect...
				#the combination of the above two informations, i.e., the number of individuals with a given number of unique values of the covariate
					#and a given number of factor levels for which at least two cases exist
		#data: data frame containing all relevant variables (i.e., response, fixed and random effects as well as those indicated in other.vars (e.g., offset terms)
			#also includes columns for dummy variables coding the levels (except the reference level) of all factors
		#pot.terms: length one vector comprising the model wrt the random slopes (for all combinations of fixed and random effects, i.e., most likely some will need to be omitted)
			#note that this comprises only the random slopes but not the correlations between random slopes and intercepts not the random intercept itself
		#pot.terms.with.corr: length one vector comprising the model wrt the random slopes (for all combinations of fixed and random effects, i.e., most likely some will need to be omitted)
			#note that this comprises the random slopes and intercepts  and also all correlations among them
	treat.covs.as.factors=treat.covs.as.factors[1]
	if(sum(grepl(x=re, pattern="", fixed=T))>0 & length(re)==1){#if random effects are handed over as formula
		re=gsub(x=re, pattern="(1|", replacement="", fixed=T)
		re=gsub(x=re, pattern=")", replacement="", fixed=T)
		re=gsub(x=re, pattern=" ", replacement="", fixed=T)
		re=unlist(strsplit(re, split="+", fixed=T))
	}
	fe.model=gsub(x=fe.model, pattern=" ", replacement="", fixed=T)#remove spaces
	model.terms=attr(terms(as.formula(fe.model)), "term.labels")#get individual terms from fixed effects model
	fe.me=model.terms[!grepl(x=model.terms, pattern=":", fixed=T)]#remove interactions
	fe.me=fe.me[!grepl(x=fe.me, pattern="^", fixed=T)]#remove squares terms
	resp=unlist(strsplit(fe.model, split="~", fixed=T))[1]#determine response
	if(substr(resp, start=1, stop=6)=="cbind("){
		resp=gsub(x=resp, pattern="cbind(", replacement="", fixed=T)
		resp=gsub(x=resp, pattern=")", replacement="", fixed=T)
		resp=gsub(x=resp, pattern=" ", replacement="", fixed=T)
		resp=unlist(strsplit(resp, split=",", fixed=T))
	}
	xx=setdiff(c(resp, fe.me, re, other.vars), names(data))
	if(length(xx)>0){
		stop(paste(c("error: preditor(s) missing in the data is/are ", paste(xx, collapse=", ")), collapse=""))
	}
	data=droplevels(as.data.frame(na.omit(data[, c(resp, fe.me, re, other.vars)])))#keep complete data wrt all relevant variables
	model.terms=model.terms[!grepl(x=model.terms, pattern="^", fixed=T)]#remove squares terms
	modes=rep(NA, length(model.terms))#initialize vector storing whether PVs are factors or not
	effect=c(rep("main", length(fe.me)), rep("int", length(model.terms)-length(fe.me)))#create vector telling for each model term whether it is a main effect of not
	for(i in 1:ncol(data)){#for all columns in data
		if(is.character(data[, i])){data[, i]=as.factor(data[, i])}#turn character column into factor
		modes[i]=class(data[, i])#and determine its class
	}
	names(modes)=names(data)#name 'm
	to.do=data.frame(expand.grid(re=re, fe=fe.me))#create data frame with one column for each combination of fixed main and random effect
	to.do$re=as.character(to.do$re)#reformat to character (for later addressing by name)
	to.do$fe=as.character(to.do$fe)#reformat to character (for later addressing by name)
	res.detailed=lapply(1:nrow(to.do), function(xrow){#create detailed results by lapply-ing over the rows of to.do
		table(data[,to.do$re[xrow]], data[,to.do$fe[xrow]])#cross tabulate the respective fixed and random effect
	})
	names(res.detailed)=paste(to.do$fe, to.do$re, sep="_within_")#name it
	modes=modes[as.character(to.do$fe)]
	modes=gsub(x=modes, pattern="integer", replacement="numeric")
	to.do$modes.cons=modes
	to.do$modes=modes
	if(is.na(treat.covs.as.factors) & any(to.do[, "modes.cons"]=="numeric")){
		xx=to.do[to.do[, "modes.cons"]=="numeric", ]
		xx["modes.cons"]="factor"
		to.do=rbind(to.do, xx)
		to.do=to.do[order(to.do[, "re"], to.do[, "fe"]), ]
	}else if(!is.na(treat.covs.as.factors) & treat.covs.as.factors){
		to.do$modes.cons="factor"
	}
	to.do=to.do[order(to.do$fe), ]
	res.summary=lapply(1:nrow(to.do), function(xrow){#begin with creating summary results by lapply-ing over the rows of to.do
		#browser()
		xwhat=paste(to.do[xrow, c("fe", "re")], collapse="_within_")
		if(to.do$modes.cons[xrow]!="factor"){#if fixed effect is not a factor
			ires=table(apply(res.detailed[[xwhat]]>0, 1, sum))#determine number of levels of the random effect per number of unique cases of the fixed effect
		}else{#if fixed effect is a factor
			ires=table(apply(res.detailed[[xwhat]]>1, 1, sum))#determine number of levels of the random effect per number of unique cases of the fixed effect with at least two observations
		}
		ires=c(ires, tot=nrow(res.detailed[[xwhat]]))
	})
	xx=to.do[, "modes.cons"]
	xx[xx=="numeric"]="covariate"
	xx[xx=="integer"]="covariate"
	xx[to.do[, "modes.cons"]!=to.do[, "modes"]]="covariate as factor"
	xx=paste("(", xx, ")", sep="")
	names(res.summary)=paste(paste(to.do$fe, to.do$re, sep="_within_"), xx, sep=" ")#name it
	#append 'factor' or 'covariate' to the names:
	#names(res.summary)[modes[to.do$fe]=="factor"]=paste(names(res.summary)[modes[to.do$fe]=="factor"], "(factor)", sep=" ")
	#names(res.summary)[!modes[to.do$fe]=="factor"]=paste(names(res.summary)[!modes[to.do$fe]=="factor"], "(covariate)", sep=" ")
	to.do=data.frame(expand.grid(re=re, int=setdiff(model.terms, fe.me)))#create data frame with one row for each combination of interaction and random effect
	if(nrow(to.do)>0){
		#add two columns denoting original mode:
		to.do$mode=unlist(lapply(strsplit(as.character(to.do$int), split=":", fixed=T), function(x){
			paste(modes[x], collapse=":")
		}))
		to.do$mode.cons=to.do$mode
		if(is.na(treat.covs.as.factors)){#if is.na(treat.covs.as.factors)
			xx=to.do[grepl(x=to.do$mode, pattern="numeric"), ]#add rows duplicating interactions involving covariates (to treat them as covariate and factor
			xx$mode.cons=gsub(x=xx$mode.cons, pattern="numeric", replacement="factor")
			to.do=rbind(to.do, xx)
			to.do=to.do[order(to.do$re, to.do$int), ]
		}else if(treat.covs.as.factors){#if treat.covs.as.factors
			to.do$mode.cons=gsub(x=to.do$mode.cons, pattern="numeric", replacement="factor")#change mode to be considered
		}
		to.do[, "int"]=as.character(to.do[, "int"])#reformat to character (for later addressing by name)
		to.do[, "re"]=as.character(to.do[, "re"])#reformat to character (for later addressing by name)
		##this needs to be rewritten!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		to.add=lapply(1:nrow(to.do), function(xrow){#treat combinations of fixed and random effect by lapply-ing over the rows of to.do
			#browser()
			iterms=unlist(strsplit(to.do[xrow, "int"], split=":", fixed=T))#determine main effects involved in the interaction...
			imodes=unlist(strsplit(to.do[xrow, "mode.cons"], split=":", fixed=T))#... and their classes
			imodes.orig=unlist(strsplit(to.do[xrow, "mode"], split=":", fixed=T))#... and their classes
			if(all(imodes.orig=="factor")){
				comb.fac=lapply(iterms[imodes=="factor"], function(x){paste(x, data[, x], sep="")})
				comb.fac=matrix(unlist(comb.fac), ncol=length(comb.fac), byrow=F)
				comb.fac=apply(comb.fac, 1, paste, collapse="__")
				ires=table(data[,to.do$re[xrow]], comb.fac)
				class(ires)="matrix"
				ires.summary=aggregate(1:nrow(ires), data.frame(ires>1), length)
				xx=as.matrix(ires.summary[, -ncol(ires.summary)])
				xx[xx]=">1"
				xx[xx!=">1"]="!>1"
				ires.summary=data.frame(xx, freq=ires.summary$x)
				ires=list(detailed=ires, summary=ires.summary)
			}else if(any(imodes.orig=="factor")){
				comb.fac=lapply(iterms[imodes.orig=="factor"], function(x){paste(x, data[, x], sep="")})
				comb.fac=matrix(unlist(comb.fac), ncol=length(comb.fac), byrow=F)
				comb.fac=apply(comb.fac, 1, paste, collapse="__")
				covs=lapply(iterms[imodes.orig=="numeric"], function(x){data[, x]})
				if(all(imodes=="factor")){
					ires=lapply(covs, function(x){
						x=tapply(x, list(data[, to.do[xrow, "re"]], comb.fac), function(x){
							sum(table(x)>1)
						})
						x[is.na(x)]=0
						return(x)
					})
				}else{
					ires=lapply(covs, function(x){
						x=tapply(x, list(data[, to.do[xrow, "re"]], comb.fac), function(x){
							length(unique(x))
						})
						x[is.na(x)]=0
						return(x)
					})
				}
				xx=matrix(unlist(lapply(ires, c)), ncol=length(ires), byrow=F)
				xx=apply(xx, 1, paste, collapse="/")
				xx=matrix(xx, ncol=ncol(ires[[1]]), byrow=F)
				colnames(xx)=colnames(ires[[1]])
				rownames(xx)=rownames(ires[[1]])
				ires.detailed=xx
				if(all(imodes=="factor")){
					xx=matrix(unlist(lapply(ires, c))>=2, ncol=length(ires), byrow=F)
					xx[xx]=">1"
					xx[xx!=">1"]="!>1"
				}else{
					xx=matrix(unlist(lapply(ires, c))>=3, ncol=length(ires), byrow=F)
					xx[xx]=">2"
					xx[xx!=">2"]="!>2"
				}
				xx=apply(xx, 1, paste, collapse="/")
				xx=matrix(xx, ncol=ncol(ires[[1]]), byrow=F)
				colnames(xx)=colnames(ires[[1]])
				rownames(xx)=rownames(ires[[1]])
				ires.summary=aggregate(1:nrow(xx), data.frame(xx), length)
				colnames(ires.summary)[ncol(ires.summary)]="freq"
				ires=list(detailed=ires.detailed, summary=ires.summary)
			}else{
				xx=lapply(iterms, function(x){
					apply(res.detailed[[paste(c(x, to.do[xrow, "re"]), collapse="_within_")]]>0, 1, sum)
				})
				yy=matrix(unlist(xx), ncol=length(xx), byrow=F)
				colnames(yy)=iterms
				rownames(yy)=names(xx[[1]])
				ires=yy
				ires.summary=aggregate(1:nrow(ires), data.frame(ires), length)
				colnames(ires.summary)=c(colnames(ires), "freq")#[ncol(ires.summary)]
				ires=list(detailed=ires, summary=ires.summary)
				#browser()
			}
			return(ires)
		})
		xx=paste(to.do$int, to.do$re, sep="_within_")#gsub(x=paste(to.do$int, to.do$re, sep="_within_"), pattern=":", replacement="_", fixed=T)#name it
		mode.mat=to.do[, "mode"]
		mode.mat=gsub(x=mode.mat, pattern="numeric", replacement="COV")
		mode.mat=gsub(x=mode.mat, pattern="factor", replacement="FAC")
		mode.mat=strsplit(mode.mat, split=":", fixed=T)#), nrow=nrow(to.do), byrow=T)
		mode.cons.mat=to.do[, "mode.cons"]
		mode.cons.mat=gsub(x=mode.cons.mat, pattern="numeric", replacement="COV")
		mode.cons.mat=gsub(x=mode.cons.mat, pattern="factor", replacement="FAC")
		mode.cons.mat=strsplit(mode.cons.mat, split=":", fixed=T)#), nrow=nrow(to.do), byrow=T)
		for(i in 1:length(mode.mat)){
			mode.mat[[i]][mode.mat[[i]]!=mode.cons.mat[[i]]]="COVasFAC"
		}
		#mode.mat[mode.mat=="numeric"]="COV"
		#mode.mat[mode.mat=="factor"]="FAC"
		#mode.mat[mode.mat!=mode.cons.mat]="COVasFAC"
		xx=paste(xx, paste("(", unlist(lapply(mode.mat, paste, collapse=":")), ")", sep=""), sep=" ")
		to.add.2=lapply(to.add, "[[", "summary")
		names(to.add.2)=xx
		res.summary=c(res.summary, to.add.2)#and append to.add to res.summary
		to.add.2=lapply(to.add, "[[", "detailed")
		names(to.add.2)=xx
		res.detailed=c(res.detailed, to.add.2)
	}
	#add columns with the dummy coded factor levels to data:
	to.code=fe.me[modes[fe.me]=="factor"]#determine factors among the fixed effects:
	if(length(to.code)>0){
		coded=lapply(to.code, function(xc){#for each factor
			lapply(levels(data[, xc])[-1], function(xl){#for all levels except the reference level
				as.numeric(data[, xc]==xl)#code it
			})
		})
		coded=matrix(unlist(coded), nrow=nrow(data), byrow=F)#reformat to matrix
		xnames=unlist(lapply(to.code, function(xc){#determine column names to be given to the matrix...
			paste(xc, levels(data[, xc])[-1], sep=".")
		}))
		colnames(coded)=xnames#... and use 'm
		data=data.frame(data, coded)#and append 'm to data
	}else{
		xnames=""
	}
	#now the model expression wrt random slopes; first no correllations between random slopes and intercepts:
	pot.terms=c(xnames, fe.me[modes[fe.me]!="factor"])#create vector with names of dummy variables and names of fixed main effects not being factors
	pot.terms=outer(pot.terms, re, Vectorize(function(x, y){#create matrix with separate model term for each combination of fixed and random effect
		paste(c("(0+", x, "|", y, ")"), collapse="")
	}))
	pot.terms=paste(c(t(pot.terms)), collapse="+")#put them all in a single entry
	#now the random slopes part of the model including random intercepts and slopes and their correlation:
	pot.terms.with.corr=paste(c(xnames, fe.me[modes[fe.me]!="factor"]), collapse="+")#get all fixed effects terms together...
	pot.terms.with.corr=paste(paste("(1+", pot.terms.with.corr, "|", re, ")", sep=""), collapse="+")#... and paste random effects, brackets and all that 
		#(and everything in a single entry)
	return(list(detailed=res.detailed, summary=res.summary, data=data, pot.terms=pot.terms, pot.terms.with.corr=pot.terms.with.corr))
}
###########################################################################################################

###########################################################################################################
# diagnostics.plot - function
###########################################################################################################
diagnostics.plot<-function(mod.res, col=grey(level=0.25, alpha=0.5), size.fac=1, weights=NULL){
	if(is.null(weights)){
		if(class(mod.res)[[1]]=="lm" | class(mod.res)[[1]]=="glm"){
			weights=mod.res$model$"(weights)"
			if(is.null(weights)){weights=rep(1, length(residuals(mod.res)))}
		}else if(class(mod.res)[[1]]=="lmerMod"){
			weights=mod.res@resp$weights
		}else{
			weights=rep(1, length(residuals(mod.res)))
		}
	}
	weights=size.fac*weights/max(weights)
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 1, 0.5))
  hist(residuals(mod.res), probability=T, xlab="", ylab="", main="")
  mtext(text="histogram of residuals", side=3, line=0)
  x=seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  qqnorm(residuals(mod.res), main="", pch=19, col=col, cex=sqrt(weights))
  qqline(residuals(mod.res))
  mtext(text="qq-plot of residuals", side=3, line=0)
  plot(fitted(mod.res), residuals(mod.res), pch=19, col=col, cex=sqrt(weights))
  abline(h=0, lty=2)
  mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}

lev.thresh<-function(model.res){
	k=length(coefficients(model.res))
	n=length(residuals(model.res))
 return(2*(k+1)/n)
}
###########################################################################################################

###########################################################################################################
# ranef.diagn.plot - function
###########################################################################################################
ranef.diagn.plot<-function(model.res, QQ=F, col=grey(0.5, alpha=0.75)){
  old.par = par(no.readonly = TRUE)
	if(class(model.res)[[1]]!="glmmTMB"){
		n.plots=sum(unlist(lapply(ranef(model.res), length)))
	}else{
		n.plots=sum(unlist(lapply(ranef(model.res)[["cond"]], ncol)))+sum(unlist(lapply(ranef(model.res)[["zi"]], ncol)))
	}
	x=ifelse(n.plots%%2==1,n.plots+1,n.plots)
	xmat=outer(1:x, 1:(1+x/2), "*")-n.plots
	colnames(xmat)=1:(1+x/2)
	rownames(xmat)=1:x
	xmat=as.data.frame(as.table(xmat))
	xmat=subset(xmat, as.numeric(as.character(xmat$Var1))>=as.numeric(as.character(xmat$Var2)) & xmat$Freq>=0)
	sum.diff=as.numeric(as.character(xmat$Var1))-as.numeric(as.character(xmat$Var2))+xmat$Freq
	xmat=xmat[sum.diff==min(sum.diff),]
	xmat=xmat[which.min(xmat$Freq),]
	par(mfrow=c(xmat$Var2, xmat$Var1))
	par(mar=c(rep(2, 3), 1))
	par(mgp=c(1, 0.5, 0))
	xnames=lapply(ranef(model.res), names)
	xnames=lapply(xnames, gsub, pattern="(Intercept)", replacement="icpt", fixed=T)
	if(class(model.res)[[1]]!="glmmTMB"){
		for(i in 1:length(ranef(model.res))){
			to.plot=ranef(model.res)[[i]]
			for(k in 1:ncol(to.plot)){
				if(QQ){
					qqnorm(to.plot[,k], main="", tcl=-0.25, xlab="", ylab="", pch=19, col=col)
					qqline(to.plot[,k])
				}else{
					hist(to.plot[,k], main="", tcl=-0.25, xlab="", ylab="")
				}
				mtext(text=paste(c(names(ranef(model.res)[i]), xnames[[i]][k]), collapse="_"), side=3)
			}
		}
	}else{
		for(i in c("cond", "zi")){
			if(length(ranef(model.res)[[i]])>0){
				for(j in 1:length(ranef(model.res)[[i]])){
					to.plot=ranef(model.res)[[i]][[j]]
					for(k in 1:ncol(to.plot)){
						if(QQ){
							qqnorm(to.plot[,k], main="", tcl=-0.25, xlab="", ylab="", pch=19, col=col)
							qqline(to.plot[,k])
						}else{
							hist(to.plot[,k], main="", tcl=-0.25, xlab="", ylab="")
						}
						mtext(text=paste(c(names(ranef(model.res)[i]), names(ranef(model.res)[[i]][j]), colnames(to.plot)[k]), collapse=", "), side=3)
					}
				}
			}
		}
	}
	par(old.par)
}
###########################################################################################################

###########################################################################################################
# glmm.model.stab - function
###########################################################################################################
glmm.model.stab<-function(model.res, contr=NULL, ind.cases=F, para=F, data=NULL, use=NULL, n.cores=c("all-1", "all"), save.path=NULL, load.lib=T, lib.loc=.libPaths()){
	print("please carefully evaluate whether the result makes sense, and if not, please contact me")
	#function determining stability of GLMMs (run using lmer or glmer) by excluding levels of random effects, one at a time;
	#supports
		#weights, offset terms and random slopes;
	#does not support
		#correlations between random intercepts and random slopes
		#any terms calculated in the model formula (e.g., log, function I, etc.); but interactions do work
	#latest additions/modifications:
		#copes also with data sets where a level of a (fixed effects) factor is entirely dropped from the data
		#new way of catching warnings (are contained in the detailed output table)
		#includes sample size in the detailed output table
		#catches errors
		#use: an argument taking the names of the random effects for which model stability should be evaluated
			#(useful for models with a random effect having a unique case for each level of the data)
	#written by Roger Mundry
	#modified April 2016 (added dealing with glmer.nb; fixed a bug in the output of the random effects of tthe original model)
	#last modified Mar 2017 (added dealing with fixed effects factor levels dropped when dropping levels of random effects)
	#last modified June 14 2017 (fixed dealing with random effects comprising correlation parameters)
	#last modified July 06 2017 (fixed small bug happening when model revealed an error)
	if(load.lib){library(lme4, lib.loc=lib.loc)}
	n.cores=n.cores[1]
  model.eq=as.formula(as.character(model.res@call)[2])
	if(class(model.res)[1]=="lmerMod"){
    REML=model.res@resp$REML==1
    xfam="gaussian"
  }else{
		xfam=model.res@resp$family$family
		if(model.res@resp$family$family=="Gamma"){
			xfam=model.res@resp$family
		}
  }
  if(grepl(x=xfam[[1]], pattern="Negative Binomial", fixed=T)){
		xfam="neg.bin"
	}
  weights=model.res@resp$weights
  if(length(data)==0){
		ii.data=model.res@frame
		offs.col=grep(x=names(ii.data), pattern="offset(", fixed=T)
		if(length(offs.col)>0){
			for(i in offs.col){
				#ii.data[,i]=exp(ii.data[,i])
				names(ii.data)[i]=substr(x=names(ii.data)[i], start=8, stop=nchar(names(ii.data)[i]) -1)
			}
		}
		wght.col=grep(x=names(ii.data), pattern="(weights)", fixed=T)
		if(length(wght.col)>0){
			names(ii.data)[wght.col]=gsub(x=names(ii.data)[wght.col], pattern="(", replacement="", fixed=T)
			names(ii.data)[wght.col]=gsub(x=names(ii.data)[wght.col], pattern=")", replacement="", fixed=T)
		}else{
			ii.data$weights=weights
		}
	}else{
		ii.data=data.frame(data, weights)
	}
	if(substr(as.character(model.eq)[2], start=1, stop=6)=="cbind("){
		ii.data$weights=1
	}
  #if(xfam=="binomial"){
  #   response=as.character(model.eq)[2]
  #   if(sum(names(ii.data)==response)==0)
  #}
  ranefs=names(ranef(model.res))
	if(length(use)==0){use=ranefs}
	ranefs=ranefs[ranefs%in%use]
  xlevels=lapply(ranefs, function(x){return(as.vector(unique(ii.data[ ,x])))})
  ranefs=rep(ranefs, unlist(lapply(xlevels, length)))
  to.do=cbind(ranefs, unlist(xlevels))
  if(ind.cases){
    ii.data=data.frame(ii.data, ic=as.factor(1:nrow(ii.data)))
    to.do=rbind(to.do, cbind("ic", levels(ii.data$ic)))
  }
	keepWarnings <- function(expr) {
		localWarnings <- list()
		value <- withCallingHandlers(expr,
			warning = function(w) {
				localWarnings[[length(localWarnings)+1]] <<- w
				invokeRestart("muffleWarning")
			}
		)
		list(value=value, warnings=localWarnings)
	}
	get.ranef<-function(x){#function returning sd associated with random effect in a nicely named vector
		x=as.data.frame(x)
		y=lapply(strsplit(as.character(x$grp), split=""), function(y){
			yy=unlist(lapply(strsplit(names(ranef(model.res)), split=""), function(z){
				if(length(z)<=length(y)){
					y=y[1:length(z)]
					return(paste(y[y==z], collapse=""))
				}else{
					return("")
				}
			}))
			yy=intersect(yy, names(ranef(model.res)))
			if(length(yy)==0){yy=paste(y, collapse="")}
			return(yy)
		})
		x$grp=y
		ires=x$sdcor
		#x$var2[is.na(x$var2)]=""
		names(ires)=apply(x[, 1:3], 1, paste, collapse="@")
		names(ires)[x$grp=="Residual"]="Residual"
		return(ires)
	}
	if(length(contr)==0){
		if(class(model.res)[1]=="lmerMod"){
			contr=lmerControl()
		}else{	
			contr=glmerControl()
		}
	}
  ifun=function(x, model.res, to.do, ii.data, contr, get.ranef){
    #sel.ii.data=subset(ii.data, ii.data[,to.do[x, 1]]!=to.do[x, 2])
    sel.ii.data=ii.data[ii.data[,to.do[x, 1]]!=to.do[x, 2], ]
		if(class(model.res)[1]=="lmerMod"){
      sel.ii.res=try(lmer(model.eq, data=sel.ii.data, weights=weights, REML=REML, control=contr), silent=T)
    }else if(xfam[[1]]!="neg.bin"){#(xfam=="binomial" | xfam=="poisson"){
      sel.ii.res=try(glmer(model.eq, data=sel.ii.data, family=xfam, control=contr), silent=T)
    }else{
			sel.ii.res=try(glmer.nb(model.eq, data=sel.ii.data, control=contr), silent=T)
		}
		if(length(save.path)>0){
			est.fixed.effects=fixef(sel.ii.res)
			est.random.effects=as.data.frame(summary(sel.ii.res)$varcor)
			model.warnings=sel.ii.res$warnings
			n=length(residuals(sel.ii.res))
			what=to.do[x, ]
			save(file=paste(c(paste(c(paste(c(save.path, "m"), collapse="/"), x), collapse="_"), ".RData"), collapse=""), list=c("what", "est.fixed.effects", "est.random.effects", "model.warnings", "n"))
		}
    if(class(sel.ii.res)!="try-error"){
			#xx=try(list(fere=c(fixef(sel.ii.res$value), get.ranef(summary(sel.ii.res$value)$varcor)), N=length(residuals(sel.ii.res$value)), warnings=paste(unlist(sel.ii.res$warnings)$message, collapse="/")), silent=T)
			#if(class(xx)=="try-error"){browser()}
      return(list(fere=c(fixef(sel.ii.res), get.ranef(summary(sel.ii.res)$varcor)), N=nrow(sel.ii.res@frame), 
				warnings=
					list(
						lme4=summary(sel.ii.res)$optinfo$conv$lme4$messages,
						opt=summary(sel.ii.res)$optinfo$warnings,
						conv.code=summary(sel.ii.res)$optinfo$conv$opt
					)))
    }else{
      return(list(fere=rep(NA, length(fixef(model.res))+nrow(as.data.frame(summary(model.res)$varcor))), N=nrow(sel.ii.data), 
				warnings=list(lme4=NA, opt=NA, conv.code=NA)))
    }
  }
  if(para){
		#on.exit(expr = parLapply(cl=cl, X=1:length(cl), fun=function(x){rm(list=ls())}), add = FALSE)
		#on.exit(expr = stopCluster(cl), add = T)
    require(parallel)
    cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores!="all"){
			if(n.cores=="all-1"){n.cores=length(cl)-1}
			if(n.cores<length(cl)){
				cl=cl[1:n.cores]
			}
		}
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      library(lme4, lib.loc=lib.loc)
      return(invisible(""))
    })
    all.coeffs=parLapply(cl=cl, X=1:nrow(to.do), fun=ifun, model.res=model.res, to.do=to.do, ii.data=ii.data, contr=contr, get.ranef=get.ranef)
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      return(rm(list=ls()))
    })
    stopCluster(cl=cl)
  }else{
		#browser()
    all.coeffs=lapply(1:nrow(to.do), ifun, model.res=model.res, to.do=to.do, ii.data=ii.data, contr=contr, get.ranef=get.ranef)
  }
  all.n=unlist(lapply(all.coeffs, function(x){x$N}))
  all.warnings=lapply(all.coeffs, function(x){
		x=sapply(x$warnings, paste, collapse="")
	})
  all.warnings=as.data.frame(matrix(unlist(all.warnings), ncol=3, byrow=T))
  colnames(all.warnings)=c("lme4.w", "opt.w", "conv.code")
	############################################################################################################################
	############################################################################################################################
	############################################################################################################################
	############################################################################################################################
	xnames=unique(unlist(lapply(all.coeffs, function(x){names(x$fere)})))
	all.coeff.mat=matrix(NA, ncol=length(xnames), nrow=length(all.coeffs))
	colnames(all.coeff.mat)=xnames
	for(i in 1:length(all.coeffs)){
		all.coeff.mat[i, names(all.coeffs[[i]]$fere)]=all.coeffs[[i]]$fere
	}
	#extract results for original model:
	orig=c(fixef(model.res), get.ranef(summary(model.res)$varcor))
	xsum=apply(all.coeff.mat, 2, range, na.rm=T)
  xsum=data.frame(what=colnames(all.coeff.mat), orig=orig[colnames(all.coeff.mat)], t(xsum))
	rownames(xsum)=as.character(xsum$what)
  colnames(to.do)=c("ranef", "level")
  xx=apply(is.na(all.coeff.mat), 1, sum)
  if(sum(xx>0 & xx<ncol(all.coeff.mat))>0){
		warning(paste(c("for", sum(xx>0 & xx<ncol(all.coeff.mat)), "subset(s) the full model could not be fitted because of fixed effects factor levels dropped from the data"), collapse=" "))
	}
	colnames(all.warnings)[1:2]=c("lme4.warnings", "opt.warnings")
  all.coeff.mat=data.frame(to.do, N=all.n, all.coeff.mat, all.warnings)
  names(xsum)[3:4]=c("min", "max")
  return(list(detailed=all.coeff.mat, summary=xsum))
}
###########################################################################################################

###########################################################################################################
# m.stab.plot - function
###########################################################################################################
m.stab.plot<-function(est, lower=NULL, upper=NULL, xnames=NULL, col="black", center.at.null=F, reset.par=T, pch=18){
	#version from Aug 12 2018
	#function to plot model stability or CIs (not really supposed to be nice, but to give a quick overview/rapid diagnostic)
	#input:
		#either a three columns data frame or matrix (with rownames) handed over to argument est (first column needs to comprise the original estimate)
		#or several vectors handed over as follows:
			#est: numeric; estimated coefficients of the model
			#lower: numeric; lower limits of the estimates (either from model stability of from bootstrap)
			#upper: numeric; upper limits of the estimates (either from model stability of from bootstrap)
			#xnames: character; names to be depicted besides the error bars
		#col: character; name of the color with which data should be depicted
	if(ncol(est)==3){
		lower=est[, 2]
		upper=est[, 3]
		xnames=rownames(est)
		est=est[, 1]
		x.at=est
	}
  old.par = par(no.readonly = TRUE)
	par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1, 0.4, 0), tcl=-0.2)
	plot(x=est, y=1:length(est), pch=pch, xlab="estimate", ylab="", yaxt="n", xlim=range(c(0, lower, upper)), type="n", ylim=c(1, length(est)+1))
	abline(v=0, lty=3)
	if(center.at.null){
		x.at=rep(0, length(est))
		text(labels=xnames, x=x.at, y=(1:length(est))+0.3, cex=0.8)
	}else{
		xx=range(est)#c(lower, upper))
		for(i in 1:length(est)){
			text(labels=xnames[i], x=x.at[i], y=i+0.3, cex=0.8, adj=c((est[i]-min(xx))/diff(xx), 0))#pos=c(2, 4)[1+as.numeric(est<0)])
		}
	}
	points(x=est, y=1:length(est), pch=pch, col=col)
	segments(x0=lower, x1=upper, y0=1:length(est), y1=1:length(est), col=col)
	if(reset.par){par(old.par)}
}
###########################################################################################################

###########################################################################################################
# boot.glmm.pred - function
###########################################################################################################
boot.glmm.pred<-function(model.res, excl.warnings=F, nboots=1000, para=F, resol=1000, level=0.95, use=NULL, circ.var.name=NULL, circ.var=NULL, use.u=F, 
	n.cores=c("all-1", "all"), save.path=NULL, load.lib=T, lib.loc=.libPaths(), set.all.effects.2.zero=F){
	if(load.lib){library(lme4, lib.loc=lib.loc)}
	n.cores=n.cores[1]
	keepWarnings<-function(expr){
		localWarnings <- list()
		value <- withCallingHandlers(expr,
			warning = function(w) {
				localWarnings[[length(localWarnings)+1]] <<- w
				invokeRestart("muffleWarning")
			}
		)
		list(value=value, warnings=localWarnings)
	}
	##define function extracting all estimated coefficients (fixed and random effects) and also the model summary wrt random effects:
	extract.all<-function(mres){
		##extract random effects model summary:
		vc.mat=as.data.frame(summary(mres)$varcor)##... and prepare/extract variance covariance matrix from the model handed over
		xx=lapply(summary(mres)$varcor, function(x){attr(x, "stddev")})##append residual variance
		##create vector with names of the terms in the model
		xnames=c(names(fixef(mres)), paste(rep(names(xx), unlist(lapply(xx, length))), unlist(lapply(xx, names)), sep="@"))
		if(class(mres)[1]=="lmerMod"){xnames=c(xnames, "Residual")}##append "Residual" to xnames in case of Gaussian model
		if(class(mres)[1]=="lmerMod"){res.sd=vc.mat[vc.mat$grp=="Residual", "sdcor"]}##extract residual sd i case of Gaussian model
		#if(vc.mat$grp[nrow(vc.mat)]=="Residual"){vc.mat=vc.mat[-nrow(vc.mat), ]}##and drop residuals-row from vc.mat in case of Gaussisn model
		##deal with names in vc.mat which aren't exactly the name of the resp. random effect
		r.icpt.names=names(ranef(mres))##extract names of random intercepts...
		not.r.icpt.names=setdiff(vc.mat$grp, r.icpt.names)##... and names of the random effects having random slopes
		not.r.icpt.names=unlist(lapply(strsplit(not.r.icpt.names, split=".", fixed=T), function(x){paste(x[1:(length(x)-1)], collapse=".")}))
		vc.mat$grp[!vc.mat$grp%in%r.icpt.names]=not.r.icpt.names
		if(vc.mat$grp[nrow(vc.mat)]=="Residual"){vc.mat$var1[nrow(vc.mat)]=""}
		xnames=paste(vc.mat$grp, vc.mat$var1, sep="@")
		re.summary=unlist(vc.mat$sdcor)
		names(re.summary)=xnames
		ranef(mres)
		##extract random effects model summary: done
		##extract random effects details:
		re.detail=ranef(mres)
		xnames=paste(
			rep(x=names(re.detail), times=unlist(lapply(re.detail, function(x){nrow(x)*ncol(x)}))),
				unlist(lapply(re.detail, function(x){rep(colnames(x), each=nrow(x))})), 
				unlist(lapply(re.detail, function(x){rep(rownames(x), times=ncol(x))})),
				sep="@")
		re.detail=unlist(lapply(re.detail, function(x){
			return(unlist(c(x)))
		}))
		#browser()
		#deal with negative binomial model to extract theta:
		xx=as.character(summary(mres)$call)
		if(any(grepl(x=xx, pattern="negative.binomial"))){
			xx=xx[grepl(x=xx, pattern="negative.binomial")]
			xx=gsub(x=xx, pattern="negative.binomial(theta = ", replacement="", fixed=T)
			xx=gsub(x=xx, pattern=")", replacement="", fixed=T)
			re.detail=c(re.detail, as.numeric(xx))
			xnames=c(xnames, "theta")
		}
		names(re.detail)=xnames
		ns=c(length(fixef(mres)), length(re.summary), length(re.detail))
		names(ns)=c("n.fixef", "n.re.summary", "n.re.detail")
		return(c(fixef(mres), re.summary, re.detail, ns))
	}	
	if(excl.warnings){
		boot.fun<-function(x, model.res., keepWarnings., use.u., save.path.){
			xdone=F
			while(!xdone){
				i.res=keepWarnings(bootMer(x=model.res., FUN=extract.all, nsim=1, use.u=use.u.)$t)
				if(length(unlist(i.res$warnings)$message)==0){
					xdone=T
				}
			}
			est.effects=i.res$value
			i.warnings=NULL
			if(length(save.path.)>0){save(file=paste(c(save.path., "/b_", x, ".RData"), collapse=""), list=c("est.effects", "i.warnings"))}
			return(i.res$value)
		}
	}else{
		boot.fun<-function(y, model.res., keepWarnings., use.u., save.path.){
			#keepWarnings.(bootMer(x=model.res., FUN=fixef, nsim=1)$t)
			i.res=keepWarnings(bootMer(x=model.res., FUN=extract.all, nsim=1, use.u=use.u.))
			if(length(save.path.)>0){
				est.effects=i.res$value$t
				i.warnings=i.res$warnings
				save(file=paste(c(save.path., "/b_", y, ".RData"), collapse=""), list=c("est.effects", "i.warnings"))
			}
			return(list(ests=i.res$value$t, warns=unlist(i.res$warnings)))
		}
	}
	if(para){
		on.exit(expr = parLapply(cl=cl, X=1:length(cl), fun=function(x){rm(list=ls())}), add = FALSE)
		on.exit(expr = stopCluster(cl), add = T)
		library(parallel)
		cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores!="all"){
			if(n.cores=="all-1"){n.cores=length(cl)-1}
			if(n.cores<length(cl)){
				cl=cl[1:n.cores]
			}
		}
		parLapply(cl=cl, 1:length(cl), fun=function(x, .lib.loc=lib.loc){
		  library(lme4, lib.loc=.lib.loc)
		  return(invisible(""))
		})
		all.res=parLapply(cl=cl, X=1:nboots, fun=boot.fun, model.res.=model.res, keepWarnings.=keepWarnings, use.u.=use.u, save.path.=save.path)
	}else{
    all.res=lapply(X=1:nboots, FUN=boot.fun, model.res.=model.res, use.u.=use.u, save.path.=save.path)#, keepWarnings=excl.warnings)
	}
	if(!excl.warnings){
		all.warns=lapply(all.res, function(x){
			xxx=unlist(x$warns)
			if(length(xxx)==0){xxx=""}
			return(xxx)
		})
		all.res=lapply(all.res, function(x){x$ests})
	}else{
		all.warns=NULL
	}
	if(length(use)>0){
		#extract fixed effects terms from the model:
		xcall=as.character(model.res@call)[2]
		model.terms=attr(terms(as.formula(xcall)), "term.labels")
		REs=names(ranef(model.res))
		#for(i in 1:length(REs)){
			model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
		#}
		model=paste(model.terms, collapse="+")
		#build model wrt to the fixed effects:
		model.terms=unique(unlist(strsplit(x=model.terms, split=":", fixed=T)))
		model.terms=model.terms[!grepl(x=model.terms, pattern="I(", fixed=T)]
		model.terms=model.terms[!grepl(x=model.terms, pattern="^2)", fixed=T)]
		#exclude interactions and squared terms from model.terms:
		#create new data to be used to determine fitted values:
		ii.data=model.res@frame
		
		if(length(circ.var.name)==1){
			set.circ.var.to.zero=sum(circ.var.name%in%use)==0
		}else{
			set.circ.var.to.zero=F
		}
		
		if(length(use)==0){use=model.terms}
		new.data=vector("list", length(model.terms))
		usel=model.terms%in%use
		#if(length(use)>0)
		for(i in 1:length(model.terms)){
			if(is.factor(ii.data[, model.terms[i]])){
				new.data[[i]]=levels(ii.data[, model.terms[i]])
			}else if(!is.factor(ii.data[, model.terms[i]]) & usel[i] & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
				new.data[[i]]=seq(from=min(ii.data[, model.terms[i]]), to=max(ii.data[, model.terms[i]]), length.out=resol)
			}else  if(!is.factor(ii.data[, model.terms[i]]) & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
				new.data[[i]]=mean(ii.data[, model.terms[i]])
			}
		}
		names(new.data)=model.terms
		if(length(circ.var.name)==1){
			new.data=new.data[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
			if(sum(grepl(pattern=circ.var.name, x=use))>0){
				new.data=c(new.data, list(seq(min(circ.var, na.rm=T), max(circ.var, na.rm=T), length.out=resol)))
				names(new.data)[length(new.data)]=circ.var.name
			}else{
				new.data=c(new.data, list(0))
			}
			model.terms=model.terms[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
		}
		xnames=names(new.data)
		#browser()
		new.data=data.frame(expand.grid(new.data))
		names(new.data)=xnames
		#names(new.data)[1:length(model.terms)]=model.terms
		#browser()
		if(length(circ.var.name)==1){
			names(new.data)[ncol(new.data)]=circ.var.name
		}
		#create predictors matrix:
# 		if(length(circ.var.name)>0 & length(intersect(circ.var.name, names(new.data)))>0){
# 			new.data=cbind(new.data, sin(new.data[, circ.var.name]), cos(new.data[, circ.var.name]))
# 			names(new.data)[(ncol(new.data)-1):ncol(new.data)]=paste(c("sin(", "cos("), circ.var.name, ")", sep="")
# 			new.data=new.data[, !names(new.data)==circ.var.name]
# 		}
# 		browser()
		r=runif(nrow(new.data))
		if(set.all.effects.2.zero){
			for(iterm in setdiff(colnames(new.data), c("(Intercept)", use))){
				new.data[, iterm]=0
			}
		}
		m.mat=model.matrix(object=as.formula(paste(c("r", model), collapse="~")), data=new.data)
		if(set.circ.var.to.zero){
			m.mat[,paste(c("sin(", circ.var.name, ")"), collapse="")]=0
			m.mat[,paste(c("cos(", circ.var.name, ")"), collapse="")]=0
		}
		#m.mat=t(m.mat)
		#get the CIs for the fitted values:
		ci=lapply(all.res, function(x){
			#return(apply(m.mat[names(fixef(model.res)), , drop=F]*as.vector(x[, names(fixef(model.res))]), 2, sum))
			return(m.mat[, names(fixef(model.res)), drop=F]%*%as.vector(x[, names(fixef(model.res))]))
		})
		ci=matrix(unlist(ci), ncol=nboots, byrow=F)
		ci=t(apply(ci, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
		colnames(ci)=c("lower.cl", "upper.cl")
		fv=m.mat[, names(fixef(model.res))]%*%fixef(model.res)
		if(class(model.res)[[1]]!="lmerMod"){
			if(model.res@resp$family$family=="binomial"){
				ci=exp(ci)/(1+exp(ci))
				fv=exp(fv)/(1+exp(fv))
			}else if(model.res@resp$family$family=="poisson" | (model.res@resp$family$family=="Gamma" & model.res@resp$family$link=="log") | substr(x=model.res@resp$family$family, start=1, stop=17)=="Negative Binomial"){
				ci=exp(ci)
				fv=exp(fv)
			}
		}
		result=data.frame(new.data, fitted=fv, ci)
	}else{
		result=NULL
	}
	all.boots=matrix(unlist(all.res), nrow=nboots, byrow=T)
	colnames(all.boots)=colnames(all.res[[1]])
	ci.est=apply(all.boots, 2, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T)
	if(length(fixef(model.res))>1){
		ci.est=data.frame(orig=fixef(model.res), t(ci.est)[1:length(fixef(model.res)), ])
	}else{
		ci.est=data.frame(orig=fixef(model.res), t(t(ci.est)[1:length(fixef(model.res)), ]))
	}
	return(list(ci.predicted=result, ci.estimates=ci.est, all.warns=all.warns, all.boots=all.boots))
}
###########################################################################################################

###########################################################################################################
# overdisp.test - function
###########################################################################################################
overdisp.test<-function(x){
	##last updated Nov 26 2019
	##additions/changes: bugfix for zeroinfl; deals with glmmTMB (beta and nbinom family)
	##function needed for Gamma family:
	gamma.par.transf<-function(mean=NULL, var=NULL, shape=NULL, scale=NULL){
		##parameterizations are
			##shape and scale
			##shape and mean
			##var and mean
		if(is.null(mean) & is.null(var) & !is.null(shape) & !is.null(scale)){
			mean=shape*scale
			var=shape*scale^2
		}else if(is.null(shape) & is.null(scale) & !is.null(mean) & !is.null(var)){
			scale=var/mean
			shape=mean/scale
		}else if(is.null(scale) & is.null(var) & !is.null(shape) & !is.null(mean)){
			scale=mean/shape
			var=shape*scale^2
		}else{
			Stop("unsupported parameter commbination")
		}
		return(list(shape=as.vector(shape), scale=as.vector(scale),  mean=as.vector(mean), var=as.vector(var)))
	}
	
  if(class(x)[[1]]!="glmmTMB"){
		pr=residuals(x, type ="pearson")
		sum.dp=sum(pr^2)
	}
  n.model.terms=NULL
  if(class(x)[[1]]=="glmerMod"){
    n.model.terms=length(fixef(x))+nrow(as.data.frame(summary(x)$varcor))
    if(grepl(x=as.character(x@call)[4], pattern="negative.binomial") | grepl(x=as.character(x@call)[4], pattern="Gamma")){n.model.terms=n.model.terms+1}
	}else if(class(x)[[1]]=="glm" | class(x)[[1]]=="negbin"){
    n.model.terms=length(x$coefficients)
    n.model.terms=n.model.terms+length(summary(x)$theta)
  }else if(class(x)[[1]]=="zeroinfl"){
		n.model.terms=sum(unlist(lapply(x$coefficients, length)))
		n.model.terms=n.model.terms+length(summary(x)$theta)
  }else if(class(x)[[1]]=="betareg"){
		n.model.terms=sum(unlist(lapply(x$coefficients, length)))
  }else if(class(x)[[1]]=="glmmTMB"){
		print("caution: support for models of class glmmTMB is experimental and embryonic")
		print("and don't worry if it takes some time")
		#browser()
		#n.model.terms=sum(unlist(lapply(fixef(x), length)))+nrow(extract.ranef.from.glmmTB(x))
		n.model.terms=sum(unlist(lapply(fixef(x), length)))##fixed effects
		n.model.terms=n.model.terms+sum(unlist(lapply(summary(x)$varcor, function(y){
			if(!is.null(y)){
				lapply(y, function(yy){
					xx=dim(attr(yy, "correlation"))
					return((xx[1]^2-xx[1])/2+length(attr(yy, "stddev")))
				})
			}else{
				0
			}
		})))
		if((grepl(x=x$modelInfo$family$family, pattern="nbinom") | x$modelInfo$family$family=="beta") &
			length(attr(terms(x$call$dispformula), "term.labels"))==0){
			#n.model.terms=n.model.terms+1
			##but maybe it needs to be done for negbin??? need to check this
		}
		pr=residuals(x, type ="response")
		if(x$modelInfo$family$family=="poisson"){
			fitted.var = fitted(x)
		}else{
			#model.terms=attr(terms(x$call$dispformula), "term.labels")#attr(terms(x$call$dispformula), "intercept")
			sigma.f=exp(as.vector(model.matrix(object=x$call$dispformula, data=x$frame)[, names(fixef(x)$disp), drop=F]%*%fixef(x)$disp))
			if(x$modelInfo$family$family=="nbinom2"){
				fitted.var = fitted(x)*(1+fitted(x)/sigma.f)#fitted(x) + fitted(x)^2/sigma(x)
			}else if(x$modelInfo$family$family=="nbinom1"){
				fitted.var = fitted(x)*(1+sigma.f)#fitted(x) + fitted(x)/sigma(x)
			}else if(x$modelInfo$family$family=="beta"){
				xfitted=fitted(x)
				shape1 = xfitted * sigma.f
				shape2 = (1 - xfitted) * sigma.f
				#xfitted(1-xfitted)/(1+sigma(x))#according to the glmmTMB help page for sigma
				fitted.var = shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))
			}else if(x$modelInfo$family$family=="Gamma"){
				#x.disp=exp(as.vector(model.matrix(object=x$call$dispformula, data=x$frame)%*%log(sigma(x))))
				xx=gamma.par.transf(mean=fitted(x), var=NULL, shape=1/sigma.f, scale=NULL)
				fitted.var = xx$var
				#xdf=length(residuals(x))-(length(coef(x))+1)
				#n.model.terms=sum(unlist(lapply(fixef(x), length)))+sum(unlist(lapply(summary(x)$varcor, function(y){length(as.data.frame(y))})))
				xdf=length(residuals(x))-n.model.terms
			}
		}
		#n.model.terms=sum(unlist(lapply(fixef(x), length)))+sum(unlist(lapply(summary(x)$varcor, function(y){length(as.data.frame(y))})))
		#xdf=length(pr)-n.model.terms
		sum.dp=sum((pr/sqrt(fitted.var))^2)
		#if(as.numeric(substr(sessionInfo()[["otherPkgs"]][["glmmTMB"]]$Version, start=1, stop=1))>0){
			#sum.dp=sum(residuals(x, type ="pearson")^2)
		#}
  }
  if(!is.null(n.model.terms)){
		xdf=length(pr)-n.model.terms
		return(data.frame(chisq=sum.dp, df=xdf, P=1-pchisq(sum.dp, xdf), dispersion.parameter=sum.dp/xdf))
	}else{
		print("model isn't of any of the currently supported classes (glm, negbin, zeroinfl, glmerMod, glmmTMB)")
	}
}
###########################################################################################################

###########################################################################################################
# factor.int.plot - function
###########################################################################################################
factor.int.plot<-function(
  plot.data, factors, coefs=NULL, response, link=c("identity", "logit", "log"),
  conf.int=NULL, yaxt.labels=NULL, yaxt.at=NULL, ylab, ylim=NULL, size.fac=1,
	factor.seqs=NULL, factor.labels=NULL, to.show=c("bubbles", "data", "bars", "boxes"), weights=NULL, which.q=c(5, 3),
	pch=1, pt.col=par("fg"), rect.col=NULL, border.col=NULL, est.ci.col=par("fg"), est.ci.lty=1, quiet=T, average.response=F, log.y=F, 
	bg="white", median.col="black", median.lwd=2, fitted.lwd=2, quant.col=NULL,
	cex.lab=1, cex.lab.x=1, cex.axis=1, reset.par=T, percentile.col=NULL, percentile.lwd=NULL, my.par=NULL, cex.axis.lab=1, xlim=NULL,
	add.range=F, range.pch=4){
	#last updated: 2015 Nov 26
	#arguments:
		#plot.data: data frame comprising all variables involved
		#factors: character, vector with the names of the factors to be depicted (the hierarchical order in which they are displayed is the reverse of factors, i.e.,
			#the first factor will be the upper most at the x-axis and the last the lowermost)
		#coefs: named vector with the model coefficients to be shown (default NULL in which case no model coefficients are shown)
		#response: vector with the values of the response
		#link: character, one of "identity" (the default), "logit", or "log"; potentially leads to the model as inferred from coefs to be transformed
		#conf.int: numeric data frame or matrix with columns to be named "lwr" and "upr", respectively, and also columns for the resp. factors 
			#(precisely, the dummmy coded levels, named as in model output) indicating the CIs; default is NULL in which no CIs are depicted
		#yaxt.labels: numeric, labels to be displayed along the y-axis (default, NULL in which labels are determined from the response)
		#yaxt.at: numeric, values along the y-axis where ticks are to be placed (default, NULL in the values labels are determined from the response)NULL
		#ylab: as  usual
		#ylim: numeric, as usual (default: NULL, in which case ylim is the range of tthe response)
		#size.fac: numeric, determines the size of the points depicting the data in case to.show comprises "bubbles"
			#(default 1; if points are too large make it smaller, otherwise larger)
		#factor.seqs: character list, determines the order with which the levels of the factors are depicted along the x-axis
			#(default is NULL in which it corresponds to their levels)
		#factor.labels: character list, determines the labels which are displayed for the levels of the factors along the x-axis
			#(default is NULL in which it corresponds to their levels)
		#to.show: character vector, can be one of "bubbles", "bars", or "boxes" or a combination of he former with one of the latter two;
			#determines whether the data are shown (to.show comprises "bubbles") and whether a summary in the form of bars depicting medians ((to.show comprises "bars") or
			#boxes etc. (to.show comprises "boxes") are depicted (see also argument 'which.q')
		#which.q: numeric, one of 5 or 3; if 5 (the default) medians, quartiles and quantiles (2.5 and 97.5%) are shown, if 3 medians and quartiles are shown
		#pch, numeric, one of 1:25 (default is 1); determines the 'point character' (has the same interpretation as usual)
		#pt.col: character, color with which points/bubbles are depicted (default: par("fg"))
		#rect.col: character, color with which bars or boxes are depicted (default: NA  in which case they are shown with a border in the default foreground color and no filling)
		#est.ci.col: character, color with which CIs are depicted (default: par("fg"))
		#quiet: logical, determines whether the function returns the number of observations per bubble (only relevant if to.show comprieses "bubbles"
	old.par = par(no.readonly = TRUE)
	#if response is cbind:
	if(length(dim(response))==2){
		response=response[, 1]/apply(response, 1, sum)
	}
	which.q=which.q[1]
	#get rid of NAs:
	plot.data=droplevels(as.data.frame(na.omit(data.frame(plot.data[, factors, drop=F], response))))
	link=link[1]#get default link (in case not specified)
	#create lists factor.seqs and factor.seqs in case they are empty:
	if(length(factor.seqs)==0){
		factor.seqs=lapply(factors, function(x){levels(plot.data[, x])})
		names(factor.seqs)=factors
	}else{
		lapply(1:length(factors), function(x){
			if(length(levels(plot.data[, factors[x]]))!=length(factor.seqs[[x]])){
				stop("error: order of factors doesn't match in factors and factor.seqs")
			}else{
				if(sum(sort(levels(plot.data[, factors[x]]))!=sort(factor.seqs[[x]]))>0){
					stop("error: factor levels don't match in factors and factor.seqs")
				}
			}
		})
	}
	names(factor.seqs)=factors
	if(length(factor.labels)==0){factor.labels=factor.seqs}
	#set ylim in case it is empty:
	if(length(ylim)==0){
		ylim=range(plot.data[, "response"])
	}
	if(length(yaxt.labels)==0){
		yaxt.labels=pretty(ylim)
		yaxt.at=yaxt.labels
	}
	#create data matrix (needed for something taking place later (guess it has to do with the order of the factor levels)):
	pred.data=data.frame(expand.grid(lapply(factors, function(x){levels(plot.data[, x])})))
	names(pred.data)=factors
	pred.data$all=apply(pred.data[, factors, drop=F], 1, paste, collapse="@@@")
	pred.data=pred.data[match(apply(data.frame(expand.grid(factor.seqs)), 1, paste, collapse="@@@"), pred.data$all), ]
	##create corresponding matrix for x-axis labels
	labels.mat=data.frame(expand.grid(factor.labels))
	###subset such that only level combinations existing in the data are considered:
	to.keep=unique(apply(plot.data[, factors, drop=F], 1, paste, collapse="@@@"))
	labels.mat=subset(labels.mat, pred.data$all%in%to.keep)
	pred.data=subset(pred.data, all%in%to.keep)
	if(length(coefs)>0){
		#extract coefficients needed:
		#paste factors and their levels together:
		to.search=lapply(factors, function(x){levels(plot.data[, x])[-1]})
		to.search=paste(rep(factors, times=unlist(lapply(to.search, length))), unlist(to.search), sep="")
		tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
			length(cf)==sum(unlist(lapply(to.search, function(v){
				cf==v
			})))
		}))
		if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
		coefs=coefs[tk]
		#determine fitted values:
		#browser()
		m.mat=model.matrix(as.formula(paste(c("~", paste(factors, collapse="*")), collapse="")), data=pred.data)
		#get fitted values:
		colnames(m.mat)=unlist(lapply(strsplit(colnames(m.mat), split=":", fixed=T), function(x){paste(sort(x), collapse=":")}))
		names(coefs)=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(x){paste(sort(x), collapse=":")}))
		xfitted=apply(t(t(m.mat[, names(coefs)])*coefs), 1, sum)
		if(link=="log"){
			xfitted=exp(xfitted)
		}else if(link=="logit"){
			xfitted=exp(xfitted)/(1+exp(xfitted))
		}
		###########################################################################################
		##probably needs to be adjusted to order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		###########################################################################################
	}
	if(log.y & min(plot.data[, "response"])>0){
		plot.data[, "response"]=log(plot.data[, "response"])
	}else if(log.y & min(plot.data[, "response"])==0){
		plot.data[, "response"]=log(plot.data[, "response"]+1)
	}
	if(sum("bubbles"%in%to.show)==1 & length(ylim)==0){
		ylim=range(plot.data[, "response"])
	}
	if(sum(c("bars", "boxes")%in%to.show)==1){
		plot.data$all=apply(plot.data[, factors, drop=F], 1, paste, collapse="@@@")
		q.mat=tapply(plot.data[, "response"], plot.data$all, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
		xnames=names(q.mat)
		q.mat=as.data.frame(matrix(unlist(q.mat), nrow=length(q.mat), byrow=T))
		colnames(q.mat)=c("x0", "x025", "x25", "x5", "x75", "x975", "x1")
		rownames(q.mat)=xnames
		q.mat$all=rownames(q.mat)
		q.mat=q.mat[match(pred.data$all, q.mat$all), ]
	}
	
	hbw=0.25
	u.all.fac=pred.data[, "all"]
	#begin with plot:
  if(is.null(my.par)){
	  par(mar=c(length(factors)+0.2*length(factors)+0.6, 3.4, 0.2, 0.2), mgp=c(1.8, 0.35, 0), tcl=-0.2, bg=bg, mex=cex.lab*0.7, las=1)
  }else{
    par(my.par)
  }
	if(sum("bubbles"%in%to.show)==1 & !average.response){
		obs.mean=as.data.frame(table(plot.data[, c("response", factors)]))
		obs.mean$response=as.numeric(as.character(obs.mean$response))
		all.fac=apply(obs.mean[, factors, drop=F], 1, paste, collapse="@@@")
		obs.mean$all.fac=all.fac
		obs.mean=subset(obs.mean, Freq>0)
		if(is.null(xlim)){xlim=c(0.5, length(u.all.fac)+0.5)}
		plot(x=match(obs.mean$all.fac, u.all.fac), y=obs.mean$response, cex=size.fac*sqrt(obs.mean$Freq), yaxt="n", 
			xlim=xlim, xaxs="i", xaxt="n", xlab="", ylab=ylab, ylim=ylim, pch=pch, col=pt.col, type="n", cex.lab=cex.axis.lab)
	}else if(sum("bubbles"%in%to.show)==1 & average.response){
		obs.mean=aggregate(plot.data$response, plot.data[, factors], mean)
		names(obs.mean)[ncol(obs.mean)]="response"
		obs.mean$Freq=aggregate(plot.data$response, plot.data[, factors], length)$x
		all.fac=apply(obs.mean[, factors, drop=F], 1, paste, collapse="@@@")
		obs.mean$all.fac=all.fac
		if(is.null(xlim)){xlim=c(0.5, length(u.all.fac)+0.5)}
		plot(x=match(obs.mean$all.fac, u.all.fac), y=obs.mean$response, cex=size.fac*sqrt(obs.mean$Freq), yaxt="n", 
			xlim=xlim, xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, pch=pch, col=pt.col, type="n", cex.lab=cex.axis.lab)
	}else if(sum("data"%in%to.show)==1){
		plot.data$all.fac=apply(plot.data[, factors, drop=F], 1, paste, collapse="@@@")
		if(is.null(weights)){weights=rep(1, nrow(plot.data))}
		if(is.null(xlim)){xlim=c(0.5, length(u.all.fac)+0.5)}
		plot(x=match(plot.data$all.fac, u.all.fac), y=plot.data$response, cex=size.fac*sqrt(weights), yaxt="n", 
			xlim=xlim, xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, pch=pch, col=pt.col, type="n", cex.lab=cex.axis.lab)
		obs.mean=data.frame(plot.data, Freq=weights)
	}else if(sum("bars"%in%to.show)==1){
		if(is.null(xlim)){xlim=xlim}
		plot(x=1, y=1, xlim=c(0.5, nrow(pred.data)+0.5), xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, type="n", yaxs="i", 
			yaxt="n", cex.lab=cex.axis.lab)
	}else{
		if(is.null(xlim)){xlim=c(0.5, nrow(pred.data)+0.5)}
		plot(x=1, y=1, xlim=xlim, xaxs="i", xaxt="n", las=1, xlab="", ylab=ylab, ylim=ylim, type="n", yaxt="n", 
      cex.lab=cex.axis.lab)
	}
	axis(side=2, las=1, at=yaxt.at, labels=yaxt.labels, cex.axis=cex.axis, cex.lab=cex.lab, las=my.par[["las"]])
	#add points in case specified;
	if(sum("bubbles"%in%to.show)==1){
		points(x=match(obs.mean$all.fac, u.all.fac), y=obs.mean$response, cex=size.fac*sqrt(obs.mean$Freq), 
			pch=pch, col=pt.col)
	}else if(sum("data"%in%to.show)==1){
		points(x=match(plot.data$all.fac, u.all.fac), y=plot.data$response, cex=size.fac*sqrt(weights), pch=pch, col=pt.col)
	}
	#add bars or boxes etc. in case specified;
	if(is.null(rect.col)){rect.col=par("bg")}
	if(is.null(border.col)){border.col=par("fg")}
	if(length(quant.col)==0){quant.col=median.col}
	if(sum("bars"%in%to.show)==1){
		if(average.response){
			obs.mean.2=aggregate(plot.data$response, plot.data[, factors], mean)
			obs.mean.2$all=apply(obs.mean.2[, factors, drop=F], 1, paste, collapse="@@@")
			obs.mean.2=obs.mean.2[match(obs.mean.2$all, pred.data$all), "x"]
		}else{
			obs.mean.2=q.mat[match(q.mat$all, pred.data$all), "x5"]
		}
		rect(xleft=(1:nrow(pred.data))-hbw, xright=(1:nrow(pred.data))+hbw, ybottom=min(ylim), ytop=obs.mean.2, border=border.col, col=rect.col)
	}else if(sum("boxes"%in%to.show)==1){
		rect(xleft=(1:nrow(pred.data))-hbw, xright=(1:nrow(pred.data))+hbw, ybottom=q.mat[match(q.mat$all, pred.data$all), "x25"], ytop=q.mat[match(q.mat$all, pred.data$all), "x75"], col=rect.col, border=border.col)
		segments(x0=(1:nrow(pred.data))-hbw, x1=(1:nrow(pred.data))+hbw, y0=q.mat[match(q.mat$all, pred.data$all), "x5"], y1=q.mat[match(q.mat$all, pred.data$all), "x5"], lwd=median.lwd, col=median.col, lend=1)
		if(which.q==5){
			segments(x0=(1:nrow(pred.data)), x1=(1:nrow(pred.data)), y0=q.mat[match(q.mat$all, pred.data$all), "x75"], y1=q.mat[match(q.mat$all, pred.data$all), "x975"], col=quant.col)
			segments(x0=(1:nrow(pred.data)), x1=(1:nrow(pred.data)), y0=q.mat[match(q.mat$all, pred.data$all), "x25"], y1=q.mat[match(q.mat$all, pred.data$all), "x025"], col=quant.col)
		}
		if(add.range){
			points((1:nrow(pred.data)), y=q.mat[match(q.mat$all, pred.data$all), "x0"], pch=range.pch)
			points((1:nrow(pred.data)), y=q.mat[match(q.mat$all, pred.data$all), "x1"], pch=range.pch)
		}
	}
	#add text at x-axis:
	#u.all.fac=u.all.fac[match(u.all.fac, pred.data$all)]
	#u.all.fac.mat=matrix(unlist(strsplit(u.all.fac, split="@@@", fixed=T)), ncol=length(factors), byrow=T)
	#u.all.fac.mat=u.all.fac.mat[, ncol(u.all.fac.mat):1, drop=F]
	mtext(text=labels.mat[, 1], side=1, line=1-1+0.3*1, 
		at=1:nrow(labels.mat), cex=cex.lab.x)
	if(ncol(labels.mat)>1){
		for(i in 2:ncol(labels.mat)){
			to.add.1=c(T, labels.mat[-nrow(labels.mat), i]!=labels.mat[-1, i])
			to.add.2=cumsum(to.add.1)
			to.add.2=tapply(1:nrow(labels.mat), to.add.2, mean)
			mtext(text=labels.mat[to.add.1, i], side=1, line=i-1+0.3*i, 
				at=to.add.2,
				#at=sort(tapply(1:nrow(u.all.fac.mat), apply(u.all.fac.mat[, 1:(length(factors)+1-i), drop=F], 1, paste, collapse="@@@"), mean)),
				cex=cex.lab.x)
		}
	}
	#add model coefficients:
	if(length(coefs)>0){
		if(log.y & min(plot.data[, "response"])>0){
			xfitted=log(xfitted)
		}else if(log.y & min(plot.data[, "response"])==0){
			xfitted=log(xfitted+1)
		}
		hll=0.15
		segments(x0=(1:nrow(pred.data))-hll, x1=(1:nrow(pred.data))+hll, y0=xfitted, y1=xfitted, lwd=fitted.lwd, col=est.ci.col, lty=est.ci.lty, lend=1)
	}
	#add CIs:
	if(length(conf.int)>0){
		conf.int$all=apply(conf.int[, factors, drop=F], 1, paste, collapse="@@@")
		conf.int=conf.int[match(u.all.fac, conf.int$all), ]
		arrows(x0=(1:nrow(pred.data)), x1=(1:nrow(pred.data)), y0=conf.int[, "lwr"], y1=conf.int[, "upr"], lwd=fitted.lwd, length=0.05, code=3, angle=90, col=est.ci.col, lty=est.ci.lty)
	}
  if(reset.par){par(old.par)}
	#if data are to be shown and function shouldn't be quiet
	if((sum(to.show%in%"bubbles")==1 | sum(to.show%in%"data")==1) & !quiet){
		obs.mean=obs.mean[, names(obs.mean)!="all.fac"]
		return(obs.mean)#return object telling sample size per bubble
	}
}

###########################################################################################################

###########################################################################################################
# glmmTMB.stability - function
###########################################################################################################

glmmTMB.stab<-function(model.res, contr=NULL, ind.cases=F, para=F, data=NULL, use=NULL, n.cores=c("all-1", "all"), save.path=NULL, load.lib=T, lib.loc=.libPaths()){
	print("please carefully evaluate whether the result makes sense, and if not, please contact me")
	#function determining stability of GLMMs (run using lmer or glmer) by excluding levels of random effects, one at a time;
	#supports
		#weights, offset terms and random slopes;
	#does not support
		#correlations between random intercepts and random slopes
		#any terms calculated in the model formula (e.g., log, function I, etc.); but interactions do work
	#latest additions/modifications:
		#copes also with data sets where a level of a (fixed effects) factor is entirely dropped from the data
		#new way of catching warnings (are contained in the detailed output table)
		#includes sample size in the detailed output table
		#catches errors
		#use: an argument taking the names of the random effects for which model stability should be evaluated
			#(useful for models with a random effect having a unique case for each level of the data)
	#written by Roger Mundry
	#modified April 2016 (added dealing with glmer.nb; fixed a bug in the output of the random effects of tthe original model)
	#last modified Mar 2017 (added dealing with fixed effects factor levels dropped when dropping levels of random effects)
	#last modified June 14 2017 (fixed dealing with random effects comprising correlation parameters)
	#last modified July 06 2017 (fixed small bug happening when model revealed an error)
	#last modified Nov 26 2019 (copes with models uncluding no randm effects)
	if(load.lib){library(glmmTMB, lib.loc=lib.loc)}
	n.cores=n.cores[1]
	xx=as.character(model.res$call)
	names(xx)=names(model.res$call)
  model.fe.re=as.formula(xx["formula"])
  model.zi=as.formula(xx["ziformula"])
	model.disp=as.formula(xx["dispformula"])

	xcall=as.character(model.res$call)
	names(xcall)=names(model.res$call)
	xcall=list(cond.form=xcall["formula"], zi.form=xcall["ziformula"], disp.form=xcall["dispformula"], xfam=xcall["family"])
	if(grepl(xcall[["xfam"]], pattern="(", fixed=T)){
		xfam=xcall[["xfam"]]
		xfam=unlist(strsplit(xfam, split="(", fixed=T))
		xfam[2]=gsub(x=xfam[2], pattern=")", replacement="", fixed=T)
		xfam[2]=gsub(x=xfam[2], pattern="link = ", replacement="", fixed=T)
		xfam[2]=gsub(x=xfam[2], pattern="\"", replacement="", fixed=T)
		if(substr(xcall[["xfam"]], start=1, stop=5)!="Gamma"){
			xcall[["xfam"]]=get(xfam[1])(xfam[2])
		}else{
			if(xfam[2]=="log"){
				xcall[["xfam"]]=Gamma(link="log")
			}else if(xfam[2]=="inverse"){
				xcall[["xfam"]]=Gamma(link="inverse")
			}else if(xfam[2]=="identity"){
				xcall[["xfam"]]=Gamma(link="identity")
			}else{
				stop("Error: family not supported")
			}
		}
	}
	xfam=xcall[["xfam"]]
# 	
# 	xfam=xx["family"]
# 	if(grepl(xfam, pattern="(", fixed=T)){
# 		#xfam=xcall[["xfam"]]
# 		xfam=unlist(strsplit(xfam, split="(", fixed=T))
# 		xfam[2]=gsub(x=xfam[2], pattern=")", replacement="", fixed=T)
# 		xfam[2]=gsub(x=xfam[2], pattern="link = ", replacement="", fixed=T)
# 		xfam[2]=gsub(x=xfam[2], pattern="\"", replacement="", fixed=T)
# 		if(substr(xfam, start=1, stop=5)!="Gamma"){
# 			xfam=get(xfam[1])(xfam[2])
# 		}else{
# 			if(xfam[2]=="log"){
# 				xfam=Gamma(link="log")
# 			}else if(xfam[2]=="inverse"){
# 				xfam=Gamma(link="inverse")
# 			}else if(xfam[2]=="identity"){
# 				xfam=Gamma(link="identity")
# 			}else{
# 				stop("Error: family not supported")
# 			}
# 		}
# 		#xfam=get(xfam[1])(xfam[2])
# 	}
	##need to address how weights are recognized
	if(any(names(xx)=="weights")){
		data$XXXweights=data[, xx["weights"]]
	}else{
		data$XXXweights=1
	}
  ranefs=unique(unlist(lapply(ranef(model.res), names)))
	if(length(use)==0){use=ranefs}
	ranefs=intersect(ranefs, use)
  xlevels=lapply(ranefs, function(x){return(as.vector(unique(data[ ,x])))})
  ranefs=rep(ranefs, unlist(lapply(xlevels, length)))
  to.do=cbind(ranefs, unlist(xlevels))
  if(ind.cases){
    data=data.frame(data, ic=as.factor(1:nrow(data)))
    to.do=rbind(to.do, cbind("ic", levels(data$ic)))
  }
	extract.ranef.from.glmmTB<-function(m){
		#to.do=VarCorr(model.res)
		#varcor.names=names(to.do)
		to.do=summary(m)$varcor
		xx=lapply(to.do, function(x){
			lapply(x, attr, "stddev")
		})
		res=data.frame(
			what=rep(names(to.do), unlist(lapply(xx, function(x){sum(unlist(lapply(x, length)))}))),
			Groups=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){unlist(lapply(x, length))}))),
			var1=unlist(lapply(xx, function(x){unlist(lapply(x, names))})),
			var2=NA,
			sdcor=unlist(xx)
		)
		xx=lapply(to.do, function(x){
			lapply(x, attr, "correlation")
		})
		xx=xx[unlist(lapply(xx, length))>0]
		xx=data.frame(
			what=rep(names(xx), unlist(sum(unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))))),
			Groups=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))),
			var1=unlist(lapply(xx, function(x){lapply(x, function(x){rep(rownames(x), times=ncol(x))})})),
			var2=unlist(lapply(xx, function(x){lapply(x, function(x){rep(colnames(x), each=nrow(x))})})),
			sdcor=unlist(xx)
		)
		xx=xx[as.character(xx[, "var1"])<as.character(xx[, "var2"]), ]
		res=rbind(res, xx)
		res=res[order(as.character(res$what), as.character(res$Groups), as.character(res$var1)), ]
		return(res)
	}
	if(length(contr)==0){
		contr=glmmTMBControl()
	}
	n.ranef=nrow(extract.ranef.from.glmmTB(model.res))
  ifun<-function(x, model.res, to.do, data, contr, extract.ranef.from.glmmTB, n.ranef, xfam){
    sel.ii.data=data[data[,to.do[x, 1]]!=to.do[x, 2], ]
		sel.ii.res=try(glmmTMB(formula=model.fe.re, ziformula=model.zi, dispformula=model.disp, family=xfam, data=sel.ii.data, weights=XXXweights, control=contr), silent=T)
		#browser()
		if(length(save.path)>0){
			if(class(sel.ii.res)!="try-error"){
				est.fixed.effects=fixef(sel.ii.res)
				est.random.effects=summary(sel.ii.res)$varcor
				n=sel.ii.res$modelInfo$nobs
				what=to.do[x, ]
				converged=sel.ii.res$sdr$pdHess
				save(file=paste(c(paste(c(paste(c(save.path, "m"), collapse="/"), x), collapse="_"), ".RData"), collapse=""), 
					list=c("what", "est.fixed.effects", "est.random.effects", "n", "converged"))
			}else{
				ires="model revealed an error"
				save(file=paste(c(paste(c(paste(c(save.path, "m"), collapse="/"), x), collapse="_"), ".RData"), collapse=""), list=c("ires"))
			}
		}
    if(class(sel.ii.res)!="try-error"){
			fe=fixef(sel.ii.res)
			xnames=paste(rep(names(fe), times=unlist(lapply(fe, length))), unlist(lapply(fe, names)), sep="@")
			fe=unlist(fe)
			names(fe)=xnames
			RE=NULL
			if(any(!unlist(lapply(summary(sel.ii.res)$varcor, is.null)))){
				RE=extract.ranef.from.glmmTB(sel.ii.res)
				xx=apply(as.matrix(RE[, c("what", "Groups", "var1", "var2")]), 1, paste, collapse="@")
				RE=RE[, "sdcor"]
				names(RE)=xx
			}
			return(list(fere=c(fe, RE), N=sel.ii.res$modelInfo$nobs, converged=sel.ii.res$sdr$pdHess))
		}else{
      return(list(fere=rep(NA, sum(unlist(lapply(fixef(model.res), length)))+n.ranef), N=NA, converged="error"))
    }
  }
  if(para){
		#on.exit(expr = parLapply(cl=cl, X=1:length(cl), fun=function(x){rm(list=ls())}), add = FALSE)
		#on.exit(expr = stopCluster(cl), add = T)
    require(parallel)
    cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores!="all"){
			if(n.cores=="all-1"){n.cores=length(cl)-1}
			if(n.cores<length(cl)){
				cl=cl[1:n.cores]
			}
		}
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      library(glmmTMB, lib.loc=lib.loc)
      return(invisible(""))
    })
    all.coeffs=parLapply(cl=cl, X=1:nrow(to.do), fun=ifun, model.res=model.res, to.do=to.do, data=data, contr=contr, 
			extract.ranef.from.glmmTB=extract.ranef.from.glmmTB, n.ranef=n.ranef, xfam=xfam)
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      return(rm(list=ls()))
    })
    stopCluster(cl=cl)
  }else{
    all.coeffs=lapply(1:nrow(to.do), ifun, model.res=model.res, to.do=to.do, data=data, contr=contr, 
			extract.ranef.from.glmmTB=extract.ranef.from.glmmTB, n.ranef=n.ranef, xfam=xfam)
  }
	all.n=unlist(lapply(all.coeffs, function(x){x$N}))
	all.conv=unlist(lapply(all.coeffs, function(x){x$converged}))
	xnames=unique(unlist(lapply(all.coeffs, function(x){names(x$fere)})))
	all.coeff.mat=matrix(NA, ncol=length(xnames), nrow=length(all.coeffs))
	colnames(all.coeff.mat)=xnames
	for(i in 1:length(all.coeffs)){
		all.coeff.mat[i, names(all.coeffs[[i]]$fere)]=as.vector(all.coeffs[[i]]$fere)
	}
	#extract results for original model:
	orig=fixef(model.res)
	xnames=paste(rep(names(orig), times=unlist(lapply(orig, length))), unlist(lapply(orig, names)), sep="@")
	orig=unlist(orig)
	names(orig)=xnames
	if(any(!unlist(lapply(summary(model.res)$varcor, is.null)))){	
		RE=extract.ranef.from.glmmTB(model.res)
		xx=apply(as.matrix(RE[, c("what", "Groups", "var1", "var2")]), 1, paste, collapse="@")
		RE=RE[, "sdcor"]
		names(RE)=xx
		orig=c(orig, RE)
	}
	#browser()
  xsum=apply(all.coeff.mat, 2, range, na.rm=T)
  xsum=data.frame(what=colnames(all.coeff.mat), orig=orig[colnames(all.coeff.mat)], t(xsum))
	rownames(xsum)=as.character(xsum$what)
  colnames(to.do)=c("ranef", "level")
  xx=apply(is.na(all.coeff.mat), 1, sum)
  if(sum(xx>0 & xx<ncol(all.coeff.mat))>0){
		warning(paste(c("for", sum(xx>0 & xx<ncol(all.coeff.mat)), "subset(s) the full model could not be fitted because of fixed effects factor levels dropped from the data"), collapse=" "))
	}
  all.coeff.mat=data.frame(to.do, N=all.n, converged=all.conv, all.coeff.mat)
  names(xsum)[3:4]=c("min", "max")
  return(list(detailed=all.coeff.mat, summary=xsum))
}

###########################################################################################################

###########################################################################################################
# boot.glmmTMB - function
###########################################################################################################

boot.glmmTMB<-function(m, data, excl.non.conv=F, nboots=1000, para=F, resol=100, level=0.95, 
	use=NULL, contr=NULL, circ.var.name=NULL, circ.var=NULL, 
	n.cores=c("all-1", "all"), save.path=NULL, load.lib=T, lib.loc=.libPaths(), set.all.effects.2.zero=F){
	if(!summary(m)$link%in%c("log", "logit", "identity")){
		stop("link functions other than log, logit, or identity aren't supported yet; please contact Roger")
	}
	if(load.lib){library(lme4, lib.loc=lib.loc)}
  if(!is.null(use)){
		if(is.list(use)){use.list=use}else{use.list=list(use)}
	}else{
		use.list=NULL
	}
	print("doesn't account for circular variables in the zero-inflated model")
	print("doesn't consider the dispersion part in the bootstrap of fitted values")
	print("haven't tested it yet with a cbind response")
	print("weights aren't supported yet")
	n.cores=n.cores[1]
	if(is.null(contr)){
		contr=glmmTMBControl()
	}
	extract.ranef.glmmTMB<-function(m){
		to.do=summary(m)$varcor
		xx=lapply(to.do, function(x){
			lapply(x, attr, "stddev")
		})
		res=suppressWarnings(data.frame(
			part=rep(names(to.do), unlist(lapply(xx, function(x){sum(unlist(lapply(x, length)))}))),
			grp=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){unlist(lapply(x, length))}))),
			var1=unlist(lapply(xx, function(x){unlist(lapply(x, names))})),
			var2=NA,
			sdcor=unlist(xx)
		))
		xx=lapply(to.do, function(x){
			lapply(x, attr, "correlation")
		})
		xx=xx[unlist(lapply(xx, length))>0]
		xx=suppressWarnings(data.frame(
			part=rep(names(xx), unlist(sum(unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))))),
			grp=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))),
			var1=unlist(lapply(xx, function(x){lapply(x, function(x){rep(rownames(x), times=ncol(x))})})),
			var2=unlist(lapply(xx, function(x){lapply(x, function(x){rep(colnames(x), each=nrow(x))})})),
			sdcor=unlist(xx)
		))
		xx=xx[as.character(xx[, "var1"])<as.character(xx[, "var2"]), ]
		res=rbind(res, xx)
		res=res[order(as.character(res$part), as.character(res$grp), as.character(res$var1)), ]
		rownames(res)=NULL
		return(res)
	}
	extract.BLUPs.from.glmmTB<-function(m){
		to.do=lapply(c("cond", "zi"), function(x){
			if(length(ranef(m)[[x]])>0){
			lapply(ranef(m)[[x]], function(y){
				return(list(
					what=rep(x, prod(dim(as.matrix(y)))),
					blup=unlist(c(as.matrix(y))),
					level=rep(rownames(y), times=ncol(y)),
					Name=rep(colnames(y), each=nrow(y))
				))
			})
			}else{
				return(list(list(what=NULL, blup=NULL, level=NULL, Name=NULL)))
			}
		})
		xc=sapply(lapply(to.do[[1]], "[[", "what"), length)
		xz=sapply(lapply(to.do[[2]], "[[", "what"), length)
		res=data.frame(
			what=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "what")}))),
			Grp=c(rep(names(xc), times=xc), rep(names(xz), times=xz)),
			Name=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "Name")}))),
			level=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "level")}))),
			blup=as.vector(unlist(lapply(to.do, function(x){lapply(x, "[[", "blup")})))
		)
		res=res[order(res$what, res$Grp, res$Name, res$level), ]
		return(res)
	}
	extract.all<-function(mres){
		extract.ranef.glmmTMB<-function(m){
			to.do=summary(m)$varcor
			xx=lapply(to.do, function(x){
				lapply(x, attr, "stddev")
			})
			res=suppressWarnings(data.frame(
				part=rep(names(to.do), unlist(lapply(xx, function(x){sum(unlist(lapply(x, length)))}))),
				grp=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){unlist(lapply(x, length))}))),
				var1=unlist(lapply(xx, function(x){unlist(lapply(x, names))})),
				var2=NA,
				sdcor=unlist(xx)
			))
			xx=lapply(to.do, function(x){
				lapply(x, attr, "correlation")
			})
			xx=xx[unlist(lapply(xx, length))>0]
			xx=suppressWarnings(data.frame(
				part=rep(names(xx), unlist(sum(unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))))),
				grp=rep(unlist(lapply(to.do, names)), unlist(lapply(xx, function(x){lapply(x, function(x){prod(dim(x))})}))),
				var1=unlist(lapply(xx, function(x){lapply(x, function(x){rep(rownames(x), times=ncol(x))})})),
				var2=unlist(lapply(xx, function(x){lapply(x, function(x){rep(colnames(x), each=nrow(x))})})),
				sdcor=unlist(xx)
			))
			xx=xx[as.character(xx[, "var1"])<as.character(xx[, "var2"]), ]
			res=rbind(res, xx)
			res=res[order(as.character(res$part), as.character(res$grp), as.character(res$var1)), ]
			rownames(res)=NULL
			return(res)
		}
		##extract random effects model summary:
		m=extract.ranef.glmmTMB(m=mres)
		##extract random effects, BLUPs:
		BLUPs=extract.BLUPs.from.glmmTB(m=mres)
		##dispersion parameter:
		dp=summary(mres)$sigma
		##fixed effects
		xx=fixef(mres)
		fe=data.frame(
			what=rep(names(xx), times=sapply(xx, length)),
			term=unlist(lapply(xx, names)),
			est=unlist(xx)
		)
		fe=fe[order(fe$what, fe$term), ]
		return(list(fe=fe, sigma=dp, m=m, BLUPs=BLUPs))
	}
	##prepare for new calls:
	xcall=as.character(m$call)
	names(xcall)=names(m$call)
	xcall=list(cond.form=xcall["formula"], zi.form=xcall["ziformula"], disp.form=xcall["dispformula"], xfam=xcall["family"])
	if(grepl(xcall[["xfam"]], pattern="(", fixed=T)){
		xfam=xcall[["xfam"]]
		xfam=unlist(strsplit(xfam, split="(", fixed=T))
		xfam[2]=gsub(x=xfam[2], pattern=")", replacement="", fixed=T)
		xfam[2]=gsub(x=xfam[2], pattern="link = ", replacement="", fixed=T)
		xfam[2]=gsub(x=xfam[2], pattern="\"", replacement="", fixed=T)
		if(substr(xcall[["xfam"]], start=1, stop=5)!="Gamma"){
			xcall[["xfam"]]=get(xfam[1])(xfam[2])
		}else{
			if(xfam[2]=="log"){
				xcall[["xfam"]]=Gamma(link="log")
			}else if(xfam[2]=="inverse"){
				xcall[["xfam"]]=Gamma(link="inverse")
			}else if(xfam[2]=="identity"){
				xcall[["xfam"]]=Gamma(link="identity")
			}else{
				stop("Error: family not supported")
			}
		}
	}
		
	xweights=try(m$frame[, "(weights)"], silent=T)
	if(class(xweights)[[1]]=="try-error"){
		xweights=rep(1, nrow(m$frame))
	}
	data$xweights.X.=xweights
	rv.name=gsub(unlist(strsplit(xcall[["cond.form"]], split="~"))[1], pattern=" ", replacement="", fixed=T)
	##define function doing the bootstrap:
	boot.fun<-function(x, xcall., data., rv.name., m., contr., excl.non.conv., save.path., extract.all.){
		xdone=F
		while(!xdone){
			done2=F
			while(!done2){
				data.[, rv.name.]=simulate(object=m.)[, 1]
				if(xcall[["xfam"]][["family"]]!="beta" | 
					(xcall[["xfam"]][["family"]]=="beta" & min(data.[, rv.name.])>0 & max(data.[, rv.name.])<1)){
					done2=T
				}
			}	
			i.res=glmmTMB(formula=as.formula(xcall.[["cond.form"]]), ziformula=as.formula(xcall.[["zi.form"]]), dispformula=as.formula(xcall.[["disp.form"]]), 
				data=data., control=contr., family=xcall[["xfam"]], weights=xweights.X.)
				if(i.res$sdr$pdHess | !excl.non.conv.){
					xdone=T
				}
		}
		conv=i.res$sdr$pdHess
		i.res=extract.all.(mres=i.res)
		if(length(save.path.)>0){save(file=paste(c(save.path., "/b_", x, ".RData"), collapse=""), list=c("i.res", "conv"))}
		return(list(boot.res=i.res, conv=conv))
	}
	##do bootstrap:
	if(para){
		library(parallel)
		cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores!="all"){
			if(n.cores=="all-1"){n.cores=length(cl)-1}
			if(n.cores<length(cl)){
				cl=cl[1:n.cores]
			}
		}
		parLapply(cl=cl, 1:length(cl), fun=function(x, .lib.loc=lib.loc){
		  library(glmmTMB, lib.loc=.lib.loc)
		  return(invisible(""))
		})
		all.res=parLapply(cl=cl, X=1:nboots, fun=boot.fun, m.=m, xcall.=xcall, data.=data, rv.name.=rv.name, contr.=contr, 
			excl.non.conv.=excl.non.conv, save.path.=save.path, extract.all.=extract.all)
	}else{
    all.res=lapply(X=1:nboots, FUN=boot.fun, m.=m, xcall.=xcall, data.=data, rv.name.=rv.name, contr.=contr, 
			excl.non.conv.=excl.non.conv, save.path.=save.path, extract.all.=extract.all)
	}
	##extract results:
	##convergence:
	all.conv=sapply(all.res, "[[", "conv")
	##sigma:
	all.sigma=sapply(all.res, function(x){x$boot.res$sigma})
	##estimates, fixef effects:	
	all.fe=lapply(all.res, function(x){x$boot.res$fe$est})
	all.fe=matrix(unlist(all.fe), nrow=length(all.fe), byrow=T)
	colnames(all.fe)=apply(all.res[[1]]$boot.res$fe[, c("what", "term")], 1, paste, collapse="@")
	##estimates, random effects:
	all.re=lapply(all.res, function(x){x$boot.res$m$sdcor})
	all.re=matrix(unlist(all.re), nrow=length(all.re), byrow=T)
	colnames(all.re)=apply(all.res[[1]]$boot.res$m[, c("part", "grp", "var1", "var2")], 1, paste, collapse="@")
	##estimates, BLUPs (well, it doesn't make any sense  to consider them):
	all.BLUPs=lapply(all.res, function(x){x$boot.res$BLUPs$blup})
	all.BLUPs=matrix(unlist(all.BLUPs), nrow=length(all.BLUPs), byrow=T)
	colnames(all.BLUPs)=apply(all.res[[1]]$boot.res$BLUPs[, c("what", "Grp", "Name", "level")], 1, paste, collapse="@")
	##store results:
	all.ind.boots=list(fe=all.fe, re=all.re, sigma=all.sigma, conv=all.conv)
	##get confidence intervals and merge them with original values:
	##fixed effects:
		orig=extract.all(m)
		xx=orig$fe
		xx$comb=apply(xx[, c("what", "term")], 1, paste, collapse="@")
		#xx$comb=paste(rownames(xx), xx$comb, sep="@")
		ci=apply(all.fe, 2, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T)
		#times=unlist(lapply(summary(m)$coefficients, nrow))
		#comb=paste(rep(c("cond", "zi", "disp")[1:length(times)], times=times), unlist(lapply(summary(m)$coefficients, rownames)), sep="@")
		ci=ci[, match(xx$comb, colnames(ci))]
		#xx=xx[match(comb, xx$comb), ]
		ci.fe=data.frame(orig=xx$est, t(ci))
		##order according to original:
		ci.fe=ci.fe[paste(rep(c("cond", "zi", "disp"), times=unlist(lapply(fixef(m), length))), unlist(lapply(fixef(m), names)), sep="@"), ]
	##random effects:
		xx=orig$m
		xx$comb=apply(xx[, c("part", "grp", "var1", "var1")], 1, paste, collapse="@")
		ci=apply(all.re, 2, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T)
		ci.re=data.frame(orig=xx$sdcor, t(ci))
	##store results:
	ci.estimates=list(fe=ci.fe, re=ci.re)

	if(length(use.list)>0){
		ci.fitted.cond=NULL
		ci.fitted.zi=NULL
		ci.fitted.disp=NULL
		#extract fixed effects terms from the model:
		model.terms=c(attr(terms(as.formula(xcall[["cond.form"]])), "term.labels"), 
			attr(terms(as.formula(xcall[["zi.form"]])), "term.labels"), 
			attr(terms(as.formula(xcall[["disp.form"]])), "term.labels"))
		#exclude random effects:
		model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
		#exclude random effects, interactions and squared terms from model.terms:
		model.terms=model.terms[!grepl(x=model.terms, pattern="I(", fixed=T)]
		model.terms=model.terms[!grepl(x=model.terms, pattern=":", fixed=T)]
		model.terms=unique(model.terms)
		
    all.fitted=lapply(use.list, function(use){
			#create new data to be used to determine fitted values:
			if(length(use)==0){use=model.terms}
			new.data=vector("list", length(model.terms))
			if(length(circ.var.name)==1){
				set.circ.var.to.zero=sum(circ.var.name%in%use)==0
			}else{
				set.circ.var.to.zero=F
			}
			usel=model.terms%in%use
			#if(length(use)>0)
			for(i in 1:length(model.terms)){
				if(is.factor(data[, model.terms[i]])){
					new.data[[i]]=levels(data[, model.terms[i]])
				}else if(!is.factor(data[, model.terms[i]]) & usel[i] & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
					new.data[[i]]=seq(from=min(data[, model.terms[i]]), to=max(data[, model.terms[i]]), length.out=resol)
				}else  if(!is.factor(data[, model.terms[i]]) & ifelse(length(circ.var.name)==0, T, !grepl(x=model.terms[i], pattern=circ.var.name))){
					new.data[[i]]=mean(data[, model.terms[i]])
				}
			}
			names(new.data)=model.terms

			if(length(circ.var.name)==1){
				new.data=new.data[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
				if(sum(grepl(pattern=circ.var.name, x=use))>0){
					new.data=c(new.data, list(seq(min(circ.var, na.rm=T), max(circ.var, na.rm=T), length.out=resol)))
					names(new.data)[length(new.data)]=circ.var.name
				}else{
					new.data=c(new.data, list(0))
				}
				model.terms=model.terms[!(model.terms%in%paste(c("sin(", "cos("), circ.var.name, ")", sep=""))]
			}
			xnames=names(new.data)
			#browser()
			new.data=data.frame(expand.grid(new.data))
			names(new.data)=xnames
			#browser()
			if(length(circ.var.name)==1){
				names(new.data)[ncol(new.data)]=circ.var.name
			}
			if(set.all.effects.2.zero){
				for(iterm in setdiff(colnames(new.data), c("(Intercept)", use))){
					new.data[, iterm]=0
				}
			}
			##prepare model frame for conditional model for prediction:
			model=attr(terms(as.formula(xcall[["cond.form"]])), "term.labels")
			model=model[!grepl(x=model, pattern="|", fixed=T)]
			if(length(model)==0){model="1"}
			cond.m.mat=model.matrix(object=as.formula(paste(c("~", paste(model, collapse="+")), collapse="")), data=new.data)
			if(set.circ.var.to.zero){
				cond.m.mat[,paste(c("sin(", circ.var.name, ")"), collapse="")]=0
				cond.m.mat[,paste(c("cos(", circ.var.name, ")"), collapse="")]=0
			}
			##get bootstrapped fitted values for the conditional part:
			ests=all.ind.boots$fe[, substr(colnames(all.ind.boots$fe), start=1, stop=5)=="cond@", drop=F]
			est.names=gsub(colnames(ests), pattern="cond@", replacement="")
			ci.cond=lapply(1:nrow(ests), function(x){
				return(cond.m.mat[, est.names, drop=F]%*%as.vector(ests[x, ]))
			})
			ci.cond=matrix(unlist(ci.cond), ncol=nboots, byrow=F)
			ci.cond=t(apply(ci.cond, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
			colnames(ci.cond)=c("lower.cl", "upper.cl")
			##get fitted value for the conditional part:
			ests=ci.estimates$fe
			ests=ests[substr(rownames(ests), start=1, stop=5)=="cond@", , drop=F]
			est.names=gsub(rownames(ests), pattern="cond@", replacement="")
			fv.cond=cond.m.mat[, est.names, drop=F]%*%ests$orig
			if(summary(m)$link=="log"){
				ci.fitted.cond=data.frame(fitted=exp(fv.cond), exp(ci.cond))
			}else if(summary(m)$link=="logit"){
				ci.fitted.cond=data.frame(fitted=plogis(fv.cond), plogis(ci.cond))
			}else if(summary(m)$link=="identity"){
				ci.fitted.cond=data.frame(fitted=fv.cond, ci.cond)
			}
			ci.fitted.cond=data.frame(new.data, ci.fitted.cond)
			
			##prepare model frame for ZI model for prediction:
			if(xcall[["zi.form"]]!="~0"){
				model=attr(terms(as.formula(xcall[["zi.form"]])), "term.labels")
				model=model[!grepl(x=model, pattern="|", fixed=T)]
				if(length(model)==0){model="1"}
				zi.m.mat=model.matrix(object=as.formula(paste(c("~", paste(model, collapse="+")), collapse="")), data=new.data)
				##get bootstrapped fitted values for the conditional part:
				ests=all.ind.boots$fe[, substr(colnames(all.ind.boots$fe), start=1, stop=3)=="zi@", drop=F]
				est.names=gsub(colnames(ests), pattern="zi@", replacement="")
				ci.zi=lapply(1:nrow(ests), function(x){
					return(zi.m.mat[, est.names, drop=F]%*%as.vector(ests[x, ]))
				})
				ci.zi=matrix(unlist(ci.zi), ncol=nboots, byrow=F)
				ci.zi=t(apply(ci.zi, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
				colnames(ci.zi)=c("lower.cl", "upper.cl")
				##get fitted value for the conditional part:
				ests=ci.estimates$fe
				ests=ests[substr(rownames(ests), start=1, stop=3)=="zi@", , drop=F]
				est.names=gsub(rownames(ests), pattern="zi@", replacement="")
				fv.zi=zi.m.mat[, est.names, drop=F]%*%ests$orig
				ci.fitted.zi=data.frame(fitted=plogis(fv.zi), plogis(ci.zi))
				#ci.fitted.zi=ci.fitted*(1-ci.fitted.zi)
				ci.fitted.zi=data.frame(new.data, ci.fitted.zi)
			}

			##prepare model frame for disp model for prediction:
			if(!is.null(summary(m)$coefficients$disp)){
				model=attr(terms(as.formula(xcall[["disp.form"]])), "term.labels")
				model=model[!grepl(x=model, pattern="|", fixed=T)]
				if(length(model)==0){model="1"}
				disp.m.mat=model.matrix(object=as.formula(paste(c("~", paste(model, collapse="+")), collapse="")), data=new.data)
				##get bootstrapped fitted values for the conditional part:
				ests=all.ind.boots$fe[, substr(colnames(all.ind.boots$fe), start=1, stop=5)=="disp@", drop=F]
				est.names=gsub(colnames(ests), pattern="disp@", replacement="")
				ci.disp=lapply(1:nrow(ests), function(x){
					return(disp.m.mat[, est.names, drop=F]%*%as.vector(ests[x, ]))
				})
				ci.disp=matrix(unlist(ci.disp), ncol=nboots, byrow=F)
				ci.disp=t(apply(ci.disp, 1, quantile, prob=c((1-level)/2, 1-(1-level)/2), na.rm=T))
				colnames(ci.disp)=c("lower.cl", "upper.cl")
				##get fitted value for the dispersion part:
				ests=ci.estimates$fe
				ests=ests[substr(rownames(ests), start=1, stop=5)=="disp@", , drop=F]
				est.names=gsub(rownames(ests), pattern="disp@", replacement="")
				fv.disp=disp.m.mat[, est.names, drop=F]%*%ests$orig
				ci.fitted.disp=data.frame(fitted=exp(fv.disp), exp(ci.disp))
				#ci.fitted.disp=ci.fitted*(1-ci.fitted.disp)
				ci.fitted.disp=data.frame(new.data, ci.fitted.disp)
			}
			ci.fitted=list(
				cond=ci.fitted.cond, 
				zi=ci.fitted.zi,
				disp=ci.fitted.disp
			)
			return(ci.fitted)
		})
		if(length(use.list)>1){
			names(all.fitted)=unlist(lapply(use.list, paste, collapse="@"))
		}else{
			all.fitted=all.fitted[[1]]
		}
	}else{
		all.fitted=NULL
	}
	return(list(ci.estimates=ci.estimates, ci.fitted=all.fitted, all.boots=all.ind.boots))
}