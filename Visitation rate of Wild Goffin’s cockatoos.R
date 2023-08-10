############################################################################################################
##								Visitation rates of Wild Goffin’s cockatoos								  ##
##																										  ##
## 					Mark O'Hara, Alice M. I. Auersperg, Ludwig Huber, Dewi M. Prawiradilaga				  ##
##						 					& Berenika Mioduszewska										  ##
############################################################################################################

############################################################################################################
## loading required packages ###############################################################################
############################################################################################################
library(lme4) #for GLMM modelling
library(lsmeans) #for factorlevel post hoc tests
library(ggplot2) #for nice plotting
library(scico) #colorscale for sensitive to colorblindness
library(chron) #for timescaling
library (glmmTMB) #for GLMM modelling of hurdle models
library(bbmle) # for AICtab
library(cowplot) 
library(patchwork)
############################################################################################################
## loading data and custom functions #######################################################################
############################################################################################################
setwd("~/[download_directory]")#set working directory to the download directoy or directory containing dataset
source('Functions.r') # load additional custom functions for model plotting
datafull <-read.table(file="Data.txt", header=TRUE,na.strings="na", stringsAsFactors=T) # load data for visitation rates and weather data
wdat<-read.table(file="Weather.txt", header=TRUE,na.strings="na", stringsAsFactors=T) # load observed weather data for study periode


# Datapreparation
data<-subset(datafull,daysbait<=10)

#combining Date and Time and turning them into radians
data$Datetime<-strptime(data$datetime,format='%Y/%m/%d:%H:%M:%S')
#-> per day two maxima in the pooled visits, hence to model the response with two peaks convert date to radians with dayly cycle (cycle duration = 24h) and multiply by 4x pi and offset by 7.11 to set the peaks at 6:00 and 18:00
data$time.rad=-(pi*6/6)+4*pi*as.numeric((data$timemin/60))/24
data$date<-as.Date(data$Date,format='%d.%m.%Y')
data$n.date<-as.numeric(data$date)
data$nested<-as.factor(paste(data$Location,data$Bout,data$n.date,sep=".")) # explicitly specify nested structure of the data
plot(cos(data$time.rad)~as.numeric((data$timemin/60)))


#check offset for peaks
wdat$Datetime<-strptime(wdat$Datetime,format='%Y/%m/%d:%H:%M:%S')
wdat$time.rad=-(pi*6/6)+(4*pi*(wdat$timemin/60)/24)         # 1unit=1/2h 12h...2pi
wdat$x<-cos(wdat$time.rad)
summary(wdat)
ggplot(wdat,aes(as.numeric((wdat$timemin/60)),x))+geom_point()
wdat[which.max(wdat$x),]
boxplot(wdat $Wind~round(wdat $timemin/60))
boxplot((wdat $Temp)~round(wdat $timemin/60))
boxplot((wdat $Lux)~round(wdat $timemin/60))
boxplot(wdat $UVI~round(wdat $timemin/60))
boxplot(wdat $Rain~round(wdat $timemin/60))
boxplot(wdat $Hum~round(wdat $timemin/60))
###########################################################################################################
### Analysis	############################################################################################
###########################################################################################################
##inspect data
str(data)
summary(data)

#check response variable
hist(subset(data,Manualcount>0)$Manualcount)
#check the predictors
plot(data$Location,data$Manualcount)
plot(data$Bout,data$Manualcount)
hist(data$Temp)
hist(data$Lux)
hist(data$Rain)
hist(data$Wind)
hist(data$Humidity)

#pooled Visits per hour over all days
boxplot(data$Manualcount~round(data$timemin/60))
boxplot(data$Wind~round(data$timemin/60))
boxplot(data$Temp~round(data$timemin/60))
boxplot(data$Lux~round(data$timemin/60))
boxplot(data$UVI~round(data$timemin/60))
boxplot(data$Rain~round(data$timemin/60))
boxplot(data$Humidity~round(data$timemin/60))

#pca for weather measures that are strongly correlated
num_data<-data[,11:16]
head(num_data)
norm_data<-scale(num_data)
head(norm_data)
cor_mat<-cor(norm_data)
p.mat <- cor_pmat(cor_mat)
library(ggcorrplot)
ggcorrplot(cor_mat,outline.color = "white",lab=T,p.mat=p.mat,insig="pch",pch=1,pch.cex=15)
ggsave("FigS1.Correlationmatrix_weatherdata.pdf", plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 20  , height = 20, units = "cm", dpi = 500)


PCA<-(prcomp(data[,c("Temp","Humidity","UVI","Lux")],scale=T))
summary(PCA) # Component 1 explains 88.6% of Variance -> use PC1 as placeholder for Temp, Lux and Humidity
PCA$rotation
screeplot(PCA, addlabels = TRUE)
eigen(PCA)
PCA$x
data$PC1<-PCA$x[,1]

#check which random slopes to include (function by Roger Mundry)
xx.fe.re=fe.re.tab(fe.model="Manualcount~Playback+Bout+time.rad+PC1+Rain+Wind+daysbait",
       re="(1|Location)+(1|Bout)+(1|n.date)",data=data)
xx.fe.re$summary    
# use random slopes for days since baiting and recurring foraging bouts to 
#for easier interpretability of the model output z-transform covariates
data$z.dbait=as.vector(scale(data$daysbait))
hist(data$z.dbait)
boxplot(data$Manualcount~data$daysbait) #!!!!!!
data$z.rain=as.vector(scale(data$Rain))
hist(data$z.rain)
boxplot(data$Manualcount~data$Rain)
data$z.wind=as.vector(scale(data$Wind))
hist(data$z.wind)
data$z.pc=as.vector(scale(data$PC1))
hist(data$PC1)
		
#Hurdle model with glmmTMB splitting visits and null visits and number of visits if visited (cannot handle random slopes)	                                                
fit_hnpois <- glmmTMB(Manualcount~Playback+sin(time.rad)+cos(time.rad)+PC1+z.rain+z.wind+z.dbait+(1+sin(time.rad)+cos(time.rad)+z.dbait||Location)+(1+sin(time.rad)+cos(time.rad)+z.dbait||nested),
                                    data=data,
                                   	ziformula=~.,
			                        family=truncated_poisson)



##Model Validation
#Autocorrelation-VIF
library(car)
xx=lm(Manualcount~Playback+sin(time.rad)+cos(time.rad)+PC1+z.rain+z.wind+z.dbait,data=data)
vif(xx)# acitvity bouts (modelled as sin(time.rad)) correlate with intrinsic weather data (as PC1-representing Temp, Humidity and Light/UVI)

#Test for overdispersion:
sqrt(nobs(fit_hnpois) / (1+df.residual(fit_hnpois)))#Dispersion parameter
library(performance)
check_overdispersion(fit_hnpois)#no indications for overdispersion

#Check Distribution of random effects
ranef.diagn.plot(fit_hnpois)

#Check model stability
full.stab=glmmTMB.stab(model.res=fit_hnpois, data=data)
table(full.stab$detailed$opt.warnings)
round(full.stab$summary[, -1], 3)
m.stab.plot(full.stab$summary[, -1])

#Full-null(reduced) model comparison

red<- glmmTMB(Manualcount~1+(1+sin(time.rad)+cos(time.rad)||Location)+(1+sin(time.rad)+cos(time.rad)||nested),
                                    data=data,
                                   	ziformula=~.,
			                        family=truncated_poisson)

as.data.frame(anova(red, fit_hnpois, test="Chisq"))

#overview of fixed effects
round(summary(fit_hnpois)$coefficients$cond,3) #poisson model
round(summary(fit_hnpois)$coefficients$zi,3) #poisson model
#overview of random effects
summary(fit_hnpois)$varcor
mean(data$Rain)
sd(data$Rain)
mean(data$Wind)
sd(data$Wind)
mean(data$daysbait)
sd(data$daysbait)

#interpretation
round(summary(fit_hnpois)$coefficients$zi,3)
#backtransformation of effect sizes based on z.transformed values
est.orig=fixef(fit_hnpois)$zi["z.rain"]/sd(data$Rain)

exp(fixef(fit_hnpois)$zi["(Intercept)"])/
      (1+exp(fixef(fit_hnpois)$zi["(Intercept)"]))
      
exp(fixef(fit_hnpois)$cond["(Intercept)"])

#bootstrapping confidence intervals (use with caution, takes a long time)
full.boot=boot.glmmTMB(fit_hnpois, data=data,
       excl.non.conv=F, nboots=1000, para=T, resol=100,
       level=0.95, use=NULL, contr=NULL, n.cores=c("all-1", "all"))
full.boot$ci.estimates$re

############################################################################################################    
###### Graphs ##############################################################################################
############################################################################################################
#### effects on cameratrap images w/wo individuals of goffins (binomial model)
#preparing data and centering factor playback
data$nz <- as.numeric(data$Manualcount>0) 
data$s <- sin(data$time.rad)
data$c <- cos(data$time.rad)
data$play.code=as.numeric(data$Playback==levels(data$Playback)[2])
data$play.code=data$play.code-mean(data$play.code)
#preparing plot for correct over PC1/time
xx=aggregate(x=data$nz, by=data[, c("nz", "PC1","timemin")],
        FUN=mean)
xx$N=aggregate(x=data$nz, by=data[, c("nz", "PC1","timemin")],
        FUN=length)$x  
summary(xx)
################################
### plotting temporal foraging bouts based on the model descriptive in terms of images w/wo goffins
ggplot(data,aes(x=timemin/60,color=as.factor(nz),fill=as.factor(nz)))+geom_point(data=xx,aes(timemin/60,x,size=N),alpha=.05)+scale_size_continuous(range = c(1,7),name="Number of\nimages")+geom_density(mapping = aes(y = after_stat(count)),alpha=.1,adjust=.3)+scale_x_discrete(limits=factor(0:24),name="Time of day [h]")+ylab("Total # of images")+scale_fill_discrete("Images",breaks=c(0,1),labels=c("without Goffins","with Goffins"))+scale_color_discrete("Images",breaks=c(0,1),labels=c("without Goffins","with Goffins"))+	theme_minimal()
#save the graph as pdf
ggsave("Figx1._daily.visits.pdf", plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 20  , height = 8, units = "cm", dpi = 500)

################################
### plotting effect of intrinsic weather variation


# fitting a binomial model for plotting
plot.res=glmer(nz~play.code+s+c+PC1+z.rain+z.wind+z.dbait+(1+s+c+z.dbait||Location)+(1+s+c+z.dbait||nested),data=data,family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))

bootstraping confidence intervals (takes a long time)
boot.res.plot=boot.glmm.pred(model.res=plot.res,
       excl.warnings=F, nboots=1000, para=F, resol=1000,
       level=0.95, use="PC1")

#add model ine and CI from bootstrapping in a dataframe
plot.xvals=seq(from=min(data$PC1), to=max(data$PC1),
        length.out=1000)
plot.PC1=data.frame(x=plot.xvals,	
				y=boot.res.plot$ci.predicted$fitted,
				ylwr= boot.res.plot$ci.predicted$lower.cl,
				yupr=boot.res.plot$ci.predicted$upper.cl) #make dataframe with model and lines for CI

# ploting the data and CI from the bootstrapped model
ggplot(data,aes(PC1,nz))+
geom_point(data=xx,aes(PC1,x,size=N),alpha=.1)+scale_size_continuous(range = c(1,5 ),name="Number of\nrecorded visits")+
	geom_line(data= plot.PC1,aes(x, y),lwd=1,inherit.aes=F)+
	geom_line(data= plot.PC1,aes(x, yupr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
	geom_line(data= plot.PC1,aes(x, ylwr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
	geom_ribbon(data= plot.PC1,aes(x,ymin=ylwr, ymax=yupr), alpha=0.2,inherit.aes=F)+
	labs(x="PC1",y="Probability for visiting a platform")+
	theme_minimal()
#save the graph as pdf
ggsave("Figx2_PC1visits.pdf", plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 20  , height = 10, units = "cm", dpi = 500)

############################################################################################################
####  plotting effects of time on number of individuals (conditional/poisson model)
# manually removing zeros for plotting
datac<-subset(data,Manualcount>0)
fit_hpois <- glmer(Manualcount~play.code+s+c+PC1+z.rain+z.wind+z.dbait+(1+s+c+z.dbait||Location)+(1+s+c+z.dbait||nested),data=datac,family="poisson", control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))
summary(fit_hpois)
################################
### plotting foraging bouts
#aggregating binnned number of individuals per hour
datac$time.b<-as.numeric(datac$timemin/60)
to.plot=aggregate(x=1:nrow(datac),
      by=datac[, c("time.b", "Manualcount","time.rad")], FUN=length)
summary(to.plot)

################################
### Individuals by weather measures
#solar intensity
wdat$time.b<-wdat$timemin/60
p1<-ggplot(wdat, aes(time.b,Lux/1000)) + geom_smooth(color=scico(1,begin=0.2,end=0.2,palette="lajolla"),fill=scico(1,alpha=0.2,begin=0.2,end=0.2,palette="imola")) + scale_x_continuous("Hours of day",breaks=c(0:24))+geom_point(data=to.plot,aes(x=time.b,y=Manualcount*10,size=x),alpha=.1,show.legend=T)+scale_y_continuous(position="left",name="Solar intensity [kLux]",sec.axis=sec_axis(trans=(~./10),"Individuals feeding",breaks=c(0,1,2,3,4,5,6,7,8,9,10,11)))+scale_size_continuous(range = c(2, 7),name="Number of\ntriggered images")+theme_minimal()+theme(axis.line.y.right = element_line(color = scico(1,begin=0.2,end=0.2,palette="lajolla")), 
        axis.ticks.y.left = element_line(color = scico(1,begin=0.2,end=0.2,palette="lajolla")),
        axis.text.y.left = element_text(color = scico(1,begin=0.2,end=0.2,palette="lajolla")), 
        axis.title.y.left = element_text(color = scico(1,begin=0.2,end=0.2,palette="lajolla")))

#Temp
p2<-ggplot(wdat, aes(time.b,Temp)) + geom_smooth(color=scico(1,begin=0.6,end=0.6,palette="lajolla"),fill=scico(1,alpha=.2,begin=0.6,end=0.6,palette="lajolla")) + scale_x_continuous("Hours of day",breaks=c(0:24))+geom_point(data=to.plot,aes(x=time.b,y=Manualcount*3,size=x),alpha=.1,show.legend=T)+scale_y_continuous(position="left",name="Temperature [°C]",breaks=c(5,10,15,20,25,30),sec.axis=sec_axis(trans=(~./3),"Individuals feeding",breaks=c(0,1,2,3,4,5,6,7,8,9,10,11)))+scale_size_continuous(range = c(2, 7),name="Number of\ntriggered images")+theme_minimal()+theme(axis.line.y.right = element_line(color = scico(1,begin=0.6,end=0.6,palette="lajolla")), 
        axis.ticks.y.left = element_line(color = scico(1,begin=0.6,end=0.6,palette="lajolla")),
        axis.text.y.left = element_text(color = scico(1,begin=0.6,end=0.6,palette="lajolla")), 
        axis.title.y.left = element_text(color = scico(1,begin=0.6,end=0.6,palette="lajolla")))

#Humidity
p3<-ggplot(wdat, aes(time.b,Hum)) +geom_point(data=to.plot,aes(x=time.b,y=Manualcount*10,size=x),alpha=.1,show.legend=T)+geom_smooth(data=wdat,aes(x=time.b,y=Temp*100/30),color=scico(1,begin=0.6,end=0.6,palette="lajolla"),fill=scico(1,alpha=.2,begin=0.6,end=0.6,palette="lajolla"))+geom_smooth(color=scico(1,begin=0.1,end=0.1,palette="imola"),fill=scico(1,alpha=0.2,begin=0.1,end=0.1,palette="imola"))+geom_smooth(data=wdat,aes(x=time.b,y=Lux/1000),color=scico(1,begin=0.2,end=0.2,palette="lajolla"),fill=scico(1,alpha=0.2,begin=0.2,end=0.2,palette="lajolla")) +scale_x_continuous("Time of day [h]",breaks=c(0:24))+scale_y_continuous(position="left",name="Humidity [%]",sec.axis=sec_axis(trans=(~./10),"Individuals feeding",breaks=c(0,1,2,3,4,5,6,7,8,9,10,11)))+scale_size_continuous(range = c(2, 7),name="Number of\ntriggered images")+theme_minimal()+theme(axis.line.y.left = element_line(color = scico(1,begin=0.2,end=0.2,palette="imola")), axis.ticks.y.left = element_line(color = scico(1,begin=0,end=0.2,palette="imola")),       axis.text.y.left = element_text(color = scico(1,begin=0.2,end=0.2,palette="imola")),         axis.title.y.left = element_text(color = scico(1,begin=0.2,end=0.2,palette="imola")))

#combining graph labels
Fig.3b<- wrap_elements(get_plot_component(p2, "ylab-l")) +
  wrap_elements(get_y_axis(p2)) +
  wrap_elements(get_plot_component(p1, "ylab-l")) +
  wrap_elements(get_y_axis(p1)) +
  p3 + 
  plot_layout(widths = c(1.2, 0.2, 1.2, .2, 40))

################################
### plotting effect of intrinsic weather variation
respc=aggregate(x=datac$Manualcount, by=datac[, c("PC1", "Manualcount")],
        FUN=mean)
respc$N=aggregate(x=datac$Manualcount, by=datac[, c("PC1", "Manualcount")],
        FUN=length)$x
plot.pcxvals=seq(from=min(data$PC1), to=max(data$PC1),
         length.out=1000)
 
#bootstrap confidence intervals
boot.respc.plot=boot.glmm.pred(model.res=fit_hpois,
        excl.warnings=F, nboots=1000, para=T, resol=1000,
        level=0.95, use="PC1")

#add model ine and CI from bootstrapping in a dataframe
plot.pc=data.frame(x=plot.pcxvals,	
				y=boot.respc.plot$ci.predicted$fitted,
				ylwr= boot.respc.plot$ci.predicted$lower.cl,
				yupr=boot.respc.plot$ci.predicted$upper.cl) #make dataframe with model and lines for CI

# plot data and bootstrapped model
Fig.3a<-ggplot(data,aes(PC1,Manualcount))+
geom_point(data=respc,aes(PC1,x,size=N),alpha=.4)+scale_size_continuous(range = c(2, 7),name="Number of\ntriggered images")+
	geom_line(data= plot.pc,aes(x, y),lwd=1,inherit.aes=F)+
	geom_line(data= plot.pc,aes(x, yupr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
	geom_line(data= plot.pc,aes(x, ylwr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
	geom_ribbon(data= plot.pc,aes(x,ymin=ylwr, ymax=yupr), alpha=0.2,inherit.aes=F)+
	labs(x="Daily weather variation[PC1]",y="Number of individuals feeding")+
	theme_minimal()

plot_grid(Fig.3a, Fig.3b,  
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

ggsave("Fig.x3_hourlyindividualsbyweather.pdf", plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 23  , height = 20, units = "cm", dpi = 500)

################################
### plotting effect of rain
resrain=aggregate(x=datac$Manualcount, by=datac[, c("Rain", "Manualcount")],
		        FUN=mean)
resrain$N=aggregate(x=datac$Manualcount, by=datac[, c("Rain", "Manualcount")],
		        FUN=length)$x
		   
#bootstrap confidence intervals
boot.resrain.plot=boot.glmm.pred(model.res=fit_hpois,
		        excl.warnings=F, nboots=1000, para=T, resol=1000,
		        level=0.95, use="z.rain")
		
#add model ine and CI from bootstrapping in a dataframe
plot.rainxvals=seq(from=min(datac$z.rain), to=max(datac$z.rain),
		         length.out=1000)
plot.rain=data.frame(x=plot.rainxvals,	
				y=boot.resrain.plot$ci.predicted$fitted,
				ylwr= boot.resrain.plot$ci.predicted$lower.cl,
				yupr=boot.resrain.plot$ci.predicted$upper.cl) #make dataframe with model and lines for CI
		
# plot data and bootstrapped model
Fig.4a<-ggplot(data,aes(Rain,Manualcount))+
geom_point(data=resrain,aes(Rain,x,size=N),alpha=.4)+scale_size_continuous(range = c(2, 7),name="Number of\ntriggered images")+
		geom_line(data= plot.rain,aes(x, y),lwd=1,inherit.aes=F)+
		geom_line(data= plot.rain,aes(x, yupr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
		geom_line(data= plot.rain,aes(x, ylwr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
		geom_ribbon(data= plot.rain,aes(x,ymin=ylwr, ymax=yupr), alpha=0.2,inherit.aes=F)+
		labs(x="Rainfall[mm/h]",y="Number of individuals feeding")+
		theme_minimal()
		
################################
### plotting effect of ressource depletion
resdep=aggregate(x=datac$Manualcount, by=datac[, c("daysbait", "Manualcount")],
        FUN=mean)
resdep$N=aggregate(x=datac$Manualcount, by=datac[, c("daysbait", "Manualcount")],
        FUN=length)$x

#preparing a the model for plotting
plot.resdep=glmer(Manualcount~play.code+s+c+PC1+z.rain+z.wind+z.dbait+(1+s+c+z.dbait||Location)+(1+s+c+z.dbait||nested),data=datac,family="poisson",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))
   
#bootstrap confidence intervals
boot.resdep.plot=boot.glmm.pred(model.res=plot.resdep,
        excl.warnings=F, nboots=1000, para=T, resol=1000,
        level=0.95, use="z.dbait")

#add model ine and CI from bootstrapping in a dataframe
plot.dbxvals=seq(from=min(datac$daysbait), to=max(datac$daysbait),
         length.out=1000)
plot.daysbait=data.frame(x=plot.dbxvals,	
				y=boot.resdep.plot$ci.predicted$fitted,
				ylwr= boot.resdep.plot$ci.predicted$lower.cl,
				yupr=boot.resdep.plot$ci.predicted$upper.cl) #make dataframe with model and lines for CI

# plot data and bootstrapped model
Fig.4b<-ggplot(data,aes(daysbait,Manualcount))+
geom_point(data=resdep,aes(daysbait,x,size=N),alpha=.4)+scale_size_continuous(range = c(1, 7),name="Number of\ntriggered images")+
	geom_line(data= plot.daysbait,aes(x, y),lwd=1,inherit.aes=F)+
	geom_line(data= plot.daysbait,aes(x, yupr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
	geom_line(data= plot.daysbait,aes(x, ylwr),alpha=.6,lwd=1,lty=3,inherit.aes=F)+
	geom_ribbon(data= plot.daysbait,aes(x,ymin=ylwr, ymax=yupr), alpha=0.2,inherit.aes=F)+
	labs(x="Days after baiting",y="Number of individuals feeding")+
	theme_minimal()
plot_grid(Fig.4a, Fig.4b,  
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

ggsave("Fig.x4_Rain&Daysafterbaiting.pdf", plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 23  , height = 20, units = "cm", dpi = 500)


####
# plotting visits by location
Fig.S1a<-ggplot(data,aes(x=timemin/60,color=as.factor(nz),fill=as.factor(nz)))+geom_point(data=xx,aes(timemin/60,x,size=N),alpha=.05)+scale_size_continuous(range = c(1,7),name="Number of\nimages")+geom_density(mapping = aes(y = after_stat(count)),alpha=.1,adjust=.3)+scale_x_discrete(limits=factor(0:24),breaks=c("6","9","12","15","18"))+ylab("Total # of images")+scale_fill_discrete("Images",breaks=c(0,1),labels=c("without Goffins","with Goffins"))+scale_color_discrete("Images",breaks=c(0,1),labels=c("without Goffins","with Goffins"))+	labs(x="")+theme_minimal()+facet_wrap(.~Location,nrow=2)

to.plot2=aggregate(x=1:nrow(datac),
      by=datac[, c("time.b", "Manualcount","time.rad","Location")], FUN=length)

Fig.S1b<-ggplot(to.plot2, aes(time.b,Manualcount,size=x)) +geom_point(alpha=.1,show.legend=T)+scale_x_continuous("Time of day [h]",limits=c(0,24),breaks=c(6,9,12,15,18))+scale_y_continuous("Individuals feeding",breaks=c(0,1,2,3,4,5,6,7,8,9,10))+scale_size_continuous(range = c(2, 7),name="Number of\ntriggered images")+theme_minimal()+facet_wrap(.~Location,nrow=2)
plot_grid(Fig.S1a, Fig.S1b,  
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
          
#save the graph as pdf
ggsave("FigS2._daily.visits by location.pdf", plot = last_plot(), device = NULL, path = NULL, scale = 1, width = 30  , height = 20, units = "cm", dpi = 500)
############################################################################################################
############################################################################################################
