#remove.packages("pomp")
#install.packages("pomp",repos="https://cran.rstudio.com/bin/windows/contrib/3.6/pomp_2.8.zip")
#install.packages("pomp",repos="https://kingaa.github.io/")
#library(devtools)
#devtools::install_version("pomp",version="2.8.0",repos="https://kingaa.github.io/")
##https://cran.rstudio.com/bin/windows/contrib/3.6/pomp_2.8.zip
#devtools::install_github("kingaa/pomp")
#library(installr)
#updateR()


library(pomp)
library(ggplot2)

library(rsample)   # data splitting 

library(earth)     # fit MARS models
library(caret)     # automating the tuning process
library(vip)       # variable importance
library(pdp)       # variable relationships

#################### LOADING THE DATA ############################################
covid <- read.csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv')

covid = covid[which(covid$Province_State == "Texas"),]

covid = covid[which(covid$Admin2 %in% c("Austin")),]

covid = covid[names(covid)[ 13: length(names(covid)) ]]

covid = as.data.frame(colSums(covid))

#covid_19 = as.data.frame(tail(covid,-46))
#covid_19 = as.data.frame(tail(covid,-46))

covid_19 = as.data.frame(covid[covid!=0])

names(covid_19) = c("covid")

covid_19$time = seq(1,dim(covid_19)[1])

new_cases = covid_19
new_cases$covid = c(diff(new_cases$covid),0)

#names(covid_19) = c("date","cases","time")
################### DEATHS ########################################
muertos <- read.csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv')
muertos = muertos[which(muertos$Province_State == "Texas"),]

muertos = muertos[which(muertos$Admin2 %in% c("Austin")),]

muertos = muertos[names(muertos)[ 13: length(names(muertos)) ]]

muertos = as.data.frame(colSums(muertos))

#rownames(muertos) == rownames(covid)

#covid_19 = as.data.frame(tail(covid,-55))
muertos_19 = as.data.frame(muertos[covid!=0])

#muertos_19 = as.data.frame(muertos[muertos!=0])

names(muertos_19) = c("muertos")


muertos_19$time = seq(1,dim(muertos_19)[1])

#names(covid_19) = c("date","cases","time")

ggplot(muertos_19,mapping=aes(x=time,y=muertos))+
  geom_line()+geom_point()

####################### PLOTTING DATA OF INFECTOUS ##########################################
ggplot(covid_19,mapping=aes(x=time,y=covid))+
  geom_line()+geom_point()

fit1 <- lm(log(covid)~time,data=subset(covid_19))

#x11()
plot(covid_19$time, log(covid_19$covid))
abline(fit1)

summary(fit1)
coef(fit1)

####################
# 
# hyper_grid <- expand.grid(
#   degree = c(1,2), 
#   nprune = c(3)
# )
# 
# # for reproducibiity
# set.seed(555)
# 
# # cross validated model
# tuned_mars <- caret::train(log(covid)~time,
#   data = covid_19,
#   method = "earth",
#   metric = "RMSE",
#   trControl = caret::trainControl(method = "cv", number = 10),
#   tuneGrid = hyper_grid
# )
# 
# # best model
# tuned_mars$bestTune
# 
# # plot results
# ggplot(tuned_mars)
# 
# x11()
# plot(covid_19$time, log(covid_19$covid))
# lines(covid_19$time, predict(tuned_mars,covid_19))
# 
# summary(tuned_mars)
############################################################
# Active cases = total confirmed - total recovered - total deaths

act_rec = {}

act_rec$time = covid_19$time

act_rec$ac_rec = covid_19$covid - muertos_19$muertos

act_rec = as.data.frame(act_rec)

ggplot(act_rec,mapping=aes(x=time,y=ac_rec))+
  geom_line()+geom_point()

#########################################################################
########## DETERMINISTIC MODEL TO FIND PARAMETERS #######################
pomp(
  data=subset(muertos_19), #, select = -date
  times="time",t0=0,
  skeleton=vectorfield(
    Csnippet("
      DS = -(Beta)*S*(Itn+Itp+Is+Ih+Iv+Irs+Irh)/N;
      DE0 = (Beta)*S*(Itn+Itp+Is+Ih+Iv+Irs+Irh)/N-(1-alpha)*sigma*E0-alpha*E0;
      DEtp = alpha*E0-sigma*Etp;
      DItn = (1-alpha)*sigma*E0-alpha*Itn-gamman*Itn;
      DItp = alpha*Itn+sigma*Etp-gammap*Itp-etas*Itp;
      DIs = etas*Itp-gammasi*Is-etah*Is;
      DIh = etah*Is-gammahi*Ih-etav*Ih;
      DIv = etav*Ih-gammavi*Iv-delta*Iv;
      DIrs = gammasi*Is-gammas*Irs;
      DIrh = gammahi*Ih+gammavi*Iv-gammah*Irh;
      DR = gamman*Itn+gammap*Itp+gammas*Irs+gammah*Irh;
      DD = delta*Iv;")),
  rinit=Csnippet("
      S = N-E0_0-Etp_0-Itn_0-Itp_0-Is_0-Ih_0-Iv_0-Irs_0-Irh_0-R_0-D_0;
      E0 = E0_0;
      Etp = Etp_0;
      Itn = Itn_0;
      Itp = Itp_0;
      Is = Is_0;
      Ih = Ih_0;
      Iv = Iv_0;
      Irs = Irs_0;
      Irh = Irh_0;
      R = R_0;
      D = D_0;"),
  statenames=c("S","E0","Etp","Itn","Itp","Is","Ih","Iv","Irs","Irh","R","D"),
  paramnames=
    c("Beta","sigma","alpha","gamman","gammap","etas","gammasi","etah","etav","gammahi","gammavi",
      "E0_0","Etp_0","Itn_0","Itp_0","Is_0","Ih_0","Iv_0","Irs_0","gammas","gammah",
      "Irh_0","R_0","D_0","delta","N")) -> deter.model

#NN = 5822434 # Wisconsin
#NN = 488081 #546695    # Dane 488,081
#NN = 592025   #Milwaukee

NN = 964254  # austin

sse <- function (params) {
  x <- trajectory(deter.model, params=params, format="data.frame")
  #discrep = ( (x["Etp"]+x["Is"]+x["Itp"]+x["Ih"]+x["Iv"]+x["Irs"]+x["Irh"]+x["R"])-act_rec$ac_rec )^2
  #discrep = ( (x["C"])-act_rec$ac_rec )^2
  ##################################################################################
  vect1 =  c(x["D"])$D
  vect2 = c(obs(deter.model))
  # 
  #   discrep2 = 0
  #   for (i in 1:200){
  #     #wh = sort(sample(1:length(vect2), as.integer(0.95*length(vect2)) )) #0.975
  #     wh = sort(sample(1:length(vect2), length(vect2), replace = TRUE))
  # 
  #     vect1samp = vect1[wh]
  #     vect2samp = vect2 #[wh]
  # 
  #     discrep2 = discrep2 + sum( (vect1samp-vect2samp )^2 )
  #   }
  #   discrep2 = discrep2/200
  ##################################################################################
  ##### agregar error a medicion 2 posisson + y - 
  
  discrep2 = sum( ( vect1-vect2 )^2 )
  
  #return( (sum(discrep) + sum(discrep2)) /(2*length(muertos_19$muertos)) )
  #return(  sum(discrep)/(length(obs(deter.model))) + sum(discrep2)/(length(obs(deter.model))) )
  #return(  sum(discrep2)/( length(obs(deter.model) )  ) ) 
  return( discrep2 ) 
}

###           
f2 <- function (par) {
  params <- c(Beta=par[1], sigma=par[2], alpha=par[3], gamman=par[11], gammap=par[12],
              etas=par[4] , gammasi=par[5],
              etah=par[6], etav=par[7], gammahi=par[8], gammavi=par[9],
              E0_0=par[13], Etp_0=0, Itn_0=0, Itp_0=par[16], Is_0=0, Ih_0=0,
              Iv_0=0, Irs_0=0, gammas=par[14], gammah=par[15], Irh_0=0, R_0=0,D_0=0, 
              delta=par[10], N=NN)
  sse(params)
  
  #x <- trajectory(deter.model, params=params, format="data.frame")
  #discrep <- (par[16]*x["Is"])-obs(deter.model)
  #return (sum(discrep^2)/length(obs(deter.model)))
}

#c(betas, sigmas, alphas, etas_s, gammas_si, etas_h, etas_v, gammas_hi, gammas_vi,
#  deltas, gammas_n, gammas_p, E0s_0, gammas, gammah)

fit2 = optim(fn=f2, par=c(0.20, 0.291, 0.82, 0.183, 0.50, 0.03, 0.25, 0.50, 0.50, 0.65, 
                          1/8, 1/8, 250, 0.5, 0.5, 5),
             lower=c(0.07, 1/4, 0.7, 1/6, 0.10, 0.02, 0.10, 0.10, 0.10,  0.20, 
                     1/14, 1/14, 0, 0, 0, 1),
             upper=c(0.90, 1/3, 0.95, 1/5, 0.90, 0.04, 0.40, 0.90, 0.90, 0.80,
                     1/3, 1/3, 500, 1, 1, 10),
             method="L-BFGS-B")



# fit2 = optim(fn=f2, par=c(0.20, 0.20, 0.95, 0.09, 0.50, 0.05, 0.4812, 0.50, 0.50, 0.75, 
#                           0.30, 0.30, 5/NN, 0.5, 0.5) ,
#              lower=c(0.07, 1/5.4, 0.913, 1/61, 0.10, 1/39, 0.417, 0.10, 0.10,  0.60, 
#                      0.05, 0.05, 0/NN, 0, 0),
#              upper=c(0.50, 1/4.2, 0.997, 1/19, 0.90, 1/13,0.5454, 0.90, 0.90, 0.90,
#                      0.80, 0.80, 5000/NN, 1, 1), 
#              method="L-BFGS-B")

#"L-BFGS-B"
#0.90, 1/7, 0.0166666
#fit2$par

par = fit2$par 

coef(deter.model) <- 
  c(Beta=par[1], sigma=par[2], alpha=par[3], gamman=par[11], gammap=par[12],
    etas=par[4] , gammasi=par[5],
    etah=par[6], etav=par[7], gammahi=par[8], gammavi=par[9],
    E0_0=as.integer(par[13]), Etp_0=0, Itn_0=0, Itp_0=as.integer(par[16]), Is_0=0, Ih_0=0,
    Iv_0=0, Irs_0=0, gammas=par[14], gammah=par[15], Irh_0=0, R_0=0, D_0=0, delta=par[10], N=NN)

x <- trajectory(deter.model,format="data.frame")
dat <- merge(as.data.frame(deter.model),x,by='time')

#NN = par[16]

#dat2 <- merge(as.data.frame(act_rec),x,by='time')

dat$data = dat$muertos

ggplot(dat,aes(x=time))+
  geom_line(aes(y=data, color='data'), size=2)+
  geom_line(aes(y=D, color='D'), size=2) + ggtitle("Compartment D v.s. Death data") +
  labs(y="Population", x = "Time in days after 03-26-2020") + theme_classic()

#####################
dat1 <- merge(new_cases,x,by='time')
dat1$data = dat1$covid

ggplot(dat1,aes(x=time))+
  geom_line(aes(y=data, color='data'), size=2)+
  geom_line(aes(y=Etp+Itp+Is+Ih+Iv, color='New'), size=2) + ggtitle("Compartments v.s. Death data") +
  labs(y="Population", x = "Time in days after 03-26-2020") + theme_classic()
##################

params1 <- c(Beta=par[1], sigma=par[2], alpha=par[3], gamman=par[11], gammap=par[12],
             etas=par[4] , gammasi=par[5],
             etah=par[6], etav=par[7], gammahi=par[8], gammavi=par[9],
             E0_0=par[13], Etp_0=0, Itn_0=0, Itp_0=par[16], Is_0=0, Ih_0=0,
             Iv_0=0, Irs_0=0, gammas=par[14], gammah=par[15], Irh_0=0, R_0=0,D_0=0, delta=par[10], N=NN)

params1

# ggplot(dat2,aes(x=time))+
#  geom_line(aes(y=ac_rec),color='black')+
#  geom_line(aes(y=C),color='red')


#"S","E0","Etp","Itn","Itp","Is","Ih","Iv","Irs","Irh","R","D"
ggplot(x,aes(x=time))+
  #geom_line(aes(y=S,color='S'),size=2)+
  geom_line(aes(y=E0,color='E0'),size=2)+
  geom_line(aes(y=Etp,color='Etp'),size=2)+
  geom_line(aes(y=Itn,color='Itn'),size=2)+
  geom_line(aes(y=Itp,color='Itp'),size=2)+
  geom_line(aes(y=Is,color='Is'),size=2)+
  geom_line(aes(y=Ih,color='Ih'),size=2)+
  geom_line(aes(y=Iv,color='Iv'),size=2)+
  geom_line(aes(y=Irs,color='Irs'),size=2)+
  geom_line(aes(y=Irh,color='Irh'),size=2)+
  #geom_line(aes(y=R,color='R'),size=2)+
  geom_line(aes(y=D,color='D'),size=2)+ #coord_trans(y="log10") +
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))

#####################################################################################
################## DETERMINISTIC TRAYECTORY  ###########

#NN = 5822434 # Wisconsin

closed.sir.ode <- Csnippet("
      DS = -(Beta)*S*(Itn+Itp+Is+Ih+Iv+Irs+Irh)/N;
      DE0 = (Beta)*S*(Itn+Itp+Is+Ih+Iv+Irs+Irh)/N-(1-alpha)*sigma*E0-alpha*E0;
      DEtp = alpha*E0-sigma*Etp;
      DItn = (1-alpha)*sigma*E0-alpha*Itn-gamman*Itn;
      DItp = alpha*Itn+sigma*Etp-gammap*Itp-etas*Itp;
      DIs = etas*Itp-gammasi*Is-etah*Is;
      DIh = etah*Is-gammahi*Ih-etav*Ih;
      DIv = etav*Ih-gammavi*Iv-delta*Iv;
      DIrs = gammasi*Is-gammas*Irs;
      DIrh = gammahi*Ih+gammavi*Iv-gammah*Irh;
      DR = gamman*Itn+gammap*Itp+gammas*Irs+gammah*Irh;
      DD = delta*Iv;")

init1 <- Csnippet("
      S = N-E0_0-Etp_0-Itn_0-Itp_0-Is_0-Ih_0-Iv_0-Irs_0-Irh_0-R_0-D_0;
      E0 = E0_0;
      Etp = Etp_0;
      Itn = Itn_0;
      Itp = Itp_0;
      Is = Is_0;
      Ih = Ih_0;
      Iv = Iv_0;
      Irs = Irs_0;
      Irh = Irh_0;
      R = R_0;
      D = D_0;")

pomp(data=data.frame(time=1:150,data=NA),
     times="time",t0=0,
     skeleton=vectorfield(closed.sir.ode),
     rinit=init1,
     statenames=c("S","E0","Etp","Itn","Itp","Is","Ih","Iv","Irs","Irh","R","D"),
     paramnames=
       c("Beta","sigma","alpha","gamman","gammap","etas","gammasi","etah","etav","gammahi","gammavi",
         "E0_0","Etp_0","Itn_0","Itp_0","Is_0","Ih_0","Iv_0","Irs_0","gammas","gammah",
         "Irh_0","R_0","D_0","delta","N")) -> closed.sir

params1 <- c(Beta=par[1], sigma=par[2], alpha=par[3], gamman=par[11], gammap=par[12],
             etas=par[4] , gammasi=par[5],
             etah=par[6], etav=par[7], gammahi=par[8], gammavi=par[9],
             E0_0=par[13], Etp_0=0, Itn_0=0, Itp_0=par[16], Is_0=0, Ih_0=0,
             Iv_0=0, Irs_0=0, gammas=par[14], gammah=par[15], Irh_0=0, R_0=0,D_0=0, delta=par[10], N=NN)

x <- trajectory(closed.sir,params=params1,format="data.frame")

ggplot(x,aes(x=time))+
  #geom_line(aes(y=S,color='S'),size=2)+
  #geom_line(aes(y=E0,color='E0'),size=2)+
  #geom_line(aes(y=Etp,color='Etp'),size=2)+
  geom_line(aes(y=Itn,color='Itn'),size=2)+
  #geom_line(aes(y=Itp,color='Itp'),size=2)+
  #geom_line(aes(y=Is,color='Is'),size=2)+
  geom_line(aes(y=Ih,color='Ih'),size=2)+
  geom_line(aes(y=Iv,color='Iv'),size=2)+
  #geom_line(aes(y=Irs,color='Irs'),size=2)+
  geom_line(aes(y=Irh,color='Irh'),size=2)+
  #geom_line(aes(y=R,color='R'),size=2)+
  #geom_line(aes(y=D,color='D'),size=2)+ #coord_trans(y="log10") +
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+ ggtitle("Deterministic compartments") +
  labs(y="Population", x = "Time in days after 03-26-2020")

################################################################################
##################### STOCHASTIC SIMULATIONS ###################################

#NN = 5822434 # Wisconsin

sir_step <- Csnippet("
      double dN_SE0 = rbinom(S, 1.0-exp(-Beta*((Itn+Itp+Is+Ih+Iv+Irs+Irh)/N)*dt));
      double dN_E0Etp = rbinom(E0, 1.0-exp( -alpha*dt ) );
      double dN_E0Itn = rbinom(E0, 1.0-exp( -(1-alpha)*sigma*dt ) );
      double dN_EtpItp = rbinom(Etp, 1.0-exp( -sigma*dt ));
      double dN_ItnItp = rbinom(Itn, 1.0-exp( -alpha*dt ));
      double dN_ItnR = rbinom(Itn, 1.0-exp( -gamman*dt ));
      double dN_ItpIs = rbinom(Itp, 1.0-exp( -etas*dt ));
      double dN_ItpR = rbinom(Itp, 1.0-exp( -gammap*dt ));
      double dN_IsIrs = rbinom(Is, 1.0-exp( -gammasi*dt ));
      double dN_IsIh = rbinom(Is, 1.0-exp(-etah*dt ));
      double dN_IhIrh = rbinom(Ih, 1.0-exp(-gammahi*dt ));
      double dN_IhIv = rbinom(Ih, 1.0-exp(-etav*dt) );
      double dN_IvIrh = rbinom(Iv, 1.0-exp(-gammavi*dt));
      double dN_IrhR = rbinom(Irh, 1.0-exp(-gammah*dt));
      double dN_IrsR = rbinom(Irs, 1.0-exp(-gammas*dt));
      double dN_IvD = rbinom(Iv, 1.0-exp(-delta*dt));
      S -= dN_SE0;
      E0 += dN_SE0 - dN_E0Itn - dN_E0Etp;
      Etp += dN_E0Etp - dN_EtpItp;
      Itn += dN_E0Itn - dN_ItnItp - dN_ItnR;
      Itp += dN_ItnItp + dN_EtpItp - dN_ItpR - dN_ItpIs;
      Is += dN_ItpIs - dN_IsIrs - dN_IsIh;
      Ih += dN_IsIh - dN_IhIrh - dN_IhIv;
      Irs += dN_IsIrs - dN_IrsR;
      Irh += dN_IhIrh + dN_IvIrh - dN_IrhR;
      Iv +=  dN_IhIv - dN_IvIrh - dN_IvD;
      R += dN_ItnR + dN_ItpR + dN_IrsR + dN_IrhR;
      D += dN_IvD;")

init <- Csnippet("
      S = N-E0_0-Etp_0-Itn_0-Itp_0-Is_0-Ih_0-Iv_0-Irs_0-Irh_0-R_0-D_0;
      E0 = E0_0;
      Etp = Etp_0;
      Itn = Itn_0;
      Itp = Itp_0;
      Is = Is_0;
      Ih = Ih_0;
      Iv = Iv_0;
      Irs = Irs_0;
      Irh = Irh_0;
      R = R_0;
      D = D_0;")

tiempos = 200

rmeas <- Csnippet("covid = Etp+Is+Itp+Ih+Iv;")
#rmeas <- Csnippet("covid = Is;")

pomp(data.frame(time=1:tiempos, covid=rep(0,tiempos)),
     rprocess= euler(sir_step, delta.t=1/8),  #euler
     times="time", t0=0, rinit=init,
     rmeasure=rmeas,
     statenames=c("S","E0","Etp","Itn","Itp","Is","Ih","Iv","Irs","Irh","R","D"),
     paramnames=c("N","Beta","sigma","alpha","gamman","gammap","etas","gammasi","etah",
                  "etav","gammahi","gammavi","E0_0","Etp_0","Itn_0","Itp_0","Is_0","Ih_0",
                  "Iv_0","Irs_0","gammas","gammah","Irh_0","R_0","D_0","delta")) -> closed.sir

numero_sims = 100

sims <- simulate(closed.sir, params = c(N = NN, Beta=par[1], sigma=par[2], alpha=par[3], gamman=par[11], 
                                        gammap=par[12],etas=par[4] , gammasi=par[5],etah=par[6], 
                                        etav=par[7], gammahi=par[8], gammavi=par[9], 
                                        E0_0=as.integer(par[13]), 
                                        Etp_0=0, Itn_0=0, Itp_0=as.integer(par[16]), Is_0=0, Ih_0=0,Iv_0=0, Irs_0=0, 
                                        gammas=par[14], gammah=par[15], Irh_0=0, R_0=0,D_0=0, 
                                        delta=par[10]), nsim=numero_sims, format="data.frame")



ggplot(sims, aes(x=time))+ggtitle("12 Compartments model")+
  labs(y="Population", x = "Time in days after 03-09-2020")+
  #geom_line(aes(y=S,group=.id,color='S'),size=2)+
  geom_line(aes(y=E0,group=.id,color='E0'),size=0.5)+
  geom_line(aes(y=Etp,group=.id,color='Etp'),size=0.5)+
  geom_line(aes(y=Itn,group=.id,color='Itn'),size=0.5)+
  geom_line(aes(y=Itp,group=.id,color='Itp'),size=0.5)+
  geom_line(aes(y=Is,group=.id,color='Is'),size=0.5)+
  geom_line(aes(y=Ih,group=.id,color='Ih'),size=0.5)+
  geom_line(aes(y=Iv,group=.id,color='Iv'),size=0.5)+
  geom_line(aes(y=Irs,group=.id,color='Irs'),size=0.5)+
  geom_line(aes(y=Irh,group=.id,color='Irh'),size=0.5)+
  #geom_line(aes(y=R,group=.id,color='R'),size=2)+
  geom_line(aes(y=D,group=.id,color='D'),size=0.5)+ #coord_trans(y="log10") +
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(size=2)))+ ggtitle("Stochastic compartments") +
  labs(y="Population", x = "Time in days after 03-26-2020")


#######################################################################################
####################### Confidence intervals of simulations ###########################

confidence_interval <- function(vector, interval) {
  # Dropping NA values
  vector <- vector[!is.na(vector)]
  # Standard deviation
  vec_sd <- sd(vector)
  # Size
  n <- length(vector)
  # Mean
  vec_mean <- mean(vector)
  
  # Error according to normal distribution
  z = qnorm(1-(1-interval)/2)
  error = z*vec_sd/sqrt(n)
  
  # Error according to t distribution
  #error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  
  # Confidence interval as a vector
  #result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(c(vec_mean - error, vec_mean, vec_mean + error))
}


ci_compart = data.frame(matrix(ncol = 5, nrow = 12*tiempos))

names(ci_compart) = c("time","compartment","lower","mean","upper")

ci_compart$compartment = rep(c("S","E0","Etp","Itn","Itp","Is","Ih","Iv","Irs","Irh","R","D"),tiempos)

repetir12 = function(ind){rep(ind,12)}

ci_compart$time = unlist(lapply(1:tiempos,repetir12))

for (compart in c("S","E0","Etp","Itn","Itp","Is","Ih","Iv","Irs","Irh","R","D")){
  for (tiempo in 1:tiempos){
    ci_compart[ci_compart$time==tiempo & ci_compart$compartment == compart, ][3:5] = 
      c( confidence_interval( sims[sims$time == tiempo, c(compart)] , 0.95 ))
  }
} 

infect = ci_compart[ci_compart$compartment %in% c("Itn","Itp","Is","Ih","Iv","Irs","Irh"),]
infect$compartment = NULL
infect = aggregate(. ~ time, infect, sum)

# ggplot(ci_compart[ci_compart$compartment %in% c("Itn","Itp","Is","Ih","Iv","Irs","Irh"),], aes(x=time))+
#   geom_line(aes(y=lower, group=compartment, color='lower'),size=0.5)+
#   geom_line(aes(y=mean, group=compartment, color='mean'),size=0.5)+
#   geom_line(aes(y=upper, group=compartment, color='upper'),size=0.5)

ggplot(infect, aes(x=time))+
  geom_line(aes(y=lower, color='lower'),size=0.5)+
  geom_line(aes(y=mean, color='mean'),size=0.5)+
  geom_line(aes(y=upper, color='upper'),size=0.5)+ ggtitle("All infectious compartments") +
  labs(y="Population", x = "Time in days after 03-26-2020")
#& ci_compart$compartment != "Etp" &
#                    ci_compart$compartment != "E0" & ci_compart$compartment != "Itp" &
#                  ci_compart$compartment != "Irs" & ci_compart$compartment != "Is" &
#                  ci_compart$compartment != "D" &

ci_compart[,c("lower", "mean", "upper")] = ci_compart[,c("lower", "mean", "upper")]/NN

ggplot(ci_compart[ci_compart$compartment != "S" & ci_compart$compartment != "R",], aes(x=time))+
  geom_line(aes(y=lower, group=compartment, color=compartment),size=0.5)+
  geom_line(aes(y=mean, group=compartment, color=compartment),size=0.5)+
  geom_line(aes(y=upper, group=compartment, color=compartment),size=0.5)+
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(size=2)))+ ggtitle("Stochastic compartments") +
  labs(y="Population", x = "Time in days after 03-26-2020")  + theme_classic() + ylim(0, 100)

ggplot(ci_compart[ci_compart$compartment != "S"  & ci_compart$compartment != "Etp" &
                    ci_compart$compartment != "E0" & ci_compart$compartment != "Itp" &
                    ci_compart$compartment != "Irs" & ci_compart$compartment != "Is" &
                    ci_compart$compartment != "D" & ci_compart$compartment != "R",], aes(x=time))+
  geom_line(aes(y=lower, group=compartment, color=compartment),size=0.5)+
  geom_line(aes(y=mean, group=compartment, color=compartment),size=0.5)+
  geom_line(aes(y=upper, group=compartment, color=compartment),size=0.5)+
  theme(legend.justification=c(1,0), legend.position=c(1,0.5))+
  theme(legend.title=element_text(size=12,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#C2C2C2',
                                size=0.25,
                                linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(size=2)))+ ggtitle("Stochastic compartments") +
  labs(y="Population", x = "Time in days after 03-26-2020") + theme_classic()



#############################################################################
########################### Computing R_0 ####################################################

para = as.data.frame(t(params1))

S_0 = NN-para$E0_0-para$Etp_0-para$Itn_0-para$Itp_0-para$Is_0-para$Ih_0-para$Iv_0-para$Irs_0-para$Irh_0-para$R_0-para$D_0

mat_F = matrix(c(0,0,rep(para$Beta*S_0/NN,7),rep(0,8*9)), nrow=9, byrow=TRUE)

mat_V = matrix(c( (1-para$alpha)*para$sigma+para$alpha, rep(0,8),-para$alpha,para$sigma,rep(0,7),-(1-para$alpha)*para$sigma,
                  0,para$alpha+para$gamman,rep(0,6),0,-para$sigma,-para$alpha,para$gammap+para$etas,rep(0,5),
                  rep(0,3),-para$etas,para$gammasi+para$etah,rep(0,4),rep(0,4),
                  -para$etah,para$gammahi+para$etav,rep(0,3),rep(0,5),-para$etav,para$gammavi+para$delta,rep(0,2),rep(0,4),
                  -para$gammasi,rep(0,2),para$gammas,0,rep(0,5),-para$gammahi,-para$gammavi,0,para$gammah),
               nrow=9, byrow=TRUE)

R_0 = max(abs(eigen(mat_F%*%solve(mat_V))$val))

R_0

mat_F%*%solve(mat_V)

para$Beta

params1
