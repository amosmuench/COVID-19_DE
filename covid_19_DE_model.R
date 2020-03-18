rm(list=ls())
library(deSolve)
library(phaseR)
##-----------------------------
## initial values and times
y <- c(S=8200000, I=6012, R=50)  # initial values, 80000000 are susceptible to the disease, 5000 already infected
parms <-c(a=0.000000019 , b=0.034, g=2.651468e-05, ns=2.507847e-05, r=0.032) # parameters 
#a = infection rate, not certain, great impact on systems stability
#b = illness mortality 3.4% according to WHO
#g =birth rate, taken from stat. bundesamt
#ns= normal death rate
# recovery, not certain
#birth rate:
(1575+585)/81464294 
#2.651468e-05

#normal death rate:
2043/81464294 
#2.507847e-05

times<-seq(0,365,1) 	 # 365 days simulation

##-----------------------------
## ODE SIR model Kermack, W; McKendrick, A (1991). "Contributions to the mathematical theory of epidemics-I". Bulletin of Mathematical Biology. 
covid_19_de<-function(times, y, parms){
  with(as.list(c(y,parms)), {
    dS <- (-a)*S*I+(g*S)-(ns*S) #ODE susceptible
    dI <- (a*S*I)-(b*I)-(ns*I)-(r*I) # ODE infeced 
    dR<- (r*I) #ODE immune, not sure whether this term is biologically proven yet
    list(c(dS, dI, dR))})}
##-----------------------------
## solve the model
out <-as.data.frame(ode(y, times, covid_19_de, parms))


##-----------------------------
## plot the model 
# base plot
plot(out$time, out$S, type="l", lwd=2)
lines(out$time, out$I, col="blue")
plot(out$time, out$I, col="blue")
# ggplot
library(ggplot2)
library(reshape2)
plot<- ggplot(melt(out, id.vars="time"), aes(x=time, y=value, group=variable)) +
  geom_line(aes(linetype=variable))+ 
  ylab("number of people")+
  xlab('t [d]')+
  theme(legend.direction = "horizontal", legend.position = "top", legend.title = element_blank(), 
        legend.background = element_blank()) 


plot

#Plot phase diagram, same ODE, given initial values
covid_19_de_phase<-function(times, y, parms){
  with(as.list(c(parms)), {
    S=y[1]
    I=y[2]
    R=y[3]
    dS <- (-a)*S*I+(g*S)-(ns*S) #ODE susceptible
    dI <- (a*S*I)-(b*I)-(ns*I)-(r*I) # ODE infeced 
    dR<- (r*I) #ODE immune, not sure whether this term is biologically proven yet
    return(list(c(dS, dI, dR)))})}

flowfield <- flowField(covid_19_de_phase, xlim=c(0,1), ylim=c(0,1), 
                       parameters = parms, points =15, 
                       xlab="Susceptible S1",  ylab="Infected S2",
                       system ="two.dim", add = FALSE)
grid()


#not working yet
growth.trajectory<- trajectory(covid_19_de_phase(), 
                               y0 = matrix(c(c(0,0.05,0.1),
                                             c(0.2,0,0.2)),
                                           nrow=3,ncol=2,byrow=F),
                               tlim=c(0,15), 
                               parameters = parms,  
                               system = "two.dim", 
                               col = c("blue","red","darkgreen"))

