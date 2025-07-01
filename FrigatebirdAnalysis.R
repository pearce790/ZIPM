## File and Data Loading
source("SourceCode.R")
data <- readxl::read_excel("Argusia database.xlsx",sheet = "Sheet1")

## Frigatebird Data Wrangling
data <- data[-1,] #remove empty row
UFB_data <- data %>% 
  filter(Species == "UFB",Year>=2007,`Site Name`=="NE Herald",Sector!="total",!str_detect(Sector,"BB")) %>%
  mutate(Time = as.numeric(paste0(Year,Month))) %>% 
  group_by(Time,Sector) %>%
  mutate(Sector =as.numeric(Sector)) %>%
  summarize(Total=sum(as.numeric(`Total (excl. old nests)`))) %>% 
  spread(key=Time,value=Total)
LFB_data <- data %>% 
  filter(Species == "LFB",Year>=2007,`Site Name`=="NE Herald",Sector!="total",!str_detect(Sector,"BB")) %>%
  mutate(Time = as.numeric(paste0(Year,Month))) %>% 
  group_by(Time,Sector) %>%
  mutate(Sector =as.numeric(Sector)) %>%
  summarize(Total=sum(as.numeric(`Total (excl. old nests)`))) %>% 
  spread(key=Time,value=Total)
GFB_data <- data %>% 
  filter(Species == "GFB",Year>=2007,`Site Name`=="NE Herald",Sector!="total",!str_detect(Sector,"BB")) %>%
  mutate(Time = as.numeric(paste0(Year,Month))) %>% 
  group_by(Time,Sector) %>%
  mutate(Sector =as.numeric(Sector)) %>%
  summarize(Total=sum(as.numeric(`Total (excl. old nests)`))) %>% 
  spread(key=Time,value=Total)
NonFB_data <- data %>% 
  filter(!(Species %in% c("UFB","LFB","GFB","Total FB")),Year>=2007,`Site Name`=="NE Herald",Sector!="total",
         !str_detect(Sector,"BB"),!str_detect(Sector,"FB")) %>%
  mutate(Time = as.numeric(paste0(Year,Month))) %>% 
  group_by(Time,Sector) %>%
  mutate(Sector =as.numeric(Sector)) %>%
  summarize(Total=sum(as.numeric(`Total (excl. old nests)`))) %>% 
  spread(key=Time,value=Total)
N <- as.matrix(UFB_data[,2:5])
I <- nrow(N)
J <- ncol(N)
t <- rep(1,I)



## Table 1
counts <- c(sum(LFB_data[,2:5]),sum(GFB_data[,2:5]),sum(UFB_data[,2:5]))
prev <- counts/sum(counts)
xtable::xtable(data.frame(Species=c("Lesser","Greater","Unidentified"),counts,prev))

## Table 2
plot_N <- as.data.frame(N)
plot_N$Site <- 1:11
plot_N <- plot_N[,c(5,1:4)]
names(plot_N)[2:5] <- c("August 2007","September 2008","October 2009","August 2012")
rownames(plot_N) <- NULL
print(xtable::xtable(plot_N,digits=0),include.rownames = FALSE)

## APPENDIX C: Analysis in which we update t
set.seed(1)
res <- EM_multistart(N,t,tol=1e-6,starts=100,update_t=T)
theta_hat <- res$theta
tau_hat <- get_se(N=N,I=I,J=J,pi=res$pi,epsilon=res$epsilon,mu=res$mu,nu=res$nu,t=res$t)
res_conf <- get_conf(theta_hat,tau_hat,I,J,alpha=0.05)

data.frame(pi=round(res$pi,2),
           epsilon=round(res$epsilon,2),
           mu=round(res$mu,2),
           nu=round(res$nu,2),
           theta=paste0(round(res$theta,2)," (",
                        round(res_conf[1],2),",",
                        round(res_conf[2],2),")"))
res$t

## Model Estimation
set.seed(1)
res <- EM_multistart(N,t,tol=1e-6,starts=1000,update_t=F)
theta_hat <- res$theta
tau_hat <- get_se(N=N,I=I,J=J,pi=res$pi,epsilon=res$epsilon,mu=res$mu,nu=res$nu,t=res$t)
res_conf <- get_conf(theta_hat,tau_hat,I,J,alpha=0.05)

## Table 3
results <- data.frame(pi=round(res$pi,2),
                      epsilon=round(res$epsilon,2),
                      mu=round(res$mu,2),
                      nu=round(res$nu,2),
                      theta=paste0(round(res$theta,2)," (",
                                   round(res_conf[1],2),",",
                                   round(res_conf[2],2),")"))
results
res$t

xtable::xtable(results)
pander::pander(results,caption="Maximum likelihood estimates of parameters",
               style="rmarkdown")
