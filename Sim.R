# Load source code
source("SourceCode.R")

# Set parameters and run simulation
mu <- 10
nu <- 5
results <- as.data.frame(matrix(NA,nrow=0,ncol=7))
for(Iteration in 1:1000){
  print(Iteration)
  for(I in c(20,40,80)){
    t <- rep(1,I)
    for(J in c(20,40,80)){
      for(pi in c(.1,.25,.40)){
        for(epsilon in c(.6,.7,.8)){
          set.seed(I*J*pi*epsilon*Iteration)
          N <- rZIPM(I,J,pi,epsilon,mu,nu,t)
          tryCatch({
            res <- EM_multistart(N,t,tol=1e-4,starts=20)
            theta_hat <- res$theta
            tau_hat <- get_se(N=N,I=I,J=J,pi=res$pi,epsilon=res$epsilon,mu=res$mu,nu=res$nu,t)
            results <- rbind(results,c(I,J,pi,epsilon,Iteration,theta_hat,tau_hat))}, 
            error = function(msg){print("Error!")}
          )
        }}}}
  save(results,file=paste0("Results.Rdata"))
}

# Load saved results and perform data wrangling
load("Results.Rdata")
names(results) <- c("I","J","pi","epsilon","iteration","theta_hat","tau_hat")
results <- cbind(results,t(apply(results,1,function(res){
  get_conf(res[6],res[7],res[1],res[2],alpha=0.05)
})))
names(results)[8:9] <- c("conf_lower","conf_upper")
results[which(results$theta_hat < .75),c(6,8,9)] <- 1/results[which(results$theta_hat < .75),c(6,9,8)]
results <- na.exclude(results)
results_sum <- results %>%
  group_by(I,J,pi,epsilon) %>%
  summarize(mae = mean(abs(theta_hat - 2)),
            coverage = mean(conf_lower <= 2 & conf_upper >= 2))

# Figures 1 and 2
p1<-ggplot(results_sum,aes(x=factor(pi),group=factor(epsilon),fill=factor(epsilon),y=mae))+
  geom_bar(stat="identity",position="dodge")+theme_bw(base_size=15)+
  facet_grid(vars(I),vars(J),labeller=label_both)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Greens")[-1])+
  theme(panel.grid.minor = element_blank())+
  labs(fill=expression(epsilon),x=expression(pi),y="Mean Absolute Error")
p2<-ggplot(results_sum,aes(x=factor(pi),group=factor(epsilon),fill=factor(epsilon),y=coverage))+
  theme_bw(base_size = 15)+
  geom_bar(stat="identity",position="dodge")+
  coord_cartesian(ylim=c(0.6,1))+
  scale_y_continuous(breaks=c(0.6,0.8,1.0))+
  facet_grid(vars(I),vars(J),labeller=label_both)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Greens")[-1])+
  labs(fill=expression(epsilon),x=expression(pi),y="Coverage of 95% CI")+
  geom_hline(yintercept=0.95,color="red",lty=2)

ggsave("ZIPM_Sim1.pdf",p1,units="in",width=6,height=4)
ggsave("ZIPM_Sim2.pdf",p2,units="in",width=6,height=4)

