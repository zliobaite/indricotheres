# 2024 02 15 I.Zliobaite
# density equations

data_traits <- read.csv('data/data_densities.csv', header = TRUE, sep = "\t")

#library(lmodel2)

log10MASS <- log10(data_traits[,'MASSKG'])
log10Dam <- log10(data_traits[,'denDam87'])
log10Tetra <- log10(data_traits[,'denTetra'])
log10Both <- log10(data_traits[,'denBoth'])

indBB <- which(data_traits[,'BB']==1)

mass1 <- 11000
mass2 <- 70
area_Helsinki <- 217

data_traits <- cbind(data_traits,log10MASS,log10Dam,log10Tetra,log10Both)

res_equations <- c()

fitTetra <- lm(log10Tetra~log10MASS,data=data_traits)
fit <- fitTetra
#fitTetraL <- lmodel2(log10Tetra~log10MASS,data=data_traits,"interval", "interval", 99)
sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred1 <- (10^intc)*(mass1^sl)
pred2 <- (10^intc)*(mass2^sl)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Santini','All',sl,intc,10^intc,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)


fitTetraBB <- lm(log10Tetra~log10MASS,data=data_traits[indBB,])
fit <- fitTetraBB
sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred1 <- (10^intc)*(mass1^sl)
pred2 <- (10^intc)*(mass2^sl)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Santini','Acute lophs',sl,intc,10^intc,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)


fitDam <- lm(log10Dam~log10MASS,data=data_traits)
fit <- fitDam
#fitDamL <- lmodel2(log10Dam~log10MASS,data=data_traits,"interval", "interval", 99)
sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred1 <- (10^intc)*(mass1^sl)
pred2 <- (10^intc)*(mass2^sl)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Damuth','All',sl,intc,10^intc,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)

fitDamBB <- lm(log10Dam~log10MASS,data=data_traits[indBB,])
fit <- fitDamBB
sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred1 <- (10^intc)*(mass1^sl)
pred2 <- (10^intc)*(mass2^sl)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Damuth','Acute lophs',sl,intc,10^intc,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)


fitBoth <- lm(log10Both~log10MASS,data=data_traits)
fit <- fitBoth
#fitBothL <- lmodel2(log10Both~log10MASS,data=data_traits,"interval", "interval", 99)
sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred1 <- (10^intc)*(mass1^sl)
pred2 <- (10^intc)*(mass2^sl)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Both','All',sl,intc,10^intc,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)


fitBothBB <- lm(log10Both~log10MASS,data=data_traits[indBB,])
fit <- fitBothBB
sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred1 <- (10^intc)*(mass1^sl)
pred2 <- (10^intc)*(mass2^sl)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Both','Acute lophs',sl,intc,10^intc,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)

ee <- 2.71828
ii <- 4.23
ss <- -0.75
pred1 <- (ee^ii)*(mass1^ss)
pred2 <- (ee^ii)*(mass2^ss)
Hki1 <- area_Helsinki*pred1
Hki2 <- area_Helsinki*pred2
res_now <- c('Dam81','published',ss,ii,ee^ii,pred1,pred2,Hki1,Hki2)
res_equations <- rbind(res_equations,res_now)

colnames(res_equations) <- c('Data source','Model fit','b','ii','a','Pred density indrik','Pred density hum','Pred indrik in Helsinki','Pred hum in Helsinki')


write.table(res_equations, file = "outputs/results_equations_estimates.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   




colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf('plots/fig_Dam.pdf',height = 5,width = 4.5)
plot(log10MASS,log10Dam,pch=1,col='darkgrey',lwd=1.5,xlab = 'log10 Mass, kg', ylab = 'log10 Density, ind/km2')
points(log10MASS[indBB],log10Dam[indBB],pch=1,col=colorBlindBlack8[2],lwd=1.5)
abline(fitDam,col='darkgrey',lwd=2)
abline(fitDamBB,col=colorBlindBlack8[2],lwd=2)
legend("bottomleft", legend=c("Damuth all, DA","Damuth ac. lop., DB"),col=c('darkgrey',colorBlindBlack8[2]), lty=1, cex=0.8,lwd=3,bty = "n")
dev.off()

pdf('plots/fig_Tetra.pdf',height = 5,width = 4.5)
plot(log10MASS,log10Tetra,pch=1,col='darkgrey',lwd=1.5,xlab = 'log10 Mass, kg', ylab = 'log10 Density, ind/km2')
points(log10MASS[indBB],log10Tetra[indBB],pch=1,col=colorBlindBlack8[3],lwd=1.5)
abline(fitTetra,col='darkgrey',lwd=2)
abline(fitTetraBB,col=colorBlindBlack8[3],lwd=2)
legend("bottomleft", legend=c("Santini all, SA","Santini ac. lop., SB"),col=c('darkgrey',colorBlindBlack8[3]), lty=1, cex=0.8,lwd=3,bty = "n")
dev.off()

pdf('plots/fig_Both.pdf',height = 5,width = 4.5)
plot(log10MASS,log10Both,pch=1,col='darkgrey',lwd=1.5,xlab = 'log10 Mass, kg', ylab = 'log10 Density, ind/km2')
points(log10MASS[indBB],log10Both[indBB],pch=1,col='#90ee90',lwd=1.5)
abline(fitBoth,col='darkgrey',lwd=2)
abline(fitBothBB,col='#90ee90',lwd=2)
legend("bottomleft", legend=c("Both all, SA","Both ac. lop., SB"),col=c('darkgrey','#90ee90'), lty=1, cex=0.8,lwd=3,bty = "n")
dev.off()

pdf('plots/fig_densities.pdf',height = 5,width = 4.5)
plot(log10MASS,log10Both,col='white',xlab = 'log10 Mass, kg', ylab = 'log10 Density, ind/km2')
points(log10MASS,log10Dam,pch=1,col=colorBlindBlack8[7],lwd=1.5)
points(log10MASS,log10Tetra,pch=1,col=colorBlindBlack8[6],lwd=1.5)
points(log10MASS,log10Both,pch=1,lwd=1.5,col='darkgrey')
abline(fitDam,col=colorBlindBlack8[7],lwd=2)
abline(fitTetra,col=colorBlindBlack8[6],lwd=2)
abline(log10(ee^ii),ss,col=colorBlindBlack8[4],lwd=2)
abline(fitBoth,col='darkgrey',lwd=2)
legend("bottomleft", legend=c("Santini, SA","Damuth, DA","Both, BA","D81"),col=c(colorBlindBlack8[6], colorBlindBlack8[7],'darkgrey',colorBlindBlack8[4]), lty=1, cex=0.8,lwd=3,bty = "n")
dev.off()

