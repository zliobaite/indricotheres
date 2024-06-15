# 2024 06 15 I.Zliobaite
# dmi estimation

data_dmi <- read.csv('data/data_dmi.csv', header = TRUE, sep = ",")

log10mass <- log10(data_dmi[,'Bmkg'])
log10dmi <- log10(data_dmi[,'DMIkgd'])

fit <- lm(log10dmi~log10mass,data=data_dmi)

sl <- fit$coefficients[2]
intc <- fit$coefficients[1]
pred <- (10^intc)*(data_dmi[,'Bmkg']^sl)

pdf('plots/fig_dmi_pred.pdf',height = 6,width = 6)
plot(log10(data_dmi[,'DMIkgd']),log10(pred))
abline(0,1)
dev.off()

pred_indrik <- (10^intc)*(11000^sl)

print('Indricothere DMI')
print(pred_indrik)