dir()
library(dplyr)
library(survival)
library(colorspace)
library(ggfortify)
library(ggplot2)
library(xtable)
library(readxl)
library(simPH)
library(survminer)
datos <- read_excel("dataRisk.xlsx")
ind_trt <- which(datos$trt == 1)
d_trt <- datos[ind_trt,]
###
str(d_trt)
str(d_sin_trt)

##Group with treatment##
time_trt <- d_trt$time
delta_trt <- d_trt$status
km_t <- data.frame(time = time_trt, delta = delta_trt)
t_trt <- Surv(km_t$time,km_t$delta)
xfit_trt <- survfit(t_trt~1)
par(mfrow=c(1,1))
plot(xfit_trt,conf.int=TRUE)

# Survival Table
res_t <- summary(xfit_trt)
cols_t <- lapply(c(2:6, 8:10) , function(x) res[x])

# Combine the columns into a data frame
tbl_t <- do.call(data.frame, cols)
str(tbl_t)

#Generate Survival Table for LateX
print(xtable(tbl_t))


##Group without treatment##
time_sin_trt <- d_sin_trt$time
delta_sin_trt <- d_sin_trt$status
km_sin_trt <- data.frame(time=time_sin_trt,delta=delta_sin_trt)
t_sin_trt <- Surv(km_sin_trt$time,km_sin_trt$delta)
xfit_sin_trt <- survfit(t_sin_trt ~ 1)
summary(xfit_sin_trt)
par(mfrow=c(1,1))
plot(xfit_sin_trt,conf.int=TRUE)

# Survival Table
res_sin_trt <- summary(xfit_sin_trt)
cols <- lapply(c(2:6, 8:10) , function(x) res_sin_trt[x])

# Combine the columns into a data frame
tbl_sin_trt <- do.call(data.frame, cols)
str(tbl_sin_trt)

#Generate Survival Table for LateX
print(xtable(tbl_sin_trt))


##Treatment with argon laser##
d_trt_arg <- d_trt[which(d_trt$laser == "argon"),]
time_trt_arg <- d_trt_arg$time
delta_trt_arg <- d_trt_arg$status
km_t_arg <- data.frame(time = time_trt_arg, delta = delta_trt_arg)
t_trt_arg <- Surv(km_t_arg$time,km_t_arg$delta)
xfit_trt_arg <- survfit(t_trt_arg~1)
par(mfrow=c(1,1))
plot(xfit_trt_arg,conf.int=TRUE)

# Surival table
res_arg <- summary(xfit_trt_arg)
cols <- lapply(c(2:6, 8:10) , function(x) res_arg[x])

# Combine the columns into a data frame
tbl_arg <- do.call(data.frame, cols)
str(tbl_arg)

# Generate Survival Table for LateX
print(xtable(tbl_arg))


##Treatment with xenon laser##
d_trt_xen <- d_trt[which(d_trt$laser == "xenon"),]
time_trt_xen <- d_trt_xen$time
delta_trt_xen <- d_trt_xen$status
km_t_xen <- data.frame(time = time_trt_xen, delta = delta_trt_xen)
t_trt_xen <- Surv(km_t_xen$time,km_t_xen$delta)
xfit_trt_xen <- survfit(t_trt_xen~1)
par(mfrow=c(1,1))
plot(xfit_trt_xen,conf.int=TRUE)

# Extract the columns you want
res_xen <- summary(xfit_trt_xen)
cols <- lapply(c(2:6, 8:10) , function(x) res_xen[x])

# Combine the columns into a data frame
tbl_xen <- do.call(data.frame, cols)
str(tbl_xen)

#Generate Latex Table
print(xtable(tbl_xen))


##Stratified Cox Model##
datos <- read_excel("dataRisk.xlsx")
time <- datos$time
delta <- datos$status
km <- data.frame(time = time, delta = delta)
t_datos <- Surv(km$time,km$delta)
mod <- coxph(t_datos ~ laserArgon + laserXenon + strata(riskGroup), data = datos)
summary(mod)
riskf_base <- basehaz(mod)
modfit <- survfit(mod)

##Risk Curves Graphics##
#Option 1
plot(modfit, fun = "cumhaz",
     main = "Funciones de riesgo base por Grupo de Riesgo", col = c("blue","red","green","orange","black","yellow"))
#Option 2
ggplot(riskf_base,aes(x=time, y = hazard, col = strata)) +
  geom_line(size=1.5) + 
  xlab("Tiempo(Meses)") + 
  ylab("Riesgo Acumulado") +
  labs(title = "Funciones de riesgo base") +
  theme(plot.title = element_text(hjust = 0.5))

##Survial Curves Graphics##
#Option 1
simPH::ggfitStrata(modfit, byStrata = FALSE, xlab = "Tiempo(Meses)", 
      ylab = "S(t) = P(X>t)", title = "Curvas de Supervivencia por Grupo de Riesgo",rcolour = "#FFFFFF")
#Option 2
levels(modfit$strata) <- c("GR1","GR2","GR3","GR4","GR5","GR6")
g_supervivencia <- ggsurvplot(modfit, data = datos,xlab = "Tiempo(Meses)",ylab = "Probabilidad de Supervivencia")

##Kaplan Meier Estimation##
#Patients with Treatment
str(datos)
ind_trt <- which(datos$trt == 1)
datos_trt <- datos[ind_trt, ]
datos_sin_trt <- datos[-ind_trt,]
KMtrt <- survfit(Surv(datos_trt$time, datos_trt$status)~riskGroup, data = datos_trt)
KM_sin_trt <- survfit(Surv(datos_sin_trt$time, datos_sin_trt$status)~riskGroup, data = datos_sin_trt)
d_trt <- summary(KMtrt)
cols_trt <- lapply(c(2:6,8:11), function(x) d_trt[x])
tbl_trt <- do.call(data.frame,cols_trt)
levels(tbl_trt$strata) <- c("GR1","GR2","GR3","GR4","GR5","GR6")
ggsurvplot(tbl_trt, data = datos_trt, color = "strata", xlim = c(0,60), break.time.by = 10)  + labs(
  title    = "Curvas de Supervivencia KM con Tratamiento"
)

#Patients without treatment
d_sin_trt <- summary(KM_sin_trt)
cols_sin_trt <- lapply(c(2:6,8:11), function(x) d_sin_trt[x])
tbl_sin_trt <- do.call(data.frame,cols_sin_trt)
levels(tbl_sin_trt$strata) <- c("GR1","GR2","GR3","GR4","GR5","GR6")
ggsurvplot(tbl_sin_trt, data = datos_sin_trt, color = "strata", xlim = c(0,60), break.time.by = 10) + labs(
  title    = "Curvas de Supervivencia KM sin Tratamiento"
)
