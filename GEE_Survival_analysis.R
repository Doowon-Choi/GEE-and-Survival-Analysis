
#### load dataset
dta = as.numeric(scan("menarche.txt", what = "numeric"))
dta = data.frame(matrix(dta, ncol = 5, byrow = TRUE))
colnames(dta) = c("id", "age", "m_age", "time", "fat")

####################### Building a model over correlated response using GEE ######################

#### Draw a spaghetti plot showing fat by age timecourses ####
library(grid)
library(gridExtra)
library(ggplot2)
library(ggfortify)

dta$id <- as.factor(dta$id)
ggplot(dta, aes(x=age, y=fat, group=id,color=id)) +
  geom_line() + geom_point() +  ggtitle("fat by age")+
  theme_bw()

#### Construct variable y_pre and y_post as the average of fat by menerche status  ####

n <- length(unique(dta$id))

y_pre <- y_post <- numeric(n)
for (i in 1:n){
  y_pre[i] <- with(dta, mean(fat[id==i & time < 0]))
  y_post[i] <- with(dta, mean(fat[id==i & time >=0]))
}

ii <- (1:162)[is.na(y_post)]
y_pre <- y_pre[-ii]
y_post <- y_post[-ii]
n <- n-2
dta1 <- cbind(y_pre, y_post)

plot(y_pre, y_post)
abline(0,1)
lines(lowess(y_pre, y_post), col=2)

####  calculate correlation y_pre and y_post  ####
cor(as.numeric(dta1[,1]),as.numeric(dta1[,2]))

####  conduct paired t-test between pre and post   ####
t.test(dta1[,2], dta1[,1], paired=TRUE, alternative = c("two.sided"))

#### create a categorical variable and fit GEE model tothe new dataset. Assume an exchangeable correlation structure ####
y_all <- c(rbind(dta1[,1],dta1[,2]))
pre_post <- rep(c(0,1),160)
y_id <- c(rbind(c(1:160),c(1:160)))

dta3 <- cbind(y_all,pre_post,y_id)
names(dta3) <- c("y_all","pre_post","y_id")
dta3 <- as.data.frame(dta3)
dta3$pre_post <- as.factor(dta3$pre_post)
dta3$y_id <- as.factor(dta3$y_id)

library(gee)
fit_gee = gee(dta3$y_all ~ dta3$pre_post, id=dta3$y_id, corstr="exchangeable")
summary(fit_gee)
attributes(fit_gee)

####  a wald test for testing the null whether beta1 is zero or not. Compare it with paired t-test ####
fit_gee$coefficients[2] / sqrt(fit_gee$robust.variance[2,2])

#### estimate rho in GEE ####
fit_gee$working.correlation


#### bootstrap confidence interval for correlation rho ####
B <- 3000
boot_est <- numeric(B)
for(b in 1:B) {
  dta_boot <- NULL
  dta_boot_candi<- sample(dta3$y_id,160, replace = TRUE)
  for ( i in 1:160){
    dta_sub<- dta3[dta3$y_id==dta_boot_candi[i],]
    dta_boot <- rbind(dta_boot,dta_sub)
  }
  
  fit_gee_boot = gee(dta_boot$y_all ~ dta_boot$pre_post, id=dta_boot$y_id, corstr="exchangeable")
  rho_b <- fit_gee_boot$working.correlation[2,1]
  # est <- t(lin_com) %*% coef(fit_gls)
  boot_est[b] <- rho_b
  if (b %% 200 == 0){
    print(paste("Iteration", b))
    
  }
  remove(dta_sub,fit_gee_boot)
}

quantile(boot_est, c(0.025, 0.975))
hist(boot_est)

################### Use the heart data in the survival package ########################

#### Data preparation ####
library(survival)
Y = Age = Surg = Trans = Delta = rep(NA, 103)
for(i in 1:103) {
  ii = (1:nrow(heart))[heart$id == i]
  Y[i] = rev(heart$stop[ii])[1]
  Age[i] = rev(heart$age[ii])[1]
  Surg[i] = rev(heart$surgery[ii])[1]
  Trans[i] = rev(heart$transplant[ii])[1]
  Delta[i] = rev(heart$event[ii])[1]
}
Trans = Trans - 1

### Y: Event times
### Age: age - 48
### Surg: an indicator variable whether the patient has ahd prior bypass surgery
### Tran: an indicator variable whether the patient reeived a heart transplant
### Delta: a censoring indicator

dta21 <- cbind(Y,Age,Surg,Trans,Delta)
names(dta21) <- c("y","Age","Surg","Trans","Delta")
dta21$Surg <- factor(dta21$Surg)
dta21$Trans <- factor(dta21$Trans)
dta21 <- as.data.frame(dta21)
n <- nrow(dta21)


#### Compute Kaplan Meier survival curves for Surg=0 and Surg=1
# divide two tretment groups 0,1
y_0 <- dta21[dta21$Surg==0,1] 
y_1 <- dta21[dta21$Surg==1,1]
delta_0 <-dta21[dta21$Surg==0,5] 
delta_1 <- dta21[dta21$Surg==1,5]

Y_comb = c(y_0, y_1)
Delta_comb = c(delta_0, delta_1)
GRP = factor(c(rep(0, length(y_0)),rep(1, length(y_1))))

km_fit_0 = survfit(Surv(y_0, delta_0) ~ 1)
km_fit_1 = survfit(Surv(y_1, delta_1) ~ 1)

par(mfrow=c(1,2))
plot(km_fit_0,main="Surg=0")
plot(km_fit_1,main="Surg=1")

#### report the mean survival times by Surg status ####

### surg=0 ###
time_hat = c(0, km_fit_0$time)
surv_hat = c(1, km_fit_0$surv)

mu_hat_km_0 = 0
for(i in 2:length(surv_hat)) {
  mu_hat_km_0 = mu_hat_km_0 + surv_hat[i - 1] * (time_hat[i] - time_hat[i - 1])
}
mu_hat_km_0

### surg=1 ###
time_hat = c(0, km_fit_1$time)
surv_hat = c(1, km_fit_1$surv)

mu_hat_km_1 = 0
for(i in 2:length(surv_hat)) {
  mu_hat_km_1 = mu_hat_km_1 + surv_hat[i - 1] * (time_hat[i] - time_hat[i - 1])
}
mu_hat_km_1


#### two curves in the same plot: From this curve, we can easily see that bypass surgery is an effective care for patinets due to better survival prognosis
plot(km_fit_0$time, km_fit_0$surv, xlab = "time", ylab = expression(hat(S)), 
     xlim = range(c(y_0, y_1)), ylim = c(0, 1), type = "s", col = "blue", lwd = 2, main="Compare two groups")
lines(km_fit_1$time, km_fit_1$surv, type = "s", col = "red", lwd = 2)
legend(850, 0.95, legend = c("Suug=0", "Surg=1"), lwd = c(2, 2), col = c("blue", "red"), 
       bty = "n")

### the below code shows the same output of the above.
km_fit_surg <- survfit(Surv(Y, Delta) ~ Surg, data=dta21) #### over time
png("Two KM curves.png", width=500,height=500)
autoplot(km_fit_surg)
dev.off()

#### check Trans for each group
tran_0 <- dta21[dta21$Surg==0,4] 
tran_1 <- dta21[dta21$Surg==1,4]

mean(tran_0)
mean(tran_1)

##### Perform a logrank test whether the survival curves with Surg = 0 and Surg = 1 are equal ####
survdiff(Surv(Y_comb, Delta_comb) ~ GRP, rho=0)


############## Fit the Cox proportioanl harzards model ######################
### Survival function directly descrirbes the survival experience of a study cohort, but the hazard functio is more popular for several reasons
### say, the hazard function is a measure of instantaneous potential while a survival curve is a cumulative measure over time. Also, it may be fitted
### by a specific parametric fomr like a Weibull, lognormal curve.

# check the significance of each term
fit_1 = coxph(with(dta21,Surv(Y,Delta) ~ Age + Surg + Trans + Surg*Trans))
cox_fit <- survfit(fit_1)
summary(fit_1)


###### check the assumption of proportional hazard with respect to each term ########
### small p value means that proportional hazard is reasonable.
fit_diag = cox.zph(fit_1)
fit_diag


#### fit a Cox proportioanl hazards model and plot the predicted survivla for a age=10, surg=0, trans = 1
plot(with(dta21, survfit(fit_1, newdata=data.frame(Age=10, Surg=0, Trans=1))),xlab = "time", 
     ylab = expression(hat(S)),main="Survival curve for a given individual")


### 95% confidence interval for log the harzard ratio over the particular individual above
v = c(10, 0, 1,0)
t(v) %*% fit_1$coefficients + c(-1, 1) * 1.96 * sqrt(t(v) %*% V_beta_hat %*% v)


#### 95% confidence interval for the harzard ratio
exp(t(v) %*% fit_1$coefficients + c(-1, 1) * 1.96 * sqrt(t(v) %*% V_beta_hat %*% v))

#################### compare and plot Kaplan-Meier vs cox by ggplot2 and ggfortify ######################
km_fit <- survfit(Surv(Y, Delta) ~ 1, data=dta21) #### over time
autoplot(km_fit)

kmi <- rep("KM",length(km_fit$time))
km_df <- data.frame(km_fit$time, km_fit$surv,kmi)
names(km_df) <- c("Time","Surv","Model")

coxi <- rep("Cox",length(cox_fit$time))
cox_df <- data.frame(cox_fit$time, cox_fit$surv,coxi)
names(cox_df) <- c("Time","Surv","Model")

plot_df <- rbind(km_df,cox_df)

p <- ggplot(plot_df, aes(x = Time, y = Surv, color = Model))
p + geom_line()

png("Compare_KM_Cox.png", width=500,height=500)
p <- ggplot(plot_df, aes(x = Time, y = Surv, color = Model))
p + geom_line()
dev.off()
######################################## alen model ##############################################
#### Aalen model assumes that the cumulative hazard can be expressed as as a(t) + X B(t), where a(t) is a time-dependent intercept term, 
## X is the matrix from of covariates for the subject, and B(t) is a time-dependent coefficients. The following plot shows 
## the effects of the covariates change over time

aa_fit <-aareg(Surv(Y, Delta) ~ Age + Surg + Trans, data=dta21)
png("Alen.png", width=500,height=500)
autoplot(aa_fit)
dev.off()



########################################## Weibull regression tree for survival analysis ##################
library(partykit)

### change the scale of Age variable as a raw value
dta21$Age <- dta21$Age + 48
dta21$treat <- rnorm(n,2,2)

weibull_reg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}

### In order to calculate objective functions in MOB, loglikelihood function is defined, which is not offered in package "survival"

logLik.survreg <- function(object, ...)
  + structure(object$loglik[2], df = sum(object$df), class = "logLik")

### to fit a Weibull-regression tree, predictors and split variables are specified as MOB tree does
### Predictor: Surg (Treatment), Split variables: Age and Trans
### By fitting the regression tree, we expect that the effect of treatment on survival time depends on the condition of Age and Trans
weibull_reg_tree <- mob(Surv(Y, Delta) ~ Surg | Age + Trans, data = dta21, fit = weibull_reg, control = mob_control(minsize = 30))
plot(weibull_reg_tree)

### factor
data(veteran)
cols <- c("trt","celltype","prior")
veteran[cols] <- lapply(veteran[cols], factor)

## remove the first and second highest time 
veteran<-veteran[!(veteran$time==max(veteran$time)),]
veteran<-veteran[!(veteran$time==max(veteran$time)),]

weibull_reg_tree <- mob(Surv(time, status) ~ trt + karno | celltype + diagtime + prior + age, data = veteran, fit = weibull_reg, control = mob_control(minsize = 30))

png("Weibull_reg.png", width=500,height=700)
plot(weibull_reg_tree, terminal_panel = gbsg2node, tnex = 3)
dev.off()

#### The following code is adopted from the vignette in "partykit" package to draw plot above
#### Note that the following code is run first before drawing the plot above.
gbsg2node <- function(mobobj, 
                      col = "black", linecol = "red", cex = 0.5, pch = NULL,
                      jitter = FALSE, xscale = NULL, yscale = NULL, ylines = 1.5,
                      id = TRUE, xlab = FALSE, ylab = FALSE)
{
  ## obtain dependent variable
  mf <- model.frame(mobobj)
  y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, rhs = 0L)
  if(isTRUE(ylab)) ylab <- names(y)[1L]
  if(identical(ylab, FALSE)) ylab <- ""
  if(is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 0, 2)
  y <- y[[1L]]
  
  ## plotting character and response
  if(is.null(pch)) pch <- y[,2] * 18 + 1
  y <- y[,1]
  y <- as.numeric(y)
  pch <- rep(pch, length.out = length(y))
  if(jitter) y <- jitter(y)
  
  ## obtain explanatory variables
  x <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, rhs = 1L)
  xnam <- colnames(x)
  z <- seq(from = min(x[,2]), to = max(x[,2]), length = 51)
  z <- data.frame(a = rep(sort(x[,1])[c(1, NROW(x))], c(51, 51)), b = z)
  names(z) <- names(x)
  z$x <- model.matrix(~ ., data = z)
  
  ## fitted node ids
  fitted <- mobobj$fitted[["(fitted)"]]
  
  if(is.null(xscale)) xscale <- range(x[,2]) + c(-0.05, 0.05) * diff(range(x[,2]))
  if(is.null(yscale)) yscale <- range(y) + c(-0.05, 0.05) * diff(range(y))
  
  ## panel function for scatter plots in nodes
  rval <- function(node) {
    
    ## node index
    nid <- id_node(node)
    ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)
    
    ## dependent variable
    y <- y[ix]
    
    ## predictions
    yhat <- if(is.null(node$info$object)) {
      refit.modelparty(mobobj, node = nid)
    } else {
      node$info$object
    }
    yhat <- predict(yhat, newdata = z, type = "quantile", p = 0.5)
    pch <- pch[ix]
    
    ## viewport setup
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"), 
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_scatterplot", nid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = "white", col = 0))
    
    ## main title
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    mainlab <- paste(ifelse(id, paste("Node", nid, "(n = "), ""),
                     info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
    grid.text(mainlab)
    popViewport()
    
    plot_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, xscale = xscale,
                        yscale = yscale, name = paste("node_scatterplot", nid, "plot", sep = ""))
    pushViewport(plot_vp)
    
    ## scatterplot
    grid.points(x[ix,2], y, gp = gpar(col = col, cex = cex), pch = pch)
    grid.lines(z[1:51,2], yhat[1:51], default.units = "native", gp = gpar(col = linecol))
    grid.lines(z[52:102,2], yhat[52:102], default.units = "native", gp = gpar(col = linecol, lty = 2))
    
    grid.xaxis(at = c(ceiling(xscale[1]*10), floor(xscale[2]*10))/10)
    grid.yaxis(at = c(ceiling(yscale[1]), floor(yscale[2])))
    
    if(isTRUE(xlab)) xlab <- xnam[2]
    if(!identical(xlab, FALSE)) grid.text(xlab, x = unit(0.5, "npc"), y = unit(-2, "lines"))
    if(!identical(ylab, FALSE)) grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2, "lines"), rot = 90)
    
    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()
    
    upViewport()
  }
  
  return(rval)
}
class(gbsg2node) <- "grapcon_generator"
