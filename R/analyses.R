library(interactions)
library(lme4)
library(lmerTest)
library(hier.part)
library(ggplot2)
library(pracma)
library(factoextra)
library(ggfortify)
library(scales)
library(tidyverse)
library(lattice)
library(imager)
library(survey)
library(gridExtra)
library(plyr)
library(jtools)
source("R/functions.R")


# color function #
t_col <- function(color, percent = 50, name = NULL) {
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col) }

### DATA MANAGEMENT ############################################################################################################
data <- read.csv("output/Complexity_sub.csv")

L <- 2 #2x2m box
scl <- L / c(1, 2, 3.125, 6.25, 12.5, 25, 50, 100, 200)
L0 <- min(scl) 

data2 <- data %>% mutate(R = 0.5*log10(R_theory_mean^2 - 1),
                         D = log10(L/L0)*(D_mean-3),
                         H = log10(H_mean/(sqrt(2)*L0))) 

#### FIGURES ###################################################################################################################

#### FIGURE 1 ####
# habitat examples
agr.pic <- load.image("figs/AGR.jpg")
pav.pic <- load.image("figs/PAV.jpg")
rob.pic <- load.image("figs/ROB.jpg")

# organize data
mdat <- data2[c("R","H")]
mdat$m <- 1

# best fit orthogonal plane 
A <- data.matrix(mdat)
#B <- data.matrix(-data2["D"])
A <- matrix(c (sum(mdat$R * mdat$R) , sum(mdat$R * mdat$H) , sum(mdat$R) , 
               sum(mdat$R * mdat$H) , sum(mdat$H * mdat$H) , sum(mdat$H) ,
               sum(mdat$R) , sum(mdat$H) , nrow(mdat)), nrow = 3, ncol = 3)
A.i <- inv(A)
B <- matrix(c (sum(mdat$R * data2["D"]) , 
               sum(mdat$H * data2["D"]) , 
               sum(data2["D"])), nrow = 3, ncol = 1)
fin <- A.i %*% B

f <- function(x, y) { z <- x * fin[1,1] + y * fin[2,1] + fin[3,1] }

# plane for plotting 
grid.lines = 350
R.pred <- seq(min(mdat$R), max(mdat$R), length=grid.lines)
H.pred <- seq(min(mdat$H), max(mdat$H), length=grid.lines)
mat <- expand.grid(R=seq(min(mdat$R), max(mdat$R), length=grid.lines), 
                   H=H.pred)

D.mat <- matrix(f(mat$R, mat$H), nrow = grid.lines, ncol = grid.lines)

# d from plane vs data d 
comp <- data.frame("TrueD" = data2["D"],
                   "PlaneD" = f(mdat$R, mdat$H))
comp$R <- comp$D - comp$PlaneD
R2 <- round(1 - ( sum((comp$D - comp$PlaneD)^2) / (sum((comp$D - mean(comp$D))^2))),3)

R_mean <- aggregate(R ~ habitat, data2, mean)
H_mean <- aggregate(H ~ habitat, data2, mean)
D_mean <- aggregate(D ~ habitat, data2, mean)

# AGR
dataAGR <- data2[data2$habitat == "AGR" ,]
AGRcollist <- c("#D55E00","#FFFFFF")
AGRcolors <- colorRampPalette(AGRcollist)(100)

# PAV
dataPAV <- data2[data2$habitat == "PAV" ,]
PAVcollist <- c("#0072B2","#FFFFFF")
PAVcolors <- colorRampPalette(PAVcollist)(100)

# ROB
dataROB <- data2[data2$habitat == "ROB" ,]
ROBcollist <- c("#009E73","#FFFFFF")
ROBcolors <- colorRampPalette(ROBcollist)(100)


# plot
png("figs/figure1_habitats.png", width = 1.5, height = 3.5, units = "in", res = 600)
par(mar=c(0, 2, 1, 1) +0.1,  family = "Times", ps = 10, mfrow = c(3,1)) 

plot(agr.pic, axes = FALSE)
title("Aggregate reef", family = "Times")
mtext("(a)", side = 3, at = -.5)
plot(pav.pic, axes = FALSE)
title("Pavement", family = "Times")
mtext("(b)", side = 3, at = -.5)
plot(rob.pic, axes = FALSE)
title("Rock & Boulder", family = "Times")
mtext("(c)", side = 3, at = -.5)

dev.off()

png("figs/figure1_all.png", width = 3, height = 3.5, units = "in", res = 600)
par(mar=c(1.2, 2, 1, .5) +0.1, family = "Times", ps = 10)

scatter3D(data2$R, data2$H, data2$D, pch = 20, cex = 1, 
          col= rev(hcl.colors(100, "grays", alpha = 0.3)), 
          xlim=c(min(data2$R)-0.25, max(data2$R)+0.25), 
          ylim=c(min(data2$H)-0.25, max(data2$H)+0.25), 
          zlim=c(min(data2$D)-0.05, max(data2$D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=data2$D), 
          theta=215, 
          phi=0,
          colkey=list(side = 2, length = 0.5, width = 1, line.clab = 1, dist=0),
          clab=expression(italic(D)))
points3D(R_mean$R - 0.2, H_mean$H + 0.2, D_mean$D,
         col = c("#009E73","#0072B2","#D55E00"), cex = 2, pch = 20,
         surf=NULL,
         colkey=FALSE,
         add=TRUE)
mtext("(d)", side = 3, at = -.5)
text3D(-0.6, 2.9, -1.2, labels=expression(italic(r)^2 == 0.974), surf = NULL, add = TRUE, family = "Times")

dev.off()

png("figs/figure1_ByHabitat.png", width = 1.5, height = 3.5, units = "in", res = 600)
par(mar=c(.6, .5, 1, .6) +0.1, family = "Times", ps = 10, mfrow = c(3,1))
scatter3D(dataAGR$R, dataAGR$H, dataAGR$D, pch = 20, cex = .8, 
          col = rev(AGRcolors),
          xlim=c(min(data2$R)-0.25, max(data2$R)+0.25), 
          ylim=c(min(data2$H)-0.25, max(data2$H)+0.25), 
          zlim=c(min(data2$D)-0.05, max(data2$D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=dataAGR$D), 
          theta=215, 
          phi=0,
          colkey=FALSE)
mtext("(e)", side = 3, at = -.5)
scatter3D(dataPAV$R, dataPAV$H, dataPAV$D, pch = 20, cex =.8, 
          col= rev(PAVcolors),
          xlim=c(min(data2$R)-0.25, max(data2$R)+0.25), 
          ylim=c(min(data2$H)-0.25, max(data2$H)+0.25), 
          zlim=c(min(data2$D)-0.05, max(data2$D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=dataPAV$D), 
          theta=215, 
          phi=0,
          colkey=FALSE)
mtext("(f)", side = 3, at = -.5)
scatter3D(dataROB$R, dataROB$H, dataROB$D, pch = 20, cex = .8, 
          col= rev(ROBcolors),
          xlim=c(min(data2$R)-0.25, max(data2$R)+0.25), 
          ylim=c(min(data2$H)-0.25, max(data2$H)+0.25), 
          zlim=c(min(data2$D)-0.05, max(data2$D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=dataROB$D), 
          theta=215, 
          phi=0,
          colkey=FALSE)
mtext("(g)", side = 3, at = -.5)

dev.off()



#### FIGURE 3 ####

# survey weights
sectors <- read.csv("data/Sectors-Strata-Areas.csv")
wtemp <- left_join(data2, sectors)

weight <- plyr::ddply(wtemp, .(SEC_NAME,DEPTH_BIN,NH), summarize, n = length(unique(site)))
weight$survey_weight <- weight$NH / weight$n

data3 <- left_join(data2, weight)

data3$conc <- c(paste0(data3$island, "_", data3$SEC_NAME, "_", data3$DEPTH_BIN, "_", data3$site))
des <- svydesign(id= ~1, strata = ~ conc, weights = ~survey_weight, data = data3) 

wave_mn <- mean(data3$wave_max_log10)
wave_sd <- sd(data3$wave_max_log10)
dep_mn  <- mean(data3$depth)
dep_sd  <- sd(data3$depth)
age_mn  <- mean(data3$age)
age_sd  <- sd(data3$age)
coral_mn <- mean(data3$CORAL)
coral_sd <- sd(data3$CORAL)
dxw_mean <- mean((data3$depth * data3$wave_max_log10))
dxw_sd <- sd((data3$depth * data3$wave_max_log10))

newdata.A <- data3[data3$habitat == "AGR",][c("depth_s","age_s","wave_max_log10_s","coral_s","habitat")]
newdata.A$depth_s <- mean(newdata.A$depth_s)
newdata.A$age_s <- mean(newdata.A$age_s)
newdata.A$wave_max_log10_s <- mean(newdata.A$wave_max_log10_s)
newdata.A$coral_s <- mean(newdata.A$coral_s)
newdata.P <- data3[data3$habitat == "PAV",][c("depth_s","age_s","wave_max_log10_s","coral_s","habitat")]
newdata.P$depth_s <- mean(newdata.P$depth_s)
newdata.P$age_s <- mean(newdata.P$age_s)
newdata.P$wave_max_log10_s <- mean(newdata.P$wave_max_log10_s)
newdata.P$coral_s <- mean(newdata.P$coral_s)
newdata.R <- data3[data3$habitat == "ROB",][c("depth_s","age_s","wave_max_log10_s","coral_s","habitat")]
newdata.R$depth_s <- mean(newdata.R$depth_s)
newdata.R$age_s <- mean(newdata.R$age_s)
newdata.R$wave_max_log10_s <- mean(newdata.R$wave_max_log10_s)
newdata.R$coral_s <- mean(newdata.R$coral_s)

# Height range #
mod.sw.H <- svyglm(H ~ 0 + habitat*depth_s + habitat*age_s + habitat*wave_max_log10_s + habitat*coral_s + 
                     habitat*depth_s*wave_max_log10_s, design=des)
summary(mod.sw.H)
summ(mod.sw.H)

# Height range & Depth 
HD.A.nd <- newdata.A
HD.A.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "AGR"]), 
                       max(data3$depth_s[data3$habitat == "AGR"]),
                       length = nrow(HD.A.nd))
pred_HD.a <- predict(mod.sw.H, newdata = HD.A.nd, type = "response", se.fit = TRUE)
pred_HD.a <- as.data.frame(pred_HD.a)
colnames(pred_HD.a) <- c("Predicted","SE")
HD.A.nd <- cbind(HD.A.nd, pred_HD.a)
HD.A.nd$Predict.lwr <- HD.A.nd$Predicted - 1.96 * HD.A.nd$SE # confidence interval upper bound
HD.A.nd$Predict.upr <- HD.A.nd$Predicted + 1.96 * HD.A.nd$SE # confidence interval lower bound
HD.P.nd <- newdata.P
HD.P.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "PAV"]), 
                       max(data3$depth_s[data3$habitat == "PAV"]),
                       length = nrow(HD.P.nd))
pred_HD.p <- predict(mod.sw.H, newdata = HD.P.nd, type = "response", se.fit = TRUE)
pred_HD.p <- as.data.frame(pred_HD.p)
colnames(pred_HD.p) <- c("Predicted","SE")
HD.P.nd <- cbind(HD.P.nd, pred_HD.p)
HD.P.nd$Predict.lwr <- HD.P.nd$Predicted - 1.96 * HD.P.nd$SE # confidence interval upper bound
HD.P.nd$Predict.upr <- HD.P.nd$Predicted + 1.96 * HD.P.nd$SE # confidence interval lower bound
HD.R.nd <- newdata.R
HD.R.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "ROB"]), 
                       max(data3$depth_s[data3$habitat == "ROB"]),
                       length = nrow(HD.R.nd))
pred_HD.r <- predict(mod.sw.H, newdata = HD.R.nd, type = "response", se.fit = TRUE)
pred_HD.r <- as.data.frame(pred_HD.r)
colnames(pred_HD.r) <- c("Predicted","SE")
HD.R.nd <- cbind(HD.R.nd, pred_HD.r)
HD.R.nd$Predict.lwr <- HD.R.nd$Predicted - 1.96 * HD.R.nd$SE # confidence interval upper bound
HD.R.nd$Predict.upr <- HD.R.nd$Predicted + 1.96 * HD.R.nd$SE # confidence interval lower bound

H.DEP <- ggplot() +
  geom_point(aes(x = data3$depth_s, y = data3$H), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=HD.A.nd$depth_s, 
                  ymin=HD.A.nd$Predict.lwr, 
                  ymax=HD.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=HD.A.nd$depth_s, 
                y=HD.A.nd$Predicted), size = .8, color = "#D55E00",
            linetype = "dashed") +
  geom_ribbon(aes(x=HD.P.nd$depth_s, 
                  ymin=HD.P.nd$Predict.lwr, 
                  ymax=HD.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=HD.P.nd$depth_s, 
                y=HD.P.nd$Predicted), size = .8, color = "#0072B2",
            linetype = "dashed") +
  geom_ribbon(aes(x=HD.R.nd$depth_s, 
                  ymin=HD.R.nd$Predict.lwr, 
                  ymax=HD.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=HD.R.nd$depth_s, 
                y=HD.R.nd$Predicted), size = .8, color = "#009E73",
            linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 0, 2, 2), "pt")) +
  labs(x = "Depth, m", y = "Height range") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (6-dep_mn)/dep_sd,
                                (12-dep_mn)/dep_sd, 
                                (18-dep_mn)/dep_sd, 
                                (24-dep_mn)/dep_sd), 
                     labels = c(0,6,12,18,24))

# Height range & Age 
HA.A.nd <- newdata.A
HA.A.nd$age_s <- seq(min(data3$age_s[data3$habitat == "AGR"]), 
                     max(data3$age_s[data3$habitat == "AGR"]),
                     length = nrow(HA.A.nd))
pred_HA.a <- predict(mod.sw.H, newdata = HA.A.nd, type = "response", se.fit = TRUE)
pred_HA.a <- as.data.frame(pred_HA.a)
colnames(pred_HA.a) <- c("Predicted","SE")
HA.A.nd <- cbind(HA.A.nd, pred_HA.a)
HA.A.nd$Predict.lwr <- HA.A.nd$Predicted - 1.96 * HA.A.nd$SE # confidence interval upper bound
HA.A.nd$Predict.upr <- HA.A.nd$Predicted + 1.96 * HA.A.nd$SE # confidence interval lower bound
HA.P.nd <- newdata.P
HA.P.nd$age_s <- seq(min(data3$age_s[data3$habitat == "PAV"]), 
                     max(data3$age_s[data3$habitat == "PAV"]),
                     length = nrow(HA.P.nd))
pred_HA.p <- predict(mod.sw.H, newdata = HA.P.nd, type = "response", se.fit = TRUE)
pred_HA.p <- as.data.frame(pred_HA.p)
colnames(pred_HA.p) <- c("Predicted","SE")
HA.P.nd <- cbind(HA.P.nd, pred_HA.p)
HA.P.nd$Predict.lwr <- HA.P.nd$Predicted - 1.96 * HA.P.nd$SE # confidence interval upper bound
HA.P.nd$Predict.upr <- HA.P.nd$Predicted + 1.96 * HA.P.nd$SE # confidence interval lower bound
HA.R.nd <- newdata.R
HA.R.nd$age_s <- seq(min(data3$age_s[data3$habitat == "ROB"]), 
                     max(data3$age_s[data3$habitat == "ROB"]),
                     length = nrow(HA.R.nd))
pred_HA.r <- predict(mod.sw.H, newdata = HA.R.nd, type = "response", se.fit = TRUE)
pred_HA.r <- as.data.frame(pred_HA.r)
colnames(pred_HA.r) <- c("Predicted","SE")
HA.R.nd <- cbind(HA.R.nd, pred_HA.r)
HA.R.nd$Predict.lwr <- HA.R.nd$Predicted - 1.96 * HA.R.nd$SE # confidence interval upper bound
HA.R.nd$Predict.upr <- HA.R.nd$Predicted + 1.96 * HA.R.nd$SE # confidence interval lower bound

H.AGE <- ggplot() +
  geom_point(aes(x = data3$age_s, y = data3$H), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=HA.A.nd$age_s,
                  ymin=HA.A.nd$Predict.lwr,
                  ymax=HA.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=HA.A.nd$age_s,
                y=HA.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=HA.P.nd$age_s, 
                  ymin=HA.P.nd$Predict.lwr, 
                  ymax=HA.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=HA.P.nd$age_s, 
                y=HA.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=HA.R.nd$age_s,
                  ymin=HA.R.nd$Predict.lwr,
                  ymax=HA.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=HA.R.nd$age_s,
                y=HA.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 6, 2, 0), "pt")) +
  labs(x = "Geological Reef Age, mya", y = "") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-age_mn)/age_sd, 
                                (1-age_mn)/age_sd, 
                                (2-age_mn)/age_sd, 
                                (3-age_mn)/age_sd, 
                                (4-age_mn)/age_sd, 
                                (5-age_mn)/age_sd),
                     labels = c(0,1,2,3,4,5))

# Height range & Wave 
HW.A.nd <- newdata.A
HW.A.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "AGR"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "AGR"]),
                                length = nrow(HW.A.nd))
pred_HW.a <- predict(mod.sw.H, newdata = HW.A.nd, type = "response", se.fit = TRUE)
pred_HW.a <- as.data.frame(pred_HW.a)
colnames(pred_HW.a) <- c("Predicted","SE")
HW.A.nd <- cbind(HW.A.nd, pred_HW.a)
HW.A.nd$Predict.lwr <- HW.A.nd$Predicted - 1.96 * HW.A.nd$SE # confidence interval upper bound
HW.A.nd$Predict.upr <- HW.A.nd$Predicted + 1.96 * HW.A.nd$SE # confidence interval lower bound
HW.P.nd <- newdata.P
HW.P.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "PAV"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "PAV"]),
                                length = nrow(HW.P.nd))
pred_HW.p <- predict(mod.sw.H, newdata = HW.P.nd, type = "response", se.fit = TRUE)
pred_HW.p <- as.data.frame(pred_HW.p)
colnames(pred_HW.p) <- c("Predicted","SE")
HW.P.nd <- cbind(HW.P.nd, pred_HW.p)
HW.P.nd$Predict.lwr <- HW.P.nd$Predicted - 1.96 * HW.P.nd$SE # confidence interval upper bound
HW.P.nd$Predict.upr <- HW.P.nd$Predicted + 1.96 * HW.P.nd$SE # confidence interval lower bound
HW.R.nd <- newdata.R
HW.R.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "ROB"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "ROB"]),
                                length = nrow(HW.R.nd))
pred_HW.r <- predict(mod.sw.H, newdata = HW.R.nd, type = "response", se.fit = TRUE)
pred_HW.r <- as.data.frame(pred_HW.r)
colnames(pred_HW.r) <- c("Predicted","SE")
HW.R.nd <- cbind(HW.R.nd, pred_HW.r)
HW.R.nd$Predict.lwr <- HW.R.nd$Predicted - 1.96 * HW.R.nd$SE # confidence interval upper bound
HW.R.nd$Predict.upr <- HW.R.nd$Predicted + 1.96 * HW.R.nd$SE # confidence interval lower bound


H.WAVE <- ggplot() +
  geom_point(aes(x = data3$wave_max_log10_s, y = data3$H), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=HW.A.nd$wave_max_log10_s,
                  ymin=HW.A.nd$Predict.lwr,
                  ymax=HW.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=HW.A.nd$wave_max_log10_s,
                y=HW.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=HW.P.nd$wave_max_log10_s, 
                  ymin=HW.P.nd$Predict.lwr, 
                  ymax=HW.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=HW.P.nd$wave_max_log10_s, 
                y=HW.P.nd$Predicted), size = .8, color = "#0072B2",
            linetype = "dashed") +
  geom_ribbon(aes(x=HW.R.nd$wave_max_log10_s,
                  ymin=HW.R.nd$Predict.lwr,
                  ymax=HW.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=HW.R.nd$wave_max_log10_s,
                y=HW.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 0, 2, 0), "pt")) +
  labs(x = "Exposure, kW/m (log10)", y = "") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) 

# Height range & Coral
HC.A.nd <- newdata.A
HC.A.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "AGR"]), 
                       max(data3$coral_s[data3$habitat == "AGR"]),
                       length = nrow(HC.A.nd))
pred_HC.a <- predict(mod.sw.H, newdata = HC.A.nd, type = "response", se.fit = TRUE)
pred_HC.a <- as.data.frame(pred_HC.a)
colnames(pred_HC.a) <- c("Predicted","SE")
HC.A.nd <- cbind(HC.A.nd, pred_HC.a)
HC.A.nd$Predict.lwr <- HC.A.nd$Predicted - 1.96 * HC.A.nd$SE # confidence interval upper bound
HC.A.nd$Predict.upr <- HC.A.nd$Predicted + 1.96 * HC.A.nd$SE # confidence interval lower bound
HC.P.nd <- newdata.P
HC.P.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "PAV"]), 
                       max(data3$coral_s[data3$habitat == "PAV"]),
                       length = nrow(HC.P.nd))
pred_HC.p <- predict(mod.sw.H, newdata = HC.P.nd, type = "response", se.fit = TRUE)
pred_HC.p <- as.data.frame(pred_HC.p)
colnames(pred_HC.p) <- c("Predicted","SE")
HC.P.nd <- cbind(HC.P.nd, pred_HC.p)
HC.P.nd$Predict.lwr <- HC.P.nd$Predicted - 1.96 * HC.P.nd$SE # confidence interval upper bound
HC.P.nd$Predict.upr <- HC.P.nd$Predicted + 1.96 * HC.P.nd$SE # confidence interval lower bound
HC.R.nd <- newdata.R
HC.R.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "ROB"]), 
                       max(data3$coral_s[data3$habitat == "ROB"]),
                       length = nrow(HC.R.nd))
pred_HC.r <- predict(mod.sw.H, newdata = HC.R.nd, type = "response", se.fit = TRUE)
pred_HC.r <- as.data.frame(pred_HC.r)
colnames(pred_HC.r) <- c("Predicted","SE")
HC.R.nd <- cbind(HC.R.nd, pred_HC.r)
HC.R.nd$Predict.lwr <- HC.R.nd$Predicted - 1.96 * HC.R.nd$SE # confidence interval upper bound
HC.R.nd$Predict.upr <- HC.R.nd$Predicted + 1.96 * HC.R.nd$SE # confidence interval lower bound

H.COR <- ggplot() +
  geom_point(aes(x = data3$coral_s, y = data3$H), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=HC.A.nd$coral_s,
                  ymin=HC.A.nd$Predict.lwr,
                  ymax=HC.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=HC.A.nd$coral_s,
                y=HC.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=HC.P.nd$coral_s, 
                  ymin=HC.P.nd$Predict.lwr, 
                  ymax=HC.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=HC.P.nd$coral_s, 
                y=HC.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=HC.R.nd$coral_s,
                  ymin=HC.R.nd$Predict.lwr,
                  ymax=HC.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=HC.R.nd$coral_s,
                y=HC.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 2, 2, 0), "pt")) +
  labs(x = "Coral Cover, %", y = "") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))

# Rugosity #
mod.sw.R <-svyglm(R ~ 0 + habitat*depth_s + habitat*age_s + habitat*wave_max_log10_s + habitat*coral_s + 
                    habitat*depth_s*wave_max_log10_s, design=des)
summary(mod.sw.R)
summ(mod.sw.R)

# Rugosity & Depth 
RD.A.nd <- newdata.A
RD.A.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "AGR"]), 
                       max(data3$depth_s[data3$habitat == "AGR"]),
                       length = nrow(RD.A.nd))
pred_RD.a <- predict(mod.sw.R, newdata = RD.A.nd, type = "response", se.fit = TRUE)
pred_RD.a <- as.data.frame(pred_RD.a)
colnames(pred_RD.a) <- c("Predicted","SE")
RD.A.nd <- cbind(RD.A.nd, pred_RD.a)
RD.A.nd$Predict.lwr <- RD.A.nd$Predicted - 1.96 * RD.A.nd$SE # confidence interval upper bound
RD.A.nd$Predict.upr <- RD.A.nd$Predicted + 1.96 * RD.A.nd$SE # confidence interval lower bound
RD.P.nd <- newdata.P
RD.P.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "PAV"]), 
                       max(data3$depth_s[data3$habitat == "PAV"]),
                       length = nrow(RD.P.nd))
pred_RD.p <- predict(mod.sw.R, newdata = RD.P.nd, type = "response", se.fit = TRUE)
pred_RD.p <- as.data.frame(pred_RD.p)
colnames(pred_RD.p) <- c("Predicted","SE")
RD.P.nd <- cbind(RD.P.nd, pred_RD.p)
RD.P.nd$Predict.lwr <- RD.P.nd$Predicted - 1.96 * RD.P.nd$SE # confidence interval upper bound
RD.P.nd$Predict.upr <- RD.P.nd$Predicted + 1.96 * RD.P.nd$SE # confidence interval lower bound
RD.R.nd <- newdata.R
RD.R.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "ROB"]), 
                       max(data3$depth_s[data3$habitat == "ROB"]),
                       length = nrow(RD.R.nd))
pred_RD.r <- predict(mod.sw.R, newdata = RD.R.nd, type = "response", se.fit = TRUE)
pred_RD.r <- as.data.frame(pred_RD.r)
colnames(pred_RD.r) <- c("Predicted","SE")
RD.R.nd <- cbind(RD.R.nd, pred_RD.r)
RD.R.nd$Predict.lwr <- RD.R.nd$Predicted - 1.96 * RD.R.nd$SE # confidence interval upper bound
RD.R.nd$Predict.upr <- RD.R.nd$Predicted + 1.96 * RD.R.nd$SE # confidence interval lower bound

R.DEP <- ggplot() +
  geom_point(aes(x = data3$depth_s, y = data3$R), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=RD.A.nd$depth_s,
                  ymin=RD.A.nd$Predict.lwr,
                  ymax=RD.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=RD.A.nd$depth_s,
                y=RD.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=RD.P.nd$depth_s, 
                  ymin=RD.P.nd$Predict.lwr, 
                  ymax=RD.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=RD.P.nd$depth_s, 
                y=RD.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=RD.R.nd$depth_s,
                  ymin=RD.R.nd$Predict.lwr,
                  ymax=RD.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=RD.R.nd$depth_s,
                y=RD.R.nd$Predicted), size = .8, color = "#009E73",
            linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(2, 0, 0, 2), "pt")) +
  labs(x = "", y = "Rugosity") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (6-dep_mn)/dep_sd,
                                (12-dep_mn)/dep_sd, 
                                (18-dep_mn)/dep_sd, 
                                (24-dep_mn)/dep_sd), 
                     labels = c(0,6,12,18,24))

# Rugosity & Age 
RA.A.nd <- newdata.A
RA.A.nd$age_s <- seq(min(data3$age_s[data3$habitat == "AGR"]), 
                     max(data3$age_s[data3$habitat == "AGR"]),
                     length = nrow(RA.A.nd))
pred_RA.a <- predict(mod.sw.R, newdata = RA.A.nd, type = "response", se.fit = TRUE)
pred_RA.a <- as.data.frame(pred_RA.a)
colnames(pred_RA.a) <- c("Predicted","SE")
RA.A.nd <- cbind(RA.A.nd, pred_RA.a)
RA.A.nd$Predict.lwr <- RA.A.nd$Predicted - 1.96 * RA.A.nd$SE # confidence interval upper bound
RA.A.nd$Predict.upr <- RA.A.nd$Predicted + 1.96 * RA.A.nd$SE # confidence interval lower bound
RA.P.nd <- newdata.P
RA.P.nd$age_s <- seq(min(data3$age_s[data3$habitat == "PAV"]), 
                     max(data3$age_s[data3$habitat == "PAV"]),
                     length = nrow(RA.P.nd))
pred_RA.p <- predict(mod.sw.R, newdata = RA.P.nd, type = "response", se.fit = TRUE)
pred_RA.p <- as.data.frame(pred_RA.p)
colnames(pred_RA.p) <- c("Predicted","SE")
RA.P.nd <- cbind(RA.P.nd, pred_RA.p)
RA.P.nd$Predict.lwr <- RA.P.nd$Predicted - 1.96 * RA.P.nd$SE # confidence interval upper bound
RA.P.nd$Predict.upr <- RA.P.nd$Predicted + 1.96 * RA.P.nd$SE # confidence interval lower bound
RA.R.nd <- newdata.R
RA.R.nd$age_s <- seq(min(data3$age_s[data3$habitat == "ROB"]), 
                     max(data3$age_s[data3$habitat == "ROB"]),
                     length = nrow(RA.R.nd))
pred_RA.r <- predict(mod.sw.R, newdata = RA.R.nd, type = "response", se.fit = TRUE)
pred_RA.r <- as.data.frame(pred_RA.r)
colnames(pred_RA.r) <- c("Predicted","SE")
RA.R.nd <- cbind(RA.R.nd, pred_RA.r)
RA.R.nd$Predict.lwr <- RA.R.nd$Predicted - 1.96 * RA.R.nd$SE # confidence interval upper bound
RA.R.nd$Predict.upr <- RA.R.nd$Predicted + 1.96 * RA.R.nd$SE # confidence interval lower bound

R.AGE <- ggplot() +
  geom_point(aes(x = data3$age_s, y = data3$R), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=RA.A.nd$age_s,
                  ymin=RA.A.nd$Predict.lwr,
                  ymax=RA.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=RA.A.nd$age_s,
                y=RA.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=RA.P.nd$age_s, 
                  ymin=RA.P.nd$Predict.lwr, 
                  ymax=RA.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=RA.P.nd$age_s, 
                y=RA.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=RA.R.nd$age_s,
                  ymin=RA.R.nd$Predict.lwr,
                  ymax=RA.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=RA.R.nd$age_s,
                y=RA.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(2, 6, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-age_mn)/age_sd, 
                                (1-age_mn)/age_sd, 
                                (2-age_mn)/age_sd, 
                                (3-age_mn)/age_sd, 
                                (4-age_mn)/age_sd, 
                                (5-age_mn)/age_sd),
                     labels = c(0,1,2,3,4,5))

# Rugosity & Wave 
RW.A.nd <- newdata.A
RW.A.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "AGR"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "AGR"]),
                                length = nrow(RW.A.nd))
pred_RW.a <- predict(mod.sw.R, newdata = RW.A.nd, type = "response", se.fit = TRUE)
pred_RW.a <- as.data.frame(pred_RW.a)
colnames(pred_RW.a) <- c("Predicted","SE")
RW.A.nd <- cbind(RW.A.nd, pred_RW.a)
RW.A.nd$Predict.lwr <- RW.A.nd$Predicted - 1.96 * RW.A.nd$SE # confidence interval upper bound
RW.A.nd$Predict.upr <- RW.A.nd$Predicted + 1.96 * RW.A.nd$SE # confidence interval lower bound
RW.P.nd <- newdata.P
RW.P.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "PAV"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "PAV"]),
                                length = nrow(RW.P.nd))
pred_RW.p <- predict(mod.sw.R, newdata = RW.P.nd, type = "response", se.fit = TRUE)
pred_RW.p <- as.data.frame(pred_RW.p)
colnames(pred_RW.p) <- c("Predicted","SE")
RW.P.nd <- cbind(RW.P.nd, pred_RW.p)
RW.P.nd$Predict.lwr <- RW.P.nd$Predicted - 1.96 * RW.P.nd$SE # confidence interval upper bound
RW.P.nd$Predict.upr <- RW.P.nd$Predicted + 1.96 * RW.P.nd$SE # confidence interval lower bound
RW.R.nd <- newdata.R
RW.R.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "ROB"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "ROB"]),
                                length = nrow(RW.R.nd))
pred_RW.r <- predict(mod.sw.R, newdata = RW.R.nd, type = "response", se.fit = TRUE)
pred_RW.r <- as.data.frame(pred_RW.r)
colnames(pred_RW.r) <- c("Predicted","SE")
RW.R.nd <- cbind(RW.R.nd, pred_RW.r)
RW.R.nd$Predict.lwr <- RW.R.nd$Predicted - 1.96 * RW.R.nd$SE # confidence interval upper bound
RW.R.nd$Predict.upr <- RW.R.nd$Predicted + 1.96 * RW.R.nd$SE # confidence interval lower bound

R.WAVE <- ggplot() +
  geom_point(aes(x = data3$wave_max_log10_s, y = data3$R), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=RW.A.nd$wave_max_log10_s,
                  ymin=RW.A.nd$Predict.lwr,
                  ymax=RW.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=RW.A.nd$wave_max_log10_s,
                y=RW.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=RW.P.nd$wave_max_log10_s, 
                  ymin=RW.P.nd$Predict.lwr, 
                  ymax=RW.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=RW.P.nd$wave_max_log10_s, 
                y=RW.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=RW.R.nd$wave_max_log10_s,
                  ymin=RW.R.nd$Predict.lwr,
                  ymax=RW.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=RW.R.nd$wave_max_log10_s,
                y=RW.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(2, 0, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) 

# Rugosity & Coral
RC.A.nd <- newdata.A
RC.A.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "AGR"]), 
                       max(data3$coral_s[data3$habitat == "AGR"]),
                       length = nrow(RC.A.nd))
pred_RC.a <- predict(mod.sw.R, newdata = RC.A.nd, type = "response", se.fit = TRUE)
pred_RC.a <- as.data.frame(pred_RC.a)
colnames(pred_RC.a) <- c("Predicted","SE")
RC.A.nd <- cbind(RC.A.nd, pred_RC.a)
RC.A.nd$Predict.lwr <- RC.A.nd$Predicted - 1.96 * RC.A.nd$SE # confidence interval upper bound
RC.A.nd$Predict.upr <- RC.A.nd$Predicted + 1.96 * RC.A.nd$SE # confidence interval lower bound
RC.P.nd <- newdata.P
RC.P.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "PAV"]), 
                       max(data3$coral_s[data3$habitat == "PAV"]),
                       length = nrow(RC.P.nd))
pred_RC.p <- predict(mod.sw.R, newdata = RC.P.nd, type = "response", se.fit = TRUE)
pred_RC.p <- as.data.frame(pred_RC.p)
colnames(pred_RC.p) <- c("Predicted","SE")
RC.P.nd <- cbind(RC.P.nd, pred_RC.p)
RC.P.nd$Predict.lwr <- RC.P.nd$Predicted - 1.96 * RC.P.nd$SE # confidence interval upper bound
RC.P.nd$Predict.upr <- RC.P.nd$Predicted + 1.96 * RC.P.nd$SE # confidence interval lower bound
RC.R.nd <- newdata.R
RC.R.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "ROB"]), 
                       max(data3$coral_s[data3$habitat == "ROB"]),
                       length = nrow(RC.R.nd))
pred_RC.r <- predict(mod.sw.R, newdata = RC.R.nd, type = "response", se.fit = TRUE)
pred_RC.r <- as.data.frame(pred_RC.r)
colnames(pred_RC.r) <- c("Predicted","SE")
RC.R.nd <- cbind(RC.R.nd, pred_RC.r)
RC.R.nd$Predict.lwr <- RC.R.nd$Predicted - 1.96 * RC.R.nd$SE # confidence interval upper bound
RC.R.nd$Predict.upr <- RC.R.nd$Predicted + 1.96 * RC.R.nd$SE # confidence interval lower bound

R.COR <- ggplot() +
  geom_point(aes(x = data3$coral_s, y = data3$R), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=RC.A.nd$coral_s,
                  ymin=RC.A.nd$Predict.lwr,
                  ymax=RC.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=RC.A.nd$coral_s,
                y=RC.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=RC.P.nd$coral_s, 
                  ymin=RC.P.nd$Predict.lwr, 
                  ymax=RC.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=RC.P.nd$coral_s, 
                y=RC.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=RC.R.nd$coral_s,
                  ymin=RC.R.nd$Predict.lwr,
                  ymax=RC.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=RC.R.nd$coral_s,
                y=RC.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(2, 2, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))

# Fractal Dimension #
mod.sw.D <-svyglm(D ~ 0 + habitat*depth_s + habitat*age_s + habitat*wave_max_log10_s + 
                    habitat*coral_s + habitat*depth_s*wave_max_log10_s, design=des)
summary(mod.sw.D)
summ(mod.sw.D)

# Fractal dimension & Depth 
DD.A.nd <- newdata.A
DD.A.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "AGR"]), 
                       max(data3$depth_s[data3$habitat == "AGR"]),
                       length = nrow(DD.A.nd))
pred_DD.a <- predict(mod.sw.D, newdata = DD.A.nd, type = "response", se.fit = TRUE)
pred_DD.a <- as.data.frame(pred_DD.a)
colnames(pred_DD.a) <- c("Predicted","SE")
DD.A.nd <- cbind(DD.A.nd, pred_DD.a)
DD.A.nd$Predict.lwr <- DD.A.nd$Predicted - 1.96 * DD.A.nd$SE # confidence interval upper bound
DD.A.nd$Predict.upr <- DD.A.nd$Predicted + 1.96 * DD.A.nd$SE # confidence interval lower bound
DD.P.nd <- newdata.P
DD.P.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "PAV"]), 
                       max(data3$depth_s[data3$habitat == "PAV"]),
                       length = nrow(DD.P.nd))
pred_DD.p <- predict(mod.sw.D, newdata = DD.P.nd, type = "response", se.fit = TRUE)
pred_DD.p <- as.data.frame(pred_DD.p)
colnames(pred_DD.p) <- c("Predicted","SE")
DD.P.nd <- cbind(DD.P.nd, pred_DD.p)
DD.P.nd$Predict.lwr <- DD.P.nd$Predicted - 1.96 * DD.P.nd$SE # confidence interval upper bound
DD.P.nd$Predict.upr <- DD.P.nd$Predicted + 1.96 * DD.P.nd$SE # confidence interval lower bound
DD.R.nd <- newdata.R
DD.R.nd$depth_s <- seq(min(data3$depth_s[data3$habitat == "ROB"]), 
                       max(data3$depth_s[data3$habitat == "ROB"]),
                       length = nrow(DD.R.nd))
pred_DD.r <- predict(mod.sw.D, newdata = DD.R.nd, type = "response", se.fit = TRUE)
pred_DD.r <- as.data.frame(pred_DD.r)
colnames(pred_DD.r) <- c("Predicted","SE")
DD.R.nd <- cbind(DD.R.nd, pred_DD.r)
DD.R.nd$Predict.lwr <- DD.R.nd$Predicted - 1.96 * DD.R.nd$SE # confidence interval upper bound
DD.R.nd$Predict.upr <- DD.R.nd$Predicted + 1.96 * DD.R.nd$SE # confidence interval lower bound

D.DEP <- ggplot() +
  geom_point(aes(x = data3$depth_s, y = data3$D), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=DD.A.nd$depth_s,
                  ymin=DD.A.nd$Predict.lwr,
                  ymax=DD.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=DD.A.nd$depth_s,
                y=DD.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=DD.P.nd$depth_s, 
                  ymin=DD.P.nd$Predict.lwr, 
                  ymax=DD.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=DD.P.nd$depth_s, 
                y=DD.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=DD.R.nd$depth_s,
                  ymin=DD.R.nd$Predict.lwr,
                  ymax=DD.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=DD.R.nd$depth_s,
                y=DD.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 2), "pt")) +
  labs(x = "", y = "Fractal dimension") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (6-dep_mn)/dep_sd,
                                (12-dep_mn)/dep_sd, 
                                (18-dep_mn)/dep_sd, 
                                (24-dep_mn)/dep_sd), 
                     labels = c(0,6,12,18,24))

# Fractal dimension & Age 
DA.A.nd <- newdata.A
DA.A.nd$age_s <- seq(min(data3$age_s[data3$habitat == "AGR"]), 
                     max(data3$age_s[data3$habitat == "AGR"]),
                     length = nrow(DA.A.nd))
pred_DA.a <- predict(mod.sw.D, newdata = DA.A.nd, type = "response", se.fit = TRUE)
pred_DA.a <- as.data.frame(pred_DA.a)
colnames(pred_DA.a) <- c("Predicted","SE")
DA.A.nd <- cbind(DA.A.nd, pred_DA.a)
DA.A.nd$Predict.lwr <- DA.A.nd$Predicted - 1.96 * DA.A.nd$SE # confidence interval upper bound
DA.A.nd$Predict.upr <- DA.A.nd$Predicted + 1.96 * DA.A.nd$SE # confidence interval lower bound
DA.P.nd <- newdata.P
DA.P.nd$age_s <- seq(min(data3$age_s[data3$habitat == "PAV"]), 
                     max(data3$age_s[data3$habitat == "PAV"]),
                     length = nrow(DA.P.nd))
pred_DA.p <- predict(mod.sw.D, newdata = DA.P.nd, type = "response", se.fit = TRUE)
pred_DA.p <- as.data.frame(pred_DA.p)
colnames(pred_DA.p) <- c("Predicted","SE")
DA.P.nd <- cbind(DA.P.nd, pred_DA.p)
DA.P.nd$Predict.lwr <- DA.P.nd$Predicted - 1.96 * DA.P.nd$SE # confidence interval upper bound
DA.P.nd$Predict.upr <- DA.P.nd$Predicted + 1.96 * DA.P.nd$SE # confidence interval lower bound
DA.R.nd <- newdata.R
DA.R.nd$age_s <- seq(min(data3$age_s[data3$habitat == "ROB"]), 
                     max(data3$age_s[data3$habitat == "ROB"]),
                     length = nrow(DA.R.nd))
pred_DA.r <- predict(mod.sw.D, newdata = DA.R.nd, type = "response", se.fit = TRUE)
pred_DA.r <- as.data.frame(pred_DA.r)
colnames(pred_DA.r) <- c("Predicted","SE")
DA.R.nd <- cbind(DA.R.nd, pred_DA.r)
DA.R.nd$Predict.lwr <- DA.R.nd$Predicted - 1.96 * DA.R.nd$SE # confidence interval upper bound
DA.R.nd$Predict.upr <- DA.R.nd$Predicted + 1.96 * DA.R.nd$SE # confidence interval lower bound

D.AGE <- ggplot() +
  geom_point(aes(x = data3$age_s, y = data3$D), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=DA.A.nd$age_s,
                  ymin=DA.A.nd$Predict.lwr,
                  ymax=DA.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=DA.A.nd$age_s,
                y=DA.A.nd$Predicted), size = .8, color = "#D55E00",
            linetype = "dashed") +
  geom_ribbon(aes(x=DA.P.nd$age_s, 
                  ymin=DA.P.nd$Predict.lwr, 
                  ymax=DA.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=DA.P.nd$age_s, 
                y=DA.P.nd$Predicted), size = .8, color = "#0072B2",
            linetype = "dashed") +
  geom_ribbon(aes(x=DA.R.nd$age_s,
                  ymin=DA.R.nd$Predict.lwr,
                  ymax=DA.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=DA.R.nd$age_s,
                y=DA.R.nd$Predicted), size = .8, color = "#009E73",
            linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 6, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-age_mn)/age_sd, 
                                (1-age_mn)/age_sd, 
                                (2-age_mn)/age_sd, 
                                (3-age_mn)/age_sd, 
                                (4-age_mn)/age_sd, 
                                (5-age_mn)/age_sd),
                     labels = c(0,1,2,3,4,5))

# Fractal dimension & Wave 
DW.A.nd <- newdata.A
DW.A.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "AGR"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "AGR"]),
                                length = nrow(DW.A.nd))
pred_DW.a <- predict(mod.sw.D, newdata = DW.A.nd, type = "response", se.fit = TRUE)
pred_DW.a <- as.data.frame(pred_DW.a)
colnames(pred_DW.a) <- c("Predicted","SE")
DW.A.nd <- cbind(DW.A.nd, pred_DW.a)
DW.A.nd$Predict.lwr <- DW.A.nd$Predicted - 1.96 * DW.A.nd$SE # confidence interval upper bound
DW.A.nd$Predict.upr <- DW.A.nd$Predicted + 1.96 * DW.A.nd$SE # confidence interval lower bound
DW.P.nd <- newdata.P
DW.P.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "PAV"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "PAV"]),
                                length = nrow(DW.P.nd))
pred_DW.p <- predict(mod.sw.D, newdata = DW.P.nd, type = "response", se.fit = TRUE)
pred_DW.p <- as.data.frame(pred_DW.p)
colnames(pred_DW.p) <- c("Predicted","SE")
DW.P.nd <- cbind(DW.P.nd, pred_DW.p)
DW.P.nd$Predict.lwr <- DW.P.nd$Predicted - 1.96 * DW.P.nd$SE # confidence interval upper bound
DW.P.nd$Predict.upr <- DW.P.nd$Predicted + 1.96 * DW.P.nd$SE # confidence interval lower bound
DW.R.nd <- newdata.R
DW.R.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s[data3$habitat == "ROB"]), 
                                max(data3$wave_max_log10_s[data3$habitat == "ROB"]),
                                length = nrow(DW.R.nd))
pred_DW.r <- predict(mod.sw.D, newdata = DW.R.nd, type = "response", se.fit = TRUE)
pred_DW.r <- as.data.frame(pred_DW.r)
colnames(pred_DW.r) <- c("Predicted","SE")
DW.R.nd <- cbind(DW.R.nd, pred_DW.r)
DW.R.nd$Predict.lwr <- DW.R.nd$Predicted - 1.96 * DW.R.nd$SE # confidence interval upper bound
DW.R.nd$Predict.upr <- DW.R.nd$Predicted + 1.96 * DW.R.nd$SE # confidence interval lower bound

D.WAVE <- ggplot() +
  geom_point(aes(x = data3$wave_max_log10_s, y = data3$D), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=DW.A.nd$wave_max_log10_s,
                  ymin=DW.A.nd$Predict.lwr,
                  ymax=DW.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=DW.A.nd$wave_max_log10_s,
                y=DW.A.nd$Predicted), size = .8, color = "#D55E00",
            linetype = "dashed") +
  geom_ribbon(aes(x=DW.P.nd$wave_max_log10_s, 
                  ymin=DW.P.nd$Predict.lwr, 
                  ymax=DW.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=DW.P.nd$wave_max_log10_s, 
                y=DW.P.nd$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=DW.R.nd$wave_max_log10_s,
                  ymin=DW.R.nd$Predict.lwr,
                  ymax=DW.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=DW.R.nd$wave_max_log10_s,
                y=DW.R.nd$Predicted), size = .8, color = "#009E73",
            linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) 

# Fractal dimension & Coral
DC.A.nd <- newdata.A
DC.A.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "AGR"]), 
                       max(data3$coral_s[data3$habitat == "AGR"]),
                       length = nrow(DC.A.nd))
pred_DC.a <- predict(mod.sw.D, newdata = DC.A.nd, type = "response", se.fit = TRUE)
pred_DC.a <- as.data.frame(pred_DC.a)
colnames(pred_DC.a) <- c("Predicted","SE")
DC.A.nd <- cbind(DC.A.nd, pred_DC.a)
DC.A.nd$Predict.lwr <- DC.A.nd$Predicted - 1.96 * DC.A.nd$SE # confidence interval upper bound
DC.A.nd$Predict.upr <- DC.A.nd$Predicted + 1.96 * DC.A.nd$SE # confidence interval lower bound
DC.P.nd <- newdata.P
DC.P.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "PAV"]), 
                       max(data3$coral_s[data3$habitat == "PAV"]),
                       length = nrow(DC.P.nd))
pred_DC.p <- predict(mod.sw.D, newdata = DC.P.nd, type = "response", se.fit = TRUE)
pred_DC.p <- as.data.frame(pred_DC.p)
colnames(pred_DC.p) <- c("Predicted","SE")
DC.P.nd <- cbind(DC.P.nd, pred_DC.p)
DC.P.nd$Predict.lwr <- DC.P.nd$Predicted - 1.96 * DC.P.nd$SE # confidence interval upper bound
DC.P.nd$Predict.upr <- DC.P.nd$Predicted + 1.96 * DC.P.nd$SE # confidence interval lower bound
DC.R.nd <- newdata.R
DC.R.nd$coral_s <- seq(min(data3$coral_s[data3$habitat == "ROB"]), 
                       max(data3$coral_s[data3$habitat == "ROB"]),
                       length = nrow(DC.R.nd))
pred_DC.r <- predict(mod.sw.D, newdata = DC.R.nd, type = "response", se.fit = TRUE)
pred_DC.r <- as.data.frame(pred_DC.r)
colnames(pred_DC.r) <- c("Predicted","SE")
DC.R.nd <- cbind(DC.R.nd, pred_DC.r)
DC.R.nd$Predict.lwr <- DC.R.nd$Predicted - 1.96 * DC.R.nd$SE # confidence interval upper bound
DC.R.nd$Predict.upr <- DC.R.nd$Predicted + 1.96 * DC.R.nd$SE # confidence interval lower bound

D.COR <- ggplot() +
  geom_point(aes(x = data3$coral_s, y = data3$D), col = "lightgrey", alpha = 0.3) +
  geom_ribbon(aes(x=DC.A.nd$coral_s,
                  ymin=DC.A.nd$Predict.lwr,
                  ymax=DC.A.nd$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=DC.A.nd$coral_s,
                y=DC.A.nd$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=DC.P.nd$coral_s, 
                  ymin=DC.P.nd$Predict.lwr, 
                  ymax=DC.P.nd$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=DC.P.nd$coral_s, 
                y=DC.P.nd$Predicted), size = .8, color = "#0072B2",
            linetype = "dashed") +
  geom_ribbon(aes(x=DC.R.nd$coral_s,
                  ymin=DC.R.nd$Predict.lwr,
                  ymax=DC.R.nd$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=DC.R.nd$coral_s,
                y=DC.R.nd$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(0, 2, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))

# plot
g_legend<-function(a.gplot){ # function to extract legend from a plot
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(ggplot(data3) +
                     geom_line(aes(x = CORAL, y = habitat, color = habitat)) + 
                     theme(legend.text = element_text(size = 10, family = "Times-Roman"), legend.title = element_text(size = 10,family = "Times-Roman"), 
                           legend.position = "bottom", legend.direction = "horizontal") +
                     guides(color = guide_legend(override.aes = list(size = 2))) +
                     scale_color_manual(values = c("Aggregate reef" = "#D55E00", 
                                                   "Pavement" = "#0072B2",
                                                   "Rock and Boulder" = "#009E73"), 
                                        name = "Habitat type"))
fig3 <- grid.arrange(arrangeGrob(R.DEP + ggtitle("(a)"), R.AGE + ggtitle("(b)"), R.WAVE + ggtitle("(c)"), R.COR + ggtitle("(d)"),
                                 D.DEP + ggtitle("(e)"), D.AGE + ggtitle("(f)"), D.WAVE + ggtitle("(g)"), D.COR + ggtitle("(h)"),
                                 H.DEP + ggtitle("(i)"), H.AGE + ggtitle("(j)"), H.WAVE + ggtitle("(k)"), H.COR + ggtitle("(l)"),
                                 nrow = 3, ncol = 4), 
                     mylegend, 
                     heights = c(10,1))
ggsave("figs/figure3.png", plot = fig3, width = 7, height = 5, units = c("in"), dpi = 600)











# Relationships by habitat #####################################################################################################
hab <- c("AGR","ROB","PAV")

temp1 <- data[data$habitat== hab[1],]
temp2 <- data[data$habitat== hab[2],]
temp3 <- data[data$habitat== hab[3],]

fitR1 <- lm(temp1$R_log10 ~ temp1$wave_max_log10)
fitH1 <- lm(temp1$H_log10 ~ temp1$wave_max_log10)
fitD1 <- lm(temp1$D_sum ~ temp1$wave_max_log10)
fitR2 <- lm(temp2$R_log10 ~ temp2$wave_max_log10)
fitH2 <- lm(temp2$H_log10 ~ temp2$wave_max_log10)
fitD2 <- lm(temp2$D_sum ~ temp2$wave_max_log10)
fitR3 <- lm(temp3$R_log10 ~ temp3$wave_max_log10)
fitH3 <- lm(temp3$H_log10 ~ temp3$wave_max_log10)
fitD3 <- lm(temp3$D_sum ~ temp3$wave_max_log10)

dd <- seq(min(data$depth.m), max(data$depth.m), 1)
modR1 <- lm(R_log10 ~  poly(depth.m, 2), temp1)
modR2 <- lm(R_log10 ~  poly(depth.m, 2), temp2)
modR3 <- lm(R_log10 ~  poly(depth.m, 2), temp3)
modH1 <- lm(H_log10 ~  poly(depth.m, 2), temp1)
modH2 <- lm(H_log10 ~  poly(depth.m, 2), temp2)
modH3 <- lm(H_log10 ~  poly(depth.m, 2), temp3)
modD1 <- lm(D_sum ~  poly(depth.m, 2), temp1)
modD2 <- lm(D_sum ~  poly(depth.m, 2), temp2)
modD3 <- lm(D_sum ~  poly(depth.m, 2), temp3)

imodR1 <- lm(R_log10 ~ rank, temp1)
imodR2 <- lm(R_log10 ~ rank, temp2)
imodR3 <- lm(R_log10 ~ rank, temp3)
imodH1 <- lm(H_log10 ~ rank, temp1)
imodH2 <- lm(H_log10 ~ rank, temp2)
imodH3 <- lm(H_log10 ~ rank, temp3)
imodD1 <- lm(D_sum ~ rank, temp1)
imodD2 <- lm(D_sum ~ rank, temp2)
imodD3 <- lm(D_sum ~ rank, temp3)

par(mfrow=c(3,3), mar=c(5, 6.1, 3, 3.5))

plot(R_log10 ~ wave_max_log10, data, pch = 16, cex = .5, ylab = "Rugosity (log10)", xlab = "Exposure, kW/m (log10)", col = "grey")
abline(fitR1, col = "red", lwd = 2)
abline(fitR2, col = "seagreen", lwd = 2)
abline(fitR3, col = "steelblue", lwd = 2)

plot(H_log10 ~ wave_max_log10, data, pch = 16, cex = .5, ylab = "Height range (log10)", xlab = "Exposure, kW/m (log10)", col = "grey")
abline(fitH1, col = "red", lwd = 2)
abline(fitH2, col = "seagreen", lwd = 2)
abline(fitH3, col = "steelblue", lwd = 2)

plot(D_sum ~ wave_max_log10, data, pch = 16, cex = .5, ylab = "Fractal dimension", xlab = "Exposure, kW/m (log10)", col = "grey")
abline(fitD1, col = "red", lwd = 2)
abline(fitD2, col = "seagreen", lwd = 2)
abline(fitD3, col = "steelblue", lwd = 2)

plot(R_log10 ~ depth.m, data, pch = 16, cex = .5, ylab = "Rugosity (log10)", xlab = "Depth, m", col = "grey")
lines(dd, predict(modR1, list(depth.m=dd)), lwd = 2, col = "red")
lines(dd, predict(modR2, list(depth.m=dd)), lwd = 2, col = "seagreen")
lines(dd, predict(modR3, list(depth.m=dd)), lwd = 2, col = "steelblue")

plot(H_log10 ~ depth.m, data, pch = 16, cex = .5, ylab = "Height range (log10)", xlab = "Depth, m", col = "grey")
lines(dd, predict(modH1, list(depth.m=dd)), lwd = 2, col = "red")
lines(dd, predict(modH2, list(depth.m=dd)), lwd = 2, col = "seagreen")
lines(dd, predict(modH3, list(depth.m=dd)), lwd = 2, col = "steelblue")

plot(D_sum ~ depth.m, data, pch = 16, cex = .5, ylab = "Fractal dimension", xlab = "Depth, m", col = "grey")
lines(dd, predict(modD1, list(depth.m=dd)), lwd = 2, col = "red")
lines(dd, predict(modD2, list(depth.m=dd)), lwd = 2, col = "seagreen")
lines(dd, predict(modD3, list(depth.m=dd)), lwd = 2, col = "steelblue")

plot(R_log10 ~ rank, data, pch = 16, cex = .5, ylab = "Rugosity (log10)", xlab = "Island rank", col = "grey")
abline(imodR1, col = "red", lwd = 2)
abline(imodR2, col = "seagreen", lwd = 2)
abline(imodR3, col = "steelblue", lwd = 2)

plot(H_log10 ~ rank, data, pch = 16, cex = .5, ylab = "Height range (log10)", xlab = "Island rank", col = "grey")
abline(imodH1, col = "red", lwd = 2)
abline(imodH2, col = "seagreen", lwd = 2)
abline(imodH3, col = "steelblue", lwd = 2)

plot(D_sum ~ rank, data, pch = 16, cex = .5, ylab = "Fractal dimension", xlab = "Island rank", col = "grey")
abline(imodD1, col = "red", lwd = 2)
abline(imodD2, col = "seagreen", lwd = 2)
abline(imodD3, col = "steelblue", lwd = 2)

# Contour plots ##################################################################################################################
lmodR <- lmer(R_log10 ~ wave_max_log10_s * depth_s + (1|habitat) + (1|site), data) #ranova() both significant
lmodH <- lmer(H_log10 ~  wave_max_log10_s * depth_s + (1|habitat)+ (1|site), data) #ranova() both significant
lmodD <- lmer(D_mean ~  wave_max_log10_s * depth_s + (1|habitat)+ (1|site), data) #ranova() both significant

dd <- seq(min(data$depth_s), max(data$depth_s), .09)
ww <- seq(min(data$wave_max_log10_s), max(data$wave_max_log10_s), 0.13)

ch <- chull(data2$depth_s, data2$wave_max_log10_s)
ch <- c(ch, ch[1])

matR <- matD <- matH <- matrix(NA, length(dd), length(ww))

for (d in 1:length(dd)) {
  for (w in 1:length(ww)) {
    matR[d,w] <- predict(lmodR, data.frame(wave_max_log10_s=ww[w], depth_s=dd[d]), re.form=NA)
    # matD[d,w] <- predict(lmodD, data.frame(wave_max_log10_s=ww[w], depth_s=dd[d]), re.form=NA)
    # matH[d,w] <- predict(lmodH, data.frame(wave_max_log10_s=ww[w], depth_s=dd[d]), re.form=NA)
    # 
    pp <- point.in.polygon(dd[d],ww[w],data$depth_s[ch],data$wave_max_log10_s[ch])
    if (pp == 0) {
      matR[d,w] <- NA
      # matD[d,w] <- NA
      # matH[d,w] <- NA
    }
  }
}


R <- aggregate(depth_s ~ habitat, data, mean)
R$wave <- aggregate(wave_max_log10_s ~ habitat, data, mean)$wave_max_log10_s
R$wave_sd <- aggregate(wave_max_log10_s ~ habitat, data, sd)$wave_max_log10_s
R$wave_n <- aggregate(wave_max_log10_s ~ habitat, data, length)$wave_max_log10_s
R$wave_se <- R$wave_sd / sqrt(R$wave_n)
R$depth_sd <- aggregate(depth_s ~ habitat, data, sd)$depth_s
R$depth_n <- aggregate(depth_s ~ habitat, data, length)$depth_s
R$depth_se <- R$depth_sd / sqrt(R$depth_n)

collist <- c("steelblue", "red") 
mycolors <- colorRampPalette(collist)(200)

#png(paste0("figures/", hab[c], "_contours.png"), width = 400, height = 600, res = 100)

par(mfrow=c(3,1), mar=c(5, 6.1, 3, 4)) 

image(dd, ww, matR, col=mycolors, xlab="Depth, m", ylab="Exposure, kW/m (log10)")
contour(dd, ww, 10^matR, add = TRUE, col = "black", nlevels = 4, labcex = .8)
mtext(expression(bold("Rugosity (log10)")), 2, 0, cex = 1, line = 4)
#points(wave_max_log10 ~ depth.m, data, col = "black")
#text(R$depth.m+2, R$wave, R$habitat, font = 2, cex = 1, col = "black")
#arrows(R$depth.m - R$depth_se, R$wave, R$depth.m + R$depth_se, R$wave, code=3, angle=90, length=0.025, lwd = .4)
#arrows(R$depth.m, R$wave - R$wave_se, R$depth.m, R$wave + R$wave_se, code=3, angle=90, length=0.025, lwd = .4)
gradientLegend(valRange = c(min(10^matR, na.rm = TRUE), max(10^matR, na.rm = TRUE)), color = mycolors, pos = 0.5,  length = .7, depth = .03, dec = 2, fit.margin = TRUE, n.seg = 1)

image(dd, ww, matH, col=mycolors, xlab="Depth, m", ylab="Exposure, kW/m (log10)")
contour(dd, ww, 10^matH, add = TRUE, col = "black", nlevels = 4, labcex = .8)
mtext(expression(bold("Height range (log10)")), 2, 0, cex = 1, line = 4)
#points(wave_max_log10 ~ depth.m, data, col = "black")
#text(R$depth.m+2, R$wave, R$habitat, font = 2, cex = 1, col = "black")
#arrows(R$depth.m - R$depth_se, R$wave, R$depth.m + R$depth_se, R$wave, code=3, angle=90, length=0.025, lwd = .4)
#arrows(R$depth.m, R$wave - R$wave_se, R$depth.m, R$wave + R$wave_se, code=3, angle=90, length=0.025, lwd = .4)
gradientLegend(valRange = c(min(10^matH, na.rm = TRUE), max(10^matH, na.rm = TRUE)), color = mycolors, pos = 0.5,  length = .7, depth = .03, dec = 2, fit.margin = TRUE, n.seg = 1)

image(dd, ww, matD, col=mycolors, xlab="Depth, m", ylab="Exposure, kW/m (log10)")
contour(dd, ww, matD, add = TRUE, col = "black", nlevels = 4, labcex = .8)
mtext(expression(bold("Fractal dimension")), 2, 0, cex = 1, line = 4)
#points(wave_max_log10 ~ depth.m, data, col = "black")
#text(R$depth.m+2, R$wave, R$habitat, font = 2, cex = 1, col = "black")
#arrows(R$depth.m - R$depth_se, R$wave, R$depth.m + R$depth_se, R$wave, code=3, angle=90, length=0.025, lwd = .4)
#arrows(R$depth.m, R$wave - R$wave_se, R$depth.m, R$wave + R$wave_se, code=3, angle=90, length=0.025, lwd = .4)
gradientLegend(valRange = c(min(matD, na.rm = TRUE),max(matD, na.rm = TRUE)), color = mycolors, pos = 0.5,  length = .7, depth = .03, dec = 2, fit.margin = TRUE, n.seg = 1)

#dev.off()

# I <- aggregate(depth ~ island, data, mean)
# I$wave <- aggregate(wave_max_log10 ~ island, data, mean)$wave_max_log10
# I$wave_sd <- aggregate(wave_max_log10 ~ island, data, sd)$wave_max_log10
# I$wave_n <- aggregate(wave_max_log10 ~ island, data, length)$wave_max_log10
# I$wave_se <- I$wave_sd / sqrt(I$wave_n)
# I$depth_sd <- aggregate(depth ~ island, data, sd)$depth
# I$depth_n <- aggregate(depth ~ island, data, length)$depth
# I$depth_se <- I$depth_sd / sqrt(I$depth_n)

# by island
isl <- unique(data$island)
mod <- lm(R_log10 ~ wave_max_log10 * poly(depth.m, 2) * island, data)

par(mfrow=c(2,4), mar=c(5, 6.1, 3, 3.5)) 

for (c in 1:8) {
  
  temp <- data[data$island== isl[c],]
  
  dd <- seq(min(data$depth.m), max(data$depth.m), 1)
  ww <- seq(min(data$wave_max_log10), max(data$wave_max_log10), 0.1)
  
  matR <- matD <- matH <- matrix(NA, length(dd), length(ww))
  
  for (d in 1:length(dd)) {
    for (w in 1:length(ww)) {
      matR[d,w] <- predict(mod, data.frame(wave_max_log10=ww[w], depth.m=dd[d], island=isl[c]), re.form=NA)
    }
  }
  
  collist <- c("#D81B60","#1E88E5","#FFC107") 
  mycolors <- colorRampPalette(collist)(200)
  
  image(dd, ww, matR, col=mycolors, xlab="Depth, m", ylab="Exposure, kW/m (log10)")
  points(wave_max_log10 ~ depth.m, temp)
  contour(dd, ww, 10^matR, add = TRUE, col = "black", nlevels = 4)
  mtext(expression(bold("Rugosity")), 2, 0, cex = 1, line = 4)
  gradientLegend(valRange = c(min(10^matR), max(10^matR)), color = mycolors, pos = 0.5,  length = .7, depth = .03, dec = 2, fit.margin = TRUE, n.seg = 1)
  title(main = isl[c], cex.main = 2)
  
}


# PCA ##########################################################################################################################
data <- read.csv("Complexity_all.csv")
sub <- read.csv("Complexity_sub.csv")

L <- 2 #2x2m box
scl <- L / c(1, 2, 3.125, 6.25, 12.5, 25, 50, 100, 200)
L0 <- min(scl) 

data.pca <- data %>% mutate(R = 0.5*log10(R_theory_mean^2 - 1),
                         D = log10(L/L0)*(D_mean-3),
                         H = log10(H_mean/(sqrt(2)*L0))) 

names(data.pca)[64] <- "Rugosity"
names(data.pca)[65] <- "Fractal dimension"
names(data.pca)[66] <- "Height range"
names(data.pca)[27] <- "Exposure"
names(data.pca)[28] <- "Depth"
names(data.pca)[29] <- "Geological Reef Age"
names(data.pca)[59] <- "Coral Cover"


pca.var <- prcomp(data.pca[,c(64,65,66,27,28,29,59)], center = TRUE, scale = TRUE) 
summary(pca.var)

# eigenvalues -- according to Kaiser criterion, keep those greater than 1
var <- pca.var$sdev ^ 2

varPer <- var/sum(var)*100
barplot(varPer, xlab= "PC", ylab = "Percent variance", names.arg=1:length(varPer), las=1, ylim=c(0, max(varPer)), col = "grey")
abline(h = (1/7)*100, col = "red")


# relationships of components
load <- pca.var$rotation
cut <- 1/7 #cutoff
cut2 <- -cut
sort1 <- load[order(load[,1]),1]
dotplot(sort1)

sort2 <- load[order(load[,2]),2]
dotplot(sort2)

pc1 <- as.numeric(load[,1]) * 100
pc1[8] <- varPer[1]
pc2 <- as.numeric(load[,2]) * 100
pc2[8] <- varPer[2]
pc3 <- as.numeric(load[,3]) * 100
pc3[8] <- varPer[3]

library(kableExtra)
sort <- data.frame(row.names = c("Rugosity","Fractal dimension","Height range","Exposure","Depth","Geological reef age","Coral cover","Total variance"))
sort$PC1 <- format(pc1, digits = 3)
sort$PC2 <- format(pc2, digits = 3)
sort$PC3 <- format(pc3, digits = 3)
sort$PC1[c(1,3,4,5,6,7)] <- cell_spec(sort$PC1[c(1,3,4,5,6,7)], bold =  T)
sort$PC2[c(2,3,4,5,7)] <- cell_spec(sort$PC2[c(2,3,4,5,7)], bold =  T)
sort$PC3[c(1,3,4,5)] <- cell_spec(sort$PC3[c(1,3,4,5)], bold =  T)


sort %>% 
  kbl(escape = F) %>%
  kable_classic(html_font = "Times New Roman") %>%
  row_spec(7, extra_css = "border-bottom: 2px solid;")
  
# loadings plot
g <- autoplot(pca.var, data.pca,  loadings = TRUE, loadings.col = "black", loadings.label = TRUE, col = "habitat",
              loadings.label.size = 1, loadings.label.col = "black", frame = TRUE, frame.type = "norm") +
  scale_color_manual(values = alpha(c("#D55E00","#0072B2","#009E73"), 0.4)) +
  scale_fill_manual(values = alpha(c("#D55E00","#0072B2","#009E73"), 0.4)) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10)) 
  

arrow_ends <- layer_data(g, 2)[,c(2,4)]
for (a in 1:nrow(arrow_ends)) {
  for (b in 1:2) {
    value <- arrow_ends[a,b]
    if (value > 0) {
      arrow_ends[a,b] <- value + 0.008
    }
    if (value < 0) {
      arrow_ends[a,b] <- value - 0.006
    }
  }
}


g2 <- autoplot(pca.var, data.pca, loadings = TRUE, loadings.col = "black", col = "habitat",
         loadings.label.size = 2, loadings.label.col = "black", frame = TRUE, frame.type = "norm") +
  scale_color_manual(values = alpha(c("#D55E00","#0072B2","#009E73"), 0.4)) +
  scale_fill_manual(values = alpha(c("#D55E00","#0072B2","#009E73"), 0.4)) +
  #geom_point(data = arrow_ends, aes(xend, yend), size = 3, col = "white") +
  geom_label(data = arrow_ends, aes(xend, yend), label = names(pca.var$center), size = 2, label.size = 0,family = "Times-Roman") +
  #guides(fill = guide_legend(title="Habitat Type")) +
  labs(colour = "Habitat Type") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") #legend.key.size = unit(.2, 'cm'), legend.position = "bottom", 

png("FINALDRAFT/Figure2.png", width = 5, height = 4, unit= "in", res = 600)
grid.arrange(arrangeGrob(g2), 
             mylegend, 
             heights = c(10,1))
dev.off()

# 3d scatterplots part 1 (scatterplot3D) #######################################################################################
# myteal <- t_col("darkcyan", perc = 50)

data_subR <- aggregate(R_log10 ~ longitude + latitude, data, FUN = mean)
lmR <- lm(data$R_log10 ~ data$longitude + data$latitude)
data_subH <- aggregate(H_log10 ~ longitude + latitude, data, FUN = mean)
lmH <- lm(data$H_log10 ~ data$longitude + data$latitude)
data_subD <- aggregate(D_sum ~ longitude + latitude, data, FUN = mean)
lmD <- lm(data$D_sum ~ data$longitude + data$latitude)
data_subV <- aggregate(vrm_1cm ~ longitude + latitude, data, FUN = mean)
lmV <- lm(data$vrm_1cm ~ data$longitude + data$latitude)

dev.off()
par(mfrow=c(1,3), mar=c(5, 6.1, 3, 2.1))

par(pty = "s")
sR <- scatterplot3d(data_subR$longitude, data_subR$latitude, data_subR$R_log10, pch = 19, type = "p", color = "black",
                    y.margin.add = .5, grid = TRUE, box = TRUE, xlab = expression(paste("Longitude (", degree, ")")),
                    ylab = expression(paste("Latitude (", degree, ")")), zlab = "Rugosity (log10)", angle = 30, cex.axis = .5)
sR$plane3d(lmR, draw_polygon = TRUE, draw_lines = TRUE, polygon_args = list(col = myteal))

sH <- scatterplot3d(data_subH$longitude, data_subH$latitude, data_subH$H_log10, pch = 19, type = "p", color = "black",
                    y.margin.add = .5, grid = TRUE, box = TRUE, xlab = expression(paste("Longitude (", degree, ")")),
                    ylab = expression(paste("Latitude (", degree, ")")), zlab = "Height range (log10)", angle = 55, cex.axis = .5)
sH$plane3d(lmH, draw_polygon = TRUE, draw_lines = TRUE, polygon_args = list(col = myteal))

sD <- scatterplot3d(data_subD$longitude, data_subD$latitude, data_subD$D_sum, pch = 19, type = "p", color = "black",
                    y.margin.add = .5, grid = TRUE, box = TRUE, xlab = expression(paste("Longitude (", degree, ")")),
                    ylab = expression(paste("Latitude (", degree, ")")), zlab = "Fractal Dimension", angle = 30, cex.axis = .5)
sD$plane3d(lmD, draw_polygon = TRUE, draw_lines = TRUE, polygon_args = list(col = myteal))

sV <- scatterplot3d(data_subV$latitude, data_subV$longitude, data_subV$vrm_1cm, pch = 19, type = "p", color = "black",
                    y.margin.add = .5, grid = TRUE, box = TRUE, xlab = expression(paste("Latitude (", degree, ")")),
                    ylab = expression(paste("Longitude (", degree, ")")), zlab = "VRM", angle = 30, cex.axis = .5)
sV$plane3d(lmV, draw_polygon = TRUE, draw_lines = TRUE, polygon_args = list(col = myteal))

# 3d scatter plots part 2 (scatter3D) ##########################################################################################
collist2 <- c("cyan","cyan3","black")
mycolors2 <- colorRampPalette(collist2)(200)
tblack <- t_col("black", perc = 70)
library(usmap)
hawaii <- us_map(include = c("HI"))

x <- data$longitude
y <- data$latitude
zR <- data$R_log10
zH <- data$H_log10
zD <- data$D_sum
fitR <- lm(zR ~ x + y)
fitH <- lm(zH ~ x + y)
fitD <- lm(zD ~ x + y)
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid(x = x.pred, y = y.pred)
z.predR <- matrix(predict(fitR, newdata = xy), nrow = grid.lines, ncol= grid.lines)
z.predH <- matrix(predict(fitH, newdata = xy), nrow = grid.lines, ncol= grid.lines)
z.predD <- matrix(predict(fitD, newdata = xy), nrow = grid.lines, ncol= grid.lines)
fitpointsR <- predict(fitR)
fitpointsH <- predict(fitH)
fitpointsD <- predict(fitD)

par(mfrow=c(2,2), mar=c(2,1,1,1))
plot(hawaii$x, hawaii$y, type = "p", pch = 16, ylab = "Latitude", xlab = "Longitude")
scatter3D(x, y, zR, ticktype = "detailed", xlab = "Longitude", ylab = "Latitude", zlab = "Rugosity (log10)", theta = 40, col = tblack, 
          cex = 1, phi = 25, pch = 16, surf = list(x = x.pred, y = y.pred, z = z.predR, fit = fitpointsR, col = mycolors2))
scatter3D(x, y, zH, ticktype = "detailed", xlab = "Longitude", ylab = "Latitude", zlab = "Height range (log10)", theta = 40, col = tblack, 
          cex = 1, phi = 25, pch = 16, surf = list(x = x.pred, y = y.pred, z = z.predH, fit = fitpointsH, col = mycolors2))
scatter3D(x, y, zD, ticktype = "detailed", xlab = "Longitude", ylab = "Latitude", zlab = "Fractal dimension", theta = 40, col = tblack,
          cex = 1, phi = 25, pch = 16, surf = list(x = x.pred, y = y.pred, z = z.predD, fit = fitpointsD, col = mycolors2))



# heirarchial partitioning #####################################################################################################
# by habitat
hab <- c("AGR","ROB","PAV")
for (h in (1:3)){
  
  temp <- data[data$habitat== hab[h],]
  
  R <- hier.part(temp$R_log10, temp[c("age","wave_max_log10","depth")]) #ordered in decreasing scale (island greatest, exposure operates at shoreline distances, then depth)
  H <- hier.part(temp$H_log10, temp[c("age","wave_max_log10","depth")])
  D <- hier.part(temp$D_mean, temp[c("age","wave_max_log10","depth")])
  
  #png(paste0("figures/", hab[h], "_hier.png"), width = 600, height = 200, res = 100)
  
  par(mfrow=c(1,3), mar=c(5, 4.1, 3, 2))
  
  xr <- barplot(R$I.perc$ind.exp.var, ylim = c(0,100), names.arg = c("island","exposure","depth"), ylab = "% independent effects", xlab = "Rugosity", cex.names = .8)
  text(x = xr, y = R$I.perc$ind.exp.var, label = round(R$I.perc$ind.exp.var, digits = 1), pos = 3, col = "dark red", cex = .8)
  xh <- barplot(H$I.perc$ind.exp.var, ylim = c(0,100), names.arg = c("island","exposure","depth"), ylab = "% independent effects", xlab = "Height Range", main = paste0(hab[h]), cex.main = 2, cex.names = .8)
  text(x = xh, y = H$I.perc$ind.exp.var, label = round(H$I.perc$ind.exp.var, digits = 1), pos = 3, col = "dark red", cex = .8)
  xd <- barplot(D$I.perc$ind.exp.var, ylim = c(0,100), names.arg = c("island","exposure","depth"), ylab = "% independent effects", xlab = "Fractal Dimension", cex.names = .8)
  text(x = xd, y = D$I.perc$ind.exp.var, label = round(D$I.perc$ind.exp.var, digits = 1), pos = 3, col = "dark red", cex = .8)
  
  dev.off()
}

# all data
R <- hier.part(data$R_log10, data[c("island","wave_max_log10","depth")]) #ordered in decreasing scale (island greatest, exposure operates at shoreline distances, then depth)
H <- hier.part(data$H_log10, data[c("island","wave_max_log10","depth")])
D <- hier.part(data$D_sum, data[c("island","wave_max_log10","depth")])

png(paste0("figures/", hab[h], "_hier.png"), width = 600, height = 200, res = 100)

par(mfrow=c(1,3), mar=c(5, 4.1, 3, 2))

xr <- barplot(R$I.perc$ind.exp.var, ylim = c(0,100), names.arg = c("island","exposure","depth"), ylab = "% independent effects", xlab = "Rugosity", cex.names = .8)
text(x = xr, y = R$I.perc$ind.exp.var, label = round(R$I.perc$ind.exp.var, digits = 1), pos = 3, col = "dark red", cex = .8)
xh <- barplot(H$I.perc$ind.exp.var, ylim = c(0,100), names.arg = c("island","exposure","depth"), ylab = "% independent effects", xlab = "Height Range", cex.main = 2, cex.names = .8)
text(x = xh, y = H$I.perc$ind.exp.var, label = round(H$I.perc$ind.exp.var, digits = 1), pos = 3, col = "dark red", cex = .8)
xd <- barplot(D$I.perc$ind.exp.var, ylim = c(0,100), names.arg = c("island","exposure","depth"), ylab = "% independent effects", xlab = "Fractal Dimension", cex.names = .8)
text(x = xd, y = D$I.perc$ind.exp.var, label = round(D$I.perc$ind.exp.var, digits = 1), pos = 3, col = "dark red", cex = .8)

dev.off()

# variation in intercepts ######################################################################################################
lmodR <- lmer(R_log10 ~ wave_max_log10*depth + wave_max_log10 + depth + develop + (1|island), data)
summary(lmodR)
anova(lmodR)
temp <- ranef(lmodR)$island
temp$island <- rownames(temp)
temp <- temp[order(temp$`(Intercept)`),]
barplot(temp$`(Intercept)`, names=temp$island)

# centroid of depth bins ######################################################################################################### 
plot(R_log10 ~ D_sum, data)
points(R_log10 ~ D_sum, data[data$depth_bin=="Deep",], col="red")
points(R_log10 ~ D_sum, data[data$depth_bin=="Shallow",], col="green")
points(R_log10 ~ D_sum, data[data$depth_bin=="Mid",], col="blue")
rg <- aggregate(R_log10 ~ depth_bin, data, mean)
fd <- aggregate(D_sum ~ depth_bin, data, mean)
text(fd$D_sum, rg$`R_log10`, rg$depth_bin)

# 3 descriptors plot #############################################################################################################
source("R/functions.R")


data <- read.csv("Complexity_all.csv")
L <- 2
L0 <- 0.01
data$R2_log10 <- log10(data$R_theory_mean^2 - 1)
data$HL0_log10 <- log10(data$H_mean / (L0 * sqrt(2)))

grid.lines <- 350 
R <- data$R2_log10 
H <- data$HL0_log10 
D <- data$D_mean

R.pred <- seq(min(R)-0.25, max(R)+0.25, length=grid.lines)
H.pred <- seq(min(H)-0.25, max(H)+0.25, length=grid.lines)
mat <- expand.grid(R=R.pred, H=H.pred)
D.mat <- matrix(D_func((sqrt(2) * L0)*(10^mat$H), sqrt(10^mat$R+1), L, L0), nrow = grid.lines, ncol = grid.lines)


R_mean <- aggregate(R2_log10 ~ habitat, data, mean)
H_mean <- aggregate(HL0_log10 ~ habitat, data, mean)
D_mean <- aggregate(D_mean ~ habitat, data, mean)

png("Figure1.png", width = 5, height = 5, units = "in", res = 300)
par(mar=c(1.2, 2, 1, 1) +0.1, mfrow = c(2,2), family = "Times", ps = 10)
scatter3D(R, H, D, pch = 20, cex = 1, 
          col= rev(hcl.colors(100, "grays", alpha = 0.3)), 
          xlim=c(min(R)-0.25, max(R)+0.25), 
          ylim=c(min(H)-0.25, max(H)+0.25), 
          zlim=c(2, max(D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=data$D_theory_mean), 
          theta=215, 
          phi=0,
          colkey=list(side = 2, length = 0.5, width = 1, line.clab = 1, dist=0),
          clab=expression(italic(D))
          )

points3D(R_mean$R2_log10 - 0.2, H_mean$HL0_log10 + 0.2, D_mean$D_mean,
         col = c("#009E73","#0072B2","#D55E00"), cex = 2, pch = 20,
         # xlim=c(min(R)-0.25, max(R)+0.25),
         # ylim=c(min(H)-0.25, max(H)+0.25),
         # zlim=c(2, max(D)+0.05),
         # ylab="Height range",
         # xlab="Rugosity",
         # zlab="Fractal dimension",
         # theta=215,
         # phi=0,
         surf=NULL,
         colkey=FALSE,
         add=TRUE)
mtext("(a)", side = 3, at = -.5)

ss_tot <- sum((data$D_mean - mean(data$D_theory_mean))^2)
ss_f <- sum((data$D_mean - mean(data$D_theory_mean))^2)
ss_res <- sum((data$D_mean - data$D_theory_mean)^2)
r_2 <- 1 - (ss_res / ss_tot)
print(r_2)
text3D(1.8, 1.2, 2.7, labels=expression(italic(r)^2 == 0.967), surf = NULL, add = TRUE, family = "Times")



#####
data.bayes.AGR <- data[data$habitat == "AGR" ,]
data.bayes.AGR$R2_log10 <- log10(data.bayes.AGR$R_theory_mean^2 - 1)
data.bayes.AGR$HL0_log10 <- log10(data.bayes.AGR$H_mean / (L0 * sqrt(2)))

grid.lines <- 350 
R2 <- data.bayes.AGR$R2_log10 
H2 <- data.bayes.AGR$HL0_log10 
D2 <- data.bayes.AGR$D_mean

AGRcollist <- c("#D55E00","#FFFFFF")
AGRcolors <- colorRampPalette(AGRcollist)(100)

scatter3D(R2, H2, D2, pch = 20, cex = 1, 
          col = rev(AGRcolors),
          #col= rev(hcl.colors(100, "reds")),
          xlim=c(min(R)-0.25, max(R)+0.25), 
          ylim=c(min(H)-0.25, max(H)+0.25), 
          zlim=c(2, max(D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=data.bayes.AGR$D_theory_mean), 
          theta=215, 
          phi=0,
          colkey=list(side = 4, length = 0.5, width = 1, line.clab = 1, dist=0))
title("Aggregate reef", family = "Times")
mtext("(b)", side = 3, at = -.5)
#text3D(1.8, 1.2, 2.7, labels="Aggregate reef", add = TRUE, cex=.8, family = "Times")


#####
data.bayes.PAV <- data[data$habitat == "PAV" ,]
data.bayes.PAV$R2_log10 <- log10(data.bayes.PAV$R_theory_mean^2 - 1)
data.bayes.PAV$HL0_log10 <- log10(data.bayes.PAV$H_mean / (L0 * sqrt(2)))

grid.lines <- 350 
R2 <- data.bayes.PAV$R2_log10 
H2 <- data.bayes.PAV$HL0_log10 
D2 <- data.bayes.PAV$D_mean

PAVcollist <- c("#0072B2","#FFFFFF")
PAVcolors <- colorRampPalette(PAVcollist)(100)

scatter3D(R2, H2, D2, pch = 20, cex = 1, 
          col= rev(PAVcolors),
          xlim=c(min(R)-0.25, max(R)+0.25), 
          ylim=c(min(H)-0.25, max(H)+0.25), 
          zlim=c(2, max(D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=data.bayes.PAV$D_theory_mean), 
          theta=215, 
          phi=0,
          colkey=list(side = 2, length = 0.5, width = 1, line.clab = 1, dist=0))
title("Pavement", family = "Times")
mtext("(c)", side = 3, at = -.5)
#text3D(1.8, 1.2, 2.7, labels="Pavement", add = TRUE, cex=.8, family = "Times")

#####
data.bayes.ROB <- data[data$habitat_original == "ROB" ,]
data.bayes.ROB$R2_log10 <- log10(data.bayes.ROB$R_theory_mean^2 - 1)
data.bayes.ROB$HL0_log10 <- log10(data.bayes.ROB$H_mean / (L0 * sqrt(2)))

grid.lines <- 350 
R2 <- data.bayes.ROB$R2_log10 
H2 <- data.bayes.ROB$HL0_log10 
D2 <- data.bayes.ROB$D_mean

ROBcollist <- c("#009E73","#FFFFFF")
ROBcolors <- colorRampPalette(ROBcollist)(100)

scatter3D(R2, H2, D2, pch = 20, cex = 1, 
          col= rev(ROBcolors),
          xlim=c(min(R)-0.25, max(R)+0.25), 
          ylim=c(min(H)-0.25, max(H)+0.25), 
          zlim=c(2, max(D)+0.05), 
          ylab="Height range", 
          xlab="Rugosity", 
          zlab="Fractal dimension", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=data.bayes.ROB$D_theory_mean), 
          theta=215, 
          phi=0,
          colkey=list(side = 4, length = 0.5, width = 1, line.clab = 1, dist=0))
title("Rock & Boulder", family = "Times")
mtext("(d)", side = 3, at = -.5)
#text3D(1.8, 1.2, 2.7, labels="Rock & boulder", add = TRUE, family = "Times")

dev.off()

# simple relationship plot with standard error ###################################################################################
plot(D_sum ~ depth, data, las=2)
mod <- lm(D_sum ~ depth, data)
summary(mod)
d <- unique(sort(data$depth))
pre <- predict(mod, list(depth=d), se.fit=TRUE)
lines(d, pre$fit)
polygon(c(d, rev(d)), c(pre$fit + pre$se.fit, rev(pre$fit - pre$se.fit)), border=NA, col=rgb(0,0,0,0.3))

# 2d descriptors plot ##########################################################################################################
col.func <- colorRampPalette(c(rgb(.8,.8,.8,0.7),rgb(0,0,0,0.7)), alpha = FALSE)

D <- seq(2, 2.6, length=50)
L <- 2
L0 <- 0.01
data$R2_log10 <- log10(data$R_theory_mean^2 - 1)
data$HL0_log10 <- log10(data$H_mean / (L0 * sqrt(2)))

R2_log10 <- seq(min(data$R2_log10)-0.1, max(data$R2_log10)+0.1, length=50)
mat_HL0 <- matrix(NA, length(D), length(R2_log10))

for (i in 1:length(D)) {
  for (j in 1:length(R2_log10)) {
    HL0_log10 <- HL0_func(D[i], sqrt((10^R2_log10[j]) + 1), L, L0)
    mat_HL0[i, j] <- HL0_log10
  }
}

par(mar=c(7.1, 4.1, 4.1, 5.1), xpd=TRUE, mfrow = c(1,3))
image(D, R2_log10, mat_HL0, col=NA, axes=FALSE, ylab=expression(paste("Rugosity (", italic(R)^2 - 1, ")")), xlab=expression(paste("Fractal dimension (", italic(D), ")")), xlim=c(1.8, 2.6), ylim=c(-1.5, 1.5))
axis(2, at=log10(c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 12, 30)), labels=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 12, 30))
axis(1, at=seq(1.8, 2.6, 0.2))
text(1.89, 1.2, "Greater height\nrange", cex=0.7)
text(2.53, -1.3, "Smaller height\nrange", cex=0.7)
points(R2_log10 ~ D_theory_mean, data, col = col.func(10)[as.numeric(cut(data$wave_max, breaks = 10))], pch = 20)
temp.m <- data[data$wave_max == max(data$wave_max),]
temp.n <- data[data$wave_max == min(data$wave_max),]
points(mean(R2_log10) ~ mean(D_theory_mean), temp.m, col = "red", pch = 20, cex = 2)
points(mean(R2_log10) ~ mean(D_theory_mean), temp.n, col = "green", pch = 20, cex = 2)
legend("bottom", legend = c("maximum exposure","minimum exposure"), col = c("red","green"), pch = 20, ncol = 2, bty = "n", inset = c(0,-.22))

image(D, R2_log10, mat_HL0, col=NA, axes=FALSE, ylab=expression(paste("Rugosity (", italic(R)^2 - 1, ")")), xlab=expression(paste("Fractal dimension (", italic(D), ")")), xlim=c(1.8, 2.6), ylim=c(-1.5, 1.5))
axis(2, at=log10(c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 12, 30)), labels=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 12, 30))
axis(1, at=seq(1.8, 2.6, 0.2))
text(1.89, 1.2, "Greater height\nrange", cex=0.7)
text(2.53, -1.3, "Smaller height\nrange", cex=0.7)
points(R2_log10 ~ D_theory_mean, data, col = col.func(10)[as.numeric(cut(data$depth, breaks = 10))], pch = 20)
temp.m <- data[data$depth == max(data$depth),]
temp.n <- data[data$depth == min(data$depth),]
points(mean(R2_log10) ~ mean(D_theory_mean), temp.m, col = "red", pch = 20, cex = 2)
points(mean(R2_log10) ~ mean(D_theory_mean), temp.n, col = "green", pch = 20, cex = 2)
legend("bottom", legend = c("maximum depth","minimum depth"), col = c("red","green"), pch = 20, ncol = 2, bty = "n", inset = c(0,-.22))

image(D, R2_log10, mat_HL0, col=NA, axes=FALSE, ylab=expression(paste("Rugosity (", italic(R)^2 - 1, ")")), xlab=expression(paste("Fractal dimension (", italic(D), ")")), xlim=c(1.8, 2.6), ylim=c(-1.5, 1.5))
axis(2, at=log10(c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 12, 30)), labels=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 12, 30))
axis(1, at=seq(1.8, 2.6, 0.2))
text(1.89, 1.2, "Greater height\nrange", cex=0.7)
text(2.53, -1.3, "Smaller height\nrange", cex=0.7)
points(R2_log10 ~ D_theory_mean, data, col = col.func(10)[as.numeric(cut(data$age, breaks = 10))], pch = 20)
temp.m <- data[data$age == max(data$age),]
temp.n <- data[data$age == min(data$age),]
points(mean(R2_log10) ~ mean(D_theory_mean), temp.m, col = "red", pch = 20, cex = 2)
points(mean(R2_log10) ~ mean(D_theory_mean), temp.n, col = "green", pch = 20, cex = 2)
legend("bottom", legend = c("oldest","youngest"), col = c("red","green"), pch = 20, ncol = 2, bty = "n", inset = c(0,-.22))

# temp <- data[data$habitat == "ROB",]
# pts <- chull(temp$D_theory_mean, temp$R2_log10)
# pts <- c(pts, pts[1])
# polygon(temp$D_theory_mean[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[1], border=NA)
# points(R2_log10 ~ D_theory_mean, temp, col=hcl.colors(6, alpha=1)[1], pch=20)

# R by D #######################################################################################################################

ggplot(data, aes(x = D_sum, y = R_log10, col = habitat, fill = habitat, shape = habitat)) +
  geom_point() +
  labs(col= "Habitat Type", shape = "Habitat Type") +
  scale_color_manual(values = alpha(c("red","steelblue","seagreen"), 0.8)) +
  annotate("text", x = 2.13, y = 2.5, label = "Greater height\nrange") +
  annotate("text", x = 2.6, y = 0.1, label = "Smaller height\nrange") +
  scale_y_continuous(name = "Rugosity (log10)", limits = c(0, 2.6), breaks = c(0,0.5,1,1.5,2,2.5)) + 
  scale_x_continuous(name = "Fractal dimension", limits = c(2, 2.7), breaks = c(2,2.2,2.4,2.6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data, aes(x = D_sum, y = R_log10, col = depth)) +
  geom_point() +
  labs(col= "Depth") +
  scale_color_gradient2(high = "red",mid = "steelblue",low = "green", midpoint = 15) +
  annotate("text", x = 2.13, y = 2.5, label = "Greater height\nrange") +
  annotate("text", x = 2.6, y = 0.1, label = "Smaller height\nrange") +
  scale_y_continuous(name = "Rugosity (log10)", limits = c(0, 2.6), breaks = c(0,0.5,1,1.5,2,2.5)) + 
  scale_x_continuous(name = "Fractal dimension", limits = c(2, 2.7), breaks = c(2,2.2,2.4,2.6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data, aes(x = D_sum, y = R_log10, col = age)) +
  geom_point() +
  labs(col= "Island Age") +
  scale_color_gradient2(high = "red",mid = "steelblue",low = "green", midpoint = 2) +
  annotate("text", x = 2.13, y = 2.5, label = "Greater height\nrange") +
  annotate("text", x = 2.6, y = 0.1, label = "Smaller height\nrange") +
  scale_y_continuous(name = "Rugosity (log10)", limits = c(0, 2.6), breaks = c(0,0.5,1,1.5,2,2.5)) + 
  scale_x_continuous(name = "Fractal dimension", limits = c(2, 2.7), breaks = c(2,2.2,2.4,2.6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Habitat intercept ranking ####################################################################################################
par(mfrow=c(1,3), mar=c(5, 6.1, 3, 3.5))

mod <- lmer(R_log10 ~ wave_max_log10 * poly(depth.m, 2) * island + (1|habitat), data)
temp <- ranef(mod)$habitat 
temp$habitat <- rownames(temp)
#temp <- temp[order(temp$`(Intercept)`),]
barplot(temp$`(Intercept)`, names=temp$habitat, main = "Rugosity (log10)", ylim = c(-0.15,0.1))
abline(h = 0)

mod <- lmer(H_log10 ~ wave_max_log10 * poly(depth.m, 2) * island + (1|habitat), data)
temp <- ranef(mod)$habitat 
temp$habitat <- rownames(temp)
#temp <- temp[order(temp$`(Intercept)`),]
barplot(temp$`(Intercept)`, names=temp$habitat, main = "Height range (log10)", ylim = c(-0.3, 0.2))
abline(h = 0)

mod <- lmer(D_sum ~ wave_max_log10 * poly(depth.m, 2) * island + (1|habitat), data)
temp <- ranef(mod)$habitat 
temp$habitat <- rownames(temp)
#temp <- temp[order(temp$`(Intercept)`),]
barplot(temp$`(Intercept)`, names=temp$habitat, main = "Fractal dimension", ylim = c(-0.06, 0.04))
abline(h = 0)

# Metric predictions with different exposure levels at each island #############################################################
# lines correspond with high, mid, and low wave energy of each island (max of each island will be different)

dd <- 1:30
mod <- lm(R_log10 ~ wave_max_log10 * poly(depth.m, 2) * island, data)

par(mfrow=c(2, 4), mar=c(6, 6.1, 3, 3.5))

for (ii in unique(data$island)) {
  temp <- data[data$island==ii,]
  plot(R_log10 ~ depth.m, temp, main=ii, col="grey", ylab = "Rugosity (log10)", xlab = "Depth, m", pch = 16)
  wmn <- mean(temp$wave_max_log10)
  wsd <- sd(temp$wave_max_log10)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn - wsd, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), ylim=c(0, 0.5), lty=2, col="blue", lwd = 2)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), lty=1, lwd = 2)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn + wsd, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), lty=2, col="red", lwd = 2)
}
legend("topright", legend = c("+1 SD", "mean exposure", "-1 SD"), lty = c(2,1,2), col = c("red","black","blue"), bty = "n")

mod <- lm(H_log10 ~ wave_max_log10 * poly(depth.m, 2) * island, data)
for (ii in unique(data$island)) {
  temp <- data[data$island==ii,]
  plot(H_log10 ~ depth.m, temp, main=ii, col="grey", ylab = "Height range (log10)", xlab = "Depth, m", pch = 16)
  wmn <- mean(temp$wave_max_log10)
  wsd <- sd(temp$wave_max_log10)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn - wsd, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), ylim=c(0, 0.5), lty=2, col="blue", lwd = 2)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), lty=1, lwd = 2)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn + wsd, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), lty=2, col="red", lwd = 2)
}

legend("topright", legend = c("+1 SD", "mean exposure", "-1 SD"), lty = c(2,1,2), col = c("red","black","blue"), bty = "n")


mod <- lm(D_sum ~ wave_max_log10 * poly(depth.m, 2) * island, data)

for (ii in unique(data$island)) {
  temp <- data[data$island==ii,]
  plot(D_sum ~ depth.m, temp, main=ii, col="grey", ylab = "Fractal dimension", xlab = "Depth, m", pch = 16)
  wmn <- mean(temp$wave_max_log10)
  wsd <- sd(temp$wave_max_log10)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn - wsd, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), ylim=c(0, 0.5), lty=2, col="blue", lwd = 2)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), lty=1, lwd = 2)
  lines(dd, predict(mod, list(wave_max_log10=rep(wmn + wsd, length(dd)), depth.m=dd, island=rep(ii, length(dd)))), lty=2, col="red", lwd = 2)
}

legend("topright", legend = c("+1 SD", "mean exposure", "-1 SD"), lty = c(2,1,2), col = c("red","black","blue"), bty = "n")

# Asner vs. Madin ##############################################################################################################

# map <- raster("habitat/comparison/ASU_GAO_Hawaii_SE_FineComplexity_2m.tif")
# mapt <- projectRaster(map, crs = crs("+proj=longlat +datum=WGS84 +units=m +no_defs"))
# writeRaster(mapP, "habitat/comparison/HawaiiSEFineRugo.tif")



data <- read.csv("Complexity_all.csv")

# mean of each site for comparison
cond <- aggregate(depth ~ site+longitude+latitude+island, data, mean)
sub <- cond[c(2,3)]

files <- list.files("habitat/comparison", (pattern = "\\.tif$"))
method <- data.frame()

for (f in 1:length(files)) {

  map <- raster(paste0("habitat/comparison/", files[f]))
  
  ex <- data.frame(sub$longitude,
                   sub$latitude,
                   raster::extract(map, sub, buffer = 10, fun = mean)) # 5m buffer around each coordinate
  
  names(ex)[1] <- "longitude"
  names(ex)[2] <- "latitude"
  names(ex)[3] <- "R_asn"
  
  compare <- merge(ex, cond, by = c("longitude","latitude"))
  compare <- na.omit(compare)
  
  if (nrow(compare) > 0) {
    compare$buffer <- 10
    method <- rbind(method, compare)
  }
}

write.csv(method, "habitat/comparison/MethodsCompare.csv")


method <- read.csv("habitat/comparison/MethodsCompare.csv")
SfM <- read.csv("habitat/comparison/SfM_SaPa.csv")

names(SfM)[2] <- "site"
names(SfM)[3] <- "R_sfm"
SfM <- SfM[-c(1)]

compare <- left_join(method, SfM)
compare <- na.omit(compare)

buf5 <- compare[compare$buffer == 5 ,]
buf10 <- compare[compare$buffer == 10 ,]

names(buf5)[3] <- "R_asn_5mbuf"
names(buf10)[3] <- "R_asn_10mbuf"

compare2 <- merge(buf5, buf10, by = c("longitude","latitude","R_sfm","site","island","depth"))
compare2$R_sfm_log10 <- log10(compare2$R_sfm)

compare2$R_sfm_test <- rank(compare2$R_sfm)/258
# correlation plot 
library(RColorBrewer)
library(corrplot)
par(mfrow = c(1,1))
full_pred <- compare2[,c("R_sfm_test","R_asn_5mbuf","R_asn_10mbuf")]
pred_cor <- cor(full_pred)
head(round(pred_cor,2)) # good because no coefficients > 0.6 or < -0.6
corrplot(pred_cor, method="color",col = COL2('RdBu', 10),tl.col = 'black',addCoef.col ='black',addgrid.col = 'white')

# correlation plot by depth 
compare2$depth_bin[compare2$depth <= 6] <- "shallow"
compare2$depth_bin[compare2$depth > 6 & compare2$depth <=18] <- "mid"
compare2$depth_bin[compare2$depth > 18] <- "deep"

shallow <- compare2[compare2$depth <= 6 , ]
mid <- compare2[compare2$depth > 6 & compare2$depth <=18 , ]
deep <- compare2[compare2$depth > 18 , ]

names(shallow)[3] <- "R_sfm_shallow"
names(shallow)[7] <- "R_asn_5mbuf_shallow"
names(shallow)[8] <- "R_asn_10mbuf_shallow"
names(mid)[3] <- "R_sfm_mid"
names(mid)[7] <- "R_asn_5mbuf_mid"
names(mid)[8] <- "R_asn_10mbuf_mid"
names(deep)[3] <- "R_sfm_deep"
names(deep)[7] <- "R_asn_5mbuf_deep"
names(deep)[8] <- "R_asn_10mbuf_deep"

test <- merge(shallow, mid, all = T)
compare.d <- merge(test, deep, all = T)

full_pred <- shallow[,c("R_sfm_shallow","R_asn_5mbuf_shallow","R_asn_10mbuf_shallow")]#,"R_sfm_mid","R_asn_5mbuf_mid","R_asn_10mbuf_mid","R_sfm_deep","R_asn_5mbuf_deep","R_asn_10mbuf_deep")]
pred_cor <- cor(full_pred)
corrplot(pred_cor, method="color",col = COL2('RdBu', 10),tl.col = 'black',addCoef.col ='black',addgrid.col = 'white')

full_pred <- mid[,c("R_sfm_mid","R_asn_5mbuf_mid","R_asn_10mbuf_mid")] #,"R_sfm_deep","R_asn_5mbuf_deep","R_asn_10mbuf_deep"
pred_cor <- cor(full_pred)
corrplot(pred_cor, method="color",col = COL2('RdBu', 10),tl.col = 'black',addCoef.col ='black',addgrid.col = 'white')

full_pred <- deep[,c("R_sfm_deep","R_asn_5mbuf_deep","R_asn_10mbuf_deep")]
pred_cor <- cor(full_pred)
corrplot(pred_cor, method="color",col = COL2('RdBu', 10),tl.col = 'black',addCoef.col ='black',addgrid.col = 'white')


# plot by buffer 
plot(log10(compare$R_sfm[compare$buffer == 5]) ~ compare$R_asn[compare$buffer == 5], pch =16)
points(log10(compare$R_sfm[compare$buffer == 10]) ~ compare$R_asn[compare$buffer == 10], pch =16, col = "hotpink")

# plot by depth ranges too (how does it compare in the shallow, etc.) 

par(mfrow = c(1,3))
plot(log10(compare.d$R_sfm_shallow) ~ compare.d$R_asn_5mbuf_shallow, pch = 16, 
     xlab = "Airborne Rugosity", ylab = "SfM Rugosity (log10)", col = "green", ylim = c(0,1), main = "Shallow")
points(log10(compare.d$R_sfm_shallow) ~ compare.d$R_asn_10mbuf_shallow, pch = 16, col = "purple")

plot(log10(compare.d$R_sfm_mid) ~ compare.d$R_asn_5mbuf_mid, pch = 17, 
     xlab = "Airborne Rugosity", ylab = "SfM Rugosity (log10)", col = "green", ylim = c(0,1), main = "Mid")
points(log10(compare.d$R_sfm_mid) ~ compare.d$R_asn_10mbuf_mid, pch = 17, col = "purple")

plot(log10(compare.d$R_sfm_deep) ~ compare.d$R_asn_5mbuf_deep, pch = 15, 
     xlab = "Airborne Rugosity", ylab = "SfM Rugosity (log10)", col = "green", ylim = c(0,1), main = "Deep")
points(log10(compare.d$R_sfm_deep) ~ compare.d$R_asn_10mbuf_deep, pch = 15, col = "purple")

legend("topright", legend = c("5m buffer","10m buffer"), pch = 16, col = c("green","purple"), bty = "n")

# mean per island 

ggplot(compare2) + 
  geom_violin(aes(x = as.factor(island), y = R_sfm_test)) +
  geom_violin(aes(x = as.factor(island), y = R_asn_10mbuf), adjust = 1) 

# mean per depth 
SfM$buffer <- "SfM"
method.con <- method[-c(1,2,5,6)]
names(method.con)[1] <- "Rugosity"
names(SfM)[2] <- "Rugosity"
SfM <- SfM[-c(3)]

plotting <- rbind(method.con, SfM)

ggplot(plotting, aes(x = depth_bin, y = Rugosity, fill = buffer)) +
  geom_violin()




# fixed effects plot for bayesian ##############################################################################################

data <- read.csv("Complexity_22sitesremoved_sub_CNet.csv")
L <- 2 #2x2m box
scl <- L / c(1, 2, 3.125, 6.25, 12.5, 25, 50, 100, 200)
L0 <- min(scl) 
data2 <- data %>% mutate(R = 0.5*log10(R_theory_mean^2 - 1),
                         D = log10(L/L0)*(D_mean-3),
                         H = log10(H_mean/(sqrt(2)*L0))) 

intercepts <- fit$summary(c("b0_a", "b0_b", "b0_c"))
output.H <- as.data.frame(fit$summary("beta_a") )
output.D <- as.data.frame(fit$summary("beta_b") )
output.R <- as.data.frame(fit$summary("beta_c") )

output.R$variable_plot <- output.D$variable_plot <- output.H$variable_plot <- c("habitat","habitat","habitat",
                            "Depth", "Depth", "Depth",
                            "Island Age", "Island Age", "Island Age",
                            "Exposure", "Exposure", "Exposure",
                            "Coral Cover", "Coral Cover", "Coral Cover",
                            "Depth x Exposure", "Depth x Exposure", "Depth x Exposure"  )
output.R$habitat <- output.D$habitat <- output.H$habitat <- c("AGR","PAV","ROB",
                                                              "AGR","PAV","ROB",
                                                              "AGR","PAV","ROB",
                                                              "AGR","PAV","ROB",
                                                              "AGR","PAV","ROB",
                                                              "AGR","PAV","ROB")
output.R$shape <- output.D$shape <- output.H$shape <- NA
output.R$width <- output.D$width <- output.H$width <- NA
output.R$color <- output.D$color <- output.H$color <- NA
for (o in 1:nrow(output.R)) {
  if(sign(output.R[o,"q5"]) == sign(output.R[o,"q95"])) {
    output.R$width[o] <- "TRUE"
    output.R$shape[o] <- 20}
  else {
    output.R$width[o] <- "FALSE" 
    output.R$shape[o] <- 20}
}
for (o in 1:nrow(output.H)) {
  if(sign(output.H[o,"q5"]) == sign(output.H[o,"q95"])) {
    output.H$width[o] <- "TRUE"
    output.H$shape[o] <- 20}
  else {
    output.H$width[o] <- "FALSE" 
    output.H$shape[o] <- 20}
}
for (o in 1:nrow(output.D)) {
  if(sign(output.D[o,"q5"]) == sign(output.D[o,"q95"])) {
    output.D$width[o] <- "TRUE"
    output.D$shape[o] <- 20}
  else {
    output.D$width[o] <- "FALSE" 
    output.D$shape[o] <- 20}
}

output.R$color[output.R$habitat == "AGR"] <- "#D55E00"
output.R$color[output.R$habitat == "PAV"] <- "#0072B2"
output.R$color[output.R$habitat == "ROB"] <- "#009E73"
output.H$color[output.H$habitat == "AGR"] <- "#D55E00"
output.H$color[output.H$habitat == "PAV"] <- "#0072B2"
output.H$color[output.H$habitat == "ROB"] <- "#009E73"
output.D$color[output.D$habitat == "AGR"] <- "#D55E00"
output.D$color[output.D$habitat == "PAV"] <- "#0072B2"
output.D$color[output.D$habitat == "ROB"] <- "#009E73"

ggplot(output.R, aes(x = variable_plot, y = mean, color = color, group = color)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, color = "black", position = position_dodge(width = .5)) +
  geom_point(shape = output.R$shape, size = 4, col = output.R$color, position = position_dodge(width = .5)) +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("\nParameter Estimate") +
  ggtitle("Rugosity")
ggplot(output.H, aes(x = variable_plot, y = mean, color = color, group = color)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, color = "black", position = position_dodge(width = .5)) +
  geom_point(shape = output.H$shape, size = 4, col = output.H$color, position = position_dodge(width = .5)) +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("\nParameter Estimate") +
  ggtitle("Height range")
ggplot(output.D, aes(x = variable_plot, y = mean, color = color, group = color)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, color = "black", position = position_dodge(width = .5)) +
  geom_point(shape = output.D$shape, size = 4, col = output.D$color, position = position_dodge(width = .5)) +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("\nParameter Estimate") +
  ggtitle("Fractal dimension")
  

output.H.sub <- output.H[-c(1:3) ,]
output.R.sub <- output.R[-c(1:3) ,]
output.D.sub <- output.D[-c(1:3) ,]

ggplot(output.R.sub, aes(x = variable_plot, y = mean, color = color, group = color)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, color = "white", position = position_dodge(width = .5)) +
  geom_point(shape = output.R.sub$shape, size = 4, col = output.R.sub$color, position = position_dodge(width = .5)) +
  coord_flip() +
  dark_theme_classic() +
  xlab("") +
  ylab("\nParameter Estimate") +
  ggtitle("Rugosity")
ggplot(output.H.sub, aes(x = variable_plot, y = mean, color = color, group = color)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, color = "white", position = position_dodge(width = .5)) +
  geom_point(shape = output.H.sub$shape, size = 4, col = output.H.sub$color, position = position_dodge(width = .5)) +
  coord_flip() +
  dark_theme_classic() +
  xlab("") +
  ylab("\nParameter Estimate") +
  ggtitle("Height range")
ggplot(output.D.sub, aes(x = variable_plot, y = mean, color = color, group = color)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, color = "white", position = position_dodge(width = .5)) +
  geom_point(shape = output.D.sub$shape, size = 4, col = output.D.sub$color, position = position_dodge(width = .5)) +
  coord_flip() +
  dark_theme_classic() +
  xlab("") +
  ylab("\nParameter Estimate") +
  ggtitle("Fractal dimension")

# coral vs other variables #####################################################################################################

hab <- c("AGR","PAV","ROB")

temp1 <- data[data$habitat== hab[1],]
temp2 <- data[data$habitat== hab[2],]
temp3 <- data[data$habitat== hab[3],]

cor.dep <- ggplot() +
  geom_point(aes(x = data2$depth, y = data2$CORAL), size = .2) +
  geom_smooth(aes(x = temp1$depth, y = temp1$CORAL), method = "lm", col = "#D55E00", fill = "grey") +
  geom_smooth(aes(x = temp2$depth, y = temp2$CORAL), method = "lm", col = "#0072B2", fill = "grey") +
  geom_smooth(aes(x = temp3$depth, y = temp3$CORAL), method = "lm", col = "#009E73", fill = "grey") +
  labs(x = "Depth, m", y = "Coral Cover, %") +
  coord_cartesian(ylim = c(0,40)) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 12), axis.title = element_text(size = 10)) 

cor.age <- ggplot() +
  geom_point(aes(x = data2$age, y = data2$CORAL), size = .2) +
  geom_smooth(aes(x = temp1$age, y = temp1$CORAL), method = "lm", col = "#D55E00", fill = "grey") +
  geom_smooth(aes(x = temp2$age, y = temp2$CORAL), method = "lm", col = "#0072B2", fill = "grey") +
  geom_smooth(aes(x = temp3$age, y = temp3$CORAL), method = "lm", col = "#009E73", fill = "grey") +
  labs(x = "Geological reef age, mya", y = "") +
  coord_cartesian(ylim = c(0,40)) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 12), axis.title = element_text(size = 10)) 

cor.wave <- ggplot() +
  geom_point(aes(x = data2$wave_max_log10, y = data2$CORAL), size = .2) +
  geom_smooth(aes(x = temp1$wave_max_log10, y = temp1$CORAL), method = "lm", col = "#D55E00", fill = "grey") +
  geom_smooth(aes(x = temp2$wave_max_log10, y = temp2$CORAL), method = "lm", col = "#0072B2", fill = "grey") +
  geom_smooth(aes(x = temp3$wave_max_log10, y = temp3$CORAL), method = "lm", col = "#009E73", fill = "grey") +
  labs(x = "Exposure, kW/m (log10)", y = "") +
  coord_cartesian(ylim = c(0,40)) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 12), axis.title = element_text(size = 10)) 

grid.arrange(arrangeGrob(cor.dep + ggtitle("(a)") + theme(plot.title = element_text(size = 12)), 
                         cor.age + ggtitle("(b)") + theme(plot.title = element_text(size = 12)), 
                         cor.wave + ggtitle("(c)") + theme(plot.title = element_text(size = 12)),
                         nrow = 1, ncol = 3), 
             mylegend,
             heights = c(10,1))


# Predicted complexity rasters #################################################################################################

R.pred <- readAll(raster("habitat/predicted_rugosity.tif"))
D.pred <- readAll(raster("habitat/predicted_fractaldimension.tif"))
H.pred <- readAll(raster("habitat/predicted_heightrange.tif"))

par(mfrow = c(3,1), family = "Times")
plot(R.pred, col =  rev(RColorBrewer::brewer.pal(10, "RdYlBu")),  
     legend.args = list(text = "Rugosity", family = "Times", cex = .8), ylim = c(21.28, 21.32),
     family = "Times")
text(-158, 21.325, labels = "(a)", xpd = NA, cex = 1.2)
plot(D.pred, col =  rev(RColorBrewer::brewer.pal(10, "RdYlBu")),  
     legend.args = list(text = "Fractal \ndimension", family = "Times", cex = .8), ylim = c(21.28, 21.32),
     family = "Times")
text(-158, 21.325, labels = "(b)", xpd = NA, cex = 1.2)
plot(H.pred, col =  rev(RColorBrewer::brewer.pal(10, "RdYlBu")),  
     legend.args = list(text = "Height \nrange", family = "Times", cex = .8), ylim = c(21.28, 21.32),
     family = "Times")
text(-158, 21.325, labels = "(c)", xpd = NA, cex = 1.2)

