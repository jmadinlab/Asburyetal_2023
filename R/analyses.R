

##### GEOLOGICAL AND ENVIRONMENTAL FACTORS SHAPE REEF #####
##### HABITAT STRUCTURE IN THE MAIN HAWAIIAN ISLANDS #####
################### ASBURY ET AL. 2023 ###################

library(ggplot2)
library(pracma)
library(imager)
library(survey)
library(gridExtra)
library(plyr)
library(dplyr)
library(jtools)
library(plot3D)
library(basemaps)
library(fruitr)
library(ggfortify)
source("R/functions.R")

### DATA MANAGEMENT ############################################################################################################
# load data
data <- read.csv("output/MHI2019_Complexity.csv")

# metric derivations
L <- 2 #2x2m box
scl <- L / c(1, 2, 3.125, 6.25, 12.5, 25, 50, 100, 200)
L0 <- min(scl) 
data2 <- data %>% mutate(R = 0.5*log10(rugosity_theory^2 - 1),
                         D = log10(L/L0)*(fractal-3),
                         H = log10(height/(sqrt(2)*L0))) 

# standardize & center predictors to put all on same scale 
data2$depth_s <- as.vector(scale(data2$depth))
data2$wave_max_log10_s <- as.vector(scale(log10(data2$wave_max)))
data2$age_s <- as.vector(scale(data2$age))
data2$coral_s <- as.vector(scale(data2$CORAL))

# survey weights
sectors <- read.csv("data/Sectors-Strata-Areas.csv")
wtemp <- left_join(data2, sectors)

weight <- plyr::ddply(wtemp, .(SEC_NAME,DEPTH_BIN,NH), summarize, n = length(unique(site)))
weight$survey_weight <- weight$NH / weight$n

data3 <- left_join(data2, weight)

data3$conc <- c(paste0(data3$island, "_", data3$SEC_NAME, "_", data3$DEPTH_BIN, "_", data3$site))
des <- svydesign(id= ~1, strata = ~ conc, weights = ~survey_weight, data = data3) 


#### FIGURES ###################################################################################################################

# function to extract legend from a plot
g_legend<-function(a.gplot){ 
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

#### FIGURE 1 ####
ll <- cbind(data3$longitude, data3$latitude)
merc <- lnglat_to_xy(ll)

ext <- draw_ext()
map <- ggplot() + 
  basemap_gglayer(ext, map_service = "esri", map_type = "world_imagery") +
  scale_fill_identity() + 
  coord_sf() +
  geom_point(aes(merc[,1], merc[,2]), size = 1, 
             shape = 21, color = "black", fill = "yellow") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill=NA, size = 2)) 

#### FIGURE 2 ####
# habitat examples
agr.pic <- load.image("figs/AGR.jpg")
pav.pic <- load.image("figs/PAV.jpg")
rob.pic <- load.image("figs/ROB.jpg")

# organize data
mdat <- data3[c("R","H")]
mdat$m <- 1

# best fit orthogonal plane 
A <- data.matrix(mdat)
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
#png("figs/figure2_habitats.png", width = 1.5, height = 3.5, units = "in", res = 600)
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

#png("figs/figure2_all.png", width = 3, height = 3.5, units = "in", res = 600)
par(mar=c(1.2, 2, 1, .5) +0.1, family = "Times", ps = 10)

scatter3D(data2$R, data2$H, data2$D, pch = 20, cex = 1, 
          col= rev(hcl.colors(100, "grays")), 
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
          #colkey=list(side = 2, length = 0.5, width = 1, line.clab = 1, dist=0),
          clab=expression(italic(D)),
          colkey=FALSE)
points3D(R_mean$R - 0.2, H_mean$H + 0.2, D_mean$D,
         col = c("#009E73","#0072B2","#D55E00"), cex = 2, pch = 20,
         surf=NULL,
         colkey=FALSE,
         add=TRUE)
mtext("(d)", side = 3, at = -.5)
text3D(-0.6, 2.9, -1.2, labels=expression(italic(r)^2 == 0.974), surf = NULL, add = TRUE, family = "Times")

dev.off()

#png("figs/figure2_ByHabitat.png", width = 1.5, height = 3.5, units = "in", res = 600)
par(mar=c(.6, .5, 1, .6) +0.1, family = "Times", ps = 10, mfrow = c(3,1))
scatter3D(dataAGR$R, dataAGR$H, dataAGR$D, pch = 20, cex = .8, 
          col = rev(AGRcolors),
          xlim=c(min(data2$R)-0.25, max(data2$R)+0.25), 
          ylim=c(min(data2$H)-0.25, max(data2$H)+0.25), 
          zlim=c(min(data2$D)-0.05, max(data2$D)+0.05), 
          ylab="", 
          xlab="", 
          zlab="", 
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
          ylab="", 
          xlab="", 
          zlab="", 
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
          ylab="", 
          xlab="", 
          zlab="", 
          surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.01), 
                    fitpoints=dataROB$D), 
          theta=215, 
          phi=0,
          colkey=FALSE)
mtext("(g)", side = 3, at = -.5)

dev.off()

#### FIGURE 3 ####

data.pca <- data3[c("R","D","H","wave_max_log10_s","depth_s","age_s","coral_s","habitat")]
names(data.pca)[1] <- "Rugosity"
names(data.pca)[2] <- "Fractal dimension"
names(data.pca)[3] <- "Height range"
names(data.pca)[4] <- "Exposure"
names(data.pca)[5] <- "Depth"
names(data.pca)[6] <- "Geological Reef Age"
names(data.pca)[7] <- "Coral Cover"

pca.var <- prcomp(data.pca[1:7], center = TRUE, scale = TRUE) 
summary(pca.var)

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
      arrow_ends[a,b] <- value + 0.008 }
    if (value < 0) {
      arrow_ends[a,b] <- value - 0.006 }}}


g2 <- autoplot(pca.var, data.pca, loadings = TRUE, loadings.col = "black", col = "habitat",
               loadings.label.size = 2, loadings.label.col = "black", frame = TRUE, frame.type = "norm") +
  scale_color_manual(values = alpha(c("#D55E00","#0072B2","#009E73"), 0.4)) +
  scale_fill_manual(values = alpha(c("#D55E00","#0072B2","#009E73"), 0.4)) +
  geom_label(data = arrow_ends, aes(xend, yend), label = names(pca.var$center), size = 2, label.size = 0,family = "Times-Roman") +
  labs(colour = "Habitat Type") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") 

#png("figs/figure3.png", width = 5, height = 4, unit= "in", res = 600)
grid.arrange(arrangeGrob(g2), 
             mylegend, 
             heights = c(10,1))
dev.off()

#### FIGURE 4 ####

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

# Height range 
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
fig4 <- grid.arrange(arrangeGrob(R.DEP + ggtitle("(a)"), R.AGE + ggtitle("(b)"), R.WAVE + ggtitle("(c)"), R.COR + ggtitle("(d)"),
                                 D.DEP + ggtitle("(e)"), D.AGE + ggtitle("(f)"), D.WAVE + ggtitle("(g)"), D.COR + ggtitle("(h)"),
                                 H.DEP + ggtitle("(i)"), H.AGE + ggtitle("(j)"), H.WAVE + ggtitle("(k)"), H.COR + ggtitle("(l)"),
                                 nrow = 3, ncol = 4), 
                     mylegend, 
                     heights = c(10,1))
#ggsave("figs/figure4.png", plot = fig4, width = 7, height = 5, units = c("in"), dpi = 600)


