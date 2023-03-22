library(ggplot2)
library(survey)
library(gridExtra)
library(plyr)
library(dplyr)
library(jtools)
library(kableExtra)
library(latex2exp)
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

#### FIGURES & TABLES ##########################################################################################################

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
#### TABLE S2 ####
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

# relationships of components
var <- pca.var$sdev ^ 2
varPer <- var/sum(var)*100
load <- pca.var$rotation
cut <- 1/7 
pc1 <- as.numeric(load[,1]) * 100
pc1[8] <- varPer[1]
pc2 <- as.numeric(load[,2]) * 100
pc2[8] <- varPer[2]
pc3 <- as.numeric(load[,3]) * 100
pc3[8] <- varPer[3]

sort <- data.frame(row.names = c("Rugosity","Fractal dimension","Height range","Exposure","Depth","Geological reef age","Coral cover","Total variance"))
sort$PC1 <- format(pc1, digits = 3)
sort$PC2 <- format(pc2, digits = 3)
sort$PC3 <- format(pc3, digits = 3)
sort$PC1[c(1,3,4,6,7)] <- cell_spec(sort$PC1[c(1,3,4,6,7)], bold =  T)
sort$PC2[c(2,3,4,7)] <- cell_spec(sort$PC2[c(2,3,4,7)], bold =  T)
sort$PC3[c(1,2,4,5)] <- cell_spec(sort$PC3[c(1,3,4,5)], bold =  T)

sort %>% 
  kbl(escape = F) %>%
  kable_classic(html_font = "Times New Roman") %>%
  row_spec(7, extra_css = "border-bottom: 2px solid;")

#### FIGURE S1 ####

full_pred <- data3[,c("depth_s",
                      "age_s",
                      "wave_max_log10_s",
                      "coral_s")]

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}


#png("figs/figureS1.png", width = 5, height = 4, units = "in", res = 600)
par(mar=c(1.2, 2, 1, 1) +0.1, mfrow = c(2,2), family = "Times", ps = 10)

pairs(full_pred,
      upper.panel = panel.cor,
      labels = c("Depth, m", "Geological Reef Age, mya", "Exposure, kW/m (log10)", "Coral Cover, %"),
      cex = .8)

dev.off()

#### FIGURE S2 ####
mod.sw.H2 <-svyglm(H ~ depth_s + age_s + wave_max_log10_s + coral_s + depth_s*wave_max_log10_s, design=des)
mod.sw.R2 <-svyglm(R ~ depth_s + age_s + wave_max_log10_s + coral_s + depth_s*wave_max_log10_s, design=des)
mod.sw.D2 <-svyglm(D ~ depth_s + age_s + wave_max_log10_s + coral_s + depth_s*wave_max_log10_s, design=des)

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

newdata <- data3[c("depth_s","age_s","wave_max_log10_s","coral_s")]
newdata$depth_s <- mean(newdata$depth_s)
newdata$age_s <- mean(newdata$age_s)
newdata$wave_max_log10_s <- mean(newdata$wave_max_log10_s)
newdata$coral_s <- mean(newdata$coral_s)

D.nd <- newdata
D.nd$depth_s <- seq(min(data3$depth_s), 
                    max(data3$depth_s),
                    length = nrow(newdata))
pred_HD <- predict(mod.sw.H2, newdata = D.nd, type = "response", se.fit = TRUE)
pred_RD <- predict(mod.sw.R2, newdata = D.nd, type = "response", se.fit = TRUE)
pred_DD <- predict(mod.sw.D2, newdata = D.nd, type = "response", se.fit = TRUE)
pred_HD <- as.data.frame(pred_HD)
pred_RD <- as.data.frame(pred_RD)
pred_DD <- as.data.frame(pred_DD)
colnames(pred_HD) <- colnames(pred_RD) <- colnames(pred_DD) <- c("Predicted","SE")
pred_HD <- cbind(D.nd, pred_HD)
pred_RD <- cbind(D.nd, pred_RD)
pred_DD <- cbind(D.nd, pred_DD)
pred_HD$Predict.lwr <- pred_HD$Predicted - 1.96 * pred_HD$SE # confidence interval upper bound
pred_HD$Predict.upr <- pred_HD$Predicted + 1.96 * pred_HD$SE # confidence interval lower bound
pred_RD$Predict.lwr <- pred_RD$Predicted - 1.96 * pred_RD$SE # confidence interval upper bound
pred_RD$Predict.upr <- pred_RD$Predicted + 1.96 * pred_RD$SE # confidence interval lower bound
pred_DD$Predict.lwr <- pred_DD$Predicted - 1.96 * pred_DD$SE # confidence interval upper bound
pred_DD$Predict.upr <- pred_DD$Predicted + 1.96 * pred_DD$SE # confidence interval lower bound

A.nd <- newdata
A.nd$age_s <- seq(min(data3$age_s), 
                  max(data3$age_s),
                  length = nrow(newdata))
pred_HA <- predict(mod.sw.H2, newdata = A.nd, type = "response", se.fit = TRUE)
pred_RA <- predict(mod.sw.R2, newdata = A.nd, type = "response", se.fit = TRUE)
pred_DA <- predict(mod.sw.D2, newdata = A.nd, type = "response", se.fit = TRUE)
pred_HA <- as.data.frame(pred_HA)
pred_RA <- as.data.frame(pred_RA)
pred_DA <- as.data.frame(pred_DA)
colnames(pred_HA) <- colnames(pred_RA) <- colnames(pred_DA) <- c("Predicted","SE")
pred_HA <- cbind(A.nd, pred_HA)
pred_RA <- cbind(A.nd, pred_RA)
pred_DA <- cbind(A.nd, pred_DA)
pred_HA$Predict.lwr <- pred_HA$Predicted - 1.96 * pred_HA$SE # confidence interval upper bound
pred_HA$Predict.upr <- pred_HA$Predicted + 1.96 * pred_HA$SE # confidence interval lower bound
pred_RA$Predict.lwr <- pred_RA$Predicted - 1.96 * pred_RA$SE # confidence interval upper bound
pred_RA$Predict.upr <- pred_RA$Predicted + 1.96 * pred_RA$SE # confidence interval lower bound
pred_DA$Predict.lwr <- pred_DA$Predicted - 1.96 * pred_DA$SE # confidence interval upper bound
pred_DA$Predict.upr <- pred_DA$Predicted + 1.96 * pred_DA$SE # confidence interval lower bound

W.nd <- newdata
W.nd$wave_max_log10_s <- seq(min(data3$wave_max_log10_s), 
                             max(data3$wave_max_log10_s),
                             length = nrow(newdata))
pred_HW <- predict(mod.sw.H2, newdata = W.nd, type = "response", se.fit = TRUE)
pred_RW <- predict(mod.sw.R2, newdata = W.nd, type = "response", se.fit = TRUE)
pred_DW <- predict(mod.sw.D2, newdata = W.nd, type = "response", se.fit = TRUE)
pred_HW <- as.data.frame(pred_HW)
pred_RW <- as.data.frame(pred_RW)
pred_DW <- as.data.frame(pred_DW)
colnames(pred_HW) <- colnames(pred_RW) <- colnames(pred_DW) <- c("Predicted","SE")
pred_HW <- cbind(W.nd, pred_HW)
pred_RW <- cbind(W.nd, pred_RW)
pred_DW <- cbind(W.nd, pred_DW)
pred_HW$Predict.lwr <- pred_HW$Predicted - 1.96 * pred_HW$SE # confidence interval upper bound
pred_HW$Predict.upr <- pred_HW$Predicted + 1.96 * pred_HW$SE # confidence interval lower bound
pred_RW$Predict.lwr <- pred_RW$Predicted - 1.96 * pred_RW$SE # confidence interval upper bound
pred_RW$Predict.upr <- pred_RW$Predicted + 1.96 * pred_RW$SE # confidence interval lower bound
pred_DW$Predict.lwr <- pred_DW$Predicted - 1.96 * pred_DW$SE # confidence interval upper bound
pred_DW$Predict.upr <- pred_DW$Predicted + 1.96 * pred_DW$SE # confidence interval lower bound

C.nd <- newdata
C.nd$coral_s <- seq(min(data3$coral_s), 
                    max(data3$coral_s),
                    length = nrow(newdata))
pred_HC <- predict(mod.sw.H2, newdata = C.nd, type = "response", se.fit = TRUE)
pred_RC <- predict(mod.sw.R2, newdata = C.nd, type = "response", se.fit = TRUE)
pred_DC <- predict(mod.sw.D2, newdata = C.nd, type = "response", se.fit = TRUE)
pred_HC <- as.data.frame(pred_HC)
pred_RC <- as.data.frame(pred_RC)
pred_DC <- as.data.frame(pred_DC)
colnames(pred_HC) <- colnames(pred_RC) <- colnames(pred_DC) <- c("Predicted","SE")
pred_HC <- cbind(C.nd, pred_HC)
pred_RC <- cbind(C.nd, pred_RC)
pred_DC <- cbind(C.nd, pred_DC)
pred_HC$Predict.lwr <- pred_HC$Predicted - 1.96 * pred_HC$SE # confidence interval upper bound
pred_HC$Predict.upr <- pred_HC$Predicted + 1.96 * pred_HC$SE # confidence interval lower bound
pred_RC$Predict.lwr <- pred_RC$Predicted - 1.96 * pred_RC$SE # confidence interval upper bound
pred_RC$Predict.upr <- pred_RC$Predicted + 1.96 * pred_RC$SE # confidence interval lower bound
pred_DC$Predict.lwr <- pred_DC$Predicted - 1.96 * pred_DC$SE # confidence interval upper bound
pred_DC$Predict.upr <- pred_DC$Predicted + 1.96 * pred_DC$SE # confidence interval lower bound

H.DEP2 <- ggplot() +
  geom_point(aes(x=data3$depth_s,
                 y=data3$H),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_HD$depth_s, 
                  ymin=pred_HD$Predict.lwr, 
                  ymax=pred_HD$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_HD$depth_s, 
                y=pred_HD$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Depth, m", y = "Height range") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (6-dep_mn)/dep_sd,
                                (12-dep_mn)/dep_sd, 
                                (18-dep_mn)/dep_sd, 
                                (24-dep_mn)/dep_sd), 
                     labels = c(0,6,12,18,24))
H.AGE2 <- ggplot() +
  geom_point(aes(x=data3$age_s,
                 y=data3$H),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_HA$age_s, 
                  ymin=pred_HA$Predict.lwr, 
                  ymax=pred_HA$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_HA$age_s, 
                y=pred_HA$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Geological reef age, mya", y = "") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-age_mn)/age_sd, 
                                (1-age_mn)/age_sd, 
                                (2-age_mn)/age_sd, 
                                (3-age_mn)/age_sd, 
                                (4-age_mn)/age_sd, 
                                (5-age_mn)/age_sd),
                     labels = c(0,1,2,3,4,5))
H.WAVE2 <- ggplot() +
  geom_point(aes(x=data3$wave_max_log10_s,
                 y=data3$H),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_HW$wave_max_log10_s, 
                  ymin=pred_HW$Predict.lwr, 
                  ymax=pred_HW$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_HW$wave_max_log10_s, 
                y=pred_HW$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Exposure, kW/m (log10)", y = "") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) 
H.COR2 <- ggplot() +
  geom_point(aes(x=data3$coral_s,
                 y=data3$H),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_HC$coral_s, 
                  ymin=pred_HC$Predict.lwr, 
                  ymax=pred_HC$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_HC$coral_s, 
                y=pred_HC$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Coral cover, %", y = "") +
  ylim(0.8, 2.6) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))
R.DEP2 <- ggplot() +
  geom_point(aes(x=data3$depth_s,
                 y=data3$R),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_RD$depth_s, 
                  ymin=pred_RD$Predict.lwr, 
                  ymax=pred_RD$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_RD$depth_s, 
                y=pred_RD$Predicted), size = 1, color = "black",
            linetype = "dotted") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "Rugosity") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (6-dep_mn)/dep_sd,
                                (12-dep_mn)/dep_sd, 
                                (18-dep_mn)/dep_sd, 
                                (24-dep_mn)/dep_sd), 
                     labels = c(0,6,12,18,24))
R.AGE2 <- ggplot() +
  geom_point(aes(x=data3$age_s,
                 y=data3$R),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_RA$age_s, 
                  ymin=pred_RA$Predict.lwr,
                  ymax=pred_RA$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_RA$age_s, 
                y=pred_RA$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-age_mn)/age_sd, 
                                (1-age_mn)/age_sd, 
                                (2-age_mn)/age_sd, 
                                (3-age_mn)/age_sd, 
                                (4-age_mn)/age_sd, 
                                (5-age_mn)/age_sd),
                     labels = c(0,1,2,3,4,5))
R.WAVE2 <- ggplot() +
  geom_point(aes(x=data3$wave_max_log10_s,
                 y=data3$R),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_RW$wave_max_log10_s, 
                  ymin=pred_RW$Predict.lwr, 
                  ymax=pred_RW$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_RW$wave_max_log10_s, 
                y=pred_RW$Predicted), size = 1, color = "black",
            linetype = "dotted") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) 
R.COR2 <- ggplot() +
  geom_point(aes(x=data3$coral_s,
                 y=data3$R),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_RC$coral_s, 
                  ymin=pred_RC$Predict.lwr, 
                  ymax=pred_RC$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_RC$coral_s, 
                y=pred_RC$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "") +
  ylim(-0.8,0.85) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))
D.DEP2 <- ggplot() +
  geom_point(aes(x=data3$depth_s,
                 y=data3$D),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_DD$depth_s, 
                  ymin=pred_DD$Predict.lwr, 
                  ymax=pred_DD$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_DD$depth_s, 
                y=pred_DD$Predicted), size = .8, color = "black") +
  theme_classic() +
  ylim(-2.2, -1.1) +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "Fractal dimension") +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (6-dep_mn)/dep_sd,
                                (12-dep_mn)/dep_sd, 
                                (18-dep_mn)/dep_sd, 
                                (24-dep_mn)/dep_sd), 
                     labels = c(0,6,12,18,24))
D.AGE2 <- ggplot() +
  geom_point(aes(x=data3$age_s,
                 y=data3$D),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_DA$age_s, 
                  ymin=pred_DA$Predict.lwr, 
                  ymax=pred_DA$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_DA$age_s, 
                y=pred_DA$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-age_mn)/age_sd, 
                                (1-age_mn)/age_sd, 
                                (2-age_mn)/age_sd, 
                                (3-age_mn)/age_sd, 
                                (4-age_mn)/age_sd, 
                                (5-age_mn)/age_sd),
                     labels = c(0,1,2,3,4,5))
D.WAVE2 <- ggplot() +
  geom_point(aes(x=data3$wave_max_log10_s,
                 y=data3$D),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_DW$wave_max_log10_s, 
                  ymin=pred_DW$Predict.lwr, 
                  ymax=pred_DW$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_DW$wave_max_log10_s, 
                y=pred_DW$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) 
D.COR2 <- ggplot() +
  geom_point(aes(x=data3$coral_s,
                 y=data3$D),
             color = "grey", 
             alpha = 0.3,
             shape = 1) +
  geom_ribbon(aes(x=pred_DC$coral_s, 
                  ymin=pred_DC$Predict.lwr, 
                  ymax=pred_DC$Predict.upr), fill = "grey") +
  geom_line(aes(x=pred_DC$coral_s, 
                y=pred_DC$Predicted), size = .8, color = "black") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "", y = "") +
  ylim(-2.2, -1.1) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))


sup2 <- grid.arrange(arrangeGrob(R.DEP2 + ggtitle("(a)"), R.AGE2 + ggtitle("(b)"), R.WAVE2 + ggtitle("(c)"), R.COR2 + ggtitle("(d)"),
                                 D.DEP2 + ggtitle("(e)"), D.AGE2 + ggtitle("(f)"), D.WAVE2 + ggtitle("(g)"), D.COR2 + ggtitle("(h)"),
                                 H.DEP2 + ggtitle("(i)"), H.AGE2 + ggtitle("(j)"), H.WAVE2 + ggtitle("(k)"), H.COR2 + ggtitle("(l)"),
                                 nrow = 3, ncol = 4))

#ggsave("figs/figureS2.png", plot = sup2, width = 7, height = 5, units = c("in"))

#### FIGURE S3 ####
mod.sw.H <- svyglm(H ~ 0 + habitat*depth_s + habitat*age_s + habitat*wave_max_log10_s + habitat*coral_s + 
                     habitat*depth_s*wave_max_log10_s, design=des)
mod.sw.R <- svyglm(R ~ 0 + habitat*depth_s + habitat*age_s + habitat*wave_max_log10_s + habitat*coral_s + 
                     habitat*depth_s*wave_max_log10_s, design=des)
mod.sw.D <- svyglm(D ~ 0 + habitat*depth_s + habitat*age_s + habitat*wave_max_log10_s + habitat*coral_s + 
                     habitat*depth_s*wave_max_log10_s, design=des)
v1 = "depth_s" 
v2 = "wave_max_log10_s" 
d = data3 

v1_seq <- seq(from = min(d[[v1]]), to = max(d[[v1]]), length = 100)
v2_seq <- seq(from = min(d[[v2]]), to = max(d[[v2]]), length = 100)

# group by habitat 
varsA = varsP = varsR = c("habitat","depth_s","age_s","wave_max_log10_s","coral_s")
vars2A = vars2P = vars2R = c("habitat","age_s","coral_s")

# generate new data frame to hold unique values of variable of interest and mean of other predictors
new_tempA <- setNames(data.frame(matrix(ncol = length(varsA), nrow = 100*100)), varsA)
new_tempA[[v1]] <- rep(v1_seq, each = 100)
new_tempA[[v2]] <- rep(v2_seq, 100)
new_tempP <- setNames(data.frame(matrix(ncol = length(varsP), nrow = 100*100)), varsP)
new_tempP[[v1]] <- rep(v1_seq, each = 100)
new_tempP[[v2]] <- rep(v2_seq, 100)
new_tempR <- setNames(data.frame(matrix(ncol = length(varsR), nrow = 100*100)), varsR)
new_tempR[[v1]] <- rep(v1_seq, each = 100)
new_tempR[[v2]] <- rep(v2_seq, 100)

dA <- data3[data3$habitat == "AGR",]
dP <- data3[data3$habitat == "PAV",]
dR <- data3[data3$habitat == "ROB",]
# calculate means of other values of interest
for(i in c(2:length(vars2A))){
  col.name <- vars2A[i]
  new_tempA[[col.name]] <- mean(dA[[col.name]])
  new_tempP[[col.name]] <- mean(dP[[col.name]])
  new_tempR[[col.name]] <- mean(dR[[col.name]])}

new_tempA$habitat <- "AGR"
new_tempP$habitat <- "PAV"
new_tempR$habitat <- "ROB"

new_tempA_H <- new_tempA_R <- new_tempA_D <- new_tempA
new_tempP_H <- new_tempP_R <- new_tempP_D <- new_tempP
new_tempR_H <- new_tempR_R <- new_tempR_D <- new_tempR

new_tempA_H$pred <- predict(mod.sw.H, newdata = new_tempA)   # add predicted values
new_tempP_H$pred <- predict(mod.sw.H, newdata = new_tempP)
new_tempR_H$pred <- predict(mod.sw.H, newdata = new_tempR)
new_tempA_R$pred <- predict(mod.sw.R, newdata = new_tempA)   # add predicted values
new_tempP_R$pred <- predict(mod.sw.R, newdata = new_tempP)
new_tempR_R$pred <- predict(mod.sw.R, newdata = new_tempR)
new_tempA_D$pred <- predict(mod.sw.D, newdata = new_tempA)   # add predicted values
new_tempP_D$pred <- predict(mod.sw.D, newdata = new_tempP)
new_tempR_D$pred <- predict(mod.sw.D, newdata = new_tempR)

chA <- chull(dA[[v1]], dA[[v2]])
chA <- c(chA, chA[1])
for (p in 1:nrow(new_tempA)) {
  pp <- point.in.polygon(new_tempA[[v1]][p], new_tempA[[v2]][p],dA[[v1]][chA],dA[[v2]][chA])
  if (pp == 0) {
    new_tempA_H$pred[p] <- NA
    new_tempA_R$pred[p] <- NA
    new_tempA_D$pred[p] <- NA} }

chP <- chull(dP[[v1]], dP[[v2]])
chP <- c(chP, chP[1])
for (p in 1:nrow(new_tempP)) {
  pp <- point.in.polygon(new_tempP[[v1]][p], new_tempP[[v2]][p],dP[[v1]][chP],dP[[v2]][chP])
  if (pp == 0) {
    new_tempP_H$pred[p] <- NA
    new_tempP_R$pred[p] <- NA
    new_tempP_D$pred[p] <- NA} }

chR <- chull(dR[[v1]], dR[[v2]])
chR <- c(chR, chR[1])
for (p in 1:nrow(new_tempR)) {
  pp <- point.in.polygon(new_tempR[[v1]][p], new_tempR[[v2]][p],dR[[v1]][chR],dR[[v2]][chR])
  if (pp == 0) {
    new_tempR_H$pred[p] <- NA
    new_tempR_R$pred[p] <- NA
    new_tempR_D$pred[p] <- NA} }

AH <- ggplot(new_tempA_H, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dA, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        legend.position = "none",
        plot.margin = margin(0,0,0,1,"cm")) +
  labs(y = "Exposure, kW/m (log10)", x = "Depth, m") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Height \n range")), 
                      high = "#B55000",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.160", size = 3) +
  annotate(geom="text", x = -4.2, y = -.5, label = "Height range", angle = 90, family = "Times-Roman", fontface = "bold", size = 5) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd))
PH <- ggplot(new_tempP_H, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dP, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        legend.position = "none",
        plot.margin = margin(0,0.5,0,0.5,"cm")) +
  labs(y = "", x = "Depth, m") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Height \n range")), 
                      high = "#005B8E",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.002", size = 3) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd))
RH <- ggplot(new_tempR_H, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dR, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        legend.position = "none",
        plot.margin = margin(0,1,0,0,"cm")) +
  labs(y = "", x = "Depth, m") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Height \n range")), 
                      high = "#006A4D",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.006", size = 3) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd)) 
AR <- ggplot(new_tempA_R, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dA, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0,0,0,1,"cm")) +
  labs(y = "Exposure, kW/m (log10)", x = "", title = "Aggregate reef") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Rugosity")), 
                      high = "#B55000",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.129", size = 3) +
  annotate(geom="text", x = -4.2, y = -.5, label = "Rugosity", angle = 90, family = "Times-Roman", fontface = "bold", size = 5) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd))
PR <- ggplot(new_tempP_R, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dP, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0,0.5,0,0.5,"cm")) +
  labs(y = "", x = "", title = "Pavement") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Rugosity")), 
                      high = "#005B8E",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.002", size = 3) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd)) 
RR <- ggplot(new_tempR_R, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dR, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0,1,0,0,"cm")) +
  labs(y = "", x = "", title = "Rock & Boulder") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Rugosity")), 
                      high = "#006A4D",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = TeX("$p = 9.06 \\times 10^{-9}"), size = 3) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd)) 
AD <- ggplot(new_tempA_D, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dA, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        legend.position = "none",
        plot.margin = margin(0,0,0,1,"cm")) +
  labs(y = "Exposure, kW/m (log10)", x = "") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Fractal \n dimension")), 
                      high = "#B55000",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.288", size = 3) +
  annotate(geom="text", x = -4.2, y = -.5, label = "Fractal dimension", angle = 90, family = "Times-Roman", fontface = "bold", size = 5) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd)) 
PD <- ggplot(new_tempP_D, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dP, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        legend.position = "none",
        plot.margin = margin(0,0.5,0,0.5,"cm")) +
  labs(y = "", x = "") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Fractal \n dimension")), 
                      high = "#005B8E",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.044", size = 3) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd)) 
RD <- ggplot(new_tempR_D, aes(depth_s, wave_max_log10_s)) + 
  geom_raster(aes(fill=as.numeric(pred))) +
  geom_point(data = dR, aes(x = .data[[v1]], y = .data[[v2]]),
             shape = 21, size = 1, stroke = 0.5, color = "black") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white', colour = "black"),
        text = element_text(family = "Times-Roman", size = 10),
        legend.position = "none",
        plot.margin = margin(0,1,0,0,"cm")) +
  labs(y = "", x = "") +
  scale_y_continuous(breaks = c((0-wave_mn)/wave_sd, 
                                (.5-wave_mn)/wave_sd, 
                                (1-wave_mn)/wave_sd, 
                                (1.5-wave_mn)/wave_sd, 
                                (2-wave_mn)/wave_sd, 
                                (2.5-wave_mn)/wave_sd), 
                     labels = c(0,.5,1,1.5,2,2.5)) +
  scale_x_continuous(breaks = c((0-dep_mn)/dep_sd, 
                                (5-dep_mn)/dep_sd, 
                                (10-dep_mn)/dep_sd, 
                                (15-dep_mn)/dep_sd, 
                                (20-dep_mn)/dep_sd,
                                (25-dep_mn)/dep_sd,
                                (30-dep_mn)/dep_sd), 
                     labels = c(0,5,10,15,20,25,30)) +
  scale_fill_gradient(name = expression(paste("Fractal \n dimension")), 
                      high = "#006A4D",
                      low = "linen",
                      na.value = "white") +
  annotate(geom = "text", x = 1.1, y = -3.2, family = "Times-Roman", label = "p = 0.008", size = 3) +
  coord_cartesian(clip = "off", xlim = c((0-dep_mn)/dep_sd,(30-dep_mn)/dep_sd)) 

int.figure <- grid.arrange(AR + labs(subtitle = "(a)"), PR + labs(subtitle = "(b)"), RR + labs(subtitle = "(c)"),
                           AD + labs(subtitle = "(d)"), PD + labs(subtitle = "(e)"), RD + labs(subtitle = "(f)"),
                           AH + labs(subtitle = "(g)"), PH + labs(subtitle = "(h)"), RH + labs(subtitle = "(i)"),
                           nrow = 3, ncol = 3)
interact.fig <- grid.arrange(int.figure, mylegend, heights = c(10,1))

#ggsave("figs/figureS3.png", plot = interact.fig, width = 7, height = 6, units = c("in"))  

#### FIGURE S4 ####
# models run in FIGURE S3

# Height range
temp.h <- as.data.frame(summary(mod.sw.H)$coefficients) 
temp.h <- tibble::rownames_to_column(temp.h, "VALUE")
h.conf <- as.data.frame(confint(mod.sw.H, level = 0.95))
h.conf <- tibble::rownames_to_column(h.conf, "VALUE")
temp.h <- merge(temp.h, h.conf)
temp.h$variable_plot <- c("Geological Reef Age","Coral Cover","Depth","Depth x Exposure", 
                          "Habitat","Habitat", "Geological Reef Age", "Coral Cover","Depth", 
                          "Depth x Exposure","Exposure","Habitat","Geological Reef Age", "Coral Cover",
                          "Depth", "Depth x Exposure","Exposure", "Exposure")
temp.h$habitat <- c("AGR","AGR","AGR","AGR","AGR",
                    "PAV","PAV","PAV","PAV",
                    "PAV","PAV","ROB","ROB",
                    "ROB","ROB","ROB","ROB","AGR")
for (o in 1:nrow(temp.h)) {
  if(sign(temp.h[o,6]) == sign(temp.h[o,7])) {
    temp.h$shape[o] <- 19}
  else {
    temp.h$shape[o] <- 1}}

H.estimate <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = temp.h$variable_plot,
                    ymin = temp.h$`2.5 %`,
                    ymax = temp.h$`97.5 %`,
                    group = temp.h$habitat),
                position = position_dodge(width = 1), width = .5) +
  geom_point(aes(x = temp.h$variable_plot,
                 y = temp.h$Estimate,
                 color = temp.h$habitat),
             shape = temp.h$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), color = "lightgrey", linetype = "dashed") +
  ggtitle("(c) Height range") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10),
        axis.text.y=element_blank())

temp.h2 <- temp.h[temp.h$variable_plot != "Habitat",]
H.estimate2 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = temp.h2$variable_plot,
                    ymin = temp.h2$`2.5 %`,
                    ymax = temp.h2$`97.5 %`,
                    group = temp.h2$habitat),
                position = position_dodge(width = 1), width = .5) +
  geom_point(aes(x = temp.h2$variable_plot,
                 y = temp.h2$Estimate,
                 color = temp.h2$habitat),
             shape = temp.h2$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), color = "lightgrey", linetype = "dashed") +
  ggtitle("(c) Height range") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10),
        axis.text.y=element_blank())

# Rugosity 
temp.r <- as.data.frame(summary(mod.sw.R)$coefficients)
temp.r <- tibble::rownames_to_column(temp.r, "VALUE")
r.conf <- as.data.frame(confint(mod.sw.R, level = 0.95))
r.conf <- tibble::rownames_to_column(r.conf, "VALUE")
temp.r <- merge(temp.r, r.conf)
temp.r$variable_plot <- c("Geological Reef Age","Coral Cover","Depth","Depth x Exposure", 
                          "Habitat","Habitat", "Geological Reef Age", "Coral Cover","Depth", 
                          "Depth x Exposure","Exposure","Habitat","Geological Reef Age", "Coral Cover",
                          "Depth","Depth x Exposure","Exposure", "Exposure")
temp.r$habitat <- c("AGR","AGR","AGR","AGR","AGR",
                    "PAV","PAV","PAV","PAV",
                    "PAV","PAV","ROB","ROB",
                    "ROB","ROB","ROB","ROB","AGR")
for (o in 1:nrow(temp.r)) {
  if(sign(temp.r[o,6]) == sign(temp.r[o,7])) {
    temp.r$shape[o] <- 19}
  else {
    temp.r$shape[o] <- 1}}

R.estimate <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = temp.r$variable_plot,
                    ymin = temp.r$`2.5 %`,
                    ymax = temp.r$`97.5 %`,
                    group = temp.r$habitat),
                position = position_dodge(width = 1), width = .5) +
  geom_point(aes(x = temp.r$variable_plot,
                 y = temp.r$Estimate,
                 color = temp.r$habitat),
             shape = temp.r$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), color = "lightgrey", linetype = "dashed") +
  ggtitle("(a) Rugosity") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10))

temp.r2 <- temp.r[temp.r$variable_plot != "Habitat",]
R.estimate2 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = temp.r2$variable_plot,
                    ymin = temp.r2$`2.5 %`,
                    ymax = temp.r2$`97.5 %`,
                    group = temp.r2$habitat),
                position = position_dodge(width = 1), width = .5) +
  geom_point(aes(x = temp.r2$variable_plot,
                 y = temp.r2$Estimate,
                 color = temp.r2$habitat),
             shape = temp.r2$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), color = "lightgrey", linetype = "dashed") +
  ggtitle("(a) Rugosity") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10))

# Fractal dimension 
temp.d <- as.data.frame(summary(mod.sw.D)$coefficients)
temp.d <- tibble::rownames_to_column(temp.d, "VALUE")
d.conf <- as.data.frame(confint(mod.sw.D, level = 0.95))
d.conf <- tibble::rownames_to_column(d.conf, "VALUE")
temp.d <- merge(temp.d, d.conf)
temp.d$variable_plot <- c("Geological Reef Age","Coral Cover","Depth","Depth x Exposure", 
                          "Habitat","Habitat", "Geological Reef Age", "Coral Cover","Depth", 
                          "Depth x Exposure","Exposure","Habitat","Geological Reef Age", "Coral Cover",
                          "Depth","Depth x Exposure","Exposure", "Exposure")
temp.d$habitat <- c("AGR","AGR","AGR","AGR","AGR",
                    "PAV","PAV","PAV","PAV",
                    "PAV","PAV","ROB","ROB",
                    "ROB","ROB","ROB","ROB","AGR")
for (o in 1:nrow(temp.d)) {
  if(sign(temp.d[o,6]) == sign(temp.d[o,7])) {
    temp.d$shape[o] <- 19}
  else {
    temp.d$shape[o] <- 1}}

D.estimate <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = temp.d$variable_plot,
                    ymin = temp.d$`2.5 %`,
                    ymax = temp.d$`97.5 %`,
                    group = temp.d$habitat),
                position = position_dodge(width = 1), width = .5) +
  geom_point(aes(x = temp.d$variable_plot,
                 y = temp.d$Estimate,
                 color = temp.d$habitat),
             shape = temp.d$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), color = "lightgrey", linetype = "dashed") +
  ggtitle("(b) Fractal dimension") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10),
        axis.text.y=element_blank())

temp.d2 <- temp.d[temp.d$variable_plot != "Habitat",]
D.estimate2 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = temp.d2$variable_plot,
                    ymin = temp.d2$`2.5 %`,
                    ymax = temp.d2$`97.5 %`,
                    group = temp.d2$habitat),
                position = position_dodge(width = 1), width = .5) +
  geom_point(aes(x = temp.d2$variable_plot,
                 y = temp.d2$Estimate,
                 color = temp.d2$habitat),
             shape = temp.d2$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), color = "lightgrey", linetype = "dashed") +
  ggtitle("(b) Fractal dimension") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10),
        axis.text.y=element_blank())

sup4.1 <- grid.arrange(R.estimate, D.estimate, H.estimate,
                       nrow = 1, ncol = 3, widths = c(1.5,1,1))
sup4 <- grid.arrange(sup4.1, mylegend, heights = c(10,1))
#ggsave("figs/figureS4_1.png", plot = sup4, width = 7, height = 3, units = c("in"))

sup4.1.2 <- grid.arrange(R.estimate2, D.estimate2, H.estimate2,
                         nrow = 1, ncol = 3, widths = c(1.5,1,1))
sup4.2 <- grid.arrange(sup4.1.2, mylegend, heights = c(10,1))
#ggsave("figs/figureS4_2.png", plot = sup4.2, width = 7, height = 3, units = c("in"))

#### FIGURE S5 ####
cor <- svyglm(coral_s ~ 0 + habitat + depth_s:habitat + age_s:habitat + wave_max_log10_s:habitat + depth_s:wave_max_log10_s:habitat, design = des)

summary(cor)
summ(cor)

temp.c <- as.data.frame(summary(cor)$coefficients)
temp.c <- tibble::rownames_to_column(temp.c, "VALUE")
c.conf <- as.data.frame(confint(cor, level = 0.95))
c.conf <- tibble::rownames_to_column(c.conf, "VALUE")
temp.c <- merge(temp.c, c.conf)
temp.c$variable_plot <- c("Habitat","Geological Reef Age","Depth","Depth x Exposure","Exposure",
                          "Habitat","Geological Reef Age","Depth","Depth x Exposure","Exposure",
                          "Habitat","Geological Reef Age","Depth","Depth x Exposure","Exposure")
temp.c$habitat <- c("AGR","AGR","AGR","AGR","AGR",
                    "PAV","PAV","PAV","PAV","PAV",
                    "ROB","ROB","ROB","ROB","ROB")
for (o in 1:nrow(temp.c)) {
  if(sign(temp.c[o,6]) == sign(temp.c[o,7])) {
    temp.c$shape[o] <- 19}
  else {
    temp.c$shape[o] <- 1}}

coral.relat <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_point(aes(x = temp.c$variable_plot,
                 y = temp.c$Estimate,
                 color = temp.c$habitat),
             shape = temp.c$shape,
             position = position_dodge(width = 1), size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("AGR" = "#D55E00",
                                "PAV" = "#0072B2",
                                "ROB" = "#009E73")) +
  labs(x = "", y = "\nParameter Estimate", color = "Habitat type") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), color = "lightgrey", linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10))

sup5 <- grid.arrange(coral.relat, mylegend, heights = c(10,1))
#ggsave("figs/figureS5.png", plot = sup5, width = 4.5, height = 4, units = c("in"))

#### FIGURE S6 ####
hab_depth <- ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = depth), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Depth, m", x = "") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))
hab_age <- ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = age), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Geological reef age, mya", x = "") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))
hab_wave <- ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = log10(wave_max)), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Exposure, kW/m (log10)", x = "Habitat type") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))
hab_coral <- ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = CORAL), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Coral Cover, %", x = "Habitat type") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))


sup6 <-  grid.arrange(hab_depth + ggtitle("(a)"), hab_age + ggtitle("(b)"), 
                      hab_wave + ggtitle("(c)"), hab_coral + ggtitle("(d)"), nrow = 2, ncol = 2)
#ggsave("figs/figureS6.png", plot = sup6, width = 7, height = 5, units = c("in"))

