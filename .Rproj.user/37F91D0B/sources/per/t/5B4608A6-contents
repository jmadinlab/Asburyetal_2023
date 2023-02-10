library(plyr)
library(survey)
library(jtools)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(latex2exp)
#library(RCPA3)

data <- read.csv("Complexity_sub.csv")

L <- 2 #2x2m box
scl <- L / c(1, 2, 3.125, 6.25, 12.5, 25, 50, 100, 200)
L0 <- min(scl) 

# load data
data2 <- data %>% mutate(R = 0.5*log10(R_theory_mean^2 - 1),
                         D = log10(L/L0)*(D_mean-3),
                         H = log10(H_mean/(sqrt(2)*L0))) 

# survey weights
sectors <- read.csv("Sectors-Strata-Areas.csv")
wtemp <- left_join(data2, sectors)

weight <- plyr::ddply(wtemp, .(SEC_NAME,DEPTH_BIN,NH), summarize, n = length(unique(site)))
weight$survey_weight <- weight$NH / weight$n

data3 <- left_join(data2, weight)

data3$conc <- c(paste0(data3$island, "_", data3$SEC_NAME, "_", data3$DEPTH_BIN, "_", data3$site))
des <- svydesign(id= ~1, strata = ~ conc, weights = ~survey_weight, data = data3) 

wave_mn <- mean(data2$wave_max_log10)
wave_sd <- sd(data2$wave_max_log10)
dep_mn  <- mean(data2$depth)
dep_sd  <- sd(data2$depth)
age_mn  <- mean(data2$age)
age_sd  <- sd(data2$age)
coral_mn <- mean(data2$CORAL)
coral_sd <- sd(data2$CORAL)
dxw_mean <- mean((data2$depth * data2$wave_max_log10))
dxw_sd <- sd((data2$depth * data2$wave_max_log10))

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

# Height range #################################################################################################################

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

# H effects plot
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

# Rugosity #####################################################################################################################
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

# R effects plot
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

# Fractal Dimension ############################################################################################################
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

# D effects plot 
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
ggsave("FINALDRAFT/Figure3_surveyweight_updated.png", plot = fig3, width = 7, height = 5, units = c("in"), dpi = 600)

sup5.1 <- grid.arrange(R.estimate, D.estimate, H.estimate,
                     nrow = 1, ncol = 3, widths = c(1.5,1,1))
sup5 <- grid.arrange(sup5.1, mylegend, heights = c(10,1))
ggsave("SupFigure5_surveyweight1.png", plot = sup5, width = 7, height = 3, units = c("in"))

sup5.1.2 <- grid.arrange(R.estimate2, D.estimate2, H.estimate2,
                       nrow = 1, ncol = 3, widths = c(1.5,1,1))
sup5.2 <- grid.arrange(sup5.1.2, mylegend, heights = c(10,1))
ggsave("SupFigure5_surveyweight2.png", plot = sup5.2, width = 7, height = 3, units = c("in"))


# No effect of habitat #########################################################################################################

mod.sw.H2 <-svyglm(H ~ depth_s + age_s + wave_max_log10_s + coral_s + depth_s*wave_max_log10_s, design=des)
mod.sw.R2 <-svyglm(R ~ depth_s + age_s + wave_max_log10_s + coral_s + depth_s*wave_max_log10_s, design=des)
mod.sw.D2 <-svyglm(D ~ depth_s + age_s + wave_max_log10_s + coral_s + depth_s*wave_max_log10_s, design=des)

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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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
             color = "coral4", 
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


sup4 <- grid.arrange(arrangeGrob(R.DEP2 + ggtitle("(a)"), R.AGE2 + ggtitle("(b)"), R.WAVE2 + ggtitle("(c)"), R.COR2 + ggtitle("(d)"),
                                 D.DEP2 + ggtitle("(e)"), D.AGE2 + ggtitle("(f)"), D.WAVE2 + ggtitle("(g)"), D.COR2 + ggtitle("(h)"),
                                 H.DEP2 + ggtitle("(i)"), H.AGE2 + ggtitle("(j)"), H.WAVE2 + ggtitle("(k)"), H.COR2 + ggtitle("(l)"),
                                 nrow = 3, ncol = 4))

ggsave("SupFigure4_surveyweight_updated.png", plot = sup4, width = 7, height = 5, units = c("in"))

# Coral vs Predictors ##########################################################################################################

cor <- svyglm(coral_s ~ 0 + habitat + depth_s:habitat + age_s:habitat + wave_max_log10_s:habitat + depth_s:wave_max_log10_s:habitat, design = des)

summary(cor)
plot_summs(cor)
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
  # geom_errorbar(aes(x = temp.c$variable_plot,
  #                   ymin = temp.c$`2.5 %`,
  #                   ymax = temp.c$`97.5 %`,
  #                   group = temp.c$habitat),
  #               position = position_dodge(width = 1), width = .5) +
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
  #ggtitle("Coral Cover") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        text = element_text(family = "Times-Roman", size = 10))

sup6 <- grid.arrange(coral.relat, mylegend, heights = c(10,1))
ggsave("SupFigure6_surveyweight.png", plot = sup6, width = 4.5, height = 4, units = c("in"))



# coral & depth
CvsD <- svyglm(CORAL ~ depth:habitat, design = des)
summary(CvsD)
CorD.A.nd <- data3[data3$habitat == "AGR",][c("depth","habitat")]
CorD.A.nd$depth <- seq(min(data3$depth[data3$habitat == "AGR"]), 
                       max(data3$depth[data3$habitat == "AGR"]),
                       length = nrow(CorD.A.nd))
CorD.P.nd <- data3[data3$habitat == "PAV",][c("depth","habitat")]
CorD.P.nd$depth <- seq(min(data3$depth[data3$habitat == "PAV"]), 
                       max(data3$depth[data3$habitat == "PAV"]),
                       length = nrow(CorD.P.nd))
CorD.R.nd <- data3[data3$habitat == "ROB",][c("depth","habitat")]
CorD.R.nd$depth <- seq(min(data3$depth[data3$habitat == "ROB"]), 
                        max(data3$depth[data3$habitat == "ROB"]),
                       length = nrow(CorD.R.nd))
CorD.nd <- rbind(CorD.A.nd, CorD.P.nd, CorD.R.nd)
pred_CorD <- predict(CvsD, newdata = CorD.nd, type = "response", se.fit = TRUE)
pred_CorD <- as.data.frame(pred_CorD)
colnames(pred_CorD) <- c("Predicted","SE")
CorD.nd <- cbind(CorD.nd, pred_CorD)

# coral & age
CvsA <- svyglm(CORAL ~ age:habitat, design = des)
summary(CvsA)
CorA.A.nd <- data3[data3$habitat == "AGR",][c("age","habitat")]
CorA.A.nd$age <- seq(min(data3$age[data3$habitat == "AGR"]), 
                       max(data3$age[data3$habitat == "AGR"]),
                       length = nrow(CorA.A.nd))
CorA.P.nd <- data3[data3$habitat == "PAV",][c("age","habitat")]
CorA.P.nd$age <- seq(min(data3$age[data3$habitat == "PAV"]), 
                       max(data3$age[data3$habitat == "PAV"]),
                       length = nrow(CorA.P.nd))
CorA.R.nd <- data3[data3$habitat == "ROB",][c("age","habitat")]
CorA.R.nd$age <- seq(min(data3$age[data3$habitat == "ROB"]), 
                       max(data3$age[data3$habitat == "ROB"]),
                       length = nrow(CorA.R.nd))
CorA.nd <- rbind(CorA.A.nd, CorA.P.nd, CorA.R.nd)
pred_CorA <- predict(CvsA, newdata = CorA.nd, type = "response", se.fit = TRUE)
pred_CorA <- as.data.frame(pred_CorA)
colnames(pred_CorA) <- c("Predicted","SE")
CorA.nd <- cbind(CorA.nd, pred_CorA)

# coral & wave
CvsW <- svyglm(CORAL ~ wave_max_log10:habitat, design = des)
summary(CvsW)
CorW.A.nd <- data3[data3$habitat == "AGR",][c("wave_max_log10","habitat")]
CorW.A.nd$wave_max_log10 <- seq(min(data3$wave_max_log10[data3$habitat == "AGR"]), 
                     max(data3$wave_max_log10[data3$habitat == "AGR"]),
                     length = nrow(CorW.A.nd))
CorW.P.nd <- data3[data3$habitat == "PAV",][c("wave_max_log10","habitat")]
CorW.P.nd$wave_max_log10 <- seq(min(data3$wave_max_log10[data3$habitat == "PAV"]), 
                     max(data3$wave_max_log10[data3$habitat == "PAV"]),
                     length = nrow(CorW.P.nd))
CorW.R.nd <- data3[data3$habitat == "ROB",][c("wave_max_log10","habitat")]
CorW.R.nd$wave_max_log10 <- seq(min(data3$wave_max_log10[data3$habitat == "ROB"]), 
                     max(data3$wave_max_log10[data3$habitat == "ROB"]),
                     length = nrow(CorW.R.nd))
CorW.nd <- rbind(CorW.A.nd, CorW.P.nd, CorW.R.nd)
pred_CorW <- predict(CvsW, newdata = CorW.nd, type = "response", se.fit = TRUE)
pred_CorW <- as.data.frame(pred_CorW)
colnames(pred_CorW) <- c("Predicted","SE")
CorW.nd <- cbind(CorW.nd, pred_CorW)

# plot
group.colors <- c("AGR" = "#D55E00", "PAV" = "#0072B2", "ROB" = "#009E73")
coral.dep <- ggplot() + 
  geom_point(aes(x = data3$depth, y = data3$CORAL), color = "lightgrey", alpha = 0.2) +
  geom_line(aes(x = CorD.nd$depth, 
                y = CorD.nd$Predicted,
                color = CorD.nd$habitat), size = 1) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Depth, m", y = "Coral cover, %") 
coral.age <- ggplot() + 
  geom_point(aes(x = data3$age, y = data3$CORAL), color = "lightgrey", alpha = 0.2) +
  geom_line(aes(x = CorA.nd$age, 
                y = CorA.nd$Predicted,
                color = CorA.nd$habitat), size = 1) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Geological reef age, mya", y = "") 
coral.wave <- ggplot() + 
  geom_point(aes(x = data3$wave_max_log10, y = data3$CORAL), color = "lightgrey", alpha = 0.2) +
  geom_line(aes(x = CorW.nd$wave_max_log10, 
                y = CorW.nd$Predicted,
                color = CorW.nd$habitat), size = 1) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none") +
  labs(x = "Exposure, kW/m (log10)", y = "") 

sup3 <- grid.arrange(arrangeGrob(coral.dep + ggtitle("(a)"), coral.age + ggtitle("(b)"), coral.wave + ggtitle("(c)"),
                      nrow = 1, ncol = 3), mylegend, heights = c(10,1))

ggsave("SupFigure3_surveyweight.png", plot = sup3, width = 7, height = 3, units = c("in"))

# Does it respect the geometric constraint? ####################################################################################

R_1 <- sqrt(10^(as.vector(predict(mod.sw.R))/.5) + 1)
R_temp <- as.vector(predict(mod.sw.H))+as.vector(predict(mod.sw.D))
R_2 <- sqrt(10^(R_temp/.5) + 1)

comp1 <- ggplot() +
  geom_point(aes(x = R_1, y = R_2), color = "black", fill = "grey", shape = 21) +
  geom_abline(intercept = 0) +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10),
        plot.margin = unit(c(2, 0, 2, 2), "pt")) +
  labs(x = expression('r'[1]),
       y = expression('r'[2])) 

comp2 <- ggplot() +
  geom_histogram(aes(R_1/R_2), color = "black", fill = "grey") +
  geom_vline(aes(xintercept = 1), linetype = "dashed", color = "darkred") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10),
        plot.margin = unit(c(2, 2, 2, 0), "pt")) +
  labs(x = expression(frac('r'[1], 'r'[2])), y = "Density") 

fig4 <- grid.arrange(comp1 + ggtitle("(a)"), comp2 + ggtitle("(b)"), 
                     nrow = 1, ncol = 2)
ggsave("FINALDRAFT/figure4.png", plot = fig4, width = 6, height = 3, units = c("in"), dpi = 600)


# Data table ###################################################################################################################
library(data.table)
tab.isl <- data.frame(Island = c("Hawaii","Maui","Lanai","Molokai","Oahu","Kauai","Niihau"),
                  n = c(65,33,32,41,32,18,13),
                  Habitat = c("Aggregate Reef","Pavement","Rock and Boulder",NA,NA,NA,NA),
                  n = c(94,62,78,NA,NA,NA,NA))

tab <- rbind(tab.isl,tab.hab)
list(tab.isl, tab.hab)%>% 
  kbl() %>%
  kable_classic(html_font = "Times New Roman") 
tab.isl %>% 
  kbl(escape = F) %>%
  kable_classic(html_font = "Times New Roman") 

# Interaction plots ############################################################################################################
library(sp)

  v1 = "depth_s" # set interaction variable 1
  v2 = "wave_max_log10_s" # set interaction variable 2
  d = data3 # set database of observations
  
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

ggsave("Interaction_figure.png", plot = int.figure, width = 7, height = 6, units = c("in"))  

# visualize differences in habitat #############################################################################################

ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = age), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Geological reef age, mya", x = "Habitat type") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))

ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = CORAL), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Coral Cover, %", x = "Habitat type") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))

ggplot(data3) + 
  geom_boxplot(aes(x = habitat, y = H), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(y = "Height range", x = "Habitat type") +
  theme_classic() + 
  theme(text = element_text(family = "Times-Roman", size = 12)) +
  scale_x_discrete(limits = c("AGR","PAV","ROB"),
                   labels = c("Aggregate reef",
                              "Pavement",
                              "Rock & Boulder"))

hab_R <- ggplot(data3) +
  geom_boxplot(aes(x = habitat, y = R), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(x = "", y = "Rugosity") +
  theme_classic()

hab_D <- ggplot(data3) +
  geom_boxplot(aes(x = habitat, y = D), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(x = "", y = "Fractal dimension") +
  theme_classic() 

hab_H <- ggplot(data3) +
  geom_boxplot(aes(x = habitat, y = H), fill = c("#D55E00", "#0072B2",  "#009E73")) +
  labs(x = "", y = "Height range") +
  theme_classic() 

grid.arrange(hab_R, hab_D, hab_H, nrow = 1, ncol = 3)

mod_R <- aov(R ~ habitat, data.bayes)
TukeyHSD(mod_R, conf.level = .95)

mod_D <- aov(D ~ habitat, data.bayes)
TukeyHSD(mod_D, conf.level = .95)

mod_H <- aov(H ~ habitat, data.bayes)
TukeyHSD(mod_H, conf.level = .95)


ggplot() +
  geom_ribbon(aes(x=newdata.P$coral_s,
                  ymin=newdata.P$Predict.lwr,
                  ymax=newdata.P$Predict.upr), fill = "#0072B2", alpha = 0.2) +
  geom_line(aes(x=newdata.P$coral_s,
                y=newdata.P$Predicted), size = .8, color = "#0072B2") +
  geom_ribbon(aes(x=newdata.A$coral_s, 
                  ymin=newdata.A$Predict.lwr, 
                  ymax=newdata.A$Predict.upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(aes(x=newdata.A$coral_s, 
                y=newdata.A$Predicted), size = .8, color = "#D55E00") +
  geom_ribbon(aes(x=newdata.R$coral_s,
                  ymin=newdata.R$Predict.lwr,
                  ymax=newdata.R$Predict.upr), fill = "#009E73", alpha = 0.2) +
  geom_line(aes(x=newdata.R$coral_s,
                y=newdata.R$Predicted), size = .8, color = "#009E73") +
  theme_classic() +
  theme(text = element_text(family = "Times-Roman", size = 10), legend.position = "none",
        plot.margin = unit(c(2, 2, 0, 0), "pt")) +
  labs(x = "", y = "") +
  ylim(-0.8,1.2) +
  scale_x_continuous(breaks = c((0-coral_mn)/coral_sd,
                                (20-coral_mn)/coral_sd,
                                (40-coral_mn)/coral_sd,
                                (60-coral_mn)/coral_sd,
                                (80-coral_mn)/coral_sd),
                     labels = c(0,20,40,60,80))


# plane of best fit ############################################################################################################

png("Figure1_orthogonal.png", width = 5, height = 5, units = "in", res = 600)
par(mar=c(1.2, 2, 1, 1) +0.1, mfrow = c(2,2), family = "Times", ps = 10)

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
mtext("(a)", side = 3, at = -.5)
text3D(-0.6, 2.9, -1.2, labels=expression(italic(r)^2 == 0.974), surf = NULL, add = TRUE, family = "Times")

#####
dataAGR <- data2[data2$habitat == "AGR" ,]

mdat <- dataAGR[c("R","H")]
mdat$m <- 1

# best fit orthogonal plane 
A <- data.matrix(mdat)
A <- matrix(c (sum(mdat$R * mdat$R) , sum(mdat$R * mdat$H) , sum(mdat$R) , 
               sum(mdat$R * mdat$H) , sum(mdat$H * mdat$H) , sum(mdat$H) ,
               sum(mdat$R) , sum(mdat$H) , nrow(mdat)), nrow = 3, ncol = 3)
A.i <- inv(A)
B <- matrix(c (sum(mdat$R * dataAGR["D"]) , 
               sum(mdat$H * dataAGR["D"]) , 
               sum(dataAGR["D"])), nrow = 3, ncol = 1)
fin <- A.i %*% B

f <- function(x, y) { z <- x * fin[1,1] + y * fin[2,1] + fin[3,1] }

# plane for plotting 
grid.lines = 350
R.pred <- seq(min(data2$R), max(data2$R), length=grid.lines)
H.pred <- seq(min(data2$H), max(data2$H), length=grid.lines)
mat <- expand.grid(R=seq(min(data2$R), max(data2$R), length=grid.lines), 
                   H=H.pred)

D.mat <- matrix(f(mat$R, mat$H), nrow = grid.lines, ncol = grid.lines)

R2 <- dataAGR$R
H2 <- dataAGR$H
D2 <- dataAGR$D

AGRcollist <- c("#D55E00","#FFFFFF")
AGRcolors <- colorRampPalette(AGRcollist)(100)

scatter3D(R2, H2, D2, pch = 20, cex = 1, 
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
title("Aggregate reef", family = "Times")
mtext("(b)", side = 3, at = -.5)

#####
dataPAV <- data2[data2$habitat == "PAV" ,]

mdat <- dataPAV[c("R","H")]
mdat$m <- 1

# best fit orthogonal plane 
A <- data.matrix(mdat)
A <- matrix(c (sum(mdat$R * mdat$R) , sum(mdat$R * mdat$H) , sum(mdat$R) , 
               sum(mdat$R * mdat$H) , sum(mdat$H * mdat$H) , sum(mdat$H) ,
               sum(mdat$R) , sum(mdat$H) , nrow(mdat)), nrow = 3, ncol = 3)
A.i <- inv(A)
B <- matrix(c (sum(mdat$R * dataPAV["D"]) , 
               sum(mdat$H * dataPAV["D"]) , 
               sum(dataPAV["D"])), nrow = 3, ncol = 1)
fin <- A.i %*% B

f <- function(x, y) { z <- x * fin[1,1] + y * fin[2,1] + fin[3,1] }

# plane for plotting 
grid.lines = 350
R.pred <- seq(min(data2$R), max(data2$R), length=grid.lines)
H.pred <- seq(min(data2$H), max(data2$H), length=grid.lines)
mat <- expand.grid(R=seq(min(data2$R), max(data2$R), length=grid.lines), 
                   H=H.pred)

D.mat <- matrix(f(mat$R, mat$H), nrow = grid.lines, ncol = grid.lines)

R2 <- dataPAV$R
H2 <- dataPAV$H
D2 <- dataPAV$D

PAVcollist <- c("#0072B2","#FFFFFF")
PAVcolors <- colorRampPalette(PAVcollist)(100)

scatter3D(R2, H2, D2, pch = 20, cex = 1, 
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
title("Pavement", family = "Times")
mtext("(c)", side = 3, at = -.5)

#####
dataROB <- data2[data2$habitat == "ROB" ,]

mdat <- dataROB[c("R","H")]
mdat$m <- 1

# best fit orthogonal plane 
A <- data.matrix(mdat)
A <- matrix(c (sum(mdat$R * mdat$R) , sum(mdat$R * mdat$H) , sum(mdat$R) , 
               sum(mdat$R * mdat$H) , sum(mdat$H * mdat$H) , sum(mdat$H) ,
               sum(mdat$R) , sum(mdat$H) , nrow(mdat)), nrow = 3, ncol = 3)
A.i <- inv(A)
B <- matrix(c (sum(mdat$R * dataROB["D"]) , 
               sum(mdat$H * dataROB["D"]) , 
               sum(dataROB["D"])), nrow = 3, ncol = 1)
fin <- A.i %*% B

f <- function(x, y) { z <- x * fin[1,1] + y * fin[2,1] + fin[3,1] }

# plane for plotting 
grid.lines = 350
R.pred <- seq(min(data2$R), max(data2$R), length=grid.lines)
H.pred <- seq(min(data2$H), max(data2$H), length=grid.lines)
mat <- expand.grid(R=seq(min(data2$R), max(data2$R), length=grid.lines), 
                   H=H.pred)

D.mat <- matrix(f(mat$R, mat$H), nrow = grid.lines, ncol = grid.lines)

R2 <- dataROB$R
H2 <- dataROB$H
D2 <- dataROB$D

ROBcollist <- c("#009E73","#FFFFFF")
ROBcolors <- colorRampPalette(ROBcollist)(100)

scatter3D(R2, H2, D2, pch = 20, cex = 1, 
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
title("Rock & Boulder", family = "Times")
mtext("(d)", side = 3, at = -.5)
#text3D(1.8, 1.2, 2.7, labels="Rock & boulder", add = TRUE, family = "Times")

dev.off()


################################################################################################################################

dat <- read.csv("Complexity_all.csv")
dat <- dat[c("site","latitude","longitude","depth","habitat_original","R_mean","H_mean","D_mean")]

colnames(dat) <- c("SITE","LATITUDE","LONGITUDE","DEPTH","HABITAT_CODE","RUGOSITY","HEIGHT_RANGE","FRACTAL_DIMENSION")

meta <- read.csv("BenthicCover_2010-2020_Tier1_SITE.csv")
meta <- meta[c("OBS_YEAR","MISSIONID","REGION","ISLAND",
               "SITE","SITEVISITID","REEF_ZONE","DEPTH_BIN","DATE_","LATITUDE","LONGITUDE")]
meta$SITE2 <- meta$SITE
meta$SITE <- gsub("-0", "-", meta$SITE)


tot <- left_join(dat, meta, by = "SITE")

final <- tot[, c(9,10,11,12,19,13,14,15,16,5,2,3,4,6,7,8)]
colnames(final)[2] <- "MISSION_ID"
colnames(final)[5] <- "SITE"
colnames(final)[11] <- "LATITUDE"
colnames(final)[12] <- "LONGITUDE"
colnames(final)[13] <- "MAX_DEPTH_M"

colnames(final)

write.csv(final, "INPORT_STRUCTURAL_COMPLEXITY.csv")
