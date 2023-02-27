

# load data 


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

