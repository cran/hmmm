### R code from vignette source 'TUTORIAL_hmmm.rnw'

###################################################
### code chunk number 1: TUTORIAL_hmmm.rnw:179-180
###################################################
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: TUTORIAL_hmmm.rnw:182-185
###################################################
library("hmmm")
data("accident", package = "hmmm")
accident[1:20,]


###################################################
### code chunk number 3: TUTORIAL_hmmm.rnw:196-197
###################################################
y <- getnames(accident, st = 9)


###################################################
### code chunk number 4: TUTORIAL_hmmm.rnw:199-202
###################################################
count<-cbind(row.names(y)[1:3],y[1:3])
colnames(count)<-c("cell names", "counts")
print(count,quote=F)


###################################################
### code chunk number 5: TUTORIAL_hmmm.rnw:217-219
###################################################
margin <- marg.list(c("marg-marg-b-b", "b-marg-b-b", 
"marg-b-b-b", "b-b-b-b"))


###################################################
### code chunk number 6: TUTORIAL_hmmm.rnw:224-227
###################################################
model <- hmmm.model(marg = margin, lev = c(3, 4, 3, 2),
names = c("Type", "Time", "Age", "Hour"))
model


###################################################
### code chunk number 7: TUTORIAL_hmmm.rnw:241-244
###################################################
modelB <- hmmm.model(marg = margin, lev = c(3, 4, 3, 2),
names = c("Type", "Time", "Age", "Hour"),
sel = c(12:13, 14:17))


###################################################
### code chunk number 8: TUTORIAL_hmmm.rnw:248-250
###################################################
modB <- hmmm.mlfit(y, modelB)
modB


###################################################
### code chunk number 9: TUTORIAL_hmmm.rnw:254-255 (eval = FALSE)
###################################################
## print(modB, aname = "model B", printflag = TRUE)


###################################################
### code chunk number 10: TUTORIAL_hmmm.rnw:259-260 (eval = FALSE)
###################################################
## summary(modB)


###################################################
### code chunk number 11: TUTORIAL_hmmm.rnw:278-282
###################################################
modelA <- hmmm.model(marg = margin, lev = c(3, 4, 3, 2),
names = c("Type", "Time", "Age", "Hour"), sel = c(12:13, 14:17),
formula = ~ Type * Age * Hour + Time * Age * Hour + Type : Time)
modA <- hmmm.mlfit(y, modelA)


###################################################
### code chunk number 12: TUTORIAL_hmmm.rnw:286-287
###################################################
anova(modA, modB)


###################################################
### code chunk number 13: TUTORIAL_hmmm.rnw:294-298 (eval = FALSE)
###################################################
## modellog <- loglin.model(lev = c(3, 4, 3, 2),
## formula = ~ Type * Age * Hour + Time * Age * Hour + Type : Time,
## names = c("Type", "Time", "Age", "Hour"))
## modlog <- hmmm.mlfit(y, modellog)


###################################################
### code chunk number 14: TUTORIAL_hmmm.rnw:333-335
###################################################
margin <- marg.list(c("marg-marg-l-l", "g-marg-l-l", 
"marg-g-l-l", "g-g-l-l"))


###################################################
### code chunk number 15: TUTORIAL_hmmm.rnw:346-349
###################################################
model <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4), 
names = c("In", "Sa", "Co", "Ho"))
model


###################################################
### code chunk number 16: TUTORIAL_hmmm.rnw:356-362
###################################################
model1 <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = c("In", "Sa", "Co", "Ho"), sel = c(18:23, 26:27, 34:39))
data("madsen", package = "hmmm")
y <- getnames(madsen, st = 6)
mod1 <- hmmm.mlfit(y, model1)
mod1


###################################################
### code chunk number 17: TUTORIAL_hmmm.rnw:366-370
###################################################
model2 <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = c("In", "Sa", "Co", "Ho"), sel = c(18:23, 34:39))
mod2 <- hmmm.mlfit(y, model2)
mod2


###################################################
### code chunk number 18: TUTORIAL_hmmm.rnw:379-383
###################################################
model3 <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = c("In", "Sa", "Co", "Ho"), sel = c(18:23, 34:39, 44:71))
mod3 <- hmmm.mlfit(y, model3)
mod3


###################################################
### code chunk number 19: TUTORIAL_hmmm.rnw:427-428
###################################################
marginals <- marg.list(c("r-marg", "marg-r", "r-r"))


###################################################
### code chunk number 20: TUTORIAL_hmmm.rnw:441-443
###################################################
rec1 <- matrix(c(-1, -1,  1,
               -1,  1,  0), 2, 3, byrow = TRUE)


###################################################
### code chunk number 21: TUTORIAL_hmmm.rnw:452-458
###################################################
rec2<-matrix(c(-1, -1, -1,  0,  1,  1,  1,
             -1, -1, -1,  1,  0,  0,  0,
              1, -1,  0,  0,  0,  0,  0,
              0, -1,  1,  0,  0,  0,  0,
              0,  0,  0,  0,  1, -1,  0,
              0,  0,  0,  0,  0, -1,  1), 6, 7, byrow = TRUE)


###################################################
### code chunk number 22: TUTORIAL_hmmm.rnw:466-467
###################################################
rec <- recursive(rec1, rec2)


###################################################
### code chunk number 23: TUTORIAL_hmmm.rnw:475-478
###################################################
model <- hmmm.model(marg = marginals, lev = c(3, 7), 
names = c("Rel", "Pol"), cocacontr = rec)
model


###################################################
### code chunk number 24: TUTORIAL_hmmm.rnw:511-515
###################################################
Emat <- cbind(matrix(0, 2, 4), matrix(c(1, 0, 0, 1, 0, -1, -1, 0), 2, 4), 
matrix(0, 2, 12))
modelE <- hmmm.model(marg = marginals, lev = c(3, 7), 
names = c("Rel", "Pol"), cocacontr = rec, E = Emat)


###################################################
### code chunk number 25: TUTORIAL_hmmm.rnw:520-524
###################################################
data("relpol", package = "hmmm")
y <- getnames(relpol, st = 4)
modE <- hmmm.mlfit(y, modelE)
print(modE)


###################################################
### code chunk number 26: TUTORIAL_hmmm.rnw:540-542
###################################################
data("depression", package="hmmm")
y <- getnames(depression, st = 9)


###################################################
### code chunk number 27: TUTORIAL_hmmm.rnw:566-571
###################################################
margin <- marg.list(c("marg-marg-marg-b-b","b-marg-marg-b-b",
"marg-b-marg-b-b", "marg-marg-b-b-b","b-b-b-b-b"))
name <- c("R3","R2","R1","T","D")
modelsat<-hmmm.model(marg = margin, lev = c(2,2,2,2,2), names = name)
modelsat


###################################################
### code chunk number 28: TUTORIAL_hmmm.rnw:589-601
###################################################
A1<-matrix(c(
   0,0,0,1,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,1,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,1,
   0,1,0,0,0,0,0,0,0,-1,0,0,
   0,0,0,0,0,1,0,0,0,-1,0,0,
   0,0,1,0,0,0,0,0,0,0,-1,0,
   0,0,0,0,0,0,1,0,0,0,-1,0,
   1,0,0,0,-2,0,0,0,1,0,0,0
   ),8,12,byrow=TRUE)

E1<-cbind(matrix(0,8,3), A1, matrix(0,8,16))


###################################################
### code chunk number 29: TUTORIAL_hmmm.rnw:606-610
###################################################
model1<-hmmm.model(marg = margin, lev =c(2,2,2,2,2), names = name, E = E1)

fitmod1 <- hmmm.mlfit(y, model1)
fitmod1


###################################################
### code chunk number 30: TUTORIAL_hmmm.rnw:615-626
###################################################
A2<-matrix(c(
    0,0,0,1,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,1,
    0,1,0,0,0,-2,0,0,0,1,0,0,
    0,0,1,0,0,0,0,0,0,0,-1,0,
    0,0,0,0,0,0,1,0,0,0,-1,0,
    1,0,0,0,-2,0,0,0,1,0,0,0),
    7,12,byrow=TRUE)

E2<-cbind(matrix(0,7,3), A2, matrix(0,7,16))


###################################################
### code chunk number 31: TUTORIAL_hmmm.rnw:631-635
###################################################
model2<-hmmm.model(marg = margin, lev = c(2,2,2,2,2), names = name, E = E2)

fitmod2 = hmmm.mlfit(y, model2)
fitmod2


###################################################
### code chunk number 32: TUTORIAL_hmmm.rnw:645-650
###################################################
model3<-hmmm.model(marg = margin, lev = c(2,2,2,2,2), names = name, E = E2,
                   formula=~R1*R2*T*D+R3*R2*T*D )

fitmod3 = hmmm.mlfit(y, model3)
fitmod3


###################################################
### code chunk number 33: TUTORIAL_hmmm.rnw:676-677
###################################################
marginals <- marg.list(c("b-marg", "marg-g", "b-g"))


###################################################
### code chunk number 34: TUTORIAL_hmmm.rnw:691-696
###################################################
al <- list(
Type = ~ Type * (Age + Hour),
Time = ~ Time * (Age + Hour),
Type.Time = ~ Type.Time * (Age + Hour)
)


###################################################
### code chunk number 35: TUTORIAL_hmmm.rnw:711-714
###################################################
model <- hmmm.model.X(marg = marginals, lev = c(3, 4), 
names = c("Type", "Time"), Formula = al, strata = c(3, 2), 
fnames = c("Age", "Hour"))


###################################################
### code chunk number 36: TUTORIAL_hmmm.rnw:720-724
###################################################
data("accident", package = "hmmm")
y <- getnames(accident, st = 9)
mod1 <- hmmm.mlfit(y, model)
mod1


###################################################
### code chunk number 37: TUTORIAL_hmmm.rnw:730-731 (eval = FALSE)
###################################################
## summary(mod1)


###################################################
### code chunk number 38: TUTORIAL_hmmm.rnw:743-748
###################################################
alind <- list(
Type = ~ Type * Age + Type * Hour,
Time = ~ Time * Age + Time * Hour,
Type.Time = "zero"
)


###################################################
### code chunk number 39: TUTORIAL_hmmm.rnw:760-765
###################################################
alpar <- list(
Type = ~ Type + Age + Hour,
Time = ~ Time + Age + Hour,
Type.Time = ~ Type.Time + Age + Hour
)


###################################################
### code chunk number 40: TUTORIAL_hmmm.rnw:798-802
###################################################
data("polbirth", package = "hmmm")
y <- getnames(polbirth)
marginals <- marg.list(c("g-marg", "marg-l", "g-l"))
names <- c("Politics", "Birth")


###################################################
### code chunk number 41: TUTORIAL_hmmm.rnw:808-809
###################################################
ineq <- list(marg = c(1, 2), int = list(c(1, 2)), types = c("g", "l"))


###################################################
### code chunk number 42: TUTORIAL_hmmm.rnw:815-817
###################################################
model <- hmmm.model(marg = marginals, dismarg = ineq, lev = c(7, 4), 
names = names)


###################################################
### code chunk number 43: TUTORIAL_hmmm.rnw:827-828
###################################################
mlr <- hmmm.mlfit(y, model, noineq = FALSE)


###################################################
### code chunk number 44: TUTORIAL_hmmm.rnw:834-835
###################################################
msat <- hmmm.mlfit(y, model)


###################################################
### code chunk number 45: TUTORIAL_hmmm.rnw:841-844
###################################################
model0 <- hmmm.model(marg = marginals, lev = c(7, 4), sel = c(10:27), 
names = names)
mnull <- hmmm.mlfit(y, model0)


###################################################
### code chunk number 46: TUTORIAL_hmmm.rnw:854-855
###################################################
test <- hmmm.chibar(nullfit = mnull, disfit = mlr, satfit = msat)


###################################################
### code chunk number 47: TUTORIAL_hmmm.rnw:874-875
###################################################
test


###################################################
### code chunk number 48: TUTORIAL_hmmm.rnw:896-897
###################################################
y <- matrix(c(104, 24, 65, 76, 146, 30, 50, 9, 166), 9, 1)


###################################################
### code chunk number 49: TUTORIAL_hmmm.rnw:908-910
###################################################
Zmat <- kronecker(diag(3), matrix(1, 3, 1))
ZFmat <- kronecker(diag(3), matrix(1, 3, 1))[,3]


###################################################
### code chunk number 50: TUTORIAL_hmmm.rnw:922-929
###################################################
Gini <- function(m) {
A<-matrix(m,3,3,byrow=TRUE)
 GNum<-rowSums(A^2)
 GDen<-rowSums(A)^2
 G<-GNum/GDen
 c(G[1], G[3]) - c(G[2], G[1])
 }


###################################################
### code chunk number 51: TUTORIAL_hmmm.rnw:939-940
###################################################
mod_eq <- mphineq.fit(y, Z = Zmat, ZF = ZFmat, h.fct = Gini)


###################################################
### code chunk number 52: TUTORIAL_hmmm.rnw:948-949
###################################################
mod_ineq <- mphineq.fit(y, Z = Zmat, ZF = ZFmat, d.fct = Gini)


###################################################
### code chunk number 53: TUTORIAL_hmmm.rnw:955-957
###################################################
mod_sat <- mphineq.fit(y, Z = Zmat, ZF = ZFmat)
hmmm.chibar(nullfit = mod_eq, disfit = mod_ineq, satfit = mod_sat)


