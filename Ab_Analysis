setwd("H:/plink-1.07-dos/")
library(psych)
library(GPArotation)

#Merged 3 populations: 35 overlapped abs, 13 abs with >200 cases(11+1+1)
eira49 <- read.table("EIRA_pheno1-49_newID_plink.txt")
fam <- read.table("SE-E.fam")
case <- subset(fam, fam[,6]==2)
case_eira49 <- subset(eira49, eira49[,2] %in% case[,2]) #2459
wtccc47 <- read.table("WTCCC_pheno47.txt")
fam <- read.table("UK.fam")
case <- subset(fam, fam[,6]==2)
case_wtccc47 <- subset(wtccc47, wtccc47[,2] %in% case[,2]) #2040
narac67 <- read.table("narac_pheno67_case-ctrl_plink.txt")
fam <- read.table("US.fam")
case <- subset(fam, fam[,6]==2)
case_narac67 <- subset(narac67, narac67[,2] %in% case[,2]) #1908

eira35 <- cbind(case_eira49[,1:2], case_eira49[,5:10], case_eira49[,17], case_eira49[,22:24], case_eira49[,28:30], case_eira49[,14:16], case_eira49[,20], case_eira49[,19], case_eira49[,21], case_eira49[,33:38], case_eira49[,40:44], case_eira49[,46:50])
wtccc35 <- cbind(case_wtccc47[,1:2], case_wtccc47[,3], case_wtccc47[,5], case_wtccc47[,7], case_wtccc47[,9], case_wtccc47[,11], case_wtccc47[,13], case_wtccc47[,15], case_wtccc47[,18:23], case_wtccc47[,25], case_wtccc47[,27], case_wtccc47[,29], case_wtccc47[,31:49])
narac35 <- cbind(case_narac67[,1:2], case_narac67[,3], case_narac67[,5], case_narac67[,7], case_narac67[,9], case_narac67[,11], case_narac67[,13], case_narac67[,24], case_narac67[,31:33], case_narac67[,39:41], case_narac67[,18], case_narac67[,20], case_narac67[,22], case_narac67[,28], case_narac67[,27], case_narac67[,29], case_narac67[,48:53], case_narac67[,56:60], case_narac67[,63:67])
names(eira35) <- names(case_eira49)[1:37]
names(wtccc35) <- names(case_eira49)[1:37]
names(narac35) <- names(case_eira49)[1:37]
allcase35 <- rbind(eira35, wtccc35, narac35)

mergedcase <- as.matrix(allcase35[,3:37])
fa.parallel(mergedcase, fa="PC", show.legend=F, main="Scree plot in merged 35 phenotypes") 
#Parallel analysis suggests that the number of factors =  11  and the number of components =  7 
merged_pc <- principal(mergedcase, nfactors=7, rotate="none")
merged_pc
merged_rc <- principal(mergedcase, nfactors=7, rotate="varimax")
merged_rc

# Principal Components Analysis
# Call: principal(r = mergedcase, nfactors = 7, rotate = "varimax")
# Standardized loadings (pattern matrix) based upon correlation matrix
# RC1   RC2   RC3   RC4   RC5   RC6   RC7    h2   u2
# V3   0.78  0.14 -0.01 -0.02 -0.01 -0.04  0.01 0.633 0.37
# V4   0.71 -0.03 -0.03  0.00  0.04  0.35  0.11 0.646 0.35
# V5   0.60  0.11  0.00  0.04 -0.04  0.36 -0.01 0.502 0.50
# V6   0.79  0.13  0.00 -0.01  0.04  0.02  0.06 0.652 0.35
# V7   0.81  0.14 -0.02 -0.01  0.01 -0.04  0.03 0.679 0.32
# V8   0.67  0.15 -0.01  0.01 -0.02 -0.18  0.03 0.512 0.49
# V9   0.73  0.06  0.00 -0.01  0.05  0.19  0.02 0.575 0.42
# V10  0.33  0.79  0.02  0.00  0.00 -0.02 -0.02 0.730 0.27
# V11  0.21  0.82  0.01  0.02 -0.03  0.07 -0.01 0.728 0.27
# V12 -0.04  0.69 -0.03  0.03  0.03  0.33  0.02 0.589 0.41
# V13  0.56  0.53  0.01  0.00  0.07 -0.28  0.05 0.688 0.31
# V14  0.60  0.52  0.01 -0.01  0.05 -0.30  0.03 0.721 0.28
# V15  0.29  0.53  0.04 -0.02  0.07 -0.21  0.21 0.458 0.54
# V16  0.62  0.33  0.01  0.00  0.02 -0.23  0.01 0.546 0.45
# V17  0.12  0.04  0.02 -0.02  0.03  0.64  0.11 0.436 0.56
# V18  0.63 -0.16  0.00  0.02 -0.06  0.29 -0.17 0.546 0.45
# V19 -0.01  0.05  0.19 -0.07  0.27  0.14 -0.09 0.146 0.85
# V20 -0.03  0.02  0.73 -0.02  0.05  0.01  0.04 0.532 0.47
# V21  0.06  0.01 -0.07 -0.08  0.13 -0.03  0.41 0.200 0.80
# V22  0.02  0.02  0.45  0.01 -0.11 -0.06  0.25 0.285 0.71
# V23 -0.01 -0.02  0.27 -0.02  0.47 -0.07 -0.08 0.308 0.69
# V24  0.00  0.14  0.09  0.01  0.20 -0.01  0.26 0.138 0.86
# V25  0.04 -0.04 -0.03  0.05  0.32  0.09  0.11 0.129 0.87
# V26  0.01  0.00 -0.07  0.05  0.61 -0.07  0.09 0.391 0.61
# V27  0.01  0.09 -0.02 -0.03  0.38  0.05  0.34 0.274 0.73
# V28  0.02  0.03  0.09  0.00  0.14  0.09  0.32 0.143 0.86
# V29  0.01  0.02  0.02  0.84  0.04  0.00  0.06 0.707 0.29
# V30 -0.02  0.00  0.56  0.06  0.09  0.10  0.02 0.333 0.67
# V31  0.00  0.00  0.07  0.04  0.08  0.01  0.25 0.076 0.92
# V32 -0.02  0.03 -0.02 -0.01 -0.04  0.14  0.35 0.143 0.86
# V33 -0.02  0.02  0.03  0.07 -0.18 -0.06  0.49 0.287 0.71
# V34  0.04 -0.06  0.03  0.02 -0.11 -0.08  0.43 0.209 0.79
# V35  0.01 -0.01  0.69  0.02  0.07 -0.06  0.03 0.484 0.52
# V36 -0.01  0.00  0.04  0.84  0.02 -0.03 -0.01 0.709 0.29
# V37  0.00  0.03  0.02  0.01  0.46 -0.04 -0.05 0.221 0.78
# 
# RC1  RC2  RC3  RC4  RC5  RC6  RC7
# SS loadings           5.48 2.88 1.67 1.44 1.32 1.29 1.28
# Proportion Var        0.16 0.08 0.05 0.04 0.04 0.04 0.04
# Cumulative Var        0.16 0.24 0.29 0.33 0.37 0.40 0.44
# Proportion Explained  0.36 0.19 0.11 0.09 0.09 0.08 0.08
# Cumulative Proportion 0.36 0.54 0.65 0.75 0.83 0.92 1.00
# 
# Test of the hypothesis that 7 components are sufficient.
# 
# The degrees of freedom for the null model are  595  and the objective function was  8.51
# The degrees of freedom for the model are 371  and the objective function was  2.03 
# The total number of observations was  6407  with MLE Chi Square =  12956.4  with prob <  0 
# 
# Fit based upon off diagonal values = 0.91

fam <- read.table("SE-E.fam")
ctrl <- subset(fam, fam[,6]==1)
ctrl_eira49 <- subset(eira49, eira49[,2] %in% ctrl[,2]) #1885
fam <- read.table("UK.fam")
ctrl <- subset(fam, fam[,6]==1)
ctrl_wtccc47 <- subset(wtccc47, wtccc47[,2] %in% ctrl[,2]) #8430
fam <- read.table("US.fam")
ctrl <- subset(fam, fam[,6]==1)
ctrl_narac67 <- subset(narac67, narac67[,2] %in% ctrl[,2]) #2134
allctrlid <- rbind(ctrl_eira49[,1:2], ctrl_wtccc47[,1:2], ctrl_narac67[,1:2]) #12449
allctrl <- matrix(rep(1,12449*35), nrow=12449, ncol=35)
allctrl <- cbind(allctrlid, allctrl)
names(allctrl) <- names(allcase35)
all35 <- rbind(allcase35, allctrl) #6407+12449=18856
write.table(all35, "Merged_pheno1-35_case-ctrl.txt", quote=F, row.names=F, col.names=F)
