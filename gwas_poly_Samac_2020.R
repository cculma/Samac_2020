rm(list=ls()) # remove all files in global environment

# BiocManager::install("qvalue")
# install.packages("devtools")
# devtools::install_github("jendelman/GWASpoly", build_vignettes=FALSE)
# lines 3 and 5 are only required once, then you can comment with: #
# install.packages("devtools")
# devtools::install_github("jendelman/GWASpoly", build_vignettes=FALSE)
# install.packages(c("Cairo", "hrbrthemes", "ggpubr", "sommer", "doParallel", "vcfR"))

library(vcfR)
library(parallel)
library(doParallel)
library(iterators)
library(foreach)
library(devtools)
library(sommer)
library(ggpubr)
library(ggthemes)
library(hrbrthemes)
library(Cairo)
library(ggplot2)
library(GWASpoly)
library(qvalue)
library(tidyverse)
library(updog)
library(ggpubr)
future::availableCores()
setwd("~/OneDrive - Washington State University (email.wsu.edu)/Samac_2020/")

# setwd("~/Documents/Cesar/GWASpoly/Samac_2020") # you must to change the directory where the files will be:
# setwd("~/Documents") # if your files are in Documents

pheno <- read.csv("~/Documents/Cesar/GWASpoly/Samac_2020/pheno4.csv", row.names = 1)
head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-5)]
trait1
# models_1 <- c("general", "additive")
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")
# models_1 <- c("general", "additive", "1-dom", "diplo-additive", "diplo-general")
params <- set.params(fixed=c("PC1","PC2","PC3", "PC4", "PC5"),
                     fixed.type=rep("numeric",5), n.PC = 5)

data_a <- read.GWASpoly(ploidy=4, pheno.file="pheno4.csv", 
                        geno.file="SMC_GWASpoly_28346.txt", format="numeric", n.traits=length(trait1), delim=",")
data_b <- read.GWASpoly(ploidy=4, pheno.file="pheno4.csv", 
                        geno.file="SMC_UPDOG_28732.txt", format="numeric", n.traits=length(trait1), delim=",")
data_c <- read.GWASpoly(ploidy=4, pheno.file="pheno4.csv", 
                        geno.file="samac_4.txt", format="ACGT", n.traits=length(trait1), delim=",")

data_a1 <- set.K(data_a, LOCO = F)
data_a2 <- set.K(data_a, LOCO = T)
data_b1 <- set.K(data_b, LOCO = F)
data_b2 <- set.K(data_b, LOCO = T)
data_c1 <- set.K(data_c, LOCO = F)
data_c2 <- set.K(data_c, LOCO = T)

data_a1_1 <- GWASpoly(data_a1, models = models_1, traits = trait1, params = params, n.core = 32)
data_a2_1 <- GWASpoly(data_a2, models = models_1, traits = trait1, params = params, n.core = 32)
data_b1_1 <- GWASpoly(data_b1, models = models_1, traits = trait1, params = params, n.core = 32)
data_b2_1 <- GWASpoly(data_b2, models = models_1, traits = trait1, params = params, n.core = 32)
data_c1_1 <- GWASpoly(data_c1, models = models_1, traits = trait1, params = params, n.core = 32)
data_c2_1 <- GWASpoly(data_c2, models = models_1, traits = trait1, params = params, n.core = 32)

save.image("~/Documents/Cesar/GWASpoly/Samac_2020/gwas_poly1.RData")
# SMC_GWASpoly_28346.txt

# save.image("~/Documents/Cesar/GWASpoly/Samac_2020/gwas_poly1.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/Samac_2020/gwas_poly1.RData")
######

data_a1_1a <- set.threshold(data_a1_1, method= "Bonferroni", level=0.4) 
data_a1_1b <- set.threshold(data_a1_1, method= "FDR", level=0.05)
data_a1_2a <-  get.QTL(data_a1_1a, traits=c ("avg_score", "avg_result"))
data_a1_2b <-  get.QTL(data_a1_1b, bp.window=10e6)

data_c1_1a <- set.threshold(data_c1_1, method= "Bonferroni", level=0.4) 

manhattan.plot(data = data_c1_1a, traits=c ("avg_score", "avg_result"))
manhattan.plot(data = data_c1_1a, traits=c ("avg_score", "avg_result"))
manhattan.plot(data = data_c1_1a, traits = "avg_score", chrom = "Chr7")



##########

plot1 <- manhattan.plot(data = data_c1_1a, traits=c ("avg_score", "avg_result")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("darkblue","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

plot2 <- manhattan.plot(data = data_c1_1a, traits = "avg_result", chrom = "Chr4") + theme_classic(base_family = "Arial", base_size = 12) + geom_point(color='darkblue') + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

plot3 <- manhattan.plot(data = data_c1_1a, traits = "avg_score", chrom = "Chr7") + theme_classic(base_family = "Arial", base_size = 12) + geom_point(color='azure4') + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
rm(myplot1)

myplot2 <- ggarrange(plot1, ggarrange(plot2, plot3, ncol = 2, labels = c("b", "c")), nrow = 2,   labels = "a")
ggsave(filename = "myplot2.jpg", plot = myplot2, width = 6, height = 6)

myplot1 <- LD.plot(data_c1_1a) + theme_classic(base_family = "Arial", base_size = 12) + ggtitle("LD plot")
ggsave(filename = "myplot1.jpg", plot = myplot1, width = 6, height = 6)

##########

p <- LD.plot(data_4)
p + xlim(0,30) 

manhattan.plot(data = data_a1_1a)

plot1 <- manhattan.plot(data = data_a1_1a, traits=c ("avg_score", "avg_result"))
plot2 <- manhattan.plot(data = data_a1_1b, traits=c ("avg_score"), chrom="Chr7")
?manhattan.plot
plot1.1 <- plot1  + theme_ipsum(base_family = "Arial", base_size = 12) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12)) + scale_y_continuous(breaks=seq(0, 16, 2)) + ggtitle("SMC_GWASpoly_28346.txt")

plot1  + theme_ipsum(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("darkblue","azure4")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12)) + scale_y_continuous(breaks=seq(0, 8, 2))

plot2.1 <- plot2  + theme_ipsum(base_family = "Times", base_size = 12) + ggtitle("Chr_7") 
class(plot2.1) 

plot2  + theme_ipsum(base_family = "Arial", base_size = 12) + ggtitle("Chr_7") + geom_point(color='darkblue')


# plot3.1 <- plot3  + theme_ipsum(base_family = "Times", base_size = 12) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12)) + scale_y_continuous(breaks=seq(0, 16, 2)) + ggtitle("GWASpoly_samac_4.txt") 
figure <- ggarrange(plot2.1, plot1.1,
                    labels = c("a", "b"),
                    ncol = 1, nrow = 2)


########

# model 

P2 <- read.csv("~/Documents/Cesar/GWASpoly/Samac_2020/pheno3.csv")
P3 <- read.csv("~/Documents/Cesar/GWASpoly/Samac_2020/pheno_samac.csv")
P4 <- read.csv("~/Documents/Cesar/GWASpoly/Samac_2020/pheno_samac1.csv")
PCA <- read.csv("~/Documents/Cesar/GWASpoly/Samac_2020/pca_5.csv")

head(P2)
head(P3)
head(P4)
head(PCA)
colnames(P3)[1] <- "Line"
colnames(P4)[1] <- "Line"
colnames(PCA)[1] <- "Line"

P3 <- P3[,c(1,2)]
P4 <- P4[,c(1,3)]

P6 <- inner_join(P2, P3, by = "Line") %>% inner_join(., P4, by = "Line") %>% inner_join(., PCA, by = "Line")
write.csv(P6, "~/Documents/Cesar/GWASpoly/Samac_2020/pheno4.csv", row.names = F, quote = F)


########


data_1 <- data_a1_2a
head(data_1)
lab1 <- c('Trait', 'Model', 'Marker', 'Chrom', 'Ref', 'Alt')
data_1[,lab1] <- lapply(data_1[,lab1], factor)
data_1$Score <- as.numeric(data_1$Score)
str(data_1)

summary(data_1$Model)
hist(data_1$Score)
data_2 <- subset(data_1, Score > 4)
data_4 <-data_2 %>% distinct(Marker, .keep_all = TRUE)

#########
# FDR

k <- as.data.frame(data_c1_1a@scores)
dim(k)

rownames(k) <- paste(data_c1_1a@map$Chrom, data_c1_1a@map$Position, sep = "_")
logp <- 10^-k
dim(logp)
logq <- na.omit(logp$avg_score.diplo.additive)
length(logq)
hist(logq)
qobj <- qvalue(p = logp$avg_score.additive, pfdr = T, lfdr.out = T)
summary(qobj)


##################
# Annotation

##################
col_headings <- c('chrom',	'chromStart',	'chromEnd',	'name',	'score',	'strand',	'thickStart',	'thickEnd',	'itemRgb',	'blockCount',	'blockSizes',	'blockStarts')
col_headings_1 <- c('gene_id',	'isoform', 'gene_name',	'trans_length_flag',	'blastp_match_flag',	'nmd_flag',	'frame')
col_headings_2 <- c('gene_id',	'isoform1')
##################

data_7 <-  data_a1_2a
######
# data table mrbean
head(data_7)
trait2 <- c("Trait", "Model", "Chrom")
data_7[,trait2] <- lapply(data_7[,trait2], factor)
str(data_7)
S1 <- data_7[,c(4,1,2,9)]
head(S1)
S3 <- as.data.frame(summary(S1$Trait))
S4 <- as.data.frame(summary(S1$Model))
colnames(S3) <- "summary"
colnames(S4) <- "summary"
S5 <- rbind(S3, S4)

S2 <- spread(data = S1, key = Model, value = Score)
S2 <- S1 %>% distinct(Marker, .keep_all = TRUE)
S2 <- S2[,-2]

S3 <- spread(data = S1, key = Trait, value = Score)
S3 <- S3 %>% distinct(Marker, .keep_all = TRUE)
S3 <- S3[,-2]

summary(data_7$Trait)
data_7.2 <- subset(data_7,Trait == "avg_result" | Trait == "avg_score")

data_8 <- spread(data = data_7.2, key = Trait, value = Score)
data_8 <-data_8 %>% distinct(Marker, .keep_all = TRUE)
# colSums(!is.na(data_8[,c(9:14)]))
# colnames(data_8)

data_4 <-data_a1_2a %>% distinct(Marker, .keep_all = TRUE)
data_4 <- unite(data = data_4, col = "Marker1", 5:6, sep = "_", remove = T)

S3 <- as.data.frame(summary(data_7.2$Trait))
S4 <- as.data.frame(summary(data_4$Model))
S5 <- as.data.frame(summary(data_4$Chrom))
colnames(S3) <- "summary"
colnames(S4) <- "summary"
colnames(S5) <- "summary"
S6 <- rbind(S3, S4, S5)
# write.csv(S6, "summary_2stage_all.csv", row.names = T, quote = FALSE, na = "")

colnames(data_4)
data4.1 <- unite(data = data_4, col = "SNP", 8:7, sep = "/", remove = T)
data4.2 <- unite(data = data4.1, col = "Marker1", 5:6, sep = "_", remove = T)
colnames(data4.2)
str(data4.2)
data4.2 <- data4.2[,-c(1,2,3,8)]

S4 <- inner_join(data4.2, S2, by = "Marker")
S5 <- inner_join(S4, S3, by = "Marker")
write.csv(S5, "markers_mrbean.csv", row.names = F, quote = FALSE, na = "")

##################
# GRanges
library(GenomicRanges)
library(GenomicFeatures)
library(genomation)
library(Repitools)
library(plyranges)



# file <- ("~/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.bed")
file <- ("~/OneDrive - Washington State University (email.wsu.edu)/Samac_2020/blast_corrected_shen.bed")

txdb <- readBed(file, track.line = FALSE, remove.unusual = FALSE,
                zero.based = TRUE)

gr5 <- GRanges(seqnames = data_4$Chrom,
               ranges = IRanges(data_4$Position, width = 1))

overlaps <- join_overlap_left(gr5, txdb)

df2 <- annoGR2DF(overlaps)

df2 <- unite(data = df2, col = "Marker1", 1:2, sep = "_", remove = F)
df2 <-df2 %>% distinct(Marker1, .keep_all = TRUE)

colnames(df2)
df2 <- df2[,c(1:4,7)]
df3 <- separate(df2, 5, col_headings_1, sep = ";", remove = TRUE, convert = FALSE, extra = "warn")
df4 <- separate(df3, 5, col_headings_2, sep = "\\.", remove = T, convert = FALSE, extra = "warn")
df4 <- df4[1:(length(df4)-5)]
df4 <- inner_join(data_4, df4, by = "Marker1")

write.csv(df4, "~/OneDrive - Washington State University (email.wsu.edu)/Samac_2020/anno1.csv", row.names = F, quote = F) # search in uniprot



df5 <- read.table("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/uniprot1.csv", sep = "\t", header = T)
df5 <- left_join(df4, df5, by = "uniprot")
df5 <- df5[,-c(2:4)]
write.table(df5, "~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/anno2.csv", row.names = F, quote = F, sep = "\t")
######################

