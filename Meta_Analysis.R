####Load packages ####
library(metafor)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(finalfit)
library(writexl)
library(wesanderson)
#Analyses done in metafor but dmetar used for some diagnostics
library(devtools)
#devtools::install("~/Documents/Sydney/Research/Dissertation/Experiments/Meta-analysis/R Codes/dmetar-master")
library(dmetar)

####Import and inspect data####
setwd("~/Documents/Sydney/Research/Dissertation/Experiments/Meta-analysis/Data")
meta_all <- read.csv(file.choose())
str(meta_all)
#Note: NAs in lnRR and lnRR_variance removed in excel previously (full data file saved as "Meta_Final_NAs.csv")

####Structure and Inspect Data ####
#Remove cell culture studies
meta<-subset(meta_all,level=="Organism")
str(meta)
#Remove Clairardin et al 2013 due no appropriate control group (only turtle study)
meta<- subset(meta, first_author != "Clairardin")
#Modify data for later analyses (from Nicholas Wu)
meta<- subset(meta, n_control != 1)  # remove any studies with sample size of 1
meta$lnRR_variance2 <- abs(meta$lnRR_variance) # change negative values to absolute values (for variance)
meta$lnRR_variance2[meta$lnRR_variance2 == 0] <- meta$lnRR_variance2 + 0.0001 #add 0.0001 to remove 0 values
str(meta)

###Descriptive summary of full data
nrow(meta) # number of effect sizes = 778
length(unique(meta$study_ID)) # number of studies = 52
length(unique(meta$scientific_name_OTL))# number of species = 12
length(unique(meta$molecule)) # number of molecules = 65

#Types of molecules measured
summary2<-table(meta$molecule_type)
summary2<-as.data.frame(summary2)
summary2
#Most entirely hormones (659/778 = ~84.7% comprised of hormones)

####Preliminary Investigation of Hormone Data ####
#Subset data to only look at hormones
horm <- subset(meta , molecule_type %in% "Hormone")
str(horm)
#Remove vitellogenin (VTG)
horm <-subset(horm, molecule != "Vitellogenin")
###Descriptive summary
nrow(horm) # number of effect sizes = 624
length(unique(horm$study_ID)) # number of studies = 46
length(unique(horm$scientific_name_OTL))# number of species = 11
length(unique(horm$molecule)) # number of molecules = 24

# Number of effects sizes per class/species
class_summary<-table(horm$class)
class_summary<-as.data.frame(class_summary)
ggplot(class_summary,aes(Var1, Freq, fill=Var1)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label=Freq), position = position_stack(vjust = 0.5), size=5) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) + theme_classic() +
  coord_flip()
#Individual species
species_summary<-table(horm$scientific_name_OTL, horm$study_ID)
species_summary<-as.data.frame(species_summary)
species_graph<-subset(species_summary, Freq != 0)
ggplot(species_graph,aes(Var1, Freq, fill=Var2)) + 
  geom_bar(position="stack",stat="identity",colour="black") +
  geom_text(aes(label=Freq), position = position_stack(vjust = 0.5), color="black", size=5) +
  coord_flip() + theme_classic() + theme(legend.position = "none")
#Clear bias of dataset towards mammals/rats (Rattus norvegicus) makes calculation of phylogenetic signal uninformative

#Look at distribution of hormones present across studies
molecule_summary<-table(horm$molecule, horm$study_ID)
molecule_summary<-as.data.frame(molecule_summary)
molecule_graph<-subset(molecule_summary, Freq != 0)
ggplot(molecule_graph,aes(Var1, Freq, fill=Var2)) + 
  geom_bar(position="stack",stat="identity",colour="black") +
  geom_text(aes(label=Freq), position = position_stack(vjust = 0.5), color="black", size=5) +
  coord_flip() + theme_classic() + theme(legend.position = "none")
# Will removed hormones where n studies is less than 2 and VTG

#Number and types of bisphenols measured
summary6<-table(horm$bisphenol)
summary6<-as.data.frame(summary6)
summary6
#Individual bisphenols measured present across studies
bisphenol_summary<-table(horm$bisphenol, horm$study_ID)
bisphenol_summary<-as.data.frame(bisphenol_summary)
bisphenol<-subset(bisphenol_summary, Freq != 0)
ggplot(bisphenol,aes(Var1, Freq, fill=Var2)) + 
  geom_bar(position="stack",stat="identity",colour="black") +
  geom_text(aes(label=Freq), position = position_stack(vjust = 0.5), color="black", size=5) +
  coord_flip() + theme_classic() + theme(legend.position = "none")
#remove BPB, TCBPA, and TBBPA

####Finalized Data  Summary ####
#Subset by hormones sampled in more than 2 studies
horm2<-subset(horm,abbreviation %in% c("T", "E2","P4", "FSH","LH", "T4", "T3","TSH", "CORT"))
#Subset by bisphenols sampled in more than 2 studies
level_key<- c(BPB="Derivative", TBBPA = "Derivative", TCBPA = "Derivative")
horm2$bisphenol<-dplyr::recode(horm2$bisphenol, !!!level_key)
horm2 <- subset(horm2,bisphenol != "Derivative")

#Investigation of data
nrow(horm2) # number of effect sizes = 514
length(unique(horm2$study_ID)) # number of studies = 44
length(unique(horm2$scientific_name_OTL))# number of species = 9
length(unique(horm2$molecule)) # number of molecules = 9

#Table of hormones
summary7<-table(horm2$molecule, horm2$type)
summary7<-as.data.frame(summary7)
summary7<-subset(summary7, Freq != 0)
summary7

#Table of species classifications (Note: Amphibians now removed from data)
class_summary2<-table(horm2$class)
class_summary2<-as.data.frame(class_summary2)
class_summary2

#Table of species
species_summary2<-table(horm2$scientific_name_OTL, horm2$study_ID)
species_summary2<-as.data.frame(species_summary2)
species_summary2

#Table of bisphenols
summary9<-table(horm2$bisphenol)
summary9<-as.data.frame(summary9)
summary9

###Final Analyses####
overall.horm <- rma.mv(yi = lnRR, V = lnRR_variance2,
                       random = list(~1 | scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID), 
                       method = "REML", data = horm2)
summary(overall.horm)
#create model without random effects
overall.horm2 <- rma.mv(yi = lnRR, V = lnRR_variance2,
                        method = "REML", data = horm2)
#Calculate heterogeneity (from Jackson 2012)
c(100 * (vcov(overall.horm)[1,1] - vcov(overall.horm2)[1,1]) / vcov(overall.horm)[1,1])
#heterogeneity = 99.86 (rounded from 99.86343)

####Publication Bias and Sensitivity (code from Wu and Seebacher 2020)####
#Time lag effect
ggplot(horm2, aes(x = year_published, y = lnRR, size = n_control)) + 
  geom_point(shape = 21, fill = "#4292c6", alpha = 0.5) + 
  labs(x = "Publication year", y = "Effect size", size = "N") +
  scale_size_area(max_size = 10) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  geom_hline(yintercept = 0, lty = 2) + 
  geom_smooth(method = "lm", size = 1, se = F, colour = "#084594") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)
pub.horm <- rma.mv(yi = lnRR, V = lnRR_variance2, mod = ~ year_published, 
                   random = list(~1 | scientific_name, ~1 | study_ID, ~1 | effect_size_ID), 
                   method = "REML", data = horm2)
summary(pub.horm)
# No time lag effect present

#Additional tests of publication bias
resid <- rstandard(overall.horm) # recover residuals, se, z, slab
resid.df <- do.call("rbind", resid)
resid.df2 <- t(resid.df)
resid.es <- rma(yi = resid, vi = se, data = resid.df2)
par(mfrow = c(2,1))
funnel(resid.es, yaxis = "seinv", legend = TRUE, main = "Inverse standard error")
taf <- trimfill(resid.es, estimator = "R0")
taf
funnel(taf, yaxis="seinv", legend=TRUE, main = "Trim-fill output")

#Egger's regression
regtest(resid.es, model="lm")

###Sensitivity & Influences Analyses; takes over 24 hours to run (output saved)
overall.cooks <- cooks.distance(overall.horm)
plot(overall.cooks, type = "o", pch = 19, xlab = "Observed Outcome", ylab = "Cook's Distance")
sample_size <- nrow(horm2)
abline(h = 4/sample_size, col="red")
text(x=1:length(overall.cooks)+1, y=overall.cooks, labels=ifelse(overall.cooks>4/sample_size, names(overall.cooks),""), col="red")
#all distances > 1 so should be fine
dfbetas(overall.horm)

#### Analysis of Different Subgroups and Moderators #####
####1. Class & species #####
#Large proportion of data is mammals/rats and so subset data to look specifically at those groups
mammal1<-subset(horm2,class=="Mammalia")
nrow(mammal1) # k =350
####1.1 Mammalian hormones ####
#Look at all mammals 
mammal.horm <- rma.mv(yi = lnRR, V = lnRR_variance2,
                   mods= ~molecule,
                   random = list( ~1 | study_ID, ~1 | effect_size_ID),
                   method = "REML", data = mammal1)
summary(mammal.horm)
#No effect of developmental stage, bisphenol, concentration, duration, lag, or sex (removed from models)
data.mammal <- data.frame(trait = substr(row.names(mammal.horm$b), 6, 50), estimate = mammal.horm$b, ci.lb = mammal.horm$ci.lb, ci.ub = mammal.horm$ci.ub, p=mammal.horm$pval)
all<-ggplot(data.mammal, aes(x = trait, y = estimate, colour = trait))+geom_hline(yintercept=0, linetype = "dashed") + geom_point(size = 3.5, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab(NULL) + ylim(-1,1) +
  scale_x_discrete(labels=c("pt"="CORT","uleEstradiol"="E2","uleFollicle Stimulating Hormone"="FSH",
                            "uleLuteinizing Hormone"="LH","uleProgesterone"="P4", "uleTestosterone"="T","uleThyroid Stimulating Hormone"= "TSH",
                            "uleThyroxine"="T4","uleTriiodothyronine"="T3"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18, face="bold",vjust=0.8, hjust = 0.5),
        axis.text.y = element_text(size=18)) + coord_fixed() +
  scale_color_manual(values = c("#74d5ce","#5dcec6","#46c7be", "#2fc0b6","#18baae","#15a79c","#13948b","#108279","#0e6f68","#0c5d57"))
#### 1.2 Rats ####
#Subset data to just look at rats, both overall and just hormone data
rat1<-subset(horm2,scientific_name_OTL=="Rattus norvegicus")
str(rat1)
#Inspect data
reference_table1<-aggregate(bisphenol_conc ~bisphenol+developmental_stage+exposure_method,rat1,mean)
reference_table1[order(reference_table1$bisphenol, reference_table1$developmental_stage, reference_table1$exposure_method), ]
#write_xlsx(reference_table1,'Reference Table 1.xlsx')
reference_table2<-rat1 %>% group_by(molecule, bisphenol, strain, sex) %>% summarise("# of Effect Sizes" = n())
as.data.frame(reference_table2)
#write_xlsx(reference_table2,'Reference Table 2.xlsx')

###Descriptive summary of data
nrow(rat1) # number of effect sizes = 312
length(unique(rat1$study_ID)) # number of studies = 28
length(unique(rat1$molecule)) # number of molecules = 9
#Look at table of when hormones are measured over life 
table(rat1$molecule,rat1$developmental_stage)
#Note molecule releveled to E2 because levels dropped if automatic reference (CORT) is used because of data distribution
rat1$molecule<-as.factor(rat1$molecule)
rat.horm <- rma.mv(yi = lnRR, V = lnRR_variance2,
                      mods= ~relevel(molecule,ref="Estradiol")*developmental_stage,
                      random = list( ~1 | study_ID, ~1 | effect_size_ID),
                      method = "REML", data = rat1)
summary(rat.horm)
#Graph 
data.rat <- data.frame(trait = substr(row.names(rat.horm$b),7,70), estimate = rat.horm$b, ci.lb = rat.horm$ci.lb, ci.ub = rat.horm$ci.ub, pval = rat.horm$pval)
#Export data into excel to more easily modify groupings, then re-imported for graphing (edited file is Rat_Data_Graph.csv)
##used "write_xlsx(data.rat,'Rat_Data.xlsx')" to create xlsx table
rat_graph<-read.csv(file.choose())
p<-ggplot(rat_graph, aes(x = Graph, y = estimate, colour = Graph))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(size = 3.5, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab("Effect Size (lnRR)") +
  ylim(-1.4,1.4) + facet_wrap(vars(Stage),nrow = 1,ncol = 5, scales = "free_x") + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size =1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 18),
        axis.line = element_line(colour = "black", size=.1),
        axis.text.x = element_text(size=18, face="bold", angle = 90, vjust=0.4, hjust = 0.95), 
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18, face="bold"),
        panel.margin.y = unit(0, "lines"))+
  scale_color_manual(values = c("#dbccf6","#c9b2f1","#b799ed", "#a57fe8","#8466ba","#7359a2","#634c8b","#534074","#42335d","#312646"))
#format graph to fit each facet to its relative size
gp <- ggplotGrob(p)
gtable::gtable_show_layout(gp)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
grid::grid.draw(gp)
#### 1.3 Mammals: No Rats ####
norats <- mammal1 %>% filter(scientific_name_OTL != "Rattus norvegicus")
aggregate(Bispheno_corrected_ug ~exposure_method+bisphenol,norats,mean)
nrow(norats) # number of effect sizes = 38
length(unique(norats$study_ID)) # number of studies = 4
length(unique(norats$molecule)) # number of molecules = 6
norats.horm <- rma.mv(yi = lnRR, V = lnRR_variance2,
                      mods= ~molecule,
                      random = list(~1 | scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID),
                      method = "REML", data = norats)
summary(norats.horm)
norats2 <- data.frame(trait = substr(row.names(norats.horm$b), 6, 50), estimate = norats.horm$b, ci.lb = norats.horm$ci.lb, ci.ub = norats.horm$ci.ub,  pval = norats.horm$pval)
norat<-ggplot(norats2, aes(x = trait, y = estimate, colour = trait))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(size = 3.5, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab(NULL) +
  ylim(-1,1) +
  scale_x_discrete(labels=c("pt"="E2","uleFollicle Stimulating Hormone"="FSH",
                            "uleLuteinizing Hormone"="LH","uleTestosterone"="T",
                            "uleThyroxine"="T4","uleTriiodothyronine"="T3")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18,face="bold", vjust=0.8, hjust = 0.5),
        axis.text.y = element_text(size=18)) +
  scale_color_manual(values = c("#99cc99","#66b366","#4da64d","#339933","#008000","#006600","#005a00","#003300"))
#pdf(file ="~/Documents/Sydney/Research/Dissertation/Experiments/Meta-analysis/All_NoRats.pdf", width = 8, height = 5)
fig1<-ggarrange(all,norat, nrow=2,ncol = 1, align = "v")
annotate_figure(fig1, left = text_grob("Effect Size (lnRR)", color = "black",face="bold", size=18, rot = 90))
#dev.off()
#####1.4 Fishes ####
fish1<-subset(horm2,class =="Actinopterygii")
#run both lines below to get correct units (Bispheno_corrected_ug is normalized unit it data set but need units for dietary exposure)
aggregate(bisphenol_conc ~exposure_method+bisphenol+bispenol_unit_standardized,fish1,mean)
aggregate(Bispheno_corrected_ug  ~exposure_method+bisphenol,fish1,mean)
nrow(fish1) # number of effect sizes = 152
length(unique(fish1$study_ID)) # number of studies = 11
length(unique(fish1$molecule)) # number of molecules = 7
fish.horm <- rma.mv(yi = lnRR, V = lnRR_variance2,
                       mods= ~molecule,
                       random = list(~1 | scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID),
                       method = "REML", data = fish1)
summary(fish.horm)
data.fish <- data.frame(trait = substr(row.names(fish.horm $b), 6, 50), estimate = fish.horm$b, ci.lb = fish.horm $ci.lb, ci.ub = fish.horm $ci.ub, p=fish.horm$pval)
#write_xlsx(data.fish,'Fish_Data.xlsx')
ggplot(data.fish, aes(x = trait, y = estimate,  colour = trait))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(size = 3.5, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab("Effect Size (lnRR)") +  ylim(-1.2,1.2) +
  scale_x_discrete(labels=c("pt"="E2","uleFollicle Stimulating Hormone"="FSH",
                            "uleLuteinizing Hormone"="LH","uleTestosterone"="T","uleThyroid Stimulating Hormone"= "TSH",
                            "uleThyroxine"="T4","uleTriiodothyronine"="T3")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18,face="bold", vjust=0.5, hjust = 0.5),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18, face="bold")) +
  scale_color_manual(values = c("#f8c194","#f7b57e","#f5a869","#f49c53","#f3903e","#f28429","#d97624","#c16920"))
####2 Bisphenols###
#Look at overall impacts of bisphenol type
bisph2<-rma.mv(yi = lnRR, V = lnRR_variance2,
               mods = ~bisphenol,
               random = list( ~1 | study_ID, ~1 | effect_size_ID),
               method = "REML", data = horm2)
summary(bisph2)
data.bisph <- data.frame(trait = substr(row.names(bisph2$b), 6, 50), estimate = bisph2$b, ci.lb = bisph2$ci.lb, ci.ub = bisph2$ci.ub, p=bisph2$pval)
write_xlsx(data.bisph, "Bisphenol_Type.xlsx")
#subset by bisphenol type (BPA, BPS, BPF, and BPAF)
bpa <-subset(horm2,bisphenol =="BPA")
bps <-subset(horm2,bisphenol =="BPS")
bpf <-subset(horm2,bisphenol =="BPF")
bpaf <-subset(horm2,bisphenol =="BPAF")

#### 2.1 BPA ####
nrow(bpa) # number of effect sizes = 342
length(unique(bpa$study_ID)) # number of studies = 36
length(unique(bpa$molecule)) # number of hormones = 9
bpa %>% group_by(molecule) %>% count()
bisph1 <-rma.mv(yi = lnRR, V = lnRR_variance2,
               mods = ~molecule,
               random = list(~1 |scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID),
               method = "REML", data = bpa)
summary(bisph1)
data.bpa <- data.frame(trait = substr(row.names(bisph1$b), 6, 50), estimate = bisph1$b, ci.lb = bisph1$ci.lb, ci.ub = bisph1$ci.ub, p=bisph1$pval)
#write_xlsx(data.bpa, "BPA.xlsx")
graph1<-ggplot(data.bpa, aes(x = trait, y = estimate))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(size = 3.5,colour = "#FF0000", show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),colour = "#FF0000", size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab(NULL) + ylim(-1.9,1.9) +
  scale_x_discrete(labels=c("pt"="CORT","uleEstradiol"="E2", "uleFollicle Stimulating Hormone"="FSH",
                                 "uleLuteinizing Hormone"="LH","uleProgesterone"="P4", "uleTestosterone"="T","uleThyroid Stimulating Hormone"= "TSH",
                                 "uleThyroxine"="T4","uleTriiodothyronine"="T3")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18,face="bold", vjust=0.5, hjust = 0.5),
        axis.text.y = element_text(size=18))
graph1
#### 2.2 BPS ####
nrow(bps) # number of effect sizes = 83
length(unique(bps$study_ID)) # number of studies = 5
length(unique(bps$molecule)) # number of hormones = 7
bps %>% group_by(molecule) %>% count()
bisph2<-rma.mv(yi = lnRR, V = lnRR_variance2,
               mods = ~molecule,
               random = list(~1 |scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID),
               method = "REML", data = bps)
summary(bisph2)
data.bps <- data.frame(trait = substr(row.names(bisph2$b), 6, 50), estimate = bisph2$b, ci.lb = bisph2$ci.lb, ci.ub = bisph2$ci.ub, p=bisph2$pval)
#write_xlsx(data.bps, "BPS.xlsx")
graph2<-ggplot(data.bps, aes(x = trait, y = estimate))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(color="#00A08A",size = 3.5, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),color="#00A08A", size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab(NULL) + ylim(-1.9,1.9) +
  scale_x_discrete(labels=c("pt"="E2","uleFollicle Stimulating Hormone"="FSH", "uleLuteinizing Hormone"="LH",
                            "uleProgesterone"="P4","uleTestosterone"="T","uleThyroxine"="T4",
                            "uleTriiodothyronine"="T3")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18,face="bold", vjust=0.5, hjust = 0.5),
        axis.text.y = element_text(size=14))
#### 2.3 BPF ####
nrow(bpf) # number of effect sizes = 46
length(unique(bpf$study_ID)) # number of studies = 4
length(unique(bpf$molecule)) # number of hormones = 7
bpf %>% group_by(molecule) %>% count()
bisph3<-rma.mv(yi = lnRR, V = lnRR_variance2,
               mods = ~molecule,
               random = list(~1 |scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID),
               method = "REML", data = bpf)
summary(bisph3)
data.bpf <- data.frame(trait = substr(row.names(bisph3$b), 6, 50), estimate = bisph3$b, ci.lb = bisph3$ci.lb, ci.ub = bisph3$ci.ub, p=bisph3$pval)
#write_xlsx(data.bpf, "BPF.xlsx")
graph3<-ggplot(data.bpf, aes(x = trait, y = estimate))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(size = 3.5, colour = "#5BBCD6", show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),colour = "#5BBCD6", size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab(NULL) + ylim(-1.9,1.9) +
  scale_x_discrete(labels=c("pt"="E2","uleFollicle Stimulating Hormone"="FSH", 
                            "uleLuteinizing Hormone"="LH",
                            "uleTestosterone"="T","uleThyroid Stimulating Hormone"="TSH", "uleThyroxine"="T4",
                            "uleTriiodothyronine"="T3")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18,face="bold", vjust=0.5, hjust = 0.5),
        axis.text.y = element_text(size=14))
#### 2.4 BPAF ####
nrow(bpaf) # number of effect sizes = 43
length(unique(bpaf$study_ID)) # number of studies = 4
length(unique(bpaf$molecule)) # number of hormones = 6
bpaf %>% group_by(molecule) %>% count()
bisph4<-rma.mv(yi = lnRR, V = lnRR_variance2,
               mods = ~molecule,
               random = list(~1 |scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID),
               method = "REML", data = bpaf)
summary(bisph4)
data.bpaf <- data.frame(trait = substr(row.names(bisph4$b), 6, 50), estimate = bisph4$b, ci.lb = bisph4$ci.lb, ci.ub = bisph4$ci.ub, p=bisph4$pval)
#write_xlsx(data.bpaf, "BPAF.xlsx")
graph4<-ggplot(data.bpaf, aes(x = trait, y = estimate, ))+geom_hline(yintercept=0, linetype = "dashed") +  geom_point(size = 3.5,colour = "#F98400", show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),colour = "#F98400", size = 1.5, width=0.1, show.legend = FALSE) + xlab(NULL) + ylab(NULL) + ylim(-1.9,1.9) +
  scale_x_discrete(labels=c("pt"="E2","uleFollicle Stimulating Hormone"="FSH", "uleLuteinizing Hormone"="LH","uleTestosterone"="T", "uleThyroxine"="T4","uleTriiodothyronine"="T3")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=18,face="bold", vjust=0.5, hjust = 0.5),
        axis.text.y = element_text(size=14))
bisphenol<-ggarrange(graph1, graph2,graph3,graph4, ncol=1, nrow=4)
annotate_figure(bisphenol, left = text_grob("Effect Size (lnRR)", color = "black",face="bold", size=18, rot = 90),)

