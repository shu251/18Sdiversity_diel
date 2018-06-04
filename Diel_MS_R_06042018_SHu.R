setwd("~/Desktop/Projects/HOT/Diel_18S/data_analysis_05022018/")
# Load libraries and raw OTU table
library(plyr)
library(dplyr)
library(vegan)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Last updated 06-04-2018

# Set-up ------------------------------------------------------------------
count<-read.table('V4_tagseq_diel.txt', sep="\t", header=TRUE,skip=1, comment.char = "")
colnames(count)[1]<-"OTU.ID"

# Name schematic for diel study
samplenum<-as.character(c(1:19))
TOD<-c("6PM","10PM","2AM","6AM","10AM","2PM","6PM","10PM","2AM","6AM","10AM","2PM","6PM","10PM","2AM","6AM","10AM","2PM","6PM")
timestamp<-c("6PM","10PM","2AM","6AM","10AM","2PM","6PM","10PM","2AM","6AM","10AM","2PM","6PM","10PM","2AM","6AM","10AM","2PM","6PM")

# After using the PR2 database, save this custom function to manually curate list of taxonomy results.
# Works by taking a column named "PR2"
### separating it out by levels and then creates a new
### column that 'renames' a simplified breakdown at a 
#### Taxa level 1 and 2.

pr2_rename_taxa_w2<-function(df){
  library(reshape2)
  split<-colsplit(df$taxonomy, "; ", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8"))
  split[ is.na(split) ] = "XXX"
  split[ split == "" ] = "XXX"
  split$Taxa<-"Other/unknown"
  split$Taxa[split$Level1 == "No blast hit"]="No blast hit"
  split$Taxa[split$Level1 == "Unassigned"]="Unassigned"
  split$Taxa[split$Level1 == "None"]="Unassigned"
  split$Taxa[split$Level2=="Amoebozoa"]="Amoebozoa"
  split$Taxa[split$Level2=="Apusozoa"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota_X"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota_Mikro"]="Other/unknown"
  split$Taxa[split$Level2=="Stramenopiles"]="Stramenopiles-Other"
  split$Taxa[split$Level2=="Alveolata"]="Alveolates-Other"
  split$Taxa[split$Level2=="Opisthokonta"]="Opisthokonts-Other"
  split$Taxa[split$Level2=="Archaeplastida"]="Archaeplastids-Other"
  split$Taxa[split$Level2=="Excavata"]="Excavates"
  split$Taxa[split$Level2=="Rhizaria"]="Rhizaria-Other"
  split$Taxa[split$Level2=="Hacrobia"]="Other/unknown"
  split$Taxa[split$Level3=="Haptophyta"]="Haptophytes"
  split$Taxa[split$Level3=="Fungi"]="Opisthokont-Fungi"
  split$Taxa[split$Level3=="Metazoa"]="Opisthokont-Metazoa"
  split$Taxa[split$Level3=="Foraminifera"]="Rhizaria-Foraminifera"
  split$Taxa[split$Level3=="Dinophyta"]="Alveolates-Dinoflagellates"
  split$Taxa[split$Level4=="Syndiniales"]="Alveolates-Syndiniales"
  split$Taxa[split$Level3=="Cryptophyta"]="Cryptophytes"
  split$Taxa[split$Level3=="Ciliophora"]="Alveolates-Ciliates"
  split$Taxa[split$Level3=="Chlorophyta"]="Archaeplastids-Chlorophytes"
  split$Taxa[split$Level3=="Cercozoa"]="Rhizaria-Cercozoa"
  split$Taxa[split$Level4=="Acantharea"]="Rhizaria-Acantharia"
  split$Taxa[split$Level4=="Chrysophyceae"]="Stramenopiles-Chrysophytes"
  split$Taxa[split$Level4=="Chrysophyceae-Synurophyceae"]="Stramenopiles-Chrysophytes"
  split$Taxa[split$Level4=="Pelagophyceae"]="Stramenopiles-Pelagophytes"
  split$Taxa[split$Level4=="Bacillariophyta"]="Stramenopiles-Diatoms"
  split$Taxa[split$Level4=="MAST"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="Polycystinea"]="Rhizaria-Polycystines"
  split$Taxa[split$Level4=="RAD-C"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-B"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-A"]="Rhizaria-RAD (A,B,C)"
  split$Taxa2<-"XXX"
  four<-c("Alveolates-Ciliates","Archaeplastids-Chlorophytes")
  split$Taxa2<-with(split, ifelse(Taxa %in% four, Level4, Taxa2)) #Take taxa named in "four" and place in the level4 name, this is the next level of tax resolution I would like to show.
  five<-c("Alveolates-Syndiniales","Stramenopiles-MAST")
  split$Taxa2<-with(split, ifelse(Taxa %in% five, Level5, Taxa2))
  two<-c("Alveolates-Other", "Other/unknown")
  split$Taxa2<-with(split, ifelse(Taxa %in% two, Level3, Taxa2))
  six<-c("Haptophytes", "Crytophytes")
  split$Taxa2<-with(split, ifelse(Taxa %in% six, Level6, Taxa2))
  four<-c("Opisthokont-Metazoa", "Opisthokonts-Other")
  split$Taxa2<-with(split, ifelse(Taxa %in% four, Level4, Taxa2))
  seven<-c("Alveolates-Dinoflagellates","Archaeplastids-Other","Excavates","Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-Chrysophytes","Stramenopiles-Other")
  split$Taxa2<-with(split, ifelse(Taxa %in% seven, Level7, Taxa2))
  # head(split)
  # unique(split$Taxa2)
  split$Taxa2<-gsub("_XXX", "", split$Taxa2)
  split$Taxa2<-gsub("_XX", "", split$Taxa2)
  split$Taxa2<-gsub("_X", "", split$Taxa2)
  #If taxa2 = "XXX" replace with Taxa, and other
  split$Taxa2<-gsub("XXX","Other-unclassified",split$Taxa2)
  head(split)
  return(split)
} 

# Get stats from raw OTU table. See Methods section in manuscript, protocols.io or github for how OTU table was generated

# Number of sequences per sample
colsum<-apply(count[2:39],2,sum); colsum
mean(colsum)
# Remove global singletons
rowsum<-apply(count[2:39],1,sum) 
count.no1 = count[ rowsum>1, ] 
# count.no1 = count[ rowsum>2, ] #option for doubleton OTUs
dim(count)[1] - dim(count.no1)[1] #total number of singletons removed

# Reports total number of unique OTUs
length(unique(count.no1$OTU.ID))

# Continue obtaining raw stats so we can generate supplementary figures to evaluate all samples
# OTU counts and information:
counts_only<-count.no1[2:39]
seq_total<-apply(counts_only,2,sum)
min(seq_total); max(seq_total); mean(seq_total)
OTU_count<-colSums(counts_only>0); OTU_count
min(OTU_count); max(OTU_count); mean(OTU_count)
OTU_single<-colSums(counts_only==1) #by sample
OTU_double<-colSums(counts_only==2) #by sample
OTU_true<-colSums(counts_only>2)
#
#Compile sample information
sample_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true);sample_info
#Visual representation of OTUs:
m.count<-melt(count.no1[2:39])
m.count$labels<-colsplit(m.count$variable, "\\.",c("DataType", "Material", "Num"))
sample_info$samples<-row.names(sample_info)
head(sample_info)
allM<-melt(sample_info[c(6,1:4)])
label<-colsplit(allM$samples, "\\.",c("DataType", "Material", "Num"))
allM<-cbind(allM,label); head(allM)
tapply(allM$value, list(allM$Material, allM$variable), mean)
samplenum<-as.character(c(1:19))
timestamp_actual<-c("18:57","22:36","2:51","6:49","12:52","14:54","18:46","22:42","2:57","6:56","10:52","14:46","18:57b","22:52","2:47","6:47","10:49","14:51","18:53")
allM$Time<-factor(allM$Num, levels=samplenum)

# Distribution of OTU sizes
bar_stats<- ggplot(allM, aes(x=Time, y=value, fill=variable))+
  geom_bar(stat="identity",position="stack",color="black")+
  labs(title="", x="",y="Total OTUs")+theme_bw()+
  scale_x_discrete(limits=1:19,labels=timestamp,expand = c(0, 0))+theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1,size=12, color = "black"),axis.text.y=element_text(size=12,color = "black"), strip.background = element_blank())+ theme(legend.title = element_blank())
OTUs<-c("OTU_count", "OTU_single", "OTU_double")
# OTU distribution
bar_stats %+% subset(allM, variable %in% OTUs)+ facet_wrap(~Material, ncol = 1, scales = "free_x")+scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("Singletons", "Doubletons","OTUs > 2 seqs"))

# Total sequences per sample
bar_stats %+% subset(allM, variable %in% "seq_total")+ facet_wrap(~Material, ncol = 1, scales = "free_x")+scale_fill_manual(values=c("#99d8c9"),labels = c("Total sequences"))

# Make Figure S2
library(cowplot)
figS2<-plot_grid(bar_stats %+% subset(allM, variable %in% "seq_total")+ facet_wrap(~Material, ncol = 1, scales = "free_x")+scale_fill_manual(values=c("#99d8c9"),labels = c("Total sequences")), labels=c("A", "B"),bar_stats %+% subset(allM, variable %in% OTUs)+ facet_wrap(~Material, ncol = 1, scales = "free_x")+scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("Singletons", "Doubletons","OTUs > 2 seqs")), align="h", ncol=2, nrow=1)
figS2 #svg save H:630, W:110
#
head(allM)
rna<-subset(allM, Material %in% "RNA")
# Get raw stats of table:
library(dplyr)
sequence_results<- allM[c(2:3,5)] %>% 
  group_by(Material, variable) %>%
  summarise_all(funs(mean(value), min(value), max(value))) %>%
  as.data.frame
sequence_results

# Subsampled data
head(count.no1);names(count.no1) # Used count.no1 from above
# Note OTU.IDs should be first column
# Save key separately, meaning the taxonomy IDs with OTUIDs
key<-count.no1[c(1,40)]; head(key)
# Set row names as OTU.IDs and remove this column
row.names(count.no1)<-count.no1$OTU.ID
keep<-count.no1[2:39] # This last step serves to make whole dataframe numeric.

# Total number of sequences that is lowest among all samples
sub<-min(colSums(keep));sub

# Rarefy sequence total to same number. This is the proper way to evaluate species richness
rare <- rrarefy(t(keep), sub)
subsampled<-as.data.frame(t(rare))
colSums(subsampled)
dim(subsampled); dim(count.no1) # total number of OTUs does not change
head(subsampled)
save(subsampled, file="subsampled_diel.RData")
# Updated 5-03-2018


# Fig2 --------------------------------------------------------------------
# Diversity analysis
head(subsampled[1:2,])
head(key[1:2,])

# Get OTU count information & format
OTU_count_rare<-colSums(subsampled>0)
min(OTU_count_rare); max(OTU_count_rare)
meltcount_rare<-melt(OTU_count_rare)
meltcount_rare$samples<-row.names(meltcount_rare)
var<-colsplit(meltcount_rare$samples, "\\.", c("DataType", "Material", "Num"))
count_all<-data.frame(var,meltcount_rare)

# Factor and plot
count_all$Time<-factor(count_all$Num, levels=samplenum)
head(count_all)
plot_div_rare <- ggplot(count_all[order(count_all$Time),], aes(x=Time, y=value,group=Material, color=Material)) +
    geom_line(stat = "identity", position = "identity", size=1.5)+geom_point(size=2.5,shape=21, color="white", aes(fill=Material))+
    scale_fill_manual(values=c("black", "red"))+scale_color_manual(values=c("black", "red"))+
    labs(title="",x="",y="Total OTUs") +
    theme(axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=45,hjust = 1,vjust = 1, color = "black"))+
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=4,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_rect(data=NULL,aes(xmin=7,xmax=10,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_rect(data=NULL,aes(xmin=13,xmax=16,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_rect(data=NULL,aes(xmin=19,xmax=Inf,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.text.y = element_text(angle=0,vjust=1, size=12, face="bold"))+
    scale_x_discrete(limits=c(), expand = c(0, 0),labels=timestamp)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ theme(legend.title = element_blank(),legend.background = element_rect(color="white"),legend.box.background = element_rect(colour = "black"))+ theme(legend.position = c(0.65, 0.25))
#
# Figure 2
plot_div_rare #w:700, H:315

head(count_all)

# last updated 05-03-2018


# Fig 3 -------------------------------------------------------------------
# Plot OTU richness by taxonomic group
head(subsampled)
subsampled$OTU.ID<-row.names(subsampled)
# Melt & convert to binary
submelt<-melt(subsampled)
submelt$bin<-ifelse(submelt$value > 0, 1, 0)
head(submelt)
binary<-dcast(submelt[c(1:2,4)], OTU.ID~variable,fill = 0)
# Add in taxonomic groups
bin_tax<-join(binary, key, by="OTU.ID", type="left", match="first")
length(unique(bin_tax$OTU.ID)); dim(bin_tax) #both should be equal to 3353

# Reformat for plotting
head(bin_tax)
tax_names<-pr2_rename_taxa_w2(bin_tax) # fxn from above
# write.csv(tax_names, file="taxnames.csv")
bin_renametax<-data.frame(bin_tax, tax_names)
# write.csv(bin_renametax, file="tmp.csv")
bin_renametax_melt<-melt(bin_renametax) #melt
head(bin_renametax_melt)
var<-colsplit(bin_renametax_melt$variable, "\\.", c("Type", "Material", "Num"))
bin_renametax_melt<-data.frame(bin_renametax_melt,var)
binary_OTUcount<-aggregate(bin_renametax_melt$value, by=list(Taxa=bin_renametax_melt$Taxa,Material=bin_renametax_melt$Material,Num=bin_renametax_melt$Num),sum)
#
head(binary_OTUcount)
unique(binary_OTUcount$Taxa)
# Factor dfs, order taxonomy and assign colors
tax_order<-c("Alveolates-Dinoflagellates","Alveolates-Ciliates","Alveolates-Syndiniales","Alveolates-Other","Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-MAST","Stramenopiles-Chrysophytes","Stramenopiles-Other","Archaeplastids-Chlorophytes","Archaeplastids-Other","Cryptophytes","Excavates","Haptophytes","Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
tax_color<-c("#67000d","#e31a1c","#dd3497","#fcbba1","#fed976","#fc8d59","#a63603","#addd8e","#7f2704","#238b45","#a1d99b","#081d58","#1f78b4","#a6cee3","#8c6bb1","#9e9ac8","#984ea3","#081d58","#662506","#ffffff","#969696","#525252","#000000")
names(tax_color)<-tax_order
binary_OTUcount$tax<-factor(binary_OTUcount$Taxa, levels=rev(tax_order))
binary_OTUcount$Time<-factor(binary_OTUcount$Num, levels=samplenum) # from name schematic above

# Remove unwanted groups
# data.agg.sub <- subset(data.agg, !grepl("Other", data.agg$Taxa)) #remove names with "Other"
# data.agg.sub <- subset(data.agg.sub, !grepl("Opisthokont", data.agg.sub$Taxa))
# rm<-c("Unassigned")
# data.agg.sub<-subset(data.agg.sub, !(Taxa %in% rm))
# unique(data.agg.sub$Taxa)

# bar plots
head(binary_OTUcount)
binary_OTUcount
plot_OTUcount_tax<-ggplot(binary_OTUcount[order(binary_OTUcount$tax),], aes(y=x,fill=tax,order=tax))+ 
  geom_bar(aes(fill= tax,x=as.numeric(Num)), position = "stack", stat = "identity", color="black")+
  facet_grid(Material~.)+ 
  scale_fill_manual(values=tax_color)+
  labs(title="", x="",y="Total OTUs by taxa")+
  theme(legend.position="right",legend.title = element_blank(),plot.title=element_text(hjust = 0,face='bold',size=20))+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=45,hjust = 1,vjust = 1, color = "black"), axis.text.y=element_text(color="black"))+
  scale_x_discrete(limits=1:19,labels=timestamp,expand = c(0, 0))+
  geom_segment(aes(x=0.5,xend=3.5,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=6.5,xend=9.5,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=12.5,xend=15.5,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=18.5,xend=19.5,y=-0.02,yend=-0.02),size=1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  scale_y_continuous(expand = c(0, 0),breaks = scales::pretty_breaks(n = 10))
# Figure S3 W: 970 H: 600
plot_OTUcount_tax %+% subset(binary_OTUcount, !(tax %in% NA))
# 
meta<-c("Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned", "Alveolates-Other","Archaeplastids-Other", "Stramenopiles-Other", NA)
alveolates<-c("Alveolates-Dinoflagellates","Alveolates-Ciliates","Alveolates-Syndiniales")
plot_OTUcount_tax %+% subset(binary_OTUcount, !(tax %in% alveolates) & !(tax %in% meta))+labs(title="non-Alveolates")
#
plot_OTUcount_tax %+% subset(binary_OTUcount, tax %in% alveolates & !(tax %in% meta))+labs(title="Alveolates")
#
alv<-get_legend(plot_OTUcount_tax %+% subset(binary_OTUcount, tax %in% alveolates & !(tax %in% meta))+labs(title="Alveolates"))
non_avl<-get_legend(plot_OTUcount_tax %+% subset(binary_OTUcount, !(tax %in% alveolates) & !(tax %in% meta))+labs(title="non-Alveolates"))
#
plots<-plot_grid(plot_OTUcount_tax %+% subset(binary_OTUcount, tax %in% alveolates & !(tax %in% meta))+theme(legend.position = "none"), plot_OTUcount_tax %+% subset(binary_OTUcount, !(tax %in% alveolates) & !(tax %in% meta))+theme(legend.position = "none"),alv, non_avl,nrow=2, ncol=2, labels=c("A", "B", "C", "D"))
plots # Figure 3 - W:1070 H:1200

# last updated SHu - 05-03-2018
#
# head(binary_OTUcount)
# wide<-dcast(binary_OTUcount[1:4], Taxa+Material~Num)
# head(wide)
# write.csv(wide, file="OTUcount_bytax.csv")


# Fig 4 -------------------------------------------------------------------
# Moving on to characterizing whole community diversity
load("Checkpoint1_dfs_dielMS.RData", verbose=T) #optional to pick up from here
count.subsampled<-join(subsampled, key, by="OTU.ID", type="left", match="first")
#head(count.subsampled)
#dim(count.subsampled)
binned<-pr2_rename_taxa_w2(count.subsampled) # fxn from above
# unique(binned$Taxa)
# unique(binned$Taxa2)
data_binned<-data.frame(count.subsampled, binned)

#total OTUs in each taxonomic group
unique(data_binned$Taxa)
x<-as.data.frame(table(data_binned$Taxa))
colSums(x[2]) #equals total OTUs
colnames(x)[2]<-"Total OTUs"
x # table 1 so far
# Add this table to Table 1 of main text (Total OTUs per taxonomic group)

# Re-format
# Convert to long form:
data.all.melt<-melt(data_binned) #melt
# head(data.all.melt)

# Separate by variables
var<-colsplit(data.all.melt$variable, "\\.", c("Type", "Material", "Num"))
data.m<-data.frame(data.all.melt,var)

# Aggregate by major taxonomic group, material, and time of day
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$Taxa,Material=data.m$Material,Num=data.m$Num),sum)

# Aggregate by major taxonomic group, taxa level 2, material, and time of day
data.agg.tax2<-aggregate(data.m$value, by=list(Taxa=data.m$Taxa, Taxa2=data.m$Taxa2,Material=data.m$Material,Num=data.m$Num),sum)

# Save aggregated data frames
# save(data.agg, data.agg.tax2, file="Aggregated_dfs_05032018.RData")

# Generate Table S2, summary of OTU table by sequence counts. Also lists manual curation of taxonomic groups
head(data.m[1:2,])
names(data.m)
tmp<-data.m[c(2,11:12,14,16)]
head(tmp)
dna<-subset(tmp, Material %in% "DNA")
rna<-subset(tmp, Material %in% "RNA")
head(dna)
# Sum to show all taxonomic information
all.dna<-aggregate(dna$value, by=list(PR2=dna$taxonomy,Taxa=dna$Taxa,Taxa2=dna$Taxa2), sum)
all.rna<-aggregate(rna$value, by=list(PR2=rna$taxonomy,Taxa=rna$Taxa,Taxa2=rna$Taxa2), sum)

# Characterize and save
head(all.dna[1:2,]); dim(all.dna); colnames(all.dna)[4]<-"DNA"
head(all.rna[1:2,]); dim(all.rna); colnames(all.rna)[4]<-"RNA"
# Join both
all<-join(all.dna, all.rna, type="left", match="first")
head(all[1:2,])
# Table S2
# write.table(all, file="TableS2_seqcounts.txt",quote=FALSE, sep="\t", row.names=FALSE)

# Aggregate by major taxonomic group and material
# tax_abun<-aggregate(data.m$value, by=list(Taxa=data.m$Taxa,Material=data.m$Material),sum)
# write.csv(tax_abun, file="taxabundance.csv")


# load("Aggregated_dfs_05032018.RData", verbose=T)

# Factor dfs, order taxonomy and assign colors
tax_order<-c("Alveolates-Dinoflagellates","Alveolates-Ciliates","Alveolates-Syndiniales","Alveolates-Other","Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-MAST","Stramenopiles-Chrysophytes","Stramenopiles-Other","Archaeplastids-Chlorophytes","Archaeplastids-Other","Cryptophytes","Excavates","Haptophytes","Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
tax_color<-c("#67000d","#e31a1c","#dd3497","#fcbba1","#fed976","#fc8d59","#a63603","#addd8e","#7f2704","#238b45","#a1d99b","#081d58","#1f78b4","#a6cee3","#8c6bb1","#9e9ac8","#984ea3","#081d58","#662506","#ffffff","#969696","#525252","#000000")
names(tax_color)<-tax_order
data.agg$tax<-factor(data.agg$Taxa, levels=rev(tax_order))
data.agg$Time<-factor(data.agg$Num, levels=samplenum) # from name schematic above

# Remove unwanted groups
data.agg.sub <- subset(data.agg, !grepl("Other", data.agg$Taxa)) #remove names with "Other"
data.agg.sub <- subset(data.agg.sub, !grepl("Opisthokont", data.agg.sub$Taxa))
rm<-c("Unassigned", "None", NA)
data.agg.sub<-subset(data.agg.sub, !(Taxa %in% rm))
unique(data.agg.sub$Taxa)

# Area plot
plot_area_sub<-ggplot(data.agg.sub[order(data.agg.sub$tax),], aes(y=x,fill=tax,order=tax))+ 
  geom_area(aes(fill= tax,x=as.numeric(Num), color=tax), position = "fill", stat = "identity")+
  facet_grid(Material~.)+ 
  scale_color_manual(values=tax_color)+
  scale_fill_manual(values=tax_color)+
  labs(title="", x="",y="Relative abundance of reads")+
  theme(legend.position="right",legend.title = element_blank(),plot.title=element_text(hjust = 0,face='bold',size=20))+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=45,hjust = 1,vjust = 1, color = "black"), axis.text.y=element_text(color="black"))+
  scale_x_discrete(limits=1:19,labels=timestamp,expand = c(0, 0))+
  geom_segment(aes(x=1,xend=4,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=1,xend=4,y=1.02,yend=1.02),size=1.5)+
  geom_segment(aes(x=7,xend=10,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=7,xend=10,y=1.02,yend=1.02),size=1.5)+
  geom_segment(aes(x=13,xend=16,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=13,xend=16,y=1.02,yend=1.02),size=1.5)+
  geom_segment(aes(x=19,xend=19.01,y=-0.02,yend=-0.02),size=1.5)+
  geom_segment(aes(x=19,xend=19.01,y=1.02,yend=1.02),size=1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  scale_y_continuous(expand = c(0, 0)) 
#plot_area_sub #Figure 2-all

# Figure 2 main text - only taxa with >0.1% of total
maj<-c("Rhizaria-RAD (A,B,C)","Stramenopiles-Pelagophytes","Archaeplastids-Chlorophytes","Opisthokont-Fungi","Opisthokont-Metazoa","Stramenopiles-Chrysophytes","Haptophytes","Alveolates-Syndiniales","Rhizaria-Acantharia","Stramenopiles-Diatoms","Opisthokonts-Other","Stramenopiles-Other","Stramenopiles-MAST","Alveolates-Ciliates","Alveolates-Dinoflagellates")
plot_area_sub %+% subset(data.agg.sub, tax %in% maj)
# save svg: W: 875, H: 600

# last updated 05-03-2018 SHu


# Fig 5-6 (ratio calc) ----------------------------------------------------
# Calculating RNA:DNA ratios. RNA sequences:DNA sequences per OTU. 
## If an OTU does NOT have RNA or DNA sequences for the same OTU, that OTU was removed
head(data_binned[1:2,]) # df from above
mdf2<-melt(data_binned)
var<-colsplit(mdf2$variable, "\\.", c("DataType", "Material", "Num"))
mdf2.vars<-data.frame(var,mdf2)
names(mdf2.vars)
wdf2.vars<-dcast(mdf2.vars[c(4,2:3,5,6:15,17)], OTU.ID+taxonomy+Level1+Level2+Level3+Level4+Level5+Level6+Level7+Level8+Taxa+Taxa2+Num~Material)
head(wdf2.vars[1:3,])
length(unique(wdf2.vars$OTU.ID))

# Calculate ratio
wdf2.vars$Ratio<-wdf2.vars$RNA/wdf2.vars$DNA
tmp<-do.call(data.frame,lapply(wdf2.vars, function(x) replace(x, is.infinite(x),NA)))
tmp1<-do.call(data.frame,lapply(tmp, function(x) replace(x, is.na(x),NA)))
tmp1[is.na(tmp1)] <- 0
# head(tmp1)
#If OTU ratio is zero, it means either RNA or DNA was absent - remove
df_wRatio<-subset(tmp1, Ratio != 0)
dim(df_wRatio)
# head(df_wRatio[1:2,])
# length(unique(df_wRatio$OTU.ID))
#
# save(key,df_wRatio, data_binned, wdf2.vars, file="Checkpoint2_dfs_wRatio_dielMS.RData")

#How many sequences start to end?
head(mdf2);names(mdf2)
head(df_wRatio);names(df_wRatio)
w<-colSums(mdf2[14])
x<-colSums(df_wRatio[14:15])
# Original total number of sequences
y<-x[1]+x[2];y
y/w*100 #Percentage of reads when we look for diel rhythmicity in OTUs with both RNA and DNA reads
length(unique(df_wRatio$OTU.ID)) # Total number of OTUs with both RNA and DNA sequences
names(df_wRatio)
unique(df_wRatio$Taxa)
uni<-unique(df_wRatio[c(1,11:12)])
#How many OTUs with Ratios? Add to Table 1 in text
names(df_wRatio)
head(uni)
unique(uni$Taxa)
a<-as.data.frame(table(uni[2]))
# write.table(a, file="Table1_OTUs_withRatios_0532018.txt",quote=FALSE, sep="\t", row.names=FALSE)

# load("Checkpoint2_dfs_wRatio_dielMS.RData", verbose=T)

# Factor dataframe:
names(df_wRatio)
df_wRatio2<-df_wRatio[c(1:2,11:16)] # rename/reorder
head(df_wRatio2[1:2,])

#Re-label with new times of day. Because I want to plot a single day period
df_wRatio2$Time<-factor(df_wRatio2$Num, levels=samplenum, labels = TOD) 
# Error will be produced: “duplicated levels in factors are deprecated”
## this is because we have repeated labels - Ignore.

# Change factors to characters
df_wRatio2$Time<-as.character(df_wRatio2$Time) #make into character not factor.
df_wRatio2$Taxa<-as.character(df_wRatio2$Taxa) #make into character not factor.
df_wRatio2$Num<-as.character(df_wRatio2$Num) #make into character not factor.
df_wRatio2$Taxa2<-as.character(df_wRatio2$Taxa2) #make into character not factor.
unique(df_wRatio2$Time) # these are the times of day

# Calculate mean and standard mean error
names(df_wRatio2)
library(dplyr)
meantax_er <- df_wRatio2[c(3,8:9)] %>% group_by(Taxa,Time) %>% summarise_all(funs(mean,length, sd, sem_ratio=sd(Ratio)/sqrt(length(Ratio))))
#older version had me calculating DNA and RNA sem.. but I removed this.
#
meantax_er<-as.data.frame(meantax_er)
names(meantax_er)
#Calculate low and high based on sem
meantax_er$ratio_low<-meantax_er$mean - meantax_er$sem_ratio
meantax_er$ratio_high<-meantax_er$mean + meantax_er$sem_ratio
colnames(meantax_er)[3]<-"Ratio_mean"
#
meantax_er$time_order<-factor(meantax_er$Time, levels=c("6AM","10AM","2PM","6PM","10PM","2AM"))
# Dataframe has the mean and standard mean error for each major taxonomic group.
head(meantax_er[1:2,])
# Factor colors for manuscript publication
tax_order<-c("Alveolates-Dinoflagellates","Alveolates-Ciliates","Alveolates-Syndiniales","Alveolates-Other","Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-MAST","Stramenopiles-Chrysophytes","Stramenopiles-Other","Archaeplastids-Chlorophytes","Archaeplastids-Other","Cryptophytes","Excavates","Haptophytes","Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
tax_color<-c("#67000d","#e31a1c","#dd3497","#fcbba1","#fed976","#fc8d59","#a63603","#addd8e","#7f2704","#238b45","#a1d99b","#081d58","#1f78b4","#a6cee3","#8c6bb1","#9e9ac8","#984ea3","#081d58","#662506","#ffffff","#969696","#525252","#000000")
label_panels<-c("A. Alveolates-Dinoflagellates","A. Alveolates-Ciliates","B. Alveolates-Syndiniales","Alveolates-Other","B. Stramenopiles-Diatoms","C. Stramenopiles-Pelagophytes","C. Stramenopiles-MAST","A. Stramenopiles-Chrysophytes","Stramenopiles-Other","D. Archaeplastids-Chlorophytes","Archaeplastids-Other","Cryptophytes","Excavates","E. Haptophytes","D. Rhizaria-Acantharia","B. Rhizaria-Cercozoa","C. Rhizaria-Polycystines","E. Rhizaria-RAD (A,B,C)","Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
meantax_er$taxfig<-factor(meantax_er$Taxa, levels=tax_order, labels = labels)
names(tax_color)<-label_panels

# Need to repeat 6AM at end for cyclical representation
sub6am<-subset(meantax_er, time_order %in% "6AM")
sub6am$Time<-NULL; sub6am$Time<-"6AMb"
meantax_er_6amdup<-rbind(meantax_er, sub6am)
#head(meantax_er_6amdup)
# Refactor:
#unique(meantax_er_6amdup$Time)
meantax_er_6amdup$time_order<-factor(meantax_er_6amdup$Time, levels=c("6AM","10AM","2PM","6PM","10PM","2AM", "6AMb"))
unique(meantax_er_6amdup$time_order)

#Ribbon plots of mean and sem RNA:DNA ratios:
lineOTU_rib <- ggplot(meantax_er_6amdup[order(meantax_er_6amdup$time_order),], aes(x=time_order, y=Ratio_mean, group=taxfig))+
  scale_x_discrete(limits=c(), expand = c(0, 0), labels=c("6AM","10AM","2PM","6PM","10PM","2AM", "6AM"))+
  geom_ribbon(aes(ymin=ratio_low, ymax=ratio_high, fill=taxfig),alpha=0.6, color=NA)+
  labs(title="Mean & standard error of the mean for OTUs\n in each taxonomic group",x="",y="RNA:DNA read ratio") +
  theme(axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=45,hjust = 1,vjust = 1, size=11, color="black"), axis.text.y = element_text(size=11, color="black"), strip.text = element_text(size=10, hjust=0, face="bold", color="black"))+
  scale_color_manual(values = rev(tax_color))+
  scale_fill_manual(values = rev(tax_color))+
  geom_rect(data=NULL,aes(xmin=4,xmax=7,ymin=-Inf,ymax=Inf),fill="#737373",alpha=0.1, color=NA)+ 
  geom_line(stat = "identity", position = "identity", size=2,aes(color=taxfig))+
  geom_point(size=3, shape=21, color="grey", aes(fill=taxfig))+ 
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())

# Plot by presumed trophic strategy
photo<-c("A. Alveolates-Dinoflagellates","B. Stramenopiles-Diatoms","C. Stramenopiles-Pelagophytes","D. Archaeplastids-Chlorophytes","E. Haptophytes")
nonphoto<-c("A. Alveolates-Ciliates","B. Alveolates-Syndiniales","C. Stramenopiles-MAST","D. Rhizaria-Acantharia","E. Rhizaria-RAD (A,B,C)")
other<-c("A. Stramenopiles-Chrysophytes","B. Rhizaria-Cercozoa","C. Rhizaria-Polycystines")
#
unique(meantax_er$taxfig) 
#Save below svg: W: 580, H: 640
# Figure 5
lineOTU_rib %+% subset(meantax_er_6amdup, taxfig %in% photo)+facet_wrap(~taxfig, scales="free", nrow = 3)+labs(title="")
# Figure 6
lineOTU_rib %+% subset(meantax_er_6amdup, taxfig %in% nonphoto)+facet_wrap(~taxfig, scales="free", nrow = 3)+labs(title="")
# Figure S5
lineOTU_rib %+% subset(meantax_er_6amdup, (taxfig %in% other))+facet_wrap(~taxfig, scales="free", nrow = 3)+labs(title="")
# write.csv(meantax_er_6amdup, file="tmp_checkplots.csv")
# last updated SHU 05-04-2018

# RAIN analysis -----------------------------------------------------------
## Set up RAIN analysis environment and seaprate RNA and DNA sequence libraries
?clr()
load("Checkpoint1_dfs_dielMS.RData", verbose=T) #optional to pick up from
#
library(rain)
t<-(1:19)
ft <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# Separate dataframes
subsampled$OTU.ID<-row.names(subsampled)
subsampled.m<-melt(subsampled)
head(subsampled.m)
var<-colsplit(subsampled.m$variable, "\\.", c("DataType", "Material", "Num"))
subsampled_df<-data.frame(var,subsampled.m)

# RNA
rna_df<-subset(subsampled_df, Material %in% "RNA"); rna_dfw<-dcast(rna_df[c(3:4,6)], OTU.ID~Num); row.names(rna_dfw)<-rna_dfw$OTU.ID; rna_dfw[1]<-NULL
# DNA
dna_df<-subset(subsampled_df, Material %in% "DNA"); dna_dfw<-dcast(dna_df[c(3:4,6)], OTU.ID~Num); row.names(dna_dfw)<-dna_dfw$OTU.ID; dna_dfw[1]<-NULL

#Perform clr normalizations
# RNA
head(rna_dfw)
rna_clr<-as.data.frame(clr(rna_dfw));head(rna_clr)
#DNA
head(dna_dfw)
dna_clr<-as.data.frame(clr(dna_dfw));head(dna_clr)

# Run RAIN for RNA
rna_RAIN<-rain(t(as.matrix(rna_dfw)), period=24, measure.sequence=ft, deltat=4, method="independent", na.rm = TRUE)
#head(rna_RAIN)

# Report how many OTUs had significantly diel rhythmicity
sigRNA<-subset(rna_RAIN, pVal < 0.05)

# Get list of significant OTUs from RNA library
sigRNA_list<-as.data.frame(row.names(sigRNA)); colnames(sigRNA_list)[1]<-"OTU.ID"
sigRNA_tax<-join(sigRNA_list,key); sigRNA_tax$origin<-"Sig_RNA"
head(sigRNA_tax[1:3,])
# dim(sigRNA_tax)

# Run RAIN for DNA
dna_RAIN<-rain(t(as.matrix(dna_clr)), period=24, measure.sequence=ft, deltat=4, method="independent", na.rm = TRUE)

# Report how many OTUs had significantly diel rhythmicity
sigDNA<-subset(dna_RAIN, pVal < 0.05)
# head(sigDNA)
dim(sigDNA)
# Get list of significant OTUs from DNA library
sigDNA_list<-as.data.frame(row.names(sigDNA)); colnames(sigDNA_list)[1]<-"OTU.ID"
sigDNA_tax<-join(sigDNA_list,key); sigDNA_tax$origin<-"Sig_DNA"
head(sigDNA_tax[1:3,])
dim(sigDNA_tax)
length(unique(sigDNA_tax$OTU.ID))
# Save checkpoint
# save(key,count.no1, subsampled, sigDNA_tax, dna_df, sigRNA_tax, rna_df, file="Checkpoint1_dfs_dielMS.RData")

# load("Checkpoint1_dfs_dielMS.RData", verbose=TRUE) # start here if picking up from checkpoint1
allsig<-rbind(sigDNA_tax, sigRNA_tax)
head(allsig)

# Unique OTUs that are sig for each RNA and DNA library based on OTU.ID
colSums(table(allsig[c(1,3)])) # 11 DNA OTUs, and 87 RNA OTUs
# Repeat, but save
a <-as.data.frame(table(allsig[c(1,3)]));head(a[1:2,])
b <-dcast(a, OTU.ID~origin); b$both<-paste(b$Sig_DNA + b$Sig_RNA)
head(b)

# Get taxonomic breakdown of OTUs that were significant
head(allsig[1:2,])
tmp<-pr2_rename_taxa_w2(allsig) # see function from up above
allsig_tax<-data.frame(allsig,tmp)
head(allsig_tax[1:2,])

# Now I have simplified naming schematic
unique(allsig_tax$Taxa)
unique(allsig_tax$Taxa2)
names(allsig_tax)

# Generate supplementary master table
head(allsig_tax[1:3,])
write.table(allsig_tax, file="TableS3_sigRhythmicity_06042018.txt", quote=FALSE, sep="\t", row.names=FALSE)

library(dplyr)
# Summarize and concatenate taxonomic group names for OTUs with sig rhythmicity
rna_tmp<-subset(allsig_tax, origin %in% "Sig_RNA")
dna_tmp<-subset(allsig_tax, origin %in% "Sig_DNA")

# Concatenate
rna_tmp2 <- rna_tmp %>% group_by(Taxa) %>% summarize(RNA_Taxa2 =paste(unique(Taxa2), collapse=", "),RNA_freq=length(OTU.ID)) %>% data.frame
# head(rna_tmp2)
dna_tmp2 <- dna_tmp %>% group_by(Taxa) %>% summarize(DNA_Taxa2 =paste(unique(Taxa2), collapse=", "),DNA_freq=length(OTU.ID)) %>% data.frame
# dim(dna_tmp2)
# dna_tmp2
# Join RNA & DNA to generate master table
tmp3 <-join(dna_tmp2, rna_tmp2, by="Taxa", type="full", match="all")
tmp3
# Write table - TABLE 1 in manuscript
write.table(tmp3, file="Table1_OTUs_sigRhythmicity_06042018.txt", quote=FALSE, sep="\t", row.names=FALSE)

# How many sequences do significantly diel OTUs account for?
head(rna_df)
head(sigRNA_tax)
rna_otus<-as.character(unique(sigRNA_tax$OTU.ID))
sigdielRNA_seq<-subset(rna_df, OTU.ID %in% rna_otus)
dim(rna_df); dim(sigdielRNA_seq)
sum(sigdielRNA_seq$value) / sum(rna_df$value) #% of RNA accounted for in diel
#
head(dna_df)
head(sigDNA_tax)
dna_otus<-as.character(unique(sigDNA_tax$OTU.ID))
sigdielDNA_seq<-subset(dna_df, OTU.ID %in% dna_otus)
dim(dna_df); dim(sigdielDNA_seq)
sum(sigdielDNA_seq$value) / sum(dna_df$value) #% of RNA accounted for in diel

# last updated 06-04-2018 SHu

# eLSA prep ---------------------------------------------------------------
# 
load("Checkpoint2_dfs_wRatio_dielMS.RData", verbose=T) #Option to pick up from here
head(wdf2.vars[1:2,])
names(wdf2.vars)
# Goal: Take subsampled data. Select only RNA, perform normalizations and run eLSA.

# Option 1: From the subsampled data, filter to specific OTUs for analysis and then perform the centered log ratio normalization.
# Option 2: From the subsampled data, perform centered log ratio normalization and then subsample to specific OTUs.

#Make abundance tables for RNA-derived OTUs that can be imported into LSA
RNA_cyto<-dcast(wdf2.vars[c(1:10,11:12,13,15)], OTU.ID+taxonomy+Level1+Level2+Level3+Level4+Level5+Level6+Level7+Level8+Taxa+Taxa2~Num, fill=0)
head(RNA_cyto);names(RNA_cyto)

#op1:
head(RNA_cyto[1:2,])
library(compositions)
# ?clr()
# names(RNA_cyto)
# 
# clr_op1<-clr(RNA_cyto)

## Filter OTUs:
# OTUs need to have more than 1 sequence
tmp<-subset(RNA_cyto, rowSums(RNA_cyto[13:29])>1)

# Remove metazoan sequences
tmp2<-subset(tmp, !grepl("Opisthokonts",tmp$Taxa)); dim(tmp2)
names(tmp2)
# Rename dataframe columns (better format for LSA)
colnames(tmp2)[c(1,13:31)]<-c("#","RNA_1_6PM","RNA_2_10PM","RNA_3_2AM","RNA_4_6AM","RNA_5_10AM","RNA_6_2PM","RNA_7_6PM","RNA_8_10PM","RNA_9_2AM","RNA_10_6AM","RNA_11_10AM","RNA_12_2PM","RNA_13_6PM","RNA_14_10PM","RNA_15_2AM","RNA_16_6AM","RNA_17_10AM","RNA_18_2PM","RNA_19_6PM")
LSA_rna<-tmp2[c(1,13:31)]

head(LSA_rna)
# Get key or attributes file
names(tmp)
key<-tmp[c(1:12)]
head(key)

# Sum of sequences across all samples has to be more than 10
names(LSA_rna)
rowsum<-apply(LSA_rna[2:20],1,sum)
LSA_rna_10 = LSA_rna[ rowsum>9, ]; dim(LSA_rna_10)

#Filter out OTUs so that the OTU had to appear in every single sample
df<-LSA_rna_10 #rename df
df2<-df[apply(df[2:20], MARGIN = 1, function(x) all(x > 0)), ]
head(df2);dim(df2) # Left with 201 OTUs

#Perform clr normalization
op1<-df2
row.names(op1)<-op1$`#`; op1$`#`<-NULL
head(op1)
op1_clr<-as.data.frame(clr(op1));head(op1_clr)
op1_clr$`#`<-row.names(op1_clr)
op1_clr<-op1_clr[c(20,1:19)]
head(op1_clr)
# LSA_rna10_allTimepoints<-df2 # Rename

#op2
OTUs_to_keep<-as.character(df2$`#`) # get OTUs to run.
# start from raw sequence counts to perform CLR:
names(RNA_cyto)
op2<-RNA_cyto[c(1,13:31)]
row.names(op2)<-op2$OTU.ID
op2$OTU.ID<-NULL
head(op2)
op2_clr<-as.data.frame(clr(op2))
head(op2_clr)
op2_clr$OTU.ID<-row.names(op2_clr)
op2_clr<-op2_clr[c(20,1:19)]
head(op2_clr)
colnames(op2_clr)[1:20]<-c("#","RNA_1_6PM","RNA_2_10PM","RNA_3_2AM","RNA_4_6AM","RNA_5_10AM","RNA_6_2PM","RNA_7_6PM","RNA_8_10PM","RNA_9_2AM","RNA_10_6AM","RNA_11_10AM","RNA_12_2PM","RNA_13_6PM","RNA_14_10PM","RNA_15_2AM","RNA_16_6AM","RNA_17_10AM","RNA_18_2PM","RNA_19_6PM")
#

#Import environmental data, add to end
envir<-read.delim("Envir_extrafile.txt",header=T)
head(envir)
colnames(envir)[1]<-"#"
# add environmental
head(op1_clr)
head(op2_clr)
#
LSA_RNA_wEnv_op1<-rbind(op1_clr, envir)
LSA_RNA_wEnv_op2<-rbind(op2_clr, envir)
#
# Write tables for LSA:
write.table(LSA_RNA_wEnv_op1, file="LSA_RNA_input_op1.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(LSA_RNA_wEnv_op2, file="LSA_RNA_input_op2.txt", quote=FALSE, sep="\t", row.names=FALSE)
#
#Also write attribute file:
head(key)
tmp<-as.data.frame(LSA_RNA_wEnv_op1$`#`); colnames(tmp)<-"OTU.ID"
key_filtered<-join(tmp, key)
head(key_filtered)
write.table(key_filtered, file="attributes_RNA_diel_Filtered10_all_opt1.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Export this table and run lsa_compute like so:
#
# lsa_compute LSA_RNA_inputTable_10_05032018.txt LSA_out_05032018.txt -d 7 -s 19 -r 1 -p perm -n percentileZ -b 0
#
# lsa_compute LSA_RNA_input_op1.txt LSA_out_op1.txt -d 7 -s 19 -r 1 -p perm -n percentileZ -b 0
#
# lsa_compute LSA_RNA_input_op2.txt LSA_out_op2.txt -d 7 -s 19 -r 1 -p perm -n percentileZ -b 0
#

# Parameters used:
# -d 7 - Delay of 7, means that delay correlations will be sought within a 12 hr time window (i.e. 6AM to 6PM acceptable)
# Analysis takes ~2 days to run.
#
# See cytoscape directions for how to open this file (along with attribute file) in cytoscape to look at OTU network diagrams



# Analyzing LSA output ----------------------------------------------------
# Output from LSA:
output<-read.delim("LSA_out_op1.txt",header=TRUE)
# output<-read.delim("LSA_out_05032018.txt",header=TRUE)
key<-read.delim("attributes_RNA_diel_Filtered10_all_opt1.txt")
# key<-read.delim("attributes_RNA_diel_Filtered10_all.txt")
names(key)
key<-key[c(1,11:12,2)] #reorder key file
head(output[1:2,]);dim(output)

# Re-format table to summarize LSA output
#Add in taxonomy information. 
## "Pairs" of co-occurring OTUs are matched with A and B. See labeling schematic below
colnames(output)[1:2]<-c("OTU.IDA", "OTU.IDB")
head(key[1:2,]);names(key)
key1<-key; colnames(key1)[1:4]<-c("OTU.IDA", "Taxa_A", "Taxa2_A", "PR2_A")
key2<-key; colnames(key2)[1:4]<-c("OTU.IDB", "Taxa_B", "Taxa2_B", "PR2_B")
#join attributes (key) with output information
output_tax<-join(output, key1, by="OTU.IDA", type="left")
dim(output_tax); dim(output) #same! 
output_tax2<-join(output_tax, key2, by="OTU.IDB", type="left")
dim(output_tax2)
head(output_tax2)[1:2,]

#Add environmental variables that were in place of OTU.IDs
# Change to characters first
output_tax2$OTU.IDA<-as.character(output_tax2$OTU.IDA)
output_tax2$OTU.IDB<-as.character(output_tax2$OTU.IDB)
output_tax2$Taxa_A<-as.character(output_tax2$Taxa_A)
output_tax2$Taxa_B<-as.character(output_tax2$Taxa_B)
# Replace NAs with OTU ID, which states the environmental variable
output_tax2$Taxa_A<-with(output_tax2, ifelse(Taxa_A %in% NA, OTU.IDA, Taxa_A))
output_tax2$Taxa_B<-with(output_tax2, ifelse(Taxa_B %in% NA, OTU.IDB, Taxa_B))

#Add column that lists "pair type"
output_tax2$Pair<-paste(output_tax2$Taxa_A, output_tax2$Taxa_B, sep="_")
#Based on major taxonomic groups, correlated OTUs:
length(unique(output_tax2$Pair)) #480 unique pairs
#unique(output_tax2$Pair)


# Stats and filtration
dim(output_tax2) #Starting number of interactions 22,578
names(output_tax2)
# Remove unassigned
rm<-c("Unassigned", "None")
tmp0<-subset(output_tax2, !(Taxa_A %in% rm | Taxa_B %in% rm))
dim(tmp0)

#Filter for significance:
# hist(tmp0$Q)
# hist(tmp0$P)
tmp1<-subset(tmp0, Q < 0.05 | P < 0.05) #Remove non-significant interactions 
dim(tmp1) # How many significant interactions
#
hist(tmp1$SPCC)
tmp2<-subset(tmp1, SPCC < -0.5 | SPCC > 0.5)
hist(tmp2$SPCC); dim(tmp2) # Visualize distribution
#
LSA_pairs_output<-tmp2
dim(LSA_pairs_output)

# save(output_tax2, LSA_pairs_output, file="LSA_Robjects.RData")
save(output_tax2, LSA_pairs_output, file="LSA_Robjects_060318.RData")
# Table S4
# write.table(LSA_pairs_output,file="All_sig_LSApairs.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(LSA_pairs_output,file="All_sig_LSApairs_opt1.txt", quote=FALSE, sep="\t", row.names=FALSE)
#
#
load("LSA_Robjects_060318.RData", verbose=T)
# load("LSA_Robjects.RData", verbose=T)
# Get summary of significant pairs
head(LSA_pairs_output)
names(LSA_pairs_output)
interactions<-as.data.frame(LSA_pairs_output[c(34,29,32)])
#head(interactions[1:2,])
interactions$Count<-"1"
head(interactions)
interactions$Pair2<-paste(interactions$Taxa2_A, interactions$Taxa2_B, sep="_")
#
head(interactions)
df<-interactions[c(1,4,5)]
colnames(df)[1:2]<-c("name", "value")
df$value<-as.numeric(df$value)

# Need to characterize all pairs. LSA output lists pairs regardless of order, so we need to account for this:
require(tidyverse)
require(purrrlyr)
df2<-df %>% 
  separate(col = name, sep = "_", c("A", "B")) %>% 
  by_row(.collate = "rows", 
         ..f = function(this_row) {
           paste0(sort(c(this_row$A, this_row$B)), collapse = "_")
         }) %>% 
  rename(sorted = ".out") %>%
  group_by(sorted) %>%
  summarize(sum(value))%>%
  as.data.frame
head(df2) #frequency of all pairs

#Get extra taxonomy information for each pair
library(dplyr)
df_tax2<- df %>%
group_by(name) %>% #new
  arrange(Pair2) %>% #new
  summarize(tax_detail=paste(unique(Pair2),collapse=", "),times=length(Pair2))%>%
  as.data.frame

# Filter out all pairs that were at least 1% of all co-occurring OTUs:
#head(df2) #frequency of all pairs
tmp<-subset(df2, !grepl("NA", df2$sorted))
tmp<-subset(tmp, !grepl("Unassigned", tmp$sorted))
df3<-subset(tmp, !grepl("Other/unknown", tmp$sorted))
colnames(df3)[2]<-"Freq"
#
sum_topLSA <- df3 %>%
    top_n(20) %>%
    arrange(desc(Freq)) %>%
    as.data.frame
#sum_topLSA

#Write csv files/tables for all interactions
# Table 2 in main text
write.csv(sum_topLSA, file="Top20_OTUpairs_06032018.csv")
# write.csv(df2, file="freq_allinteractions_06032018.csv")
# write.csv(df_tax2, file="alltax_combinations_06032018.csv") # Additional information regarding Table S4.


# last updated - 06-03-2018 SHu


# Select Rhizarian and Syndiniales OTU interactions -----------------------
load(file="LSA_Robjects_060318.RData", verbose=T) #Option to start here
# From these significant interactions:
# head(LSA_pairs_output[1:2,])
# Isolate Syndiniales
syn<-LSA_pairs_output %>% filter(grepl("Syndini", LSA_pairs_output$Pair))
total<-dim(syn)[1]
total # 266 total
head(syn)
# How many significant interactions are positive and negative?
x<-sum(syn$LS>0);print("Positive");x #Positive
print("Positive %"); 100*(x/total)
y<-sum(syn$LS<0);print("Negative");y #Negative
print("Negative %");100*(y/total) #% neg
# How many sig interactions are time-delayed?
t<-sum(syn$Delay!=0);print("Total time-delayed");t
t_pos<-sum(syn$Delay!=0 & syn$LS>0);print("Time delayed and positive");t_pos
t_neg<-sum(syn$Delay!=0 & syn$LS<0);print("Time delayed and negative");t_neg
#
# Ecologically with parasitic syndinales:
# What about negative OR positive time-delayed interactions?
# negative = y and time delayed and positive= t_pos
z<-t_pos + y; print("Percentage Neg or Pos w/ time delayed");z/total;z
#
# Isolate Syndiniales OTUs with positive and time delayed or negative:
syn_sub<-subset(syn, syn$LS<0 | syn$Delay!=0 & syn$LS>0)
write.csv(syn_sub, file="Syndinales_putativeParasitic_06042018.csv")

# Repeat for Rhizaria
rhiz<-LSA_pairs_output %>% filter(grepl("Rhizar", LSA_pairs_output$Pair))
total<-dim(rhiz)[1];total

# How many significant interactions are positive and negative?
x<-sum(rhiz$LS>0);print("Positive");x;print("Positive %");100*(x/total) # % pos
y<-sum(rhiz$LS<0);print("Neg");y;print("Neg %");100*(y/total) #% neg
# How many sig interactions are time-delayed?
t<-sum(rhiz$Delay!=0);print("time delayed");t
t_pos<-sum(rhiz$Delay!=0 & rhiz$LS>0);print("Positive time delayed");t_pos
t_neg<-sum(rhiz$Delay!=0 & rhiz$LS<0);print("Negative time delayed");t_neg
#
# Ecologically with rhizaria:
# Mutualism or commensalism:
print("Mutualism or Commensalism"); x #Positive
print("Mutualism or Commensalism %"); 100*(x/total) # % pos
#
# OR predator-prey (neg or positive time-delayed)
rhiz_pp<-subset(rhiz, rhiz$LS<0 | rhiz$Delay!=0 & rhiz$LS>0)
print("Predator-prey"); dim(rhiz_pp[1])
print("Predator-prey %"); dim(rhiz_pp[1])/total
#
rhiz_sub<-subset(rhiz, rhiz$LS>0)
write.csv(rhiz_sub, file="Rhizaria_putativeMutualism_06042018.csv")

# Deconstruct composition of Syndiniales interactions:
require(tidyverse)
require(purrrlyr)
#dim(syn_sub)
#names(syn_sub)
head(syn[1:2,]) # from here
syn2<-syn[c(34, 32, 29, 3, 9)]
syn2$Count<-"1"
#head(syn2)
syn2$Pair2<-paste(syn2$Taxa2_A, syn2$Taxa2_B, sep="_")
df<-syn2[c(1,6,7,4:5)]
#head(df)
colnames(df)[1:2]<-c("name", "value")
#
df$value<-as.numeric(df$value)
#Get extra taxonomy information for each pair
require(tidyverse)
require(purrrlyr)
syn_summary<-df %>% 
  separate(col = name, sep = "_", c("A", "B")) %>% 
  by_row(.collate = "rows", 
         ..f = function(this_row) {
           paste0(sort(c(this_row$A, this_row$B)), collapse = "_")
         }) %>% 
  rename(sorted = ".out") %>%
  group_by(sorted) %>%
  summarize(tax_detail=paste(unique(Pair2),collapse=", "),positive_timedelay=(sum(LS>0 & Delay!=0)), negative=sum(LS<0), total=sum(LS!=0))%>%
  as.data.frame
  # summarize(tax_detail=paste(unique(Pair2),collapse=", "),positive=sum(LS>0), negative=sum(LS<0), delay_neg=sum(Delay<0), delay_pos=sum(Delay>0), delay_zero=sum(Delay==0))%>%
  # as.data.frame
#head(syn_summary) #frequency of all pairs
syn_summary
write.csv(syn_summary, file="syndinales_freq_Table3_06042018.csv") # Table 3

# # Get taxonomic information 
# syn_tax2<- df %>%
#   group_by(name) %>% #new
#   arrange(Pair2) %>% #new
#   summarize(tax_detail=paste(unique(Pair2),collapse=", "),positive_timedelay=(sum(LS>0 & Delay!=0)), negative=sum(LS<0), total=sum(LS!=0))%>%
#   as.data.frame
# head(syn_tax2[1:3,])
# write.csv(syn_tax2, file = "syndinales_freq_taxID_05062018.csv") # use to populate Table 3 with taxonomic information

# Deconstruct composition of Rhizarian interactions:
require(tidyverse)
require(purrrlyr)
#dim(rhiz_sub)
#names(rhiz_sub)
#head(rhiz_sub[1:2,]) # from here
rhiz2<-rhiz[c(34, 32, 29, 3, 9)]
rhiz2$Count<-"1"

rhiz2$Pair2<-paste(rhiz2$Taxa2_A, rhiz2$Taxa2_B, sep="_")
df<-rhiz2[c(1,6,7,4:5)]
#head(df)
colnames(df)[1:2]<-c("name", "value")
#
df$value<-as.numeric(df$value)
#Get extra taxonomy information for each pair
require(tidyverse)
require(purrrlyr)
rhiza_summary<-df %>% 
  separate(col = name, sep = "_", c("A", "B")) %>% 
  by_row(.collate = "rows", 
         ..f = function(this_row) {
           paste0(sort(c(this_row$A, this_row$B)), collapse = "_")
         }) %>% 
  rename(sorted = ".out") %>%
  group_by(sorted) %>%
  summarize(tax_detail=paste(unique(Pair2),collapse=", "),positive=sum(LS>0), total=sum(value))%>%
  as.data.frame
head(rhiza_summary) #frequency of all pairs
write.csv(rhiza_summary, file="Rhizaria_freq_Table3_06042018.csv") # Table 3

# last updated 06-04-2018 with CLR normalization


# Figure 1B - Environmental data ------------------------------------------------
env<-read.delim("metadata_diel.txt",header = TRUE, sep="\t")
samplenum<-as.character(c(1:19))
timestamp<-c("6PM","10PM","2AM","6AM","10AM","2PM","6PM","10PM","2AM","6AM","10AM","2PM","6PM","10PM","2AM","6AM","10AM","2PM","6PM")
timestamp_actual<-c("18:57","22:36","2:51","6:49","12:52","14:54","18:46","22:42","2:57","6:56","10:52","14:46","18:57b","22:52","2:47","6:47","10:49","14:51","18:53")
env$Time<-factor(env$Sample_Name, levels=samplenum, labels = timestamp_actual)
m.env<-melt(env[1:5])
head(m.env)
cols<-c("#e41a1c", "#4daf4a", "#ff7f00","#984ea3")
line<-ggplot(m.env, aes(x=Time, y=value,group=variable))+geom_path(stat="identity",linetype="blank",aes(x=Time, y=value))+theme(axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=45,hjust = 1,vjust = 1, color="black"),axis.text.y = element_text(color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_rect(data=NULL,aes(xmin=-Inf,xmax=4,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_rect(data=NULL,aes(xmin=7,xmax=10,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_rect(data=NULL,aes(xmin=13,xmax=16,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_rect(data=NULL,aes(xmin=19,xmax=Inf,ymin=-Inf,ymax=Inf),color=NA,fill="#737373",alpha=0.01)+geom_path(stat="identity",linetype="blank",aes(x=Time, y=value))+scale_x_discrete(limits=c(), expand=c(0,0), labels=timestamp)
#
library(gridExtra)
#
grid.arrange(line %+% subset(m.env, variable %in% "Temp")+ scale_y_continuous(limits=c(25.5, 27.5))+geom_path(linetype = "solid",size=1)+geom_point(size=3, shape=21,fill="black", color="white")+labs(y = expression(paste("Temperature ",degree,"C"))), 
             line %+% subset(m.env, variable %in% "Sal")+ scale_y_continuous(limits=c(35,36))+geom_path(linetype = "dashed",size=1)+labs(y="Salinity"), 
             line %+% subset(m.env, variable %in% "DO")+ scale_y_continuous(limits=c(206,213))+geom_line(linetype = "solid",size=1)+geom_point(size=3, shape=23,fill="black", color="white")+labs(y="Diss. oxygen"), 
             line %+% subset(m.env, variable %in% "Chl")+ scale_y_continuous(limits=c(0.15,0.31))+geom_line(linetype = "dotdash",size=1)+labs(y="Chlorophyll a"), ncol=1)
#Saved svg: W:615, H:800
#
summary(env)

# MDS - bray --------------------------------------------------------------
load("Checkpoint2_dfs_wRatio_dielMS.RData", verbose=T) #Option to pick up from here
head(wdf2.vars[1:2,])

# Pre-process dataset:
#Remove unwanted taxa
unique(wdf2.vars$Taxa)
rm<-c("Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
tmp<-subset(wdf2.vars, !(Taxa %in% rm))

#Normalize by relative abundance
library(dplyr)
relAbun<-tmp %>% group_by(Num) %>% mutate(DNA_relabun=DNA/sum(DNA))
relAbun<-relAbun %>% group_by(Num) %>% mutate(RNA_relabun=RNA/sum(RNA))
relAbun<-as.data.frame(relAbun)
names(relAbun)

# Revert back to basic OTU and sample information, separate for DNA and RNA
names(relAbun)
samplenum<-as.character(c(1:19))
new_samplename<-c("1_6PM","2_10PM","3_2AM","4_6AM","5_10AM","6_2PM","7_6PM","8_10PM","9_2AM","10_6AM","11_10AM","12_2PM","13_6PM","14_10PM","15_2AM","16_6AM","17_10AM","18_2PM","19_6PM")
relAbun$Sample<-factor(relAbun$Num, levels = samplenum, labels = new_samplename)
names(relAbun)
dna_base<-dcast(relAbun[c(1,19,17)], Sample~OTU.ID, fill=0)
rna_base<-dcast(relAbun[c(1,19,18)], Sample~OTU.ID, fill=0)

# Format for PCA analysis input
row.names(dna_base)<-dna_base$Sample; dna_base$Sample<-NULL
row.names(rna_base)<-rna_base$Sample; rna_base$Sample<-NULL

# Remove all zero rows
head(dna_base[1:4])
rowsum<-apply(dna_base,1,sum);rowsum # all equal to zero
rowsum<-apply(rna_base,1,sum);rowsum # all equal to zero

# NMDS - calculate for RNA and DNA
## All based on relative abundance
NMDS_rna=metaMDS(rna_base,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
head(NMDS_rna$points) #given points for MDS plot
NMDS_rna$stress
RNA_pts <- NMDS_rna$points[1:nrow(NMDS_rna$points),]
RNA_pts<-as.data.frame(RNA_pts)
plot(RNA_pts$MDS1, RNA_pts$MDS2)
RNA_pts$Sample<-row.names(RNA_pts)


#DNA
NMDS_dna=metaMDS(dna_base,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
head(NMDS_dna$points) #given points for MDS plot
NMDS_dna$stress
DNA_pts <- NMDS_dna$points[1:nrow(NMDS_dna$points),]
DNA_pts<-as.data.frame(DNA_pts)
plot(DNA_pts$MDS1, DNA_pts$MDS2)
DNA_pts$Sample<-row.names(DNA_pts)

#Plot:
key<-read.delim("key_CCA.txt")# Manual curation of colors for sample names.
head(key)
key$Sample<-paste(key$Sample.Number, key$Time, sep="_")
#
label<-as.character(key$Sample)
col<-as.character(key$Color)
shap<-as.numeric(key$Shape)
shap
# names(col)<-Sample
colScale<-scale_color_manual(values=col)
# factor:
DNA_pts$order<-factor(DNA_pts$Sample, levels = label)
RNA_pts$order<-factor(RNA_pts$Sample, levels = label)

#DNA
DNA_MDS <- ggplot(DNA_pts, aes(x = MDS1, y = MDS2, label=order, shape=order), color="black") + geom_point(size=6,aes(fill=order)) + theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +scale_shape_manual(values = shap)+scale_fill_manual(values=col)+ geom_hline(yintercept=0, linetype="dashed", color = "#252525")+ geom_vline(xintercept=0, linetype="dashed", color = "#252525")
DNA_MDS+theme(legend.position = "none")
# 
#RNA
RNA_MDS <- ggplot(RNA_pts, aes(x = MDS1, y = MDS2, label=order, shape=order), color="black") + geom_point(size=6,aes(fill=order)) + theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +scale_shape_manual(values = shap)+scale_fill_manual(values=col)+ geom_hline(yintercept=0, linetype="dashed", color = "#252525")+ geom_vline(xintercept=0, linetype="dashed", color = "#252525")
RNA_MDS+theme(legend.position = "none")
# Make Figure S4
library(cowplot)
figS4<-plot_grid(RNA_MDS+theme(legend.position = "none"), DNA_MDS+theme(legend.position = "none"), align="h", ncol=2, nrow=1,labels=c("A", "B"))
figS4 #svg save W:1000, H: 478
