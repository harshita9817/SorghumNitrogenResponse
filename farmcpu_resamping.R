library(rMVP)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)

MVP.Data(fileVCF = "SAP_IMPUTED_filtered_maf0.05_het0.1.vcf",
         filePhe = "Lowhnitrogen_333_filteredgenes.csv", sep.phe = ",", 
         fileKin = TRUE,
         filePC = TRUE,
         out = "mvp.vcf")

genotype <- attach.big.matrix("mvp.vcf.geno.desc")
phenotype <- read.table("mvp.vcf.phe",head=TRUE)
map <- read.table("mvp.vcf.geno.map" , head = TRUE)
Kinship <- attach.big.matrix("mvp.vcf.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.vcf.pc.desc"))

for(x in 1:100){
  phe1=phenotype # make a copy of phenotype
  nline = nrow(phe1)
  phe1[sample(c(1:nline), as.integer(nline*0.1)),2:ncol(phe1)]=NA  # randomly choose 10% phenotype to be NA
  colnames(phe1)=paste0(colnames(phenotype),x)  # rename the phenotype by attaching bootstrap number 
  for(i in 2:ncol(phe1)){
    imMVP<-MVP(phe = phe1[,c(1,i)], geno = genotype, map = map, K=Kinship, CV.FarmCPU=Covariates_PC, file.output = "pmap.signal",
               nPC.FarmCPU = 3, maxLoop = 10, method = "FarmCPU", ,threshold=0.2,p.threshold=1e-08)
  }
}
# Get all filenames ending with .FarmCPU_signals.csv
# List all files in the 'farmcpu' directory
files <- list.files( pattern = ".FarmCPU_signals.csv$", full.names = TRUE)
# Extract trait names (the part before .FarmCPU_signals.csv)
traits <- c("Chl_slope","Chl_intercept","FloweringTime_slope", "Floweringtime_intercept","GrainWeight_slope","GrainWeight_intercept", "Plantheight_slope", "Plantheight_intercept")
# Check the extracted traits
print(traits)
get.support = function(trait) {
  files = list.files(pattern = paste0(trait, ".*FarmCPU_signals.csv"))
  
  if (length(files) >= 1) {
    signals <- files %>%
      map_df(~read.csv(., skip = 1, header = F, colClasses = c("factor", "factor", "integer", "factor", "factor", "numeric", "numeric", "numeric"), check.names = TRUE)) %>%
      setNames(c("SNP", "CHROM", "POS", "REF", "ALT", "Effect", "SE", "something", "pvalue"))
    
    signals <- signals %>%
      group_by(SNP, CHROM, POS) %>%
      summarise(P = mean(pvalue), support = n()/100)
    
    write.table(signals, file = paste0("Z", trait, "signals.csv"), quote = F, row.names = F, sep = ",")
  } else {
    print(paste0("File not found for trait: ", trait))
  }
}
# Loop through all traits
for (trait in traits) {
  get.support(trait)
}
##########################################################################################################
#PLOTTING
library(tidyverse)
library(tidyverse)

# File list with trait labels
files <- list(
  "ZChl_interceptsignals.csv"= "Chl_intercept",
  #"ZLN_DaysToBloom_2020signals.csv"= "DaysToBloom_2020",
  #"ZLN_PanicleGrain_Weight_2020signals.csv"= "PanicleGrain_weight_2020", 
  "ZChl_slopesignals.csv"= "chl_slope",
  "ZGrainWeight_interceptsignals.csv"= "grainweight_intercept",
  "ZGrainWeight_slopesignals.csv" = "grainweightlope") 
  #"ZLN_Flowering_2022signals.csv" = "flowering_2022", 
  #"ZLN_PlantHeight_2022signals.csv" = "PlantHeight2022", 
  #"ZLN_Chlorophyll_2023signals.csv" = "Chl2023")
 # "ZLN_chl_2024"= "Chl_2024", 
 # "ZLN_plantheight_2024" = "PLantHeight_2024")
# Load and combine all files
df_list <- lapply(names(files), function(f) {
  df <- read.csv(f)
  df$Trait <- files[[f]]
  df
})
df_all <- bind_rows(df_list)
# Ensure CHROM is numeric or properly ordered factor
df_all$CHROM <- as.character(df_all$CHROM)
# Compute cumulative positions
data_cum <- df_all %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS), .groups = "drop") %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0))
gwas_data <- df_all %>%
  left_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add)
# For x-axis labels
axis_set <- gwas_data %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum), .groups = "drop")
library(viridis)
ggplot(gwas_data, aes(x = bp_cum, y = support, color = Trait)) +
  geom_vline(xintercept = data_cum$bp_add, color = "grey80", linetype = "dashed", lwd = 0.2) +
  geom_point(alpha = 0.8, size = 2) +
  #geom_hline(yintercept = 0.015, color = "grey50", linetype = "dashed", lwd = 0.5) +
  geom_hline(yintercept = 0.2, color = "grey50", linetype = "dashed", lwd = 0.5) +
  scale_color_viridis_d(option = "turbo") +
  scale_x_continuous(
    breaks = axis_set$center,
    labels = paste0("", axis_set$CHROM)
  ) +
  labs(x = "Chromosome", y = "Support", color = "Trait") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
