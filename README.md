# **Genome-wide analysis of Structural Variants in Parkinsons Disease**

**Written By:** Kimberley Billingsley, Cornelis Blauwendraat, Mike Nalls

**Last Updated**: June 2022 

**Quick Description**: Detect and genotype SVs in the LNG PD genomes and run a GWAS with the GATK-SV bi-allelic autosomal 1% FDR  SV calls.  

**Link to Manuscript** 

https://docs.google.com/document/d/1g9Kx_FMqsf_kEZg_EhwJ-wuEjxx0I0hTnV9WztJpoJ4/edit?usp=sharing 

### **Background:**

Parkinson’s disease (PD) is a complex neurodegenerative disorder, affecting approximately one million individuals in the USA alone. A significant proportion of risk for PD is driven by genetics. Despite this, the majority of the common genetic variation that contributes to PD is unknown, in-part because  previous genetic studies have focussed solely on the contribution of single nucleotide variants (SNVs). Structural variants (SVs) represent a significant source of genetic variation in the human genome. Yet, SVs have not been cataloged on a genome-wide scale, and their contribution to the risk of PD remains unknown. For this analysis study we 1) leveraged the GATK-SV pipeline to detect and genotype SVs in 7772 short-read sequencing (SRS) data and 2) generated a subset of matched whole-genome Oxford Nanopore (ONT) long-read sequencing (LRS) data from the PPMI cohort to allow for comprehensive SV confirmation. We detected, genotyped, and tested 3154 “high-confidence” SV loci, representing over 412 million nucleotides of non-reference genetic variation. Utilizing the LRS data we validated three SVs that may drive the association signals at known PD risk loci, including a 2kb deletion within the gene LRRN4. Further we confirm that the majority of the SVs in the human genome cannot be not detected using SRS alone, encompassing on average around 4 million nucleotides of inaccessible sequence per genome. Therefore, although these data provide the most comprehensive survey of the contribution of SVs to the genetic risk of PD to date, this study highlights the need for  large-scale long-read datasets to fully elucidate the role of SV in PD. 

### **Structure of README:**

**1.Short-Read Sequencing (SRS) GATK-SV SV GWAS** 

 1.1. Variant QC = Extract "PASS" SV 
 
 1.2. samples QC= Remove samples that are non-European, do not carry known PD mutations and that have missing pheno/age/sex
 
 1.3. Make SV specific PCs
 
 1.4. Run stepwise logistic regression to find associated PCs with SNV and SV
 
 1.5. Run  GATK-SV risk PD GWAS
 
 1.6. Plot Manhatten (on the hudson) Vs Summary stats of the Nalls 2019 SNV Meta-analysis 5 
 
 1.7. Plot QQ-plot
 

**2. Linkage Disequilibruim Analysis SNV & SV**

2.1. Calc R2 for PR risk SNPs

2.2. Pull out LD SV of interest 

**3. Long-Read Sequencing (LRS) SV analysis with Sniffles2**

3.1. Run Sniffle2 on hg38 mapped LRS bam

3.2. Run variant QC with SURVIVOR

**4. SRS SV confirmation with LRS and Truvari**

4.1. Create SRS GATK-SV vcf per sample

4.2. Run Truvari comparing SRS and LRS 

### **1.Short-Read Sequencing (SRS) GATK-SV SV GWAS** 

#### 1.1. Variant QC = Extract "PASS" SV 


```python
module load bcftools
```


```python
bcftools view -f PASS gatksv_1percent_samples_for_gwas.vcf > pass_gatksv_1percent_samples_for_gwas.vcf
```

#### 1.2. samples QC= Remove samples that are non-European, do not carry known PD mutations and that have missing pheno/age/sex


```python
bcftools view -S sample_IDs_gwas.txt --force-samples pass_gatksv_1percent_samples_for_gwas.vcf -> clean_gatksv_1percent_samples_for_gwas.vcf
```


```python
echo "How many samples are left for GWAS"
bcftools query -l  clean_gatksv_1percent_samples_for_gwas.vcf | wc -l
```

#### 1.3. Make SV specific PCs


```python
module load flashpca
plink --vcf clean_gatksv_1percent_samples_for_gwas.vcf --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --make-bed --out FILENAME_2 
plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
flashpca --bfile FILENAME_3 --suffix _SV_filter_pruned_forPCA.txt --numthreads 19
```


```python
## in R
SV.PC <-read.table(file="pcs_SV_filter_pruned_forPCA.txt", header =T)
SNP.PC <-read.table(file="pcs_filter_pruned_forPCA.txt", header =T)
old.pheno <- read.table(file="pheno_for_tr.association_tr_snp_pcs.txt", header =T)
colnames(SV.PC) <- c("FID", "IID", "SV_PC1", "SV_PC2", "SV_PC3", "SV_PC4", "SV_PC5", "SV_PC6", "SV_PC7", "SV_PC8", "SV_PC9", "SV_PC10")
SVandSNPpc <- merge(SV.PC, SNP.PC, by="FID")
old_covar_needed <- old.pheno[c("FID","IID","SEX","PHENO", "AGE_ANALYSIS")] 
new.pheno <- merge(SVandSNPpc, old_covar_needed, by="FID")
new_covar <- new.pheno[c("FID","IID","SEX","PHENO", "AGE_ANALYSIS","SV_PC1","SV_PC2", "SV_PC3", "SV_PC4", "SV_PC5", "SV_PC6", "SV_PC7", "SV_PC8", "SV_PC9", "SV_PC10", "SNV_PC1", "SNV_PC2", "SNV_PC3", "SNV_PC4", "SNV_PC5", "SNV_PC6", "SNV_PC7", "SNV_PC8", "SNV_PC9", "SNV_PC10")] 
write.table(new_covar, file="pheno_for_SV.association_SV_snp_pcs.tab",sep="\t", row.names = F, quote=F)
```

#### 1.4. Run stepwise logistic regression to find associated PCs with SNV and SV


```python
## in R 

##Run Step wise regression to indentify what covariates are associated with PD(pheno) so we know what to include in association analysis 

library(MASS)
PCs <- read.table(file="pheno_for_SV.association_SV_snp_pcs.tab", header =T)
full.model <- lm(PHENO ~ SEX + AGE_ANALYSIS + SNV_PC1 + SNV_PC2 + SNV_PC3 + SNV_PC4 + SNV_PC5 + SNV_PC6 + SNV_PC7 + SNV_PC8 + SNV_PC9 + SNV_PC10 + SV_PC1 + SV_PC2 + SV_PC3  + SV_PC4 + SV_PC5 + SV_PC6 + SV_PC7 + SV_PC8 + SV_PC9 + SV_PC10, data = PCs)
summary(full.model)
```


```python
# Stepwise regression model
step.model <- stepAIC(full.model, direction = "both", trace = FALSE)
summary(step.model)
```


```python
## make a reduced covariate file based on all the siginicant PCs shown in the list 
```


```python
for_reduced_covar <- PCs[c("FID","IID","SEX","AGE_ANALYSIS","SNV_PC1","SNV_PC2","SNV_PC4","SNV_PC5","SNV_PC6","SNV_PC7","SNV_PC8","SNV_PC9","SNV_PC10","SV_PC1","SV_PC2","SV_PC3","SV_PC4","SV_PC5","SV_PC6","SV_PC8","SV_PC9","SV_PC10")]
```


```python
write.table(for_reduced_covar, file="1percent_all_reduced_covar_sv.tab",sep="\t", row.names = F, quote=F)
```

#### 1.5. Run  GATK-SV risk PD GWAS


```python
module load plink
plink --vcf clean_gatksv_1percent_samples_for_gwas.vcf --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --make-bed --out 1percent_clean_PD_gaksv_pass_bnd_maf_0.01_hwe_5e6_geno0.01
module load plink/2.0_alpha_1_final
plink2 --covar 1percent_all_reduced_covar_sv.tab  --covar-variance-standardize --glm hide-covar --out gatksv_gwas_clean_1percent_calls  --bfile 1percent_clean_PD_gaksv_pass_bnd_maf_0.01_hwe_5e6_geno0.01  --pheno risk_pheno.tab
```


```python
plink2 --adjust-file gatksv_gwas_clean_1percent_calls.PHENO.glm.logistic test=ADD --out adjusted. gatksv_gwas_clean_1percent_calls.PHENO.glm.logistic
```

#### 1.6. Plot Manhatten (on the hudson) Vs Summary stats of the Nalls 2019 SNV Meta-analysis 5 


```python
awk '{print $3, $1, $2, $12}'gatksv_gwas_clean_1percent_calls.PHENO.glm.logistic > for_plotting_man_gatk_gwas.txt
```


```python
#manually change header to SNP	CHR	POS	pvalue
```


```python
## Do this in R/4.0
devtools::install_github('anastasia-lucas/hudson')
```


```python
library(hudson)
```


```python
meta5_plot <- read.table(file="/data/CARD/PD/GAK_SV/GWAS/plot_meta_5_hg38.tab", header =T)
gatksv_plot <- read.table(file="/data/CARD/PD/GATK-SV_22/GWAS/for_plotting_man_gatk_gwas.txt", header =T)
gangstr_plot <- read.table(file="/data/CARD/PD/GAK_SV/GWAS/plot_gangstr_hg38.tab", header =T)
```


```python
gmirror(top=gatksv_plot, bottom=meta5_plot, tline=5.34E-6, bline=5E-8,background = "white", file = "new_gatksv_vs_meta_5_no23_HudsonPlot",highlight_p = c(0.05/nrow(gatksv_plot),0.05/nrow(meta5_plot)),highlighter="green3",toptitle="GAK-SV genome-wide association analysis", bottomtitle = " Meta 5 SNP genome-wide association analysis", hgt = 14, hgtratio = 0.5, wi = 24, res = 300)
```

#### 1.7.Plot QQ-plot


```python
vignette('qqman')
```


```python
gaksv_plot <- read.table(file="for_plotting_man_gatk_gwas.txt", header =T)
```


```python
qq(gaksv_plot$pvalue, main = "Q-Q plot of GATK_SV GWAS p-values", xlim = c(0, 10), ylim = c(0, 10), pch = 18, col = "blue4", cex = 1.5, las = 1)
```


```python
## Lambda value
```


```python
data = read.table("for_plotting_man_gatk_gwas.txt",header=T)
p <- data$pvalue
n <- length(data$pvalue)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda
```

### **2. Linkage Disequilibruim Analysis SNV & SV**

#### 2.1. Calc R2 for PR risk SNPs


```python
module load plink/1.9
for f in `cat PD_risk.SNVs.txt`
    do
         plink --bfile  1percent_clean_PD_gaksv_pass_bnd_maf_0.01_hwe_5e6_geno0.01 --ld-snp ${f} --out ${f}.risk.ld.pairs --r2 inter-chr dprime --silent 
    done
```

#### 2.2. Pull out SV of interest 


```python
## Lets see how many SV are in LD
grep "PD_" *risk.ld.pairs.ld*
```

grep "PD_" *risk.ld.pairs.ld* > GATK_SV_in_ld_PD_risk_snvs_inter_chr_dprime.txt


```python
 awk -v OFS='\t' '{print $7}' GATK_SV_in_ld_PD_risk_snvs_inter_chr_dprime.txt > GATK_SV_ID_in_ld_PD.txt
```


```python

```

### **3.Long-Read Sequencing (LRS) SV analysis with Sniffles2**

#### 3.1. Run Sniffle2 on hg38 mapped LRS bam


```python
module load sniffles/2.0.3
```


```python
sniffles -i PPMI_3404_BLOOD_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI3404_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_4018_BLOOD_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI4018_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_3223_BLOOD_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI3223_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_3173-BLOOD_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI3173_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_3150_BLOOD_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI3150_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_4011_BLOOD_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed  -v PPMISI4011_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_PPMI3951_Blood_test_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI3951_sniffles2_tandemrepeats.vcf
sniffles -i PPMI_PPMI3469_Blood_test_merged_sorted.bam --tandem-repeats  human_GRCh38_no_alt_analysis_set.bed -v PPMISI3469_sniffles2_tandemrepeats.vcf
```

#### 3.2. Run variant QC with SURVIVOR


```python
### The number of SV is still high per genome because we need to do more QC to remove short reads 
```


```python
module load survivor 
```


```python
# Remove SV > 50bp with SURVIVOR 
for f in `cat ppmi_samples_lrs_blood.txt`
    do
    SURVIVOR filter ${f}_sniffles2_tandemrepeats.vcf NA 50 -1 -1 -1 ${f}_no_50bpSV.vcf
    done
```


```python
for f in `cat ppmi_samples_lrs_blood.txt`
    do
    grep -vc '#' ${f}_no_50bpSV.vcf
    done 
```


```python
 for f in `cat ppmi_samples_lrs_blood.txt`
    do
    SURVIVOR stats ${f}_no_50bpSV.vcf -1 -1 -1 ${f}_no_50bpSV.sum.stats
    done
```


```python
for f in `cat ppmi_samples_lrs_blood.txt`
    do
    bgzip -c  ${f}_no_50bpSV.vcf  > ${f}_no_50bpSV.vcf.gz
    tabix -p vcf ${f}_no_50bpSV.vcf.gz
done
```

#### **4. SRS SV confirmation with LRS and Truvari**

#### 4.1. Create SRS GATK-SV vcf per sample


```python
module load bcftools
```


```python
for f in `cat ppmi_samples_lrs_blood.txt`
    do
    bcftools view -s ${f} clean_gatksv_1percent_samples_for_gwas.vcf > ${f}_pass_1percent_sv.vcf
    done
```


```python
 for f in `cat ppmi_samples_lrs_blood.txt`
    do       
        bcftools view ${f}_pass_1percent_sv.vcf --min-ac=1 > carriers_${f}_pass_1percent_sv.vcf
        bgzip -c carriers_${f}.pass_1percent_sv.vcf > carriers_${f}.pass_1percent_sv.vcf.gz
        tabix -p vcf carriers_${f}.pass_1percent_sv.vcf.gz  
    done
```

#### 4.2. Run Truvari comparing SRS and LRS 


```python
## GATK-SV "PASS" SV
```


```python
source  /data/CARD/PD/GAK_SV/GWAS/GATK_ONT_VAPOR/vapor/conda/etc/profile.d/conda.sh
```


```python
conda activate truvari
```


```python
for f in `cat ppmi_samples_lrs_blood.txt`
    do
        truvari bench -b ${f}_no_50bpSV.vcf.gz  -c carriers_${f}.pass_1percent_sv.vcf.gz  --pctsim=0 -f /fdb/GENCODE/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa -o ${f}_truvari_repeat_filteredLRSSV
    done
```
