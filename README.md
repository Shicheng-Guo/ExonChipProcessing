ExonChipProcessing
==================

# Introduction #

This site contains the codes and resources for exome chip processing protocol.

The list of codes:

 Name        |  Step  | Called By  | Notes  
 ------------- | -----:|------:|-------:
 MergeSampleSheet.pl       | 1B |User|Merging sample sheets
runZcall.py      | 34A |User|Run zCall
 Gender.R      | 39 |User|Checking for sex mismatch
 PCAPlot.R |     43 |User|Draw scatter plot of principle Components
 PlotHWE.R | 48 |User|Plot histograms of HWE test
 PlotHeterozygosity.R | 50 |User|Compute heterozygosity and plot histograms of heterozygosity and inbreeding coefficient
ConsistencyDupSNP.sh	|51	|User	|Prepare data for checking consistency of duplicated SNPs
ConsistencyDupSNP.pl	|51	|ConsistencyDupSNP.sh	|Checking genotyping consistency of duplciated SNPs, called by ConsistencyDupSNP.sh
Consistency1000G.sh|		52|	User	|Prepare data for checking consistency with 1000G
Consistency1000GSNP.pl|52	|Consistency1000G.sh	|Checking genotyping consistency with 1000G, called by Consistency1000G.sh 
exclude.pl	|52	|Consistency1000G.sh	|Exclude bad SNPs
AlleleFreq1000G.sh	|53	|User	|Compute allele frequency of 1000G
vcf_to_ped.py	|53	|AlleleFreq1000G.sh	|Convert VCF to ped
AlleleFreqExome.sh	|55	|User	|Compute allale frequency of exome chip
MAFtoAF.py	|55	|AlleleFreqExome.sh	|Change MAF to allele frequency
1000GAlleleFreqPlot.R	|56	|User	|Plot allele frequency scatter plot between 1000G and exome chip
BatchAlleleFreqMatrix.R	|57	|User	|Plot correlation matrix between batches
filter.pl	|52, 55	|AlleleFreqExome.sh, Consistency1000G.sh	|Filter out non-overlapping SNPs


For the resources files, you need to download them from the following links and then copy them to the resources folder under exome chip processing protocol codes folder. And you need to unzip 1000G_ExomeChipOverlapVCF.zip to get G1000.vcf.

The list of resources for 12V1_A exome chip:

 Name        | Used by Command           | Called by   | Notes 
 ------------- |:-----------:|:-----------:| -----:
[PAR_SNPs.txt](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/PAR_SNPs.txt)|	13|User in GenomeStudio|This is a list of all PAR SNPs on the exome chip, can be used for filtering them out in GenomeStudio
[Aims.txt](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/AIMs.txt)|	40|User|List of all AIMs markers on exome chip
[g1k_HumanExome-12v1_A_SNPs](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/g1k_HumanExome-12v1_A_SNPs)|	52|Consistency1000G.sh|	1000G Overlapped SNP list
[g1k_HumanExome-12v1_A_SNPs.bed](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/g1k_HumanExome-12v1_A_SNPs.bed)|	52|Consistency1000G.sh|	1000G Overlapped SNP list
[g1k_HumanExome-12v1_A_SNPs.bim](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/g1k_HumanExome-12v1_A_SNPs.bim)|	52|Consistency1000G.sh|	1000G Overlapped SNP list
[g1k_HumanExome-12v1_A_SNPs.fam](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/g1k_HumanExome-12v1_A_SNPs.fam)|	52|Consistency1000G.sh|	1000G Overlapped SNP list
[dup_snp_pair](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/dup_snp_pair)	|51|ConsistencyDupSNP.sh|	Duplicated SNP list
[1000G_ExomeChipOverlapVCF.zip (G1000.vcf)](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/1000G_ExomeChipOverlapVCF.zip)	|53, 55|AlleleFreq1000G.sh, AlleleFreqExome.sh, vcf_to_ped.py|	VCF file of 1000G data which only contains SNP overlapped with exome chip
[chr23_26.txt](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/chr23_26.txt)	|44	|plink|list of SNPs from Chr X, Y and MT
[integrated_call_samples.20101123.ped](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/integrated_call_samples.20101123.ped)	|53|vcf_to_ped.py|	Downloaded from 1000G
[integrated_call_samples.20101123.ALL.panel](https://github.com/slzhao/ExonChipProcessing/releases/download/resources.12V1_A/integrated_call_samples.20101123.ALL.panel)	|52	|Consistency1000G.sh|1000 Genome sample information downloaded from 1000G

# Procedure #

the following procedure are based on plink 1.9, plink 2.0, R 3.62 and Python 3.0

1. Converting all SNPs to the forward strand
```
input="RA2020"
plink --file $input --make-bed --out $input
```
2. Checking for gender mismatch
```
plink --bfile $input --maf 0.1 --check-sex --out $input.1
wget https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/Gender.R
Rscript Gender.R $input.1.sexcheck sex.checking.jpeg
```
3. Checking for race mismatch
```
exmAims<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/exmAIMs.txt")
db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/Illumina_CoreExome_Beadchip.hg19.exm2rs.bed.txt")
rsAims<-db[na.omit(match(exmAims[,1],db$V3)),]
write.table(rsAims,file="rsAims.bed",quote=F,row.names = F,col.names = F,sep="\t")
write.table(rsAims[,6],file="rsAims.txt",quote=F,row.names = F,col.names = F,sep="\t")

plink --bfile $input.1 --extract rsAims.txt --recode --out AIMs

perl -p -i -e '{/-9/1/g}' AIMs.map

cd ~/hpd/tools/
git clone https://github.com/chrchang/eigensoft.git

vim par.PED.EIGENSTRAT

genotypename:    AIMs.ped
snpname:         AIMs.map
indivname:       AIMs.ped
outputformat:    EIGENSTRAT
genotypeoutname: AIMs.geno
snpoutname:      AIMs.snp
indivoutname:    AIMs.ind
familynames:     NO

smartpca.perl -i AIMs.geno -a AIMs.snp -b AIMs.ind -o AIMs.pca -e AIMs.eval -l AIMs.log -p AIMs -m 0

smartpca -p AIMs.pca.par >AIMs.log
ploteig -i AIMs.pca.evec -c 1:2  -p Control  -x  -y  -o AIMs.xtxt
evec2pca.perl 10 AIMs.pca.evec AIMs.ind AIMs.pca

wget https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/PCAPlot.R
Rscript PCAPlot.R AIMs.pca

```
4. Checking for relatedness
```
plink --bfile {exome} --maf 0.1 --exclude chr23_26.txt --indep-pairwise 50 5 0.2 --out {indepSNP}
plink --bfile {exome} --extract {indepSNP.prune.in} --genome --out {outfile}
```
5. Checking for Hardy-Weinberg equilibrium (HWE) outliers
```
plink --bfile {exome} --keep {Caucasian_list} --make-bed --out {Caucasian}
plink --bfile {Caucasian} --maf 0.05 --hardy --out {outfile}
Rscript PlotHWE.R {outfile}.hwe
```
6. Checking for heterozygosity outliers
```
plink --bfile {exome} --extract {indepSNP}.prune.in --het --out {outfile}
Rscript PlotHeterozygosity.R {outfile}.het
```
7. Checking consistency between exome chip genotype and 1000 Genomes Project17 or HapMap18 genotype
```
file= "{exome}" output="{outfile}.result" sh ConsistencyDupSNP.sh
file="{exome}" output="{outfile}" sh Consistency1000G.sh
```
8. Checking for minor allele frequency (MAF) consistency between exome chip and 1000 Genomes Project genotypes
```
G1kRace="{race},.., {race}" sh AlleleFreq1000G.sh
exome_dir="{exomePath}" exomeRace="{racefile}" sh AlleleFreqExome.sh
Rscript 1000GAlleleFreqPlot.R final_{race}_g1k_exm {out}.jpeg
```
9. Checking for batch effects
```
 Rscript BatchAlleleFreqMatrix.R final_exm_${race} {out}.jpeg
```
