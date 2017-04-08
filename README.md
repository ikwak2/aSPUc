## R package, aSPU

[![License](https://img.shields.io/badge/license-GPL%203-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html) [![CRAN](http://www.r-pkg.org/badges/version/aSPU)](https://CRAN.R-project.org/package=aSPU)


Il-Youp Kwak <ikwak@umn.edu>

R/aSPU is an R package for Genetic association testing methods such as aSPU, aSPUw, aSPUpath, aSPUs, aSPUsPath, GEEaSPU, MTaSPUs, GATES, GATE-Simes, HYST etc.

### Summary table 

| Function Name         | Data Type  | Description                                    |
|-----------------------|------------|------------------------------------------------|
| aSPU, aSPUw, aSPUr    | Individual | Single trait; gene-based                       |
| aSPUs                 | Summary    | Single trait; gene-based                       |
| aSPUpath              | Individual | Single trait; pathway-based                    |
| aSPUsPath             | Summary    | Single trait; pathway-based                    |
| GEEaSPU               | Individual | Multiple traits; single SNP based              |
| MTaSPUs               | Summary    | Multiple traits; single SNP based              |
| MTaSPUsSet            | Summary    | Multiple traits; gene-based                    |
| MTaSPUsSetPath        | Summary    | Multiple traits; pathway-based                 |


*Data type indicate the structure of data set. "Individual" for individual level data. "Summary" for summary statistics data (such as Z scores or p-values of each SNP) 

### Tutorials 

Some [tutorials](https://github.com/ikwak2/aSPU_tutorials) for aSPUs, aSPUsPath and MTaSPUsSet are available. 
 - [Mapping Snp to Gene](http://www.tc.umn.edu/~ikwak/tutorials/mappingSnpToGene2.html)
 - [Getting correlation estimate among SNPs from reference panel](http://www.tc.umn.edu/~ikwak/tutorials/CorrFromRef.html)
 - [R and Perl codes to perform aSPUs and MTaSPUsSet](http://www.tc.umn.edu/~ikwak/tutorials/ForMTgenes.html)
 - [Pathway data manipulation](http://www.tc.umn.edu/~ikwak/tutorials/pathwayMani.html)
 - [Speed comparison using R, awk and Perl](http://www.tc.umn.edu/~ikwak/tutorials/SpeedComp.html)
 - [Vignette for aSPUs and aSPUsPath](http://www.tc.umn.edu/~ikwak/tutorials/aSPUstat.html)
 - [Vignette for MTaSPUsSet](http://www.tc.umn.edu/~ikwak/tutorials/MTaSPUsSet.html)

### Citations

For 'aSPU'
```
Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
A powerful and adaptive association test for rare variants,
Genetics, 197(4), 1081-95
```

For 'aSPUw'
```
Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014)
Comparison of statistical tests for group differences in brain functional networks,
NeuroImage, 1;101:681-694
```

For 'aSPUr'
```
Peng Wei, Ying Cao, Yiwei Zhang, Zhiyuan Xu, Il-Youp Kwak, Eric Boerwinkle, Wei Pan (In press)
On Robust Association Testing for Quantitative Traits and Rare Variants, G3.
```

For 'aSPUpath'
```
Wei Pan, Il-Youp Kwak and Peng Wei (2015)
A Powerful and Pathway-Based Adaptive Test for Genetic Association With Common or Rare Variants,
The American Journal of Human Genetics 97, 86-98
```

For 'aSPUs' and 'aSPUsPath'
```
Il-Youp Kwak, Wei Pan (2015)
Adaptive Gene- and Pathway-Trait Association Testing with GWAS Summary Statistics,
Bioinformatics, 32(8), 1178-1184
```

For 'GEEaSPU'
```
Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014)
Testing for association with multiple traits in generalized estimation equations, with application to neuroimaging data,
Neuroimage. 96, 309-325
```

For 'MTaSPUs'
```
Junghi Kim, Yun Bai and Wei Pan (2015)
An Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics,
Genetic Epidemiology, 8:651-663
```

For 'MTaSPUsSet' and 'MTaSPUsSetPath'
```
Il-Youp Kwak, Wei Pan (2016)
Gene- and pathway-based association tests for multiple traits with GWAS summary statistics, Bioinformatics. doi:10.1093/bioinformatics/btw577
```


For 'GATES'
```
Miao-Xin Li, Hong-Sheng Gui, Johnny S.H. Kwan and Pak C. Sham (2011)
GATES: A Rapid and Powerful Gene-Based Association Test Using Extended Simes Procedure,
The American Journal of Human Genetics 88, 283-293
```

For 'GATES-Simes'
```
Hongsheng Gui, Miaoxin Li, Pak C Sham and Stacey S Cherny (2011)
Comparisons of seven algorithms for pathway analysis using the WTCCC Crohn's Disease
BMC Research Notes, 4:386
```

For 'HYST'
```
Miao-Xin Li, Johnny S.H. Kwan and Pak C. Sham (2012)
HYST: A Hybrid Set-Based Test for Genome-wide Association Studies, with Application to Protein-Protein Interaction-Based Association Analysis
The American Journal of Human Genetics, 91, 478-488.
```



### installation
From `CRAN` :
```S
install.packages("aSPU")
```

Or, with `devtools`:
```S
library(devtools)
install_github("ikwak2/aSPU")
```
