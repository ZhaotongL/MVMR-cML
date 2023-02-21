
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MVMRcML

<!-- badges: start -->
<!-- badges: end -->

The goal of MVMRcML is to conduct Multivariable Mendelian Randomization
(MVMR) analysis.

## Installation

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ZhaotongL/MVMR-cML")
```

## Example

Here is an example which shows how to apply MVMR-cML to infer the causal
relationship from cholesterol levels (triglycerides, LDL and HDL) to
CAD.

First extract GWAS summary data with package:

``` r
library(TwoSampleMR)
#> TwoSampleMR version 0.5.6 
#> [>] New: Option to use non-European LD reference panels for clumping etc
#> [>] Some studies temporarily quarantined to verify effect allele
#> [>] See news(package='TwoSampleMR') and https://gwas.mrcieu.ac.uk for further details
gwas_id = c(
      'ebi-a-GCST002216',
      'ebi-a-GCST002222',
      'ebi-a-GCST002223')
exposure_dat = mv_extract_exposures(gwas_id)
#> API: public: http://gwas-api.mrcieu.ac.uk/
#> Please look at vignettes for options on running this locally if you need to run many instances of this command.
#> Clumping 1, 201 variants, using EUR population reference
#> Removing 68 of 201 variants due to LD with other variants or absence from LD reference panel
#> Extracting data for 133 SNP(s) from 3 GWAS(s)
#> Finding proxies for 1 SNPs in outcome ebi-a-GCST002216
#> Extracting data for 1 SNP(s) from 1 GWAS(s)
#> Harmonising Triglycerides || id:ebi-a-GCST002216 (ebi-a-GCST002216) and LDL cholesterol || id:ebi-a-GCST002222 (ebi-a-GCST002222)
#> Removing the following SNPs for being palindromic with intermediate allele frequencies:
#> rs3758348, rs7534572, rs9491696
#> Harmonising Triglycerides || id:ebi-a-GCST002216 (ebi-a-GCST002216) and HDL cholesterol || id:ebi-a-GCST002223 (ebi-a-GCST002223)
#> Removing the following SNPs for being palindromic with intermediate allele frequencies:
#> rs3758348, rs7534572, rs9491696
outcome_dat <- extract_outcome_data(exposure_dat$SNP,'ebi-a-GCST005195')
#> Extracting data for 132 SNP(s) from 1 GWAS(s)
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
#> Harmonising Triglycerides || id:ebi-a-GCST002216 (ebi-a-GCST002216) and Coronary artery disease || id:ebi-a-GCST005195 (ebi-a-GCST005195)
#> Removing the following SNPs for being palindromic with intermediate allele frequencies:
#> rs3758348, rs7534572, rs9491696
```

MVMR-cML analysis:

``` r
library(MVMRcML)
## rho_mat is pre-calculated by bivariate LDSC
rho_mat = matrix(c(1,0.1988,-0.3347,0,0.1988,1,-0.1114,0,-0.3347,-0.1114,1,0,0,0,0,1),ncol=4)
Sig_inv_l = invcov_mvmr(se_bx=mvdat$exposure_se,se_by=mvdat$outcome_se,rho_mat = rho_mat)

MVcML_res = MVmr_cML_DP(b_exp=mvdat$exposure_beta,
                        b_out=as.matrix(mvdat$outcome_beta),
                        se_bx=mvdat$exposure_se,
                        Sig_inv_l=Sig_inv_l,n = 188577,num_pert = 100,
                        K_vec = 0:20 # try a small range of K first to save time. Change this accordingly.
                        )
```

In above code, we first tried a small range of *K* from 0 to 20 (which
can be changed accordingly in practice), and MVMR-cML-BIC selected 8
invalid IVs.

``` r
## Invalid IVs selected by MVMR-cML-BIC
MVcML_res$BIC_invalid
#>      [,1]
#> [1,]    8
#> [2,]   12
#> [3,]   21
#> [4,]   25
#> [5,]   73
#> [6,]   98
#> [7,]  107
#> [8,]  129

## Inference by MVMR-cML-BIC (NOT ACCOUNTING FOR MODEL SELECTION UNCERTAINTY!)
MVcML_BIC_SE = MVcML_SdTheta(b_exp=mvdat$exposure_beta,
              b_out=as.matrix(mvdat$outcome_beta),
              Sig_inv_l=Sig_inv_l,
              theta=MVcML_res$BIC_theta,
              zero_ind = setdiff(1:length(mvdat$outcome_beta),MVcML_res$BIC_invalid))

MVcMLBIC_pval = pnorm(-abs(MVcML_res$BIC_theta/MVcML_BIC_SE))*2
```

MVMR-cML-DP suggested significant direct effects of TG and LDL on CAD,
and a null effect of HDL on CAD.

``` r
## Inference by MVMR-cML-DP
MVcMLDP_pval = pnorm(-abs(MVcML_res$BIC_DP_theta/MVcML_res$BIC_DP_se))*2
MVcMLDP_pval
#>              [,1]
#> [1,] 5.281724e-07
#> [2,] 5.758572e-17
#> [3,] 1.749562e-01
```
