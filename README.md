# PRSice-nf
Pipeline to compute polygenic risk scores

## Input
Type `--help` to get the full list of options. All parameters are prefixed with a double-dash like in `--help`.

| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| `bgen_list` |  `bgen_list.txt`| Txt file containing the list of bgen files to consider (one file per chromosome) |
| `base_file`   |        `base.txt`   | file containing association analysis results for SNPs on the base phenotype  |
| `pheno_file`   | `pheno.txt` | phenotype file: the first two column of the phenotype file should be the FID and the IID and the rest of the columns the phenotypes |

For more details on the files content, see [PRSice2 documentation](https://www.prsice.info/step_by_step/#input-data)

## Parameters
  * #### Mandatory

| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| `bgenix`   | `bgenix` | path to bgenix executable |
| `PRSice_path`   | `path_PRSice` | path to PRSice scripts |

  * #### Optional

| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| `a1_col`   |        `4`   | Column number associated with the effect allele (Default value is 4) |
| `a2_col`   | `5` | Column number associated with the reference allele (Default value is 5) |
| `qual_val`   |        `0.3`   | Quality threshold. Any SNPs with a value bellow this threshold will be excluded. (Default value is 0) |
| `pval_thr`   |        `0.0001,0.1,1`   | pvalue threshold to consider to select SNPs to include in the PRS |

To install:
- [bgen](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk)
- [PRSice2](https://www.prsice.info/)

## Nextflow command line
```
nextflow run  IARCbioinfo/PRSice-nf \
  --base_file base.txt \
  --pheno_file /data/gep/MR_Signatures/work/gabriela/GIT_Rstudio_project/results/PRSice_res/TCGA_samples.txt  --pval_thr 1 \
  --output_name ${trait} --output_folder /data/gcs/prs-computation/work/gabriela/PRS_R_project/results/PRSice_res/${trait}/ \
  --bgenix path_to_bgenix --PRSice_path path_to_PRSice \
  --bgen_list bgen_list.txt \
  --bgen_snps_perCHR /data/gep/MR_Signatures/work/gabriela/GIT_Rstudio_project/results/PRSice_res/TCGA_bgen_SNPs/ --qual_val 0.3 

```

