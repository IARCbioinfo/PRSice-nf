#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "---------------------------------------------------------"
log.info "  Imputed snps extraction <VERSION>: <SHORT DESCRIPTION> "
log.info "---------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/PRSice-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--base_file                      FILE                      File containing all SNPs to analyse"
    log.info "--pheno_file                      FILE                      File containing phenotypes information: two columns, one for the samples IDs and one for the phenotype. The phenotype can be a binary variable (e.g case or control) or a quantitative trait"
    log.info "--input_folder                      PATH                      Path to folder containing the imputation results of the ukbiobank data (bgen files)"
    log.info ""
    log.info "Optional arguments:"
    log.info "--target_file                     FILE                      File containing the imputation data, bgen file or prefix of the plink files"
    log.info "--samples_file                     FILE                     File containing the samples IDs (one sample per row) to keep for the PRS analysis"
    log.info "--PRSice_path                     PATH                      Path to PRSice executables"
    log.info "--stat_col                     STRING                      Column name in the base file for the weights to be used in the PRS computation (default value is BETA)"
    log.info "--output_folder                     PATH                      Path to output folder"

    log.info ""

    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}


params.base_file = null
params.pheno_file = null
params.input_folder = '/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/imputation_results/'
params.samples_file = null
params.target_file = null
params.PRSice_path = '/home/gabriela/softwares/'
params.output_folder = '/data/gep/MR_Signatures/work/gabriela/GIT_Rstudio_project/results/PRSice_test/'
params.output_name="res"

params.snp_col=1
params.chr_col=2
params.pos_col=3
params.a1_col=4
params.a2_col=5
params.stat_col=6
params.se_col=7
params.pval_col=8

// add a if to check if a samples file is provided if yes add the option --keep
// add a if to deal with plink files
// add a parameter (string) allowing the user to add multiple PRSice flag options

process get_chr {

  output:
  file 'chr_list.txt' into chr_list
  file 'base_file.txt' into base_file_updated

  shell:
  '''
  tail -n +2 !{params.base_file} | awk '{print $2}' | sort | uniq  > chr_list.txt

  awk '{print $snp,$chr,$pos,$a1,$a2,$stat,$se,$pval}' snp=!{params.snp_col} chr=!{params.chr_col} pos=!{params.pos_col} a1=!{params.a1_col} a2=!{params.a2_col} stat=!{params.stat_col} se=!{params.se_col} pval=!{params.pval_col} !{params.base_file} > base_file.txt
  '''
}

process compute_PRS_per_chr {

  input:
  val chr from chr_list.splitText(by:1)
  file base from base_file_updated

  output:
  file '*.all.score' into scores
  file '*.snp' into snps_included

  shell:
  '''
  chr_val=!{chr}
  echo $chr_val
  /home/gabriela/softwares/gavinband-bgen-0b7a2803adb5/build/apps/bgenix -g !{params.input_folder}chr_${chr_val}/imputation/concatenated_chr${chr_val}.bgen -list > chr_${chr_val}_SNPs.txt
  awk '{if($2==c)print $0}' c=${chr_val} !{base} > chr_${chr_val}_base.txt
  if [[ $(wc -l <chr_${chr_val}_base.txt) -ge 1 ]]
  then
    grep -Ff <(cat chr_${chr_val}_base.txt | cut -d ' ' -f 2,3  | sed 's/ /:/g') chr_${chr_val}_SNPs.txt > chr_${chr_val}_SNPs_subset.txt

    /home/gabriela/R3.5.1/bin/Rscript !{baseDir}/match_IDs.r chr_${chr_val}_base.txt chr_${chr_val}_SNPs_subset.txt ${chr_val}

    /home/gabriela/R3.5.1/bin/Rscript !{params.PRSice_path}PRSice.R --prsice !{params.PRSice_path}PRSice_linux \
    --A1 A1 --A2 A2 --chr CHR --bp BP --beta --pvalue PVAL --snp SNP --stat STAT \
    --bar-levels 0.4,0.6,1 --fastscore \
    --base PRSice_input_chr_${chr_val}.txt --pheno-file !{params.pheno_file} \
    --id-delim "_" --target !{params.input_folder}chr_${chr_val}/imputation/concatenated_chr${chr_val} --quantile 20 --all-score  --print-snp --score sum --binary-target T --type bgen  --no-clump \
    --no-regress --out chr${chr_val}
  fi

  '''

}
//--ignore-fid
// beta vs or, binary-target type of outcome, quantile

process concatenation {
  publishDir params.output_folder, mode: 'copy'

  input:
  file all_scores from scores.collect()
  file all_snps from snps_included.collect()

  output:
  file '*_combined.txt' into prs_combined

  shell:
  '''
  /home/gabriela/R3.5.1/bin/Rscript !{baseDir}/combine_PRSice_res.r !{params.pheno_file} !{params.output_name}
  '''

}
