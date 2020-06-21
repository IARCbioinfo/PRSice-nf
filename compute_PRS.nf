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
    log.info "--bgen_list                      FILE                      File containing the list of bgen files to consider"
    log.info "--PRSice_path                     PATH                      Path to PRSice executables"
    log.info "--output_folder                     PATH                      Path to output folder"
    log.info ""
    log.info "Optional arguments:"
    //log.info "--samples_file                     FILE                     File containing the samples IDs (one sample per row) to keep for the PRS analysis"
    log.info "--snp_col                     STRING                      Column name in the base file for SNP id (default value is the 1st column)"
    log.info "--chr_col                     STRING                      Column name in the base file for chromosome info (default value is the 2nd column)"
    log.info "--pos_col                     STRING                      Column name in the base file for position info (default value is the 3rd column)"
    log.info "--a1_col                     STRING                      Column name in the base file for the effect allele (default value is the 4th column)"
    log.info "--a2_col                     STRING                      Column name in the base file for the reference allele (default value is the 5th column)"
    log.info "--stat_col                     STRING                      Column name in the base file for the weights to be used in the PRS computation (default value is the 6th column)"
    log.info "--se_col                     STRING                      Column name in the base file for the standard error info (default value is the 7th column)"
    log.info "--pval_col                     STRING                      Column name in the base file for the pvalue info (default value is the 8th column)"
    log.info "--output_name                     STRING                      Suffix name to give to the outputs"

    log.info ""

    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}


params.base_file = null
params.pheno_file = null
params.bgen_list=null
//params.samples_file = null
params.PRSice_path = null
params.output_folder = null
params.output_name="res"
params.bgenix="bgenix"
params.Rscript="Rscript"
params.bgen_pattern=".bgen"
params.pval_thr="0.001,0.01,1"
params.bgen_snps_perCHR=""
params.qual_val = 0

params.snp_col=1
params.chr_col=2
params.pos_col=3
params.a1_col=4
params.a2_col=5
params.stat_col=6
params.se_col=7
params.pval_col=8
params.qual_col=9

// add a if to check if a samples file is provided if yes add the option --keep
// add a if to deal with plink files?
// add a parameter (string) allowing the user to add multiple PRSice flag options

process get_chr {

  output:
  file 'chr_list.txt' into chr_list
  file 'base_file.txt' into base_file_updated

  shell:
  '''
  awk '{print $snp,$chr,$pos,$a1,$a2,$stat,$se,$pval,$qual}' snp=!{params.snp_col} chr=!{params.chr_col} pos=!{params.pos_col} a1=!{params.a1_col} a2=!{params.a2_col} stat=!{params.stat_col} se=!{params.se_col} pval=!{params.pval_col} qual=!{params.qual_col} !{params.base_file} > base_file.txt
  tail -n +2 !{params.base_file} | awk '{print $2}' | sort | uniq  > chr_list.txt
  '''
}

process compute_PRS_per_chr {

  input:
  val chr from chr_list.splitText(by:1)
  file base from base_file_updated

  output:
  file '*.all.score' optional true into scores
  file '*.snp' optional true into snps_included

  shell:
  '''
  chr_val=!{chr}
  echo $chr_val
  bgen_file=$(grep chr${chr_val}!{params.bgen_pattern} !{params.bgen_list} | sed 's/.bgen//g')
  if [ ! -f "!{params.bgen_snps_perCHR}chr_${chr_val}_SNPs.txt" ]; then
    !{params.bgenix} -g ${bgen_file}.bgen -list > !{params.bgen_snps_perCHR}chr_${chr_val}_SNPs.txt
  fi
  awk '{if($2==c)print $0}' c=${chr_val} !{base} > chr_${chr_val}_base.txt

  n_snps=$(grep -Ff <(cat chr_${chr_val}_base.txt | cut -d ' ' -f 2,3  | sed 's/ /:/g') !{params.bgen_snps_perCHR}chr_${chr_val}_SNPs.txt | wc -l)
  if [[ ${n_snps} -ge 1 ]]
  then
    grep -Ff <(cat chr_${chr_val}_base.txt | cut -d ' ' -f 2,3  | sed 's/ /:/g') !{params.bgen_snps_perCHR}chr_${chr_val}_SNPs.txt > chr_${chr_val}_SNPs_subset.txt
    !{params.Rscript} !{baseDir}/bin/match_IDs.r chr_${chr_val}_base.txt chr_${chr_val}_SNPs_subset.txt ${chr_val}

    n_snps2=$(tail -n +2 PRSice_input_chr_${chr_val}.txt | awk '{ if (($4=="T" && $5=="A")||($4=="A" && $5=="T")||($4=="C" && $5=="G")||($4=="G" && $5=="C")) print $2, "ambig" ; else print $2 ;}' | grep -v ambig | wc -l)
    n_snps3=$(tail -n +2 PRSice_input_chr_${chr_val}.txt | awk '{ if ($9>=qual_thr) print $2 }' qual_thr=!{params.qual_val} | wc -l)
    if [[ ${n_snps2} -ge 1 && ${n_snps3} -ge 1 ]]
    then
      !{params.Rscript} !{params.PRSice_path}PRSice.R --prsice !{params.PRSice_path}PRSice_linux \
      --A1 A1 --A2 A2 --chr CHR --bp BP --beta --pvalue PVAL --snp SNP --stat STAT \
      --bar-levels !{params.pval_thr} --fastscore \
      --base PRSice_input_chr_${chr_val}.txt --pheno-file !{params.pheno_file} \
      --id-delim "_" --target ${bgen_file} \
      --quantile 20 --all-score  --print-snp --score sum --binary-target T --type bgen  --no-clump \
      --no-regress --out chr${chr_val} --base-info QUAL,!{params.qual_val}
    fi
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
  file 'prs_scores*.txt' into prs_combined
  file 'snp_file*.txt' into snps_combined

  shell:
  '''
  !{params.Rscript} !{baseDir}/bin/combine_PRSice_res.r !{params.pheno_file} !{params.output_name}
  '''

}
