# Requires trimotatic, bwa, samtools, VarScane, and variantCombiner.py

import os
from datetime import datetime
import time
import re

configfile: "/home/nate/PycharmProjects/iVar/snakemake_config.yaml"

in_dir = config["in_dir"]
ext="fastq.gz"
ref_sequence = config["ref_sequence"]
gz = True if ".gz" in ext else False
#s = datetime.now().strftime("%Y.%m.%d."+str(int(time.time())))
#out_dir="/gpfs/home/natem/analysis/"+s
out_dir = config["out_dir"]
s_locations = config["script_locations"]

SAMPLES = []
for dirname, dirnames, filenames in os.walk(in_dir):
     filenames.sort()
     for f in range(0, len(filenames)-1, 2):
          if ext in filenames[f] and ext in filenames[f+1]:
               new_name1 = re.sub('_S\d{1}_L\d{3}_R', '_R', filenames[f])
               new_name1 = re.sub('_\d{3}', '', new_name1)
               os.rename(os.path.join(dirname,filenames[f]), os.path.join(dirname, new_name1))
               new_name2 = re.sub('_S\d{1}_L\d{3}_R', '_R', filenames[f+1])
               new_name2 = re.sub('_\d{3}', '', new_name2)
               os.rename(os.path.join(dirname,filenames[f+1]), os.path.join(dirname, new_name2))
               generalName = new_name1.replace("."+ext, "").replace("_R1", "")[:-1]
               if generalName not in SAMPLES:
                   SAMPLES.append( new_name1.replace("."+ext, "").replace("_R1", "")[:-1] )


def get_read2(wildcards):
    print(wildcards)
     
rule all:
    input:
        expand("{out_dir}/_output/{sample}.csv", out_dir = out_dir, sample = SAMPLES )
 
rule variant_combining:
    input:
        expand("{{out_dir}}/_variants/{{sample}}{rep}.csv", rep=config["replicates"])
    output:
        "{out_dir}/_output/{sample}.csv"
    shell:
        "python3 {s_locations}/variantCombiner.py -r {config[virus_map]} -o {output} -t {config[vco_test]} -p {config[vco_pvalue]} -f {config[vco_minFreq]} -i {input[0]} {input[1]}"

if config["vca_compareWithReference"]:
    rule variant_calling:
        input:
            "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam"
        output:
            "{out_dir}/_variants/{sample}.csv"
        shell:
            #"module load samtools &&"
            "samtools faidx {ref_sequence} &&"
            "samtools mpileup -f {ref_sequence} {input[0]} | java -jar {s_locations}/VarScan.v2.3.9.jar pileup2snp --min-coverage {config[vca_minCoverage]} --strand-filter 0 --min-var-freq {config[vca_minVarFreq]} --min-avg-qual {config[vca_minAvgQual]} --p-value {config[vca_pvalue]} --output-vcf 1 > {output}"
else:
    rule variant_calling: 
        input:
            "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam",
            "{out_dir}/_consensus/{sample}.consensus.fasta"
        output:
            "{out_dir}/_variants/{sample}.csv"
        shell:
            #"module load samtools &&"
            "samtools faidx {ref_sequence} &&"
            "samtools mpileup -f {input[1]} {input[0]} | java -jar {s_locations}/VarScan.v2.3.9.jar pileup2snp --min-coverage {config[vca_minCoverage]} --strand-filter 0 --min-var-freq {config[vca_minVarFreq]} --min-avg-qual {config[vca_minAvgQual]} --p-value {config[vca_pvalue]} --output-vcf 1 > {output}"

rule generate_consensus:
    input:
        "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam"
    output:
        "{out_dir}/_consensus/{sample}.consensus.fasta"
    shell:
        #"module load samtools &&"
        #"module load tabix &&"
    	"samtools mpileup -uf {ref_sequence} {input} | bcftools call -mv -Oz -o {input}.vcf.gz &&"
    	"tabix {input}.vcf.gz &&"
        "cat {ref_sequence} | bcftools consensus {input}.vcf.gz > {output}"

rule trim_reads:
    input:
        "{out_dir}/_aligned_bams/{sample}.aligned.sorted.bam",
        "{bed}".format( bed = config["tr_path_to_bed"] )
    output:
        "{out_dir}/_aligned_bams/{sample}.trimmed.aligned.sorted.bam",
        "{out_dir}/_coverage/{sample}.coverage.csv",
        "{out_dir}/_output/{sample}.coverage.png"
    params:
        rep="{out_dir}/_reports/{sample}.alignreport.txt",
        tmp="{out_dir}/_tmp/"
    shell:
        "mkdir -p $(dirname {params[0]}) &&"
        "mkdir -p {params[1]} &&"
        "python {s_locations}/align_trim.py {input[1]} --report {params[0]} < {input[0]} | samtools view -F 4 -Sb - | samtools sort -T {params[1]} - -o {output[0]} &&"
        "samtools depth -a {output[0]} > {output[1]} &&"
        "python3 {s_locations}/coverageGraph.py {output[1]} {output[2]}"

rule align_reads:
    input:
        "{out_dir}/_reads/{sample}_R1.fastq",
        "{out_dir}/_reads/{sample}_R2.fastq"
    output:
        "{out_dir}/_aligned_bams/{sample}.aligned.sorted.bam"
    params:
        tmp = "{out_dir}/_tmp/"
    shell:
        #"mkdir -p {params} &&"
        "bwa mem {ref_sequence} {input[0]} {input[1]} | samtools view -F 4 -Sb | samtools sort -T {params} -o {output}"

rule extract_read_1:
    input:
        "{in_dir}".format(in_dir = in_dir)+"/{sample}_R1.fastq.gz"
    output:
        "{out_dir}/_reads/{sample}_R1.fastq"
    shell:
        "if [[ {input} == *.gz ]]; then gunzip -c {input} > {output}; else cp {input} {output};fi;"

rule extract_read_2:
    input:
        "{in_dir}".format(in_dir = in_dir)+"/{sample}_R2.fastq.gz"
    output:
        "{out_dir}/_reads/{sample}_R2.fastq"
    shell:
        "if [[ {input} == *.gz ]]; then gunzip -c {input} > {output}; else cp {input} {output}; fi"
