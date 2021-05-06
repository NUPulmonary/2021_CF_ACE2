import pandas as pd


DATA_DIR = "data"
PROJECT = "CF_bulk"


SAMPLES = pd.read_csv("samples.csv")


rule all:
    input:
        expand("{data_dir}/fastqc", data_dir=DATA_DIR),
        expand("{data_dir}/fastqc_trimmed", data_dir=DATA_DIR),
        expand("{data_dir}/star/report.csv", data_dir=DATA_DIR),
        expand("{data_dir}/star/report.txt", data_dir=DATA_DIR),
        expand("{data_dir}/star/{sample}_Aligned.sortedByCoord.out.bam.bai", data_dir=DATA_DIR, sample=SAMPLES.Name),
        "short_ace2.txt",


rule short_ace2:
    input:
        expand("{data_dir}/star/{sample}_Aligned.sortedByCoord.out.bam", data_dir=DATA_DIR, sample=SAMPLES.Name)
    output:
        "short_ace2.txt"
    params:
        slurm__skip=True
    shell:
        """
        rm -f {output}

        echo -e "Sample\tExon9a\tExon9" >> {output}
        dir=`dirname {input[0]}`
        for i in `ls $dir/*.bam`; do 
            sample=`basename $i`
            echo -en "${{sample%_001_Aligned.sortedByCoord.out.bam}}\t" >> {output}; 
            samtools view -c -h $i chrX:15580200-15580400 | tr '\n' '\t' >> {output}; 
            samtools view -c -h $i chrX:15578089-15578315 >> {output}; 
        done
        """


rule align_report:
    input:
        expand("{{data_dir}}/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES.Name)
    output:
        "{data_dir}/star/report.csv",
        "{data_dir}/star/report.txt"
    params:
        slurm__skip = True
    shell:
        """
        star_dir=`dirname {input[0]}`

        echo "sample
# total
# uni map
% uni map
# multi map
% multi map
# too many
% too many
% unmap mism
% unmap short
% unmap other
# chimeric
% chimeric" | perl -0777 -pe 's/\\n/,/g' \
            | sed 's/,$/\\n/' > $star_dir/report.csv

        for i in `ls $star_dir/*Log.final.out`; do
            sample=`basename $i`
            sample=${{sample%_Log.final.out}}
            out=`grep "mapped\|input\|chimeric" $i \
                    | grep -v "length" | cut -d"|" -f 2 \
                    | perl -0777 -pe 's/\s+/,/g' | sed 's/,$//'`
            echo $sample$out >> $star_dir/report.csv
        done
        column -s, -t $star_dir/report.csv > $star_dir/report.txt
        """


rule index:
    input: "{data_dir}/star/{sample}_Aligned.sortedByCoord.out.bam"
    output: "{data_dir}/star/{sample}_Aligned.sortedByCoord.out.bam.bai"
    params:
        slurm__hours = 1,
        slurm__cores = 6,
        slurm__mem = 40
    shell:
        """
        module purge all
        module load samtools

        samtools index {input}
        """


rule align:
    input: "{data_dir}/fastq_trimmed/{sample}_trimmed.fastq.gz"
    output: "{data_dir}/star/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        slurm__hours = 1,
        slurm__cores = 12,
        slurm__mem = 64
    shell:
        """
        module purge all
        module load STAR/2.7.5

        if [[ ! -d `dirname {output}` ]]; then
            mkdir `dirname {output}`
        fi
        STAR --runMode alignReads \
            --genomeDir /projects/b1038/Pulmonary/nmarkov/GRCh38-gencode33/index \
            --runThreadN 12 \
            --readFilesCommand 'gunzip -c' \
            --readFilesIn {input} \
            --outFileNamePrefix `dirname {output}`/{wildcards.sample}_ \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000
        """

rule fastqc_trimmed:
    input: expand("{{data_dir}}/fastq_trimmed/{sample}_trimmed.fastq.gz", sample=SAMPLES.Name)
    output: directory("{data_dir}/fastqc_trimmed")
    params:
        slurm__hours = 1,
        slurm__cores = 12,
        slurm__mem = 36
    shell:
        """
        module purge all
        module load fastqc
        module load multiqc

        mkdir {output}
        fastqc `dirname {input[0]}`/*.fastq.gz \
            -o {output} \
            -t 12

        multiqc {output}/* -o {output}
        """

rule trim:
    input:
        "{data_dir}/fastq/{sample}.fastq.gz"
    output: "{data_dir}/fastq_trimmed/{sample}_trimmed.fastq.gz"
    params:
        slurm__hours = 1,
        slurm__cores = 12,
        slurm__mem = 12
    shell:
        """
        module purge all
        module load java

        if [[ ! -d `dirname {output}` ]]; then
            mkdir `dirname {output}`
        fi
        java -jar /projects/b1038/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE \
            {input[0]} \
            {output} \
            -threads 12 \
            TRAILING:30 MINLEN:20
        """

rule fastqc:
    input: expand("{{data_dir}}/fastq/{sample}.fastq.gz", sample=SAMPLES.Name)
    output: directory("{data_dir}/fastqc")
    params:
        slurm__hours = 1,
        slurm__cores = 12,
        slurm__mem = 36
    shell:
        """
        module purge all
        module load fastqc
        module load multiqc

        mkdir {output}
        fastqc `dirname {input[0]}`/*.fastq.gz \
            -o {output} \
            -t 12

        multiqc {output}/* -o {output}
        """
