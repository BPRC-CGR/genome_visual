import pandas as pd

FAs, = glob_wildcards("input/{fa}.fasta")
BEDs, = glob_wildcards("input/{bed}.bed")

SAMPLES = [sample for sample in FAs if sample in BEDs]

print (SAMPLES)


rule final:
    input:
#        expand("output/{sample}_genes.fasta", sample = SAMPLES),
        "output/allGenes_clean.txt",
        expand("output/{sample}_genes_logic.txt", sample = SAMPLES),


rule format_files:
    input:
        fa = "input/{sample}.fasta",
        bed = "input/{sample}.bed"
    output:
        temp(touch("output/{sample}.converted")),
    conda:
        "env.yaml"
    shell:
        """
        dos2unix {input.fa} 
        dos2unix {input.bed} 
        """
    
rule format_bed:
    input:
        bed = "input/{sample}.bed"
    output:
        "output/{sample}_genes.bed",
    shell:
        """
        sed '/track/d' {input} > {output}
        """

rule fasta_extract:
    input:
        fa = "input/{sample}.fasta",
        bed = "output/{sample}_genes.bed"
    output:
        fa = "output/{sample}_genes.fasta",
        fai = temp("input/{sample}.fasta.fai"),
    conda:
        "env.yaml"
    shell:
        """
        bedtools getfasta -name -fi {input.fa} -bed {input.bed} -fo {output.fa} 
        """

rule fasta_combine:
    input:
        expand("output/{sample}_genes.fasta", sample = SAMPLES)
    output:
        temp("output/all_genes.fasta")
    shell:
        """
        cat {input} | sed 's/ /_/g' > {output}
        """

rule blast_db:
    input:
        "output/all_genes.fasta",
    output:
        temp(multiext("output/all_genes.fasta.", "ndb","nhr","nin","njs","not","nsq","ntf","nto"))
    conda:
        "env.yaml"
    shell:
        """
        makeblastdb -in {input} -dbtype nucl
        """

## Identify closely related fragment (core) to put in center.
## Calculate genes per fasta
## Identify top hit and other matche
## Calculate "matched" average
## Should do all vs all or sample vs all?

rule blast_all:
    input:
        samp = "output/all_genes.fasta",
        db = multiext("output/all_genes.fasta.", "ndb","nhr","nin","njs","not","nsq","ntf","nto")
    output:
        "output/allGenes.txt"
    conda:
        "env.yaml"
    shell:
        """
        blastn -query {input.samp} -db {input.samp} -out {output} -outfmt 6
        """

rule no_self_align:
    input:
        "output/allGenes.txt"
    output:
        "output/allGenes_clean.txt"
    shell:
        """
        # Split to get values from blast
        awk -F'\t' 'BEGIN {{OFS="\t"}} {{
            split($1, field1, "::")
            split(field1[2], value1, "[:-]")
            split($2, field2, "::")
            split(field2[2], value2, "[:-]")
            print $0, field1[1], value1[1], value1[2], value1[3], field2[1], value2[1], value2[2], value2[3]
        }}' {input} | awk -F"\t" '$14!=$18' | awk -F'\t' 'BEGIN {{OFS="\t"}} {{
            for (i=13; i<=20; i++) {{
                printf "%s%s", $i, (i==20 ? "\\t" : OFS)
            }}
            for (i=3; i<=12; i++) {{
                printf "%s%s", $i, (i==12 ? "\\n" : OFS)
            }}
        }}' | awk '{{
              key = $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8; #use the first 8 columns to create a key
              if (key == prev_key) {{
                  if ($9 > max9) max9 = $9;
                  min13 = ($13 < min13) ? $13 : min13;
                  max14 = ($14 > max14) ? $14 : max14;
              }} else {{
                  if (NR > 1) print prev_key, max9, min13, max14;
                  prev_key = key;
                  max9 = $9;
                  min13 = $13;
                  max14 = $14;
              }}
          }}
          END {{
              if (NR > 1) print prev_key, max9, min13, max14;
          }}' > {output}
        """

rule genes_extract:
    input:
        allGenes = "output/allGenes_clean.txt",
        bed = "output/{sample}_genes.bed"
    output:
        "output/{sample}_genes_logic.txt",
    params:
        "output/{sample}_genes_logic.temp",
    shell:
        """
        touch {params}
        newFile=$(echo $RANDOM | tr '[0-9]' '[a-z]')
        awk -v smp="{wildcards.sample}" '$2==smp {{print $0, $8 - $7, $11 - $10}}' {input.allGenes} | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13/$12*100}}' | sort -k2,2 -k1,1 -k6,6 -k5,5 | awk '$10>50' > output/{wildcards.sample}"_"$newFile".csv"
        sed 's/_ / /g;s/_gene//g' output/{wildcards.sample}"_"$newFile".csv" > {output}
#        awk -F"\t" '{{print $1}}' output/{wildcards.sample}"_"$newFile".csv" | awk '!s[$0]++' | for i in $(cat); do
#            awk -F"\t" -v var=$i '$1==var {{print $6}}' output/{wildcards.sample}"_"$newFile".csv" | sort | uniq -c | awk '$1==1' | for z in $(cat); do
#                awk -F"\t" -v var=$i -v var2=$z '$1==var && $6==var2' output/{wildcards.sample}"_"$newFile".csv" >> {params}
#            done
#            awk -F"\t" -v var=$i '$1==var {{print $6}}' output/{wildcards.sample}"_"$newFile".csv" | sort | uniq -c | awk '$1>1 {{print $2}}' | for y in $(cat); do
#                awk -F"\t" -v var=$i -v var2=$y '$1==var && $6==var2' output/{wildcards.sample}"_"$newFile".csv" 
#            done
#        done  
#        awk '!s[$0]++' {params} > {output}
        rm output/{wildcards.sample}"_"$newFile".csv" {params} 
        """
'''
# Map out what is best

        # sort based on name, score and evalue
        # group the the first hit based on column 1
        awk '$1!=$2' {input} | sort -k1,1 -k12,12nr -k11,11n | sort --merge -u  -k1,1 > {output}
#| sed 's/\t/@/1;s/::.*@/\t/g;s/\t/@/2;s/::.*@/\t/g' | awk '$4>1000' | cut -f1-4 > {output}

grep " gene" HLA-DRB-regio__GrCh38.bed | sed 's/ gene//g' > genes.bed 
grep -v "Exon" h1tg000011l_EAW.bed | sed '2d;/track/d' >> genes.bed 

while read -r g1 g2 perc len; do
	echo $(awk -v var=${g1} '$4==var {print $1,$2,$3}' genes.bed) $(awk -v var=${g2} '$4==var {print $1,$2,$3}' genes.bed) $g1 $g2 $perc
done < matched_genes  | tr " " "\t" > track_link.csv 

 awk -F'\t' 'BEGIN {{OFS="\t"}} {{
            split($1, field1, "::")
            split(field1[2], value1, "[:-]")
            split($2, field2, "::")
            split(field2[2], value2, "[:-]")
            print $0, field1[1], value1[1], value1[2], value1[3], field2[1], value2[1], value2[2], value2[3]
        }}'

awk '{{
              key = $1 " " $2 " " $5 " " $6 " " $7 " " $8; #use the first 8 columns to create a key
              if (key == prev_key) {{
                  if ($9 > max9) max9 = $9;
                  min13 = ($13 < min13) ? $13 : min13;
                  min15 = ($15 < min15) ? $15 : min15;
                  max14 = ($14 > max14) ? $14 : max14;
                  max16 = ($16 > max16) ? $16 : max16;
              }} else {{
                  if (NR > 1) print prev_key, max9, min13, max14, min15, max16;
                  prev_key = key;
                  max9 = $9;
                  min13 = $13;
                  min15 = $15;
                  max14 = $14;
                  max16 = $16;
              }}
          }}
          END {{
              if (NR > 1) print prev_key, max9, min13, max14, min15, max16;
          }}' 
'''
