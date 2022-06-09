# Commands Used

## Info Table

| Program    | Version |
|:----------:|:-------:|
| Centrifuge | 1.0.4   |
| CLARK      | 1.2.6.1 |
| Kaiju      | 1.7.4   |
| Kraken2    | 2.1.2   |
| MetaCache  | 2.0.0   |

## CAMI Dataset

Download samples 0 to 9.

```bash
for i in {0..9}; do
    echo "Began $i"
    mkdir s$i
    java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT s$i -p 2017.12.29_11.37.26_sample_$i/reads/anonymous_reads.fq.gz
    mv s$i/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_$i/reads/anonymous_reads.fq.gz s$i/s$i.fq.gz
    
    echo "Splitting $i reads"
    zcat s$i/s$i.fq.gz | tee >(awk 'c-->0;$0~s{if(b)for(c=b+1;c>1;c--)print r[(NR-c+1)%b];print;c=a}b{r[NR%b]=$0}' b=0 a=3 s="/1" > s$i/s${i}_R1.fq) | awk 'c-->0;$0~s{if(b)for(c=b+1;c>1;c   --)print r[(NR-c+1)%b];print;c=a}b{r[NR%b]=$0}' b=0 a=3 s="/2" > s$i/s${i}_R2.fq
    java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT s$i -p 2017.12.29_11.37.26_sample_$i/reads/reads_mapping.tsv.gz
    mv s$i/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_$i/reads/reads_mapping.tsv.gz s$i/s${i}_reads_mapping.tsv.gz
    
    echo "Sorting $i mapping"
    zcat s$i/s${i}_reads_mapping.tsv.gz | awk -F "\t" '$1 ~ /1$/' | awk -F "\t" '{gsub(/\/1/,"",$1)}1' OFS="\t" | sort -t$'\t' -n -k1.4 > s$i/s${i}_reads_mapping.tsv
    rm -r s$i/19122017_mousegut_scaffolds s$i/s${i}_reads_mapping.tsv.gz
done
```

## Command for RAM and execution time monitoring.

```bash
/usr/bin/time -v <COMMAND> 2>&1 | tee <LOG_FILE>
```

## Database Criation

### Genome and proteins download

assembly_summary.txt downloaded in 01/07/2021

There were 379 - Achaeas and  21887 - Bacterias

```bash
genome_downloader.sh -f assembly_summary_archaea_01-07.txt -d archaea_fasta
genome_downloader.sh -f assembly_summary_archaea_01-07.txt -d archaea_gbff -p
genome_downloader.sh -f assembly_summary_bac_01-07.txt -d bacteria_fasta
genome_downloader.sh -f assembly_summary_bac_01-07.txt -d bacteria_gbff -p
```

- Order of execution:
    1. Centrifuge
    2. Kraken2
    3. CLARK
    4. Kaiju
    5. MetaCache

### Centrifuge

Lines 349 to 359 of "centrifuge-download" supressed to avoid re-download of the genomes.

```bash
# Not timed
centrifuge-download -o taxonomy taxonomy
centrifuge-download -P 5 -o library -m -d "archaea,bacteria" refseq > seqid2taxid.map
cat library/*/*.fna > input-sequences.fna

# Timed
centrifuge-build -p 5 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna ab # LOG: centrifuge/centrifuge-build.log
```

### Kraken2

```bash
# Not timed
kraken2-build --download-taxonomy --db .  --threads 5
find ../genomes/ -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db .

#Timed
kraken2-build --build --db . --threads 5 # LOG: kraken2/kraken2-build.log
```

### MetaCache

MetaCache has no option to choose how many cores to be used.

```bash
# Not timed
download-ncbi-taxonomy ncbi_taxonomy

#Timed
metacache build mydatabase genomes/ -taxonomy ncbi_taxonomy # LOG: metacache/metacache_build.log
```

### Kaiju

```bash
#Timed
kaiju-makedb -s refseq -t 5 --no-download # LOG: kaiju/kaiju-makedb.log
```

### CLARK - CANCELLED

Used an empty fastq file to avoid classification and only measure the database creation.

14/07/2021 - Server crashed after about 40h of execution. Log says it was in the database creation phase with about 394GB of RAM -> LOG:.clark/classify_metagenome_db_1.log
15/07/2021 - Server crashed after about 25h of execution. Log says it was still processing the genomes with about 350GB of RAM -> LOG:.clark/classify_metagenome_db.log

```bash
# Not timed
download_taxondata.sh taxonomy

#Timed
set_targets.sh . custom # LOG: clark/set_targets.log
classify_metagenome.sh -O fake.fq -R . -m 0 -n 5 # LOG: clark/classify_metagenome_db.log
```

## Classification

### CAMI

#### Centrifuge

```bash
centrifuge -x ab -p 5 --sample-sheet ../../classified/centrifuge/centrifuge_sheet.tsv -q # LOG: centrifuge/centrifuge.log
```

#### Kraken2

```bash
bash run_kraken2.sh # LOG: kraken2/kraken2.log
```

##### run_kraken2.sh

```bash
for i in {0..9}; do
    /usr/bin/time -v kraken2 --threads 5 --db . --paired ../../cami_reads/s${i}/s${i}_R{1,2}.fq --report ../../classified/kraken2/${i}/${i}.report.krk --output ../../classified/kraken2/${i}/${i}.out.krk 2>&1 | tee ../../logs/kraken2/kraken2_s${i}.log
done
```

#### MetaCache

```bash
metacache query mydatabase ../../cami_reads/s?/s{0..9}_R1.fq ../../cami_reads/s?/s{0..9}_R2.fq -pairfiles -split-out ../../classified/metacache/metacache -taxids-only -omit-ranks -threads 5 -lowest species -separator '\t' # LOG: metacache/metacache_query.log
```

#### Kaiju

##### run_kaiju.sh

```bash
R1="../../cami_reads/s0/s0_R1.fq"
R2="../../cami_reads/s0/s0_R2.fq"
out="../../classified/metacache/s0.kaiju"

for i in {1..9}; do
    R1="${R1},../../cami_reads/s${i}/s${i}_R1.fq"
    R2="${R2},../../cami_reads/s${i}/s${i}_R2.fq"
    out="${out},../../classified/metacache/s${i}.kaiju"
done

kaiju-multi -z 5 -t nodes.dmp -f refseq/kaiju_db_refseq.fmi -i $R1  -j $R2 -o $out 2>&1 # LOG: kaiju/kaiju-multi.log
```


### REAL

#### Centrifuge

```bash
centrifuge -x ab -p 5 --sample-sheet ../../classified/centrifuge/centrifuge_sheet.tsv -q # LOG: centrifuge/centrifuge.log
```

#### Kraken2

```bash
bash run_kraken2.sh # LOG: kraken2/kraken2.log
```

##### run_kraken2.sh

```bash
for i in {37..48}; do
    /usr/bin/time -v kraken2 --threads 5 --db . --paired ../../trimmed/S${i}/S${i}_R{1,2}_P.fq --report ../../classified/kraken2/S${i}/S${i}.report.krk --output ../../classified/kraken2/S${i}/S${i}.out.krk 2>&1 | tee ../../logs_real/kraken2/S${i}.krklog
done
```

#### MetaCache

```bash
metacache query mydatabase ../../cami_reads/s?/s{37..48}_R1.fq ../../cami_reads/s?/s{37..48}_R2.fq -pairfiles -split-out ../../classified/metacache/metacache -taxids-only -omit-ranks -threads 5 -lowest species -separator '\t' # LOG: metacache/metacache_query.log
```

#### Kaiju

```bash
R1="../../trimmed/S37/S37_R1_P.fq"
R2="../../trimmed/S37/S37_R2_P.fq"
out="../../classified/kaiju/S37.kaiju"
for i in {38..48}; do
    R1="${R1},../../trimmed/S${i}/S${i}_R1_P.fq"
    R2="${R2},../../trimmed/S${i}/S${i}_R2_P.fq"
    out="${out},../../classified/kaiju/S${i}.kaiju"
done

kaiju/bin/kaiju-multi -z 5 -t nodes.dmp -f refseq/kaiju_db_refseq.fmi -i $R1 -j $R2 -o $out # LOG: kaiju/kaiju.log
```
