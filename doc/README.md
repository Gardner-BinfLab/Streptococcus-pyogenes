## Create an profile HMM database of *S. pyogenes* from EggNOGs 5.0:

### Download bacNOG models
~~~
cd /Volumes/scratch/gardnerlab/avoidance1.0/winter_projects/spyogenes/ref
wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_hmms.tar
tar xvf 2_hmms.tar
find . -name '*.gz' > ../2.txt
sed -i 's|\./||' 2.txt
find ./2 -name '*.gz' -exec cat {} + | gzip > 2.hmm.gz
~~~

### Download S. pyogenes proteomes from NCBI
#### Get GCA accession numbers from Table S2 in https://doi.org/10.1038/s41588-019-0417-8
#### acc.txt contains the accession numbers
~~~
mkdir xml
while read i <&3; do \
    esearch -db assembly -query $i | efetch -db assembly -format docsum > xml/$i.xml
    xmllint --xpath /DocumentSummarySet/DocumentSummary/AssemblyAccession xml/$i.xml
done 3< acc.txt \
| sed 's|<[A-Za-z]\+>||g;s|</[A-Za-z]\+>|\n|g' > gcf.txt
cd /Volumes/scratch/gardnerlab/avoidance1.0/winter_projects/spyogenes/doc
../bin/datasets download genome accession --inputfile gcf.txt --exclude-seq --exclude-gff3 --exclude-rna #This produces ncbi_dataset.zip
unzip ncbi_dataset.zip
~~~

### Reduce sequence redundancy
~~~
for i in `cat gcf.txt`; do \
    awk -v var=$i 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS var "|" $0}' ncbi_dataset/data/$i/protein.faa \
    | awk '{print $1 "\n" $NF}'
done > protein.faa
cd-hit -M 6400 -c 0.99 -i protein.faa -o protein.99.faa
~~~

### Get S pyrogenes specific profile hmm models from bacNOG models
~~~
nice -10 hmmsearch --cpu 100 --tblout protein.99.2.txt ../ref/2.hmm.gz protein.99.faa
~~~

### Parse hmmsearch tblout to get the top eggNOG model for each protein
~~~
awk 'BEGIN{OFS="\t"} !/#/ {print $1,$3,$5,$6}' protein.99.2.txt \
| sort -k4,4gr | awk '!seen[$1]++ {print $2}' | sort -u \
| sed 's/\./\.hmm.gz\t/' | cut -f1 > protein.99.2.eggnog.txt
~~~

### Concatenate Bacilli models
~~~
cd ..
for i in `cat doc/protein.99.2.eggnog.txt`; do \
    cat ref/2/$i; 
done | gzip > ref/spyogenes.hmm.gz
~~~

### Rescan spyogenes proteomes
~~~
cd doc
nice -10 hmmsearch --cpu 100 --tblout protein.txt ../ref/spyogenes.hmm.gz protein.faa
awk 'BEGIN{OFS="\t"} !/#/ {print $1,$3,$5,$6}' protein.txt \
| sort -k4,4gr | awk '!seen[$1]++' \
| sed 's/\.[0-9]\+//;s/GCF/GCA/;s/|/\t/' \
| gzip > ../Jupyter/protein.spyogenes.txt.gz
~~~

### Find orphan S. pygenes proteins with no corresponding EggNOG family, and build hmms of these (if not too many)
~~~
cd ..
cat  doc/protein.99.eggnog.txt | cut -f 2 | sort | uniq > doc/protein.99.eggnog.uniq-ids.txt
hmmfetch --index 91061.hmm &&  hmmfetch -f -o 91061-spyogenes-eggnogs.hmm ./91061.hmm  ../doc/protein.99.eggnog.uniq-ids.txt

grep ^">" protein.99.faa | perl -lane 'if(/^>(\S+)/){print $1}' | sort -d > protein.99.faa-ids.txt
cat protein.99.2.txt | awk '{print $1}' | sort -d | uniq > protein.99.2.eggnog-sorted-protein-ids.txt
#damned fast way of printing lines not shared by 2 files!!!
diff --new-line-format="" --unchanged-line-format="" <(sort -d protein.99.faa-ids.txt) <(sort -d protein.99.2.eggnog-sorted-protein-ids.txt) > orphan-spyogenes-protein-ids.2.txt

esl-sfetch --index protein.99.faa
esl-sfetch -o orphan-spyogenes-protein-ids.2.faa -f protein.99.faa  orphan-spyogenes-protein-ids.2.txt
~~~

### Cluster orphan sequences into families using jackhmmer. Clusters with >=5 sequences are saved in 99084.hmm & 99084.stk
~~~
../script/fasta2families.pl -v -o orphan-fams -t 0.00001 -N 5 -m 5 -i orphan-spyogenes-protein-ids.2.faa
~~~

