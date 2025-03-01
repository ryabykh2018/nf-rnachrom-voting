Example command: 
`bash voting.sh -g hg38 -d 100 -i ./data/contacts.bed -o ./output/`

Where
 - hg38 - genome version, genes annotation should be in ./genes/hg38.genes.bed; tab-delimited with columns chr, start, end, name, score(not used), strand.
 - -d 100 - clustering distance, it should be larger, than max read length
 - -i contacts.bed - bed file with all contacts, tab-delimited file with columns: chr, start, end, id, score(not used), strand.
 - -o ./output/ - output dir
    
Several tools should be installed:
bedmap version 2.4.39
bedtools version 2.30.0
and python packages from requirements.txt