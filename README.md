# nf-rnachrom-voting

"Smart" annotation of RNA parts of contacts

bash voting.sh -q /path/to/genes_file_bedrc/dir -g RNA_Chrom_genes_hg38 -d 0 -e /path/to/contacts_data/dir -c contacts_data.tab.rc -o /path/to/output/dir

example 1:

bash voting.sh -q /home/snap/projects/lncRNA_app/rnachrom-new-voting-pipeline-master/genes -g RNA_Chrom_genes_hg38 -d 0 -e /home/snap/projects/lncRNA_app/rnachrom-new-voting-pipeline-master/data -c filtered_SRR17331253_K562_CTCF_RedChIP_rep2_Unique_RNA_true_strand.tab.rc -o /home/snap/projects/lncRNA_app/rnachrom-new-voting-pipeline-master/output_filtered_SRR17331253_K562_CTCF_RedChIP_rep2_Unique_RNA_true_strand

example 2:

bash voting.sh -q /home/snap/projects/lncRNA_app/rnachrom-new-voting-pipeline-master/genes -g test_genes -d 0 -e /home/snap/projects/lncRNA_app/rnachrom-new-voting-pipeline-master/data -c test_header.tab.rc -o /home/snap/projects/lncRNA_app/rnachrom-new-voting-pipeline-master/output_test_header
