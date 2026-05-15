awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="CHROMOSOME") chr=i
        else if ($i=="GENOME_START") start=i
        else if ($i=="GENOME_STOP") stop=i
        else if ($i=="GENE_SYMBOL") gene=i
    }
    next
}
{
    c=$chr
    if (c=="MT") c="M"
    if (c !~ /^chr/) c="chr" c

    print c, $start, $stop, $gene
}' Cosmic_CancerGeneCensus_v103_GRCh38.tsv \
> cosmic_genes.tsv
