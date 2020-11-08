# This data is download form https://www.ebi.ac.uk/ena/browser/view/SRX3415289
# Reference: DOI: 10.1038/s41467-018-03369-8
# Møller HD, Mohiyuddin M, Prada-Luengo I, Sailani MR, Halling JF,
# Plomgaard P, Maretty L, Hansen AJ, Snyder MP, Pilegaard H,
# Lam HYK, Regenberg B. Circular DNA elements of chromosomal
# origin are common in healthy human somatic tissue.
# Nat Commun. 2018 Mar 14;9(1):1069. doi: 10.1038/s41467-018-03369-8.
# PMID: 29540679; PMCID: PMC5852086.

echo 'Begin to download test_data. '
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/008/SRR6315418/SRR6315418_1.fastq.gz
zcat SRR6315418_1.fastq.gz | head -n 80000 > test_1.fastq
rm SRR6315418_1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR631/008/SRR6315418/SRR6315418_2.fastq.gz
zcat SRR6315418_2.fastq.gz | head -n 80000 > test_2.fastq
rm SRR6315418_2.fastq.gz
echo 'Download test_data successfully. '
