grep "ID=gene:" ../Homo_sapiens.GRCh38.107.gff3 | grep "biotype=protein_coding;" | grep -o "ID=gene:\w*" | cut -d ":" -f 2 > Ensembl_all_geneID.txt