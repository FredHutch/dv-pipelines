## using the gencode ref: https://www.gencodegenes.org/human/, *ALL* regions

## ------------------------------
## Generate a transcript to gene mapping from a GTF
zcat gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq | sed 's/\"//g' > gencode.v33.chr_patch_hapl_scaff.annotation.t2g

## Generate a dictionary of chromosome lengths
gunzip GRCh*.fa.gz
faidx GRCh*.fa -i chromsizes > GRCh38.p13.chromsizes
gzip GRCh*.fa