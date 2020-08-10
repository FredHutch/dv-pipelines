#!/bin/bash

# get the data
# gene count matrix is downloaded by link: https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/gene_count.txt
curl https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/gene_count.txt -o gene_count.txt -retry 10
# cell annotation file is downloaded by link: https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/cell_annotate_20200119.csv

curl https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/cell_annotate_20200119.csv -o cell_annotate.csv
# gene annotation file is downloaded by link: https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/gene_annotate.csv

curl https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/gene_annotate.csv -o gene_annotate.csv
