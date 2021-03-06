# Cellranger for Dataset 1

2020-06-16

# Cellranger

## Start up an EC2 instance

Ubuntu, h1.2xlarge (8 CPUs & 32GB RAM needed, plus a ton of storage)
Not a spot instance

```


ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-34-222-121-180.us-west-2.compute.amazonaws.com



```



```
# list the NVMe devices
lspci

# Show existing drives + partitions
fdisk -l

# make a file system
mkfs -t ext4 /dev/xvdb

mkdir /scratch
mount /dev/xvdb /scratch
df -h

```

Install Python and Nextflow

```
sudo -i
apt-get -y update
apt-get -y install unzip awscli python3-pip csvkit
python3 -m pip install boto3 pandas

apt install -y bzip2 wget openjdk-8-jdk docker.io

export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

cd /root
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /root/miniconda
rm Miniconda3-latest-Linux-x86_64.sh
/root/miniconda/bin/conda install -c conda-forge -y awscli
PATH="/root/miniconda/bin:${PATH}"
echo $PATH

mkdir -p /opt/inst
cd /opt/inst
wget -qO- https://get.nextflow.io | bash
mv nextflow /usr/local/bin

PATH="/root/miniconda/bin:${PATH}"
```


Install Cellranger

```
cd /scratch
aws s3 cp s3://hutch.cloud.dev/cellranger/install/ . --recursive

tar -zxvf cellranger-3.1.0.tar.gz

export PATH=/scratch/cellranger-3.1.0:$PATH
```


## Get the data

```
cd /scratch
## Get the data

mkdir -p data
cd data

aws s3 cp s3://test-nextflow-data/3a27dce3-08b5-4957-a9ef-4b273a43dc2f . --recursive
```





### Figure out the Cellranger Call

READ OVER https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input#combine

AND https://www.reddit.com/r/bioinformatics/comments/ha2uag/cellranger_invalid_pathprefix_combination/


Needed variables:
--id (unique ID run string)
--transcriptome
   > probably genome/hg38_merkerCellPolyomavirus_CellRanger_v2
--fastqs
   > dataset1/patient1/outs/fastq_path/HKMNCBCXY/*/<files>/

CHECK INTO:
--sample ? <- look for sample sheet given to cellranger mkfastq
   > patient1/outs/input_samplesheet.csv


TBD variables:
--chemistry <- should be automatic. Make sure it found a 3' result

```
screen


cellranger count --id=$ID \
                   --transcriptome=$GENOME \
                   --sample=$SAMPLE \
                   --fastqs=$FAST_PATH



cd /scratch





```

Commandgen.sh:

```
ID="test1"
GENOME="/scratch/data/genome/hg38_merkerCellPolyomavirus_CellRanger_v2"
FAST_PATH="/scratch/data/patient1/outs/fastq_path/HKMNCBCXY"
SAMPLE="/scratch/data/patient1/outs/input_samplesheet.csv"


while IFS=, read -r col1 col2 col3
do
   command="cellranger count --id=$ID_$col2 \
   --transcriptome=$GENOME \
   --sample=$col2 \
   --fastqs=$FAST_PATH/$col2"
    echo "$command"
done < <(tail -n +2 $SAMPLE)

```

Input samplesheet:

```
Lane,Sample,Index,,,,
1,GS_Tumor_9_2013_2559,SI-GA-A2,,,,
1,GS_Tumor_10_2016_6525,SI-GA-A3,,,,
2,GS_PBMC_2_2015_GRN0304,SI-GA-A4,,,,
2,GS_PBMC_2_2016_GRN0535,SI-GA-A5,,,,
2,GS_PBMC_10_2016_GRN0760,SI-GA-A6,,,,
2,GS_PBMC_8_2013_MCCB1050,SI-GA-A7,,,,
```

Are these the samples? Looks like it

GS_PBMC_10_2016_GRN0760
GS_PBMC_2_2016_GRN0535
GS_Tumor_10_2016_6525
GS_PBMC_2_2015_GRN0304
GS_PBMC_8_2013_MCCB1050
GS_Tumor_9_2013_2559







## Save off the results

```
cd /scratch/data

aws s3 cp . s3://test-nextflow-data/70d766be-5da2-4cfc-9515-9bf416ce964b/patient1/counts/ --exclude "*" --include "GS*" --recursive

```


## Do the missing cellranger count call

```
screen

export PATH=/scratch/cellranger-3.1.0:$PATH
cellranger count --id=GS_PBMC_8_2013_MCCB1050    --transcriptome=/scratch/data/genome/hg38_merkerCellPolyomavirus_CellRanger_v2    --sample=GS_PBMC_8_2013_MCCB1050    --fastqs=/scratch/data/patient1/outs/fastq_path/HKMNCBCXY/GS_PBMC_8_2013_MCCB1050

```

Save off these results as well:

```
cd /scratch/data/GS_PBMC_8_2013_MCCB1050

aws s3 cp . s3://test-nextflow-data/70d766be-5da2-4cfc-9515-9bf416ce964b/patient1/counts/GS_PBMC_8_2013_MCCB1050/ --recursive

```


## Restore the other results

```
screen
cd /scratch/data

aws s3 cp s3://test-nextflow-data/70d766be-5da2-4cfc-9515-9bf416ce964b/patient1/counts/ . --recursive --dryrun

```



## Make a cellranger aggr command


Do cellranger agg if absolutely necessary, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
   Look over Rob's code. Do I have to do that? 

Make a CSV file:

```
library_id,molecule_h5
GS_PBMC_10_2016_GRN0760,/scratch/data/GS_PBMC_10_2016_GRN0760/outs/molecule_info.h5
GS_PBMC_2_2016_GRN0535,/scratch/data/GS_PBMC_2_2016_GRN0535/outs/molecule_info.h5
GS_Tumor_10_2016_6525,/scratch/data/GS_Tumor_10_2016_6525/outs/molecule_info.h5
GS_PBMC_2_2015_GRN0304,/scratch/data/GS_PBMC_2_2015_GRN0304/outs/molecule_info.h5
GS_PBMC_8_2013_MCCB1050,/scratch/data/GS_PBMC_8_2013_MCCB1050/outs/molecule_info.h5
GS_Tumor_9_2013_2559,/scratch/data/GS_Tumor_9_2013_2559/outs/molecule_info.h5


GS_PBMC_10_2016_GRN0760  
GS_Tumor_10_2016_6525
GS_PBMC_2_2015_GRN0304  
 GS_Tumor_9_2013_2559
GS_PBMC_2_2016_GRN0535
```

Setup for the call

```
mkdir -p /scratch/data/agg
cd /scratch/data/agg
nano #copy in the aggregate CSV to 'aggregate.csv'
```

```
screen


export PATH=/scratch/cellranger-3.1.0:$PATH
ID="patient1run1"
CSV="/scratch/data/agg/aggregate.csv"
cellranger aggr \
  --id=$ID \
  --csv=$CSV

```


### Save the results


```
cd /scratch/data/agg

aws s3 cp . s3://test-nextflow-data/70d766be-5da2-4cfc-9515-9bf416ce964b/patient1/agg/ --exclude "*" --include "aggregate.csv" --recursive


cd /scratch/data/agg/patient1run1

aws s3 cp . s3://test-nextflow-data/70d766be-5da2-4cfc-9515-9bf416ce964b/patient1/agg/ --recursive

```













