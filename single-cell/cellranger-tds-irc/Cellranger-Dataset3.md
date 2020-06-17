# Cellranger for Dataset 3

2020-06-16

# Cellranger

## Start up an EC2 instance

Ubuntu, h1.2xlarge (8 CPUs & 32GB RAM needed, plus a ton of storage)
Not a spot instance

```

ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-54-185-14-198.us-west-2.compute.amazonaws.com


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

Install Python, Nextflow, and Cellranger

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

# Install Cellranger
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

aws s3 cp s3://test-nextflow-data/ec8e5540-1132-4eba-900c-c546d8199aac . --recursive

```





### Figure out the Cellranger Call

READ OVER https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input#combine

AND https://www.reddit.com/r/bioinformatics/comments/ha2uag/cellranger_invalid_pathprefix_combination/

Do cellranger agg if absolutely necessary, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
   Look over Rob's code. Do I have to do that? 



Setup for cellranger calls:

```
mkdir -p /scratch/data/counts
cd /scratch/data/counts
mkdir -p data1
mkdir -p data2

```

Commandgen.sh for data 1:

```
ID="data1"
GENOME="/scratch/data/genome/GRCh38_CAR"
FAST_PATH="/scratch/data/data1/outs/fastq_path/CCG3JANXX"
SAMPLE="/scratch/data/data1/outs/input_samplesheet.csv"

cd /scratch/data/counts
echo "cd /scratch/data/counts/data1
" >> commandexec.sh
while IFS=, read -r col1 col2 col3
do
   command="cellranger count --id=$ID_$col2 \
   --transcriptome=$GENOME \
   --sample=$col2 \
   --fastqs=$FAST_PATH/$col2 \
   "
    echo "$command"
    echo $command >> commandexec.sh
done < <(tail -n +2 $SAMPLE)

```


Commandgen.sh for data 2

```
cd /scratch/data/counts
echo "cd /scratch/data/counts/data2" >> commandexec.sh

ID="SSC-PBMC-4-2016"
GENOME="/scratch/data/genome/GRCh38_CAR"
FAST_PATH="/scratch/data/data2"

command="cellranger count --id=$ID \
   --transcriptome=$GENOME \
   --fastqs=$FAST_PATH \
"
echo "$command"
echo $command >> commandexec.sh



chmod +x commandexec.sh
```


## Run the cellranger count command

```
screen 

export PATH=/scratch/cellranger-3.1.0:$PATH
cd /scratch/data/counts
. commandexec.sh


```



## Save off the results

```
cd /scratch
aws s3 cp comparerun1/ s3://test-nextflow-data/40d55f14-5ae8-4e7a-af9b-85d1122a0aaa/alignment/cellranger/ --recursive --dryrun
```


















