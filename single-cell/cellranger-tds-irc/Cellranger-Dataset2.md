# Cellranger for Dataset 2

2020-06-16

# Cellranger

## Start up an EC2 instance

Ubuntu, h1.2xlarge (8 CPUs & 32GB RAM needed, plus a ton of storage)
Not a spot instance

```



ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-35-162-131-10.us-west-2.compute.amazonaws.com


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

aws s3 cp s3://test-nextflow-data/ade9b62f-ae3f-4817-81b4-ff53e1d73de1 . --recursive

```





### Figure out the Cellranger Call

READ OVER https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input#combine

AND https://www.reddit.com/r/bioinformatics/comments/ha2uag/cellranger_invalid_pathprefix_combination/

Do cellranger agg if absolutely necessary, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
   Look over Rob's code. Do I have to do that? 


Commandgen.sh:

```
ID="test1"
GENOME="/scratch/data/genome"
FAST_PATH="/scratch/data/patient10"

cellranger count --id=$ID \
  --transcriptome=$GENOME \
  --fastqs=$FAST_PATH



ID="patient18"
GENOME="/scratch/data/genome"
FAST_PATH="/scratch/data/patient18"

cellranger count --id=$ID \
  --transcriptome=$GENOME \
  --fastqs=$FAST_PATH

```








## Save off the results


For patient 10, IN PROGRESS:
```
/scratch/data/genome/test1

aws s3 cp . s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a/patient10/counts/ --recursive


```

For patient 18, TBD:

```
cd /scratch/data/patient18/patient18

aws s3 cp . s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a/patient18/counts/ --recursive


```















