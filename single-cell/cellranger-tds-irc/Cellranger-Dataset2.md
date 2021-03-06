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





# Cellranger Aggr, 2020-06-30


Ubuntu, i3.2xlarge (8 CPUs & 61GB RAM needed, plus a ton of storage). Spot instance

```
ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-18-237-87-43.us-west-2.compute.amazonaws.com


```




```
# list the NVMe devices
lspci

# Show existing drives + partitions
fdisk -l

# make a file system
mkfs -t ext4 /dev/nvme0n1

mkdir /scratch
mount /dev/nvme0n1 /scratch
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

aws s3 cp s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a . --recursive

```



## Make a cellranger aggr command


Do cellranger agg if absolutely necessary, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
   Look over Rob's code. Do I have to do that? 

Make a CSV file:

```
library_id,molecule_h5
patient10,/scratch/data/patient10/counts/outs/molecule_info.h5
patient18,/scratch/data/patient18/counts/outs/molecule_info.h5
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
ID="run1"
CSV="/scratch/data/agg/aggregate.csv"
cellranger aggr \
  --id=$ID \
  --csv=$CSV

```



### Copy the results


```
cd run1
aws s3 cp . s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a/agg --recursive
```



## Run the pubweb code against the aggr results

```bash

ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-54-218-27-88.us-west-2.compute.amazonaws.com

ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-34-222-188-113.us-west-2.compute.amazonaws.com

ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-44-234-59-250.us-west-2.compute.amazonaws.com




```

**Setup**

```bash
# list the NVMe devices
lspci

# Show existing drives + partitions
fdisk -l

# make a file system
mkfs -t ext4 /dev/nvme0n1

mkdir /scratch
mount /dev/nvme0n1 /scratch
df -h

```

**Install the necessary packages & libraries**

```bash
add-apt-repository -y ppa:deadsnakes/ppa

apt -y update
apt-get -y install python3.8 awscli python3-pip csvkit unzip

python3.8 --version

update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1
update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1
update-alternatives --config python
update-alternatives --config python3

cd /usr/lib/python3/dist-packages
ln -s apt_pkg.cpython-{36m,38m}-x86_64-linux-gnu.so
ln -s apt_pkg.cpython-36m-x86_64-linux-gnu.so apt_pkg.so
ls | grep apt_pkg

cd $HOME
python -m pip install boto3 pandas scanpy tables h5py mygene shapely monet

```

Get the pubweb code

```bash
cd $HOME
mkdir -p /opt/pubweb
cd /opt/pubweb

# copy in the code
aws s3 cp s3://dv-code-dev/pubweb/ . --recursive

nano #copy in invoke-cellranger.py

# install the library
python -m pip install .


```

Get the data

```bash
mkdir -p /scratch
cd /scratch
mkdir -p data
cd data

aws s3 ls s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a/agg/

aws s3 cp s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a/agg/ . --exclude "SC_RNA_AGGREGATOR_CS*" --recursive

# 350 MB
```

Call pubweb

```bash
mkdir -p /scratch/output

python /opt/pubweb/invoke-cellranger.py \
      --input '/scratch/data/outs' \
      --output '/scratch/output' \
      --name 'tds-2agg' \
      --species 'human'

```

Save the results

```bash

cd /scratch/output
aws s3 cp . s3://test-nextflow-data/1c0f2946-df95-49df-ad20-b852ce1cd57a/pubweb/  --recursive

# 2.8GB
```


