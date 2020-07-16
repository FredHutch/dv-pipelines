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
cd /scratch/data/counts
aws s3 cp . s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data1/counts/ --recursive



cd /scratch/data/counts2
aws s3 cp . s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data2/counts/ --recursive

```







# Cellranger Aggr, 2020-07-01


Ubuntu, i3.2xlarge (8 CPUs & 61GB RAM needed, plus a ton of storage). Spot instance

```
ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-35-162-117-170.us-west-2.compute.amazonaws.com


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

Copy this to 'download.sh'
```
cd /scratch
## Get the data

mkdir -p data1
mkdir -p data2

cd /scratch/data1
aws s3 cp s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data1/ . --recursive

cd /scratch/data2
aws s3 cp s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data2/ . --recursive

```

```
nano #copy to download.sh
chmod +x download.sh
screen
. download.sh

Ctrl-a+d
```





## Make a cellranger aggr command


Do cellranger agg if absolutely necessary, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
   Look over Rob's code. Do I have to do that? 

Make a CSV file:


```
aws s3 ls s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data1/counts/ --profile sttrdevmfa

aws s3 ls s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data1/counts/CLJ_x424_IP_sampleGEX/ --profile sttrdevmfa

aws s3 ls s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data2/counts/ --profile sttrdevmfa

aws s3 ls s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/data2/counts/SSC-PBMC-4-2016/ --profile sttrdevmfa

```

```
library_id,molecule_h5
CLJ_x424_IP_sampleGEX,/scratch/data1/counts/CLJ_x424_IP_sampleGEX/outs/molecule_info.h5
CLJ_x424_d102_sampleGEX,/scratch/data1/counts/CLJ_x424_d102_sampleGEX/outs/molecule_info.h5
CLJ_x424_d12_sampleGEX,/scratch/data1/counts/CLJ_x424_d12_sampleGEX/outs/molecule_info.h5
CLJ_x424_d29_sampleGEX,/scratch/data1/counts/CLJ_x424_d29_sampleGEX/outs/molecule_info.h5
JAT_x483_IP_sampleGEX,/scratch/data1/counts/JAT_x483_IP_sampleGEX/outs/molecule_info.h5
JAT_x483_d12_sampleGEX,/scratch/data1/counts/JAT_x483_d12_sampleGEX/outs/molecule_info.h5
JAT_x483_d28_sampleGEX,/scratch/data1/counts/JAT_x483_d28_sampleGEX/outs/molecule_info.h5
REB_x254_IP_sampleGEX,/scratch/data1/counts/REB_x254_IP_sampleGEX/outs/molecule_info.h5
REB_x254_d112_sampleGEX,/scratch/data1/counts/REB_x254_d112_sampleGEX/outs/molecule_info.h5
REB_x254_d21_sampleGEX,/scratch/data1/counts/REB_x254_d21_sampleGEX/outs/molecule_info.h5
REB_x254_d38_sampleGEX,/scratch/data1/counts/REB_x254_d38_sampleGEX/outs/molecule_info.h5
TAB_x311_IP_sampleGEX,/scratch/data1/counts/TAB_x311_IP_sampleGEX/outs/molecule_info.h5
TAB_x311_d12_sampleGEX,/scratch/data1/counts/TAB_x311_d12_sampleGEX/outs/molecule_info.h5
TAB_x311_d29_sampleGEX,/scratch/data1/counts/TAB_x311_d29_sampleGEX/outs/molecule_info.h5
TAB_x311_d83_sampleGEX,/scratch/data1/counts/TAB_x311_d83_sampleGEX/outs/molecule_info.h5
SSC-PBMC-4-2016,/scratch/data2/counts/SSC-PBMC-4-2016/outs/molecule_info.h5
```

Setup for the call

```
mkdir -p /scratch
cd /scratch/agg
nano #copy in the aggregate CSV to 'aggregate.csv' in /scratch
```

```
cd /scratch
export PATH=/scratch/cellranger-3.1.0:$PATH
ID="run1"
CSV="/scratch/aggregate.csv"
cellranger aggr \
  --id=$ID \
  --csv=$CSV

cd run1
aws s3 cp . s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/agg/ --recursive

poweroff

```


```
nano #save the aggr script to 'aggregate.sh'
chmod +x aggregate.sh

screen

. aggregate.sh
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

aws s3 ls s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/agg/

aws s3 cp s3://test-nextflow-data/fa820f31-c24d-4dd1-8a06-957abe9d84eb/agg/ . --exclude "SC_RNA_AGGREGATOR_CS*" --recursive

# 354 MB
```

Call pubweb

```bash
mkdir -p /scratch/output

python /opt/pubweb/invoke-cellranger.py \
      --input '/scratch/data/outs' \
      --output '/scratch/output' \
      --name 'tds-3agg' \
      --species 'human'

```

Error encountered:

```python
Traceback (most recent call last):
  File "/opt/pubweb/invoke-cellranger.py", line 14, in <module>
    CellRanger(
  File "/opt/pubweb/pubweb/singlecell.py", line 49, in CellRanger
    Hdf5.load(outputFile, "a") | \
  File "/opt/pubweb/pubweb/commands/export/spatial.py", line 64, in __ror__
    hull = ConvexHull(points)
  File "qhull.pyx", line 2431, in scipy.spatial.qhull.ConvexHull.__init__
  File "qhull.pyx", line 356, in scipy.spatial.qhull._Qhull.__init__
scipy.spatial.qhull.QhullError: QH6214 qhull input error: not enough points(2) to construct initial simplex (need 3)

While executing:  | qhull i Qt
Options selected for Qhull 2019.1.r 2019/06/21:
  run-id 2026654714  incidence  Qtriangulate  _pre-merge  _zero-centrum
  _maxoutside  0
```

