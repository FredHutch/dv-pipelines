# Pubweb + NF


Copy the code up to S3 in development, to the ```dv-code-dev``` bucket:

```bash
cd GitHub/data-vizualization-center/hdf5 # wherever this is in your computer

aws s3 rm s3://dv-code-dev/pubweb --recursive
aws s3 cp pubweb/ s3://dv-code-dev/pubweb/ --recursive

```


## Start up an EC2 instance

Start up a VM. I've been using Ubuntu 18.04 as the OS since it's stable. Any machine size will do (I did a t3.small)

```bash

ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-44-233-150-38.us-west-2.compute.amazonaws.com

```


### Install Docker

```bash

apt-get -y install docker.io awscli

```

### Make the container

copy in the Dockerfile

```bash
cd $HOME
mkdir -p container
cd container
nano #copy in the Dockerfile
nano #copy in startwhile.sh
```

Build the container

```bash
docker build -t hdf5 .
docker images | grep hdf5

screen

docker run hdf5

screen -a+D

docker exec -ti  0fbcd69b9b1b /bin/bash
```

Create the source data

```bash
cd $HOME
mkdir -p source
cd source


echo '1' > 1.txt
echo '2' > 2.txt
echo 'foo' > foo.bar

cd $HOME
cd source
tar -czvf ../input.tar.gz .

# test the decompress
cd $HOME
tar -xzf input.tar.gz -C source2


# copy the results
aws s3 cp input.tar.gz s3://test-nextflow-data/e202b949-283e-435e-a14c-35293517760b/20200827/input/  --dryrun
```

Copy in the nextflow script

```bash
cd $HOME
mkdir -p wf
cd wf
nano #main.nf
nano #pubwebconf.json

```

Save this as ```cd $HOME/wf/rerun.sh```

```bash

TARGET_DIR=s3://test-nextflow-data/e202b949-283e-435e-a14c-35293517760b/20200827/
SOURCE_DIR=s3://test-nextflow-data/e202b949-283e-435e-a14c-35293517760b/20200827/input
NEXTFLOW_PARAMS="--s3target $TARGET_DIR --s3source $SOURCE_DIR"
echo "Nextflow params are $NEXTFLOW_PARAMS"

nextflow run main.nf $NEXTFLOW_PARAMS --wfconfig pubwebfiles.json --C nextflow.config

```

Run the script

```bash
chmod +x rerun.sh


. rerun.sh

```













