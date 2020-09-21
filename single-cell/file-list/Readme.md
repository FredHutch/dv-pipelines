# File List workflow

This is a simple workflow that downloads the files from the source directory, lists them, and saves the list to a `filelist.txt` file, which is then saved to the target directory.

In doing so it tests file loading, publishing, and Nextflow execution. It's a quick sanity check.





## Example Setup (local)

```bash
cd /data/input
touch 1.txt 2.txt 3.foo 4.cat

mkdir -p /data/input/input
touch foo.bar foo.bar2 foo.bar3

```


## Run locally


Save this as ```cd /data/wf/rerun.sh```

```bash

TARGET_DIR=/data/output
SOURCE_DIR=/data/input
NEXTFLOW_PARAMS="--s3target $TARGET_DIR --s3source $SOURCE_DIR"
echo "Nextflow params are $NEXTFLOW_PARAMS"

nextflow run main.nf $NEXTFLOW_PARAMS --C local.config

```

Run the script

```bash
chmod +x rerun.sh

. rerun.sh

```


## Run on awscli

Copy in the `nextflow.config` file.

Save this as ```cd /data/wf/rerun.sh```

```bash

TARGET_DIR=s3://test-nextflow-data/5a30f3af-2d54-4c92-a0b7-dcbb94c145e7/filelist
SOURCE_DIR=s3://test-nextflow-data/5a30f3af-2d54-4c92-a0b7-dcbb94c145e7/hdf5/output
NEXTFLOW_PARAMS="--s3target $TARGET_DIR --s3source $SOURCE_DIR"
echo "Nextflow params are $NEXTFLOW_PARAMS"

nextflow run main.nf $NEXTFLOW_PARAMS --C nextflow.config

```

Run the script

```bash
chmod +x rerun.sh

. rerun.sh

```










