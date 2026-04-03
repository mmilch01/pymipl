# XNAT AI workflow demo -- full tutorial

This tutorial illustrates the workflow to run lung mask segmentation on structural CT scans in an XNAT project. 
It is based on https://github.com/JoHof/lungmask codebase.

First, we assume that the XNAT workflow is installed into $PYMIPL_DIR:
```
export PYMIPL_DIR="~/src/pymipl"
```
We start by making environment setup dir on local machine: 
```
mkdir -p adapt_demo
cd adapt_demo
ADAPT_DEMO_DIR=`pwd`
```
Then (optionally) clone repository to this dir:
```
git clone https://github.com/JoHof/lungmask.git
``
The repository is saved under ./lungmask dir. Now we need to create the matching python environment, and install the user package.

```
#Create the environment (few seconds)
micromamba create -p ./lungmask-env python=3.12.3
micromamba activate -p ./lungmask-env

#Install the user codebase (in this case, they supply a Python package, so it's automatic)
pip install lungmask

```

Now we test the local installation on 'DICOM' directory 
```
cd <TEST_DATA_DIR> 
dcm2nii -o ./ -f test_in.nii DICOM
time lungmask DICOM test_mask.nii.gz
fsleyes test_in.nii test_mask.nii.gz

```
Then we need to create the dedicated Docker image configuration file:
```
cd $ADAPT_DEMO_DIR
cp $PYMIPL_DIR/xnat_env_bootstrap/example-custom-image.conf lungmask-image.conf
gedit lungmask-image.conf
```
If using docker.io, the only part we need to set is the local environment location:
```
REMOTE_IMAGE=mmilch01/xnat-lungmask:latest
LOCAL_ENV_FOLDER="/home/mmilchenko/src/adapt_demo/lungmask-env"
```

Now ready to build the Docker image:
```
$PYMIPL_DIR/build_deploy_custom_image.sh ./lungmask-image.conf
```

