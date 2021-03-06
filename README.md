# gbrSnake
Snakemake pipeline for automating a parallel gbrs workflow on a large HPC node, https://github.com/churchill-lab/gbrs. 
Includes instructions for building a singularity image from a gbrs docker pull. 

gbrSnake takes the workflow and containerized environment of gbrs and gives it the power and flexibility of snakemake. If provided a folder of fastq files and a samples folder with a csv describing the inputs, gbrSnake produces the gbrs outputs from the emase compression to the diploid quantification step, automatically pooling technical replicates and assigning the correct transition probability file for the genome reconstruction step.  

![gbrSnake Sample DAG](https://github.com/exsquire/gbrSnake/blob/master/img/gbrSnake_img2.PNG)
**Main command:**
Runs pipeline using a local gbrs singularity image (--use-singularity), continues with independent files if an error occurs (-k), limits the number of parallel align_fastq rules (--resources load = number of parallel alignment jobs), and redirects the stdout and stderr to a master log file. (&> master.log).

```bash
conda activate snakemake
snakemake --use-singularity -k --resources load=80 &> master.log
```

**Building a gbrs singularity image:**
Root privilleges are impicit to working with docker and users of shared HPC resources are typically barred from these privilleges for the well-being the ecosystem. 

Recently, singularity (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459) has become the solution for containerization in HPC environments. The following instructions are a basic workflow for converting docker pulls to singularity images. 

```bash
#Install docker and test installation
sudo apt install docker.io
docker --version

#For ease-of-use, give yourself root access to docker commands
sudo groupadd docker
sudo gpasswd -a $USER docker
newgrp docker
docker ps #test your new privilleges

#Pull gbrs docker
docker pull kbchoi/gbrs
CID=$(docker run -dt kbchoi/gbrs)
echo $CID > CID.txt
docker cp CID.txt $CID:./
docker exec -it $CID /bin/bash

#Ensure image works as expected. Can optionally modify the image before committing it, but beware of what you might break in the process
ctrl+d

#Roll up a local docker registry and commit your image
docker run -d -p 5000:5000 --restart=always --name registry registry:2
docker rename $CID gbrsnake
docker commit $CID gbrsnake
docker tag gbrsnake localhost:5000/gbrsnake
docker push localhost:5000/gbrsnake

#Create a definition file
vi def
# Create def file for singularity
Bootstrap: docker
Registry: http://localhost:5000
Namespace:
From: gbrsnake:latest
:wq

#Install singularity - instructions are for Ubuntu 18.04
sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
sudo apt-get update
sudo apt-get install -y singularity-container
singularity --version

# Build singularity container
sudo SINGULARITY_NOHTTPS=1 singularity build gbrsnake.simg def

#Confirm changes 
singularity shell gbrsnake.simg

```

