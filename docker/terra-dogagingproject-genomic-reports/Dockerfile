FROM ubuntu
#FROM marketplace.gcr.io/google/ubuntu1804
WORKDIR /usr

#######################
# Prerequisites
#######################

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y --no-install-recommends \
    curl \
    unzip \
    nano \
    vim \
    sudo \
    wget \
    git \
    perl \
    autoconf \
    automake \
    make \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libncurses-dev \
    gfortran \
    r-base \
 && apt-get install -y locales \
 && rm -rf /var/lib/apt/lists/* \
 && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8

ENV LANG en_US.utf8

#######################
# Python
#######################

RUN apt update \
 && apt install -yq software-properties-common \
 && add-apt-repository ppa:deadsnakes/ppa \
 && apt install -yq python3.7-dev \
 && update-alternatives --install /usr/local/bin/python3 python3 /usr/bin/python3.7 1 \
 && apt-get -yq install python3-pip \
 && /usr/local/bin/python3 -m pip install --upgrade pip

#ENV PIP_USER=true

RUN pip3 -V \
 && pip3 install --upgrade pip \
 && pip3 install numpy \
 && pip3 install pandas \
 && pip3 install pandas-gbq \
 && pip3 install pandas-profiling \
 && pip3 install h5py \
 && pip3 install scipy \
 && pip3 install scikit-learn \
 && pip3 install scikit-allel \
 && pip3 install google-cloud-bigquery \
 && pip3 install google-api-core \
 && pip3 install google-cloud-bigquery-datatransfer \
 && pip3 install google-cloud-datastore \
 && pip3 install google-cloud-resource-manager \
 && pip3 install google-cloud-storage \
 && pip3 install gencove \
 && pip3 install firecloud

#######################
# Google Cloud SDK
#######################
# Downloading gcloud package
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz

# Installing the package
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh

# Adding the package path to local
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin

#######################
# Software
#######################

# install R
#  From r-base:latest
  RUN R -e 'install.packages("caret")'
  RUN R -e 'install.packages("randomForest")'
  RUN R -e 'install.packages("scales")'

# install htslib
 RUN cd /usr/bin \
  && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
  && tar -vxjf htslib-1.9.tar.bz2 \
  && cd htslib-1.9 \
  && ./configure \
  && make

# install samtools
  RUN cd /usr/bin \
   && wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
   && tar -vxjf samtools-1.11.tar.bz2 \
   && cd samtools-1.11 \
   && ./configure \
   && make

# install bcftools
 RUN cd /usr/bin \
  && wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
  && tar -vxjf bcftools-1.9.tar.bz2 \
  && cd bcftools-1.9 \
  && ./configure \
  && make

# install plink and plink2
RUN mkdir /usr/bin/plink \
 && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200921.zip \
 && unzip plink_linux_x86_64_20200921.zip -d /usr/bin/plink \
 && rm plink_linux_x86_64_20200921.zip

RUN mkdir /usr/bin/plink2 \
 && wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20211125.zip \
 && unzip plink2_linux_x86_64_20211125.zip -d /usr/bin/plink2 \
 && rm plink2_linux_x86_64_20211125.zip

# install EAGLE
RUN mkdir /usr/bin/eagle \
 && wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz \
 && tar -xvf Eagle_v2.4.1.tar.gz -C /usr/bin/eagle \
 && rm Eagle_v2.4.1.tar.gz \
 && ln -s /usr/bin/eagle/Eagle_v2.4.1/eagle /usr/bin/eagle

# install admixture
RUN mkdir /usr/bin/admixture \
 && wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz \
 && tar -xvf admixture_linux-1.3.0.tar.gz -C /usr/bin/admixture \
 && rm admixture_linux-1.3.0.tar.gz \
 && ln -s /usr/bin/admixture/dist/admixture_linux-1.3.0/admixture /usr/bin/admixture

# install rfmix
RUN cd /usr/bin \
 && git clone https://github.com/slowkoni/rfmix.git \
 && cd rfmix \
 && autoreconf --force --install \
 && ./configure \
 && make
