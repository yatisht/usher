FROM ubuntu:20.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git wget \
    ca-certificates \
    sudo python3 python3-pip
WORKDIR /HOME/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize
RUN chmod 775 *
WORKDIR /HOME
RUN git clone https://github.com/yatisht/usher.git 
WORKDIR usher
## Checkout latest release
#RUN git checkout $(git describe --tags `git rev-list --tags --max-count=1`)
RUN ./install/installUbuntu.sh 
# faSomeRecords and faSize are needed for the UShER WDL workflow 
## set the path
ENV PATH="/HOME/usher/build:/HOME/kentsource:${PATH}"
