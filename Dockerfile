FROM ubuntu:20.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
ENV DEBIAN_FRONTEND=noninteractive
USER root
WORKDIR HOME
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git \
    ca-certificates \
    sudo 
RUN git clone https://github.com/yatisht/usher.git 
WORKDIR usher
## Checkout latest release
#RUN git checkout $(git describe --tags `git rev-list --tags --max-count=1`)
RUN ./installUbuntu.sh 
# faSomeRecords and faSize are needed for the UShER WDL workflow 
WORKDIR /HOME/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize
RUN chmod 775 *
## set the path
ENV PATH="/HOME/usher/build:/HOME/kentsource:${PATH}"
