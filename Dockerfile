FROM ubuntu:18.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
USER root
WORKDIR HOME
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git \
    ca-certificates \
    sudo 
RUN git clone https://github.com/yatisht/usher.git 
WORKDIR usher
RUN ./installUbuntu.sh 
## set the path
ENV PATH="/HOME/usher/build/:${PATH}"
