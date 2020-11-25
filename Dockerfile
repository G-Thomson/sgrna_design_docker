FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    build-essential \
    r-base \
    python3.8 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    mafft \
    libxml2-dev \
    libcairo2-dev \
    libsqlite-dev \
    libmariadbd-dev \
    libmariadbclient-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libsasl2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    && Rscript -e "install.packages('tidyverse', dependencies=TRUE)" \
    && Rscript -e "install.packages('splitstackshape', dependencies=TRUE)" \
    && ln -s /usr/bin/python3 /usr/bin/python \ 
    && ln -s /usr/bin/pip3 /usr/bin/pip \
    ## Clean up from R source install
    && cd / \
    && rm -rf /tmp/* \
    && apt-get remove --purge -y $BUILDDEPS \
    && apt-get autoremove -y \
    && apt-get autoclean -y \
    && rm -rf /var/lib/apt/lists/* 

RUN pip3 install pandas "biopython==1.76" gffutils regex

WORKDIR /home

ADD ./CRISPys/ $HOME/home/CRISPys
ADD ./Pipeline/ $HOME/home/Pipeline

ENV PATH="/home/CRISPys/:${PATH}"

ENTRYPOINT ["/bin/bash"]