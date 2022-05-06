# Use an official Python runtime as a parent image
FROM python:3.7.4-slim

# Install freebayes dependencies
RUN apt-get update -yqq && \
    apt-get install -yqq \
        make build-essential libssl-dev zlib1g-dev libbz2-dev \
        libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
        xz-utils tk-dev libffi-dev liblzma-dev python-openssl git tabix && \
    apt-get clean

# Install freebayes
RUN \
    # Clone and compile freebayes
    git config --global url.https://github.com/.insteadOf git://github.com/ && \
    cd /opt && \
    git clone --recursive --branch v1.2.0 git://github.com/ekg/freebayes.git && \
    cd freebayes && \
    make && \
    make install && \
    \
    # Install htslib
    cd /tmp && \
    wget "https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2" && \
    tar -vxjf htslib-1.8.tar.bz2 && \
    cd htslib-1.8 && \
    ./configure --disable-bz2 --disable-lzma && \
    make && \
    make install && \
    rm -rf /tmp/htslib-1.8 && \
    \
    # Install bcftools
    cd /tmp && \
    wget "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2" && \
    tar -vxjf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    ./configure --disable-bz2 --disable-lzma && \
    make && \
    make install && \
    rm -rf /tmp/bcftools-1.9

ENV PATH $PATH:/opt/freebayes/scripts/:/opt/freebayes/vcflib/scripts/:/opt

# Install gsort
RUN \
    cd /usr/local/bin && \
    curl -s https://api.github.com/repos/brentp/gsort/releases/latest \
      | grep browser_download_url \
      | grep -i $(uname) \
      | cut -d '"' -f 4 \
      | wget -O gsort -qi - && \
    chmod +x gsort

# Install mity from dev server (but first install the previous version to get the dependencies from pypi)
#RUN pip install mitywgs==0.2.2
#RUN pip install -i https://test.pypi.org/simple/ mitywgs==0.2.2rc2

# Install mity
RUN pip install mitywgs==0.4.0

WORKDIR /home

ENTRYPOINT [ "mity" ]
