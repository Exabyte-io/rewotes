FROM ubuntu:latest

RUN apt-get update \
  && apt-get install -y \
    build-essential \ 
    gfortran \
    wget \
    python3 \
    git

WORKDIR /app
RUN git clone https://github.com/QEF/q-e.git
RUN cd q-e && \
    git checkout 6eaef358375f923a38805e5c8560088df1e9ab57 && \
    ./configure && \
    make pw
