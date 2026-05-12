FROM ubuntu:20.04

LABEL maintainer="wzhou@broadinstitute.org"
LABEL org.opencontainers.image.title="SAIGE-RareEffect"
LABEL org.opencontainers.image.description="SAIGE 1.5.1 with RareEffect for rare-variant effect-size estimation and PRS"
LABEL org.opencontainers.image.version="1.5.1"
LABEL org.opencontainers.image.source="https://github.com/namks/SAIGE_RareEffect"
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

RUN apt-get update -qq && \
    apt-get -y install curl time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && \
    pixi install && \
    rm -rf /root/.cache && \
    pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'


RUN curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx && \
    mv plink-ng-2.0.0-a.6.16 plink-ng && \
    pixi run x86_64-conda-linux-gnu-cc -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o libplink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate && \
    mv libplink2_includes.a .pixi/envs/default/lib

RUN pixi run R CMD INSTALL .

RUN mv extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/ && \
   chmod a+x /usr/local/bin/step1_fitNULLGLMM.R /usr/local/bin/step2_SPAtests.R /usr/local/bin/step3_LDmat.R /usr/local/bin/createSparseGRM.R && \
   for f in run_RareEffect.R calculate_RareEffect_PRS.R; do \
     mv extdata/$f /usr/local/bin/$f && \
     sed -i '1i #!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript' /usr/local/bin/$f && \
     chmod a+x /usr/local/bin/$f; \
   done && \
   createSparseGRM.R --help && \
   step1_fitNULLGLMM.R --help && \
   step2_SPAtests.R --help && \
   step3_LDmat.R --help && \
   run_RareEffect.R --help && \
   calculate_RareEffect_PRS.R --help
