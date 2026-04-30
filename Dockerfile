FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for OpenFOAM (without pulling the .deb through apt's redirect restriction)
RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    gnupg \
    g++ \
    libreadline-dev \
    flex \
    make \
    binutils-dev \
    libopenmpi-dev \
    libopenmpi3 \
    openmpi-bin \
    libxt-dev \
    zlib1g-dev \
    gnuplot \
    && rm -rf /var/lib/apt/lists/*

# Add OpenFOAM repo key and list so apt knows the package metadata,
# then download the .deb manually with curl (which follows https->http redirects)
# and install with dpkg
RUN curl -s https://dl.openfoam.org/gpg.key | gpg --dearmor -o /usr/share/keyrings/openfoam.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/openfoam.gpg] https://dl.openfoam.org/ubuntu jammy main" \
       > /etc/apt/sources.list.d/openfoam.list \
    && apt-get update \
    && DEB_URL=$(apt-cache show openfoam13 | grep ^Filename: | head -1 | awk '{print "https://dl.openfoam.org/ubuntu/" $2}') \
    && curl -fSL -o /tmp/openfoam13.deb "$DEB_URL" \
    && dpkg -i /tmp/openfoam13.deb || true \
    && apt-get install -f -y \
    && rm /tmp/openfoam13.deb \
    && rm -rf /var/lib/apt/lists/*

RUN echo ". /opt/openfoam13/etc/bashrc" >> /root/.bashrc

WORKDIR /case

ENTRYPOINT ["/bin/bash", "-c", "source /opt/openfoam13/etc/bashrc && \"$@\"", "--"]
CMD ["/bin/bash"]
