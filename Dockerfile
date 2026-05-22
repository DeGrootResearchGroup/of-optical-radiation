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
RUN set -eux \
    && curl -s https://dl.openfoam.org/gpg.key | gpg --dearmor -o /usr/share/keyrings/openfoam.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/openfoam.gpg] https://dl.openfoam.org/ubuntu jammy main" \
       > /etc/apt/sources.list.d/openfoam.list \
    && apt-get update \
    && DEB_URL=$(apt-cache show openfoam13 | grep ^Filename: | head -1 | awk '{print "https://dl.openfoam.org/ubuntu/" $2}') \
    && curl -fSL -o /tmp/openfoam13.deb "$DEB_URL" \
    # dpkg -i can return non-zero on unmet deps; the "apt-get install -f -y"
    # line below resolves them. We don't mask the failure with "|| true"
    # any more -- if dpkg fails for a reason apt-get -f can't fix
    # (corrupt .deb, missing repo, etc.), fail the build now rather
    # than letting a broken image masquerade as healthy until a user
    # hits a missing-symbol error inside a container.
    && (dpkg -i /tmp/openfoam13.deb || apt-get install -f -y) \
    && dpkg -l openfoam13 | grep -q '^ii' \
    && rm /tmp/openfoam13.deb \
    && rm -rf /var/lib/apt/lists/*

RUN echo ". /opt/openfoam13/etc/bashrc" >> /root/.bashrc

# Mesh-tooling layer: python3-gmsh (apt; PyPI gmsh wheels are x86_64-only),
# pip + git for installing blockmeshbuilder from upstream, and the in-tree
# `uvmesh` helper. blockmeshbuilder ships only via git -- there is no PyPI
# release. uvmesh is installed in non-editable mode from /code (bind-
# mounted at runtime); the install step here would fail before the source
# is present, so we install the *dependencies* here and `pip install
# /code/tools/uvMesh` is left to per-case Allruns (or CI's outer wrapper).
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-gmsh \
    git \
    && rm -rf /var/lib/apt/lists/*
RUN pip3 install --no-cache-dir \
        numpy \
        git+https://github.com/NauticalMile64/blockmeshbuilder.git

WORKDIR /case

ENTRYPOINT ["/bin/bash", "-c", "source /opt/openfoam13/etc/bashrc && \"$@\"", "--"]
CMD ["/bin/bash"]
