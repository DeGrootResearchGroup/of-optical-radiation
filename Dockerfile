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

# Mesh-tooling layer: gmsh built from source, pip + git for installing
# blockmeshbuilder from upstream, and the in-tree `uvmesh` helper.
# blockmeshbuilder ships only via git -- there is no PyPI release. uvmesh
# is installed in non-editable mode from /code (bind-mounted at runtime);
# the install step here would fail before the source is present, so we
# install the *dependencies* here and `pip install /code/tools/uvMesh` is
# left to per-case Allruns (or CI's outer wrapper).
#
# gmsh is built from source rather than taken from apt (Ubuntu's
# python3-gmsh is 4.8, built without Netgen) or PyPI (whose Linux wheels
# are x86_64 only): uvmesh's `ReactorBody(optimize_netgen=True)` runs
# gmsh's Netgen tet optimizer, which gmsh bundles in its source. Built
# without its GUI, as a shared library with the Python API, against
# Ubuntu's OpenCASCADE.
ARG GMSH_VERSION=4.15.2
ARG GMSH_BUILD_JOBS=4
RUN apt-get update && apt-get install -y \
    python3-pip \
    git \
    cmake \
    libocct-foundation-dev \
    libocct-modeling-data-dev \
    libocct-modeling-algorithms-dev \
    libocct-data-exchange-dev \
    libocct-ocaf-dev \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fSL -o /tmp/gmsh.tgz \
       "https://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz" \
    && tar -xzf /tmp/gmsh.tgz -C /tmp \
    && cmake -S /tmp/gmsh-${GMSH_VERSION}-source -B /tmp/gmsh-build \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX=/usr/local \
       -DENABLE_BUILD_DYNAMIC=ON \
       -DENABLE_FLTK=OFF \
       -DENABLE_NETGEN=ON \
       -DENABLE_OCC=ON \
    && cmake --build /tmp/gmsh-build -j "${GMSH_BUILD_JOBS}" \
    && cmake --install /tmp/gmsh-build \
    && rm -rf /tmp/gmsh.tgz /tmp/gmsh-${GMSH_VERSION}-source /tmp/gmsh-build \
    && ldconfig \
    && python3 -c "import sys; sys.path.insert(0, '/usr/local/lib'); import gmsh; \
gmsh.initialize(); opts = gmsh.option.getString('General.BuildOptions'); \
assert 'Netgen' in opts and 'OpenCASCADE' in opts, opts; gmsh.finalize()"
# The Python API is installed beside the library, not in site-packages.
ENV PYTHONPATH=/usr/local/lib
RUN pip3 install --no-cache-dir \
        numpy \
        pytest \
        git+https://github.com/NauticalMile64/blockmeshbuilder.git

WORKDIR /case

ENTRYPOINT ["/bin/bash", "-c", "source /opt/openfoam13/etc/bashrc && \"$@\"", "--"]
CMD ["/bin/bash"]
