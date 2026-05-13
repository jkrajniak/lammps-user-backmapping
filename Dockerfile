# syntax=docker/dockerfile:1

ARG LAMMPS_VERSION=stable_22Jul2025_update2

# =============================================================================
# Stage 1: Build LAMMPS with the backmapping package
# =============================================================================
FROM ubuntu:22.04 AS builder

ARG LAMMPS_VERSION
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates \
        libopenmpi-dev openmpi-bin \
        python3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN git clone --depth 1 --branch "${LAMMPS_VERSION}" \
        https://github.com/lammps/lammps.git

COPY src/ /build/backmap-pkg/src/
COPY install.sh /build/backmap-pkg/install.sh

RUN cd /build/backmap-pkg && ./install.sh /build/lammps

RUN mkdir -p /build/lammps/build && cd /build/lammps/build && \
    cmake ../cmake \
        -D PKG_BACKMAP=yes \
        -D PKG_MOLECULE=yes \
        -D BUILD_MPI=on \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=/usr/local \
        -D BUILD_SHARED_LIBS=off && \
    cmake --build . -j "$(nproc)" && \
    cmake --install .

# =============================================================================
# Stage 2: Minimal runtime image
# =============================================================================
FROM ubuntu:22.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi3 libgomp1 \
        python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/lmp /usr/local/bin/lmp

COPY python/ /tmp/backmap-prep/
RUN pip3 install --no-cache-dir /tmp/backmap-prep/ && rm -rf /tmp/backmap-prep/

COPY scripts/run-backmap.sh /usr/local/bin/run-backmap.sh

WORKDIR /work

ENTRYPOINT []
CMD ["lmp", "-h"]
