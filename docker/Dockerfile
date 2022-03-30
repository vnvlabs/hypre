
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}:latest

#Install Hypre VnV
COPY ./ ./hypre 
RUN echo ${VNV_DIR}
RUN cd hypre/src && mkdir -p build && cd build && cmake -DHYPRE_ENABLE_BIGINT=1 -DHYPRE_ENABLE_SHARED=1 -DHYPRE_ENABLE_VNV=1 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_INSTALL_PREFIX=/root/software/hypre -DInjection_DIR=${VNV_DIR}/lib/cmake .. && cmake .. && make  && make install 
ENV HYPRE_DIR=/root/software/hypre

ARG COMMIT_HASH=unknown
LABEL vnv.version.hypre=${COMMIT_HASH}
ENV VNV_CONFIG=${VNV_CONFIG}:/hypre/.vnv-config.json

