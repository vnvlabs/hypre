
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}

ARG COMMIT_HASH=unknown
LABEL vnv.version.hypre=${COMMIT_HASH}

#Install Hypre VnV
COPY ./ ${SOURCE_DIR}/hypre 

RUN cd ${SOURCE_DIR}/hypre/src \
    && mkdir -p build \
    && cd build \
    && cmake -DHYPRE_ENABLE_BIGINT=1\
             -DHYPRE_ENABLE_SHARED=1\
             -DHYPRE_ENABLE_VNV=1\
             -DHYPRE_BUILD_EXAMPLES=1\
	     -DCMAKE_EXPORT_COMPILE_COMMANDS=ON\
	     -DCMAKE_INSTALL_PREFIX=${SOFTWARE_DIR}/hypre\
             -DInjection_DIR=${VNV_DIR}/lib/cmake .. \
   && cmake .. \
   && make \
   && make install 

ENV HYPRE_DIR=${SOFTWARE_DIR}/hypre

ENV VNV_CONFIG=${VNV_CONFIG}:${SOURCE_DIR}/hypre/src/build/vnv.__registration__
ENV VNV_CONFIG=${VNV_CONFIG}:${SOURCE_DIR}/hypre/docker/config/vnv.__registration__

