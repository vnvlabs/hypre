
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}

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
   && make install \
   && cp -r ${SOURCE_DIR}/hypre/vnv/examples ${SOFTWARE_DIR}/hypre/examples  \
   && cp -r ${SOURCE_DIR}/hypre/vnv/gui ${SOFTWARE_DIR}/hypre/gui \
   && rm -r ${SOURCE_DIR}/hypre

ENV HYPRE_DIR=${SOFTWARE_DIR}/hypre

RUN ${VNV_REGISTER} gui hypre ${SOFTWARE_DIR}/hypre/gui/vnv.__registration__ 

