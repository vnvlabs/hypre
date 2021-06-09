

from vnv:latest



#Install Hypre VnV
RUN cd /vnv-mount/applications/hypre/src/ && mkdir -f build && cd build && cmake -DHYPRE_ENABLE_BIGINT=1 -DHYPRE_ENABLE_SHARED=1 -DHYPRE_ENABLE_VNV=1 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_INSTALL_PREFIX=~/software/hypre -DInjection_DIR=${VNV_DIR} .. && cmake .. && make && make install 
ENV HYPRE_DIR=~/software/hypre
   
