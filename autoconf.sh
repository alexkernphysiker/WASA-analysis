#/bin/bash
(git submodule init;git submodule update)
MAGIC=$ROOTSORTERSYS/core/bin/sorter-config
WASA_CPP_FLAGS=`$MAGIC -cpp`
WASA_CXX_FLAGS=`$MAGIC -cxx`
WASA_LIBS=`$MAGIC -ld -libs-wasa -libs-wasa-ana`
echo ROOT_DIR=$ROOTSYS$'\n'WASA_DIR=$ROOTSORTERSYS$'\n'WASALIBS=$WASA_LIBS$'\n'WASACXX=$WASA_CXX_FLAGS$'\n'WASACPP=$WASA_CPP_FLAGS > settings.pri

