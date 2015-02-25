#/bin/bash
(git submodule init;git submodule update)
echo ROOT_DIR=$ROOTSYS$'\n'WASA_DIR=$ROOTSORTERSYS$'\n'LIBS+="$($ROOTSORTERSYS/core/bin/sorter-config -ld -libs-wasa -libs-wasa-ana)" > settings.pri
