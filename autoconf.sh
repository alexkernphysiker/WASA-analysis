#/bin/bash
(git submodule init;git submodule update)
WASA_LIBS=`$ROOTSORTERSYS/core/bin/sorter-config -ld -libs-wasa -libs-wasa-ana`
echo ROOT_DIR=$ROOTSYS$'\n'WASA_DIR=$ROOTSORTERSYS$'\n'LIBS+=$WASA_LIBS > settings.pri
