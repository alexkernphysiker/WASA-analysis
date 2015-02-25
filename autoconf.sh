#/bin/bash
(git submodule init;git submodule update)
echo ROOT_DIR=$ROOTSYS$'\n'WASA_DIR=$ROOTSORTERSYS$'\n' > settings.pri
