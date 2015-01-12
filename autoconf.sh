#/bin/bash
(git submodule init;git submodule update)
echo ROOTDEPENDPATH = $ROOTSYS/lib/$'\n'ROOTINCLUDEPATH = $ROOTSYS/include/$'\n'ROOTSORTERPATH = $ROOTSORTERSYS/ > qt_build.pri
