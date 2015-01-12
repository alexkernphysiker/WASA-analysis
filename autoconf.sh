#/bin/bash
(git submodule init;git submodule update)
echo ROOTDEPENDPATH = /usr/lib64/root/$'\n'ROOTINCLUDEPATH = /usr/include/root/$'\n'ROOTSORTERPATH = ~/software/RootSorter/ > qt_build.pri
