#/bin/bash
(cd AnalyseBase;qmake analysis.pro;make clean all)
if (( $? )); then exit 1; fi
make clean all
if (( $? )); then exit 1; fi
