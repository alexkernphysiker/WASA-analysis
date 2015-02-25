#/bin/bash
(cd AnalyseBase;qmake analysis.pro;make)
if (( $? )); then exit 1; fi
