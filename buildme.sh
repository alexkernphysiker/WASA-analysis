#/bin/bash
(cd Analysis_test;qmake Analysis_test.pro;make clean; make; rootcint -f DictOutput.cxx -c -llibAnalyser.so analysisjob.hh LinkDef.hh)
