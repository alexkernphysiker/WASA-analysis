include("../qt_build.pri")
QMAKE_CXXFLAGS+= -std=c++11
INCLUDEPATH += ..
HEADERS += LinkDef.hh analysisjob.hh
SOURCES += main.cc analysisjob.cc
