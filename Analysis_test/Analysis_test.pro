include("../qt_build.pri")
QMAKE_CXXFLAGS+= -std=c++11
INCLUDEPATH += $$ROOTINCLUDEPATH
DEPENDPATH += $$ROOTDEPENDPATH
LIBS += -lCore -lGui -lGpad -lGraf -lHist
INCLUDEPATH += $$ROOTSORTERPATH/core/include
DEPENDPATH += $$ROOTSORTERPATH/core/lib
LIBS += -lWasaRecFD -lWasaRecSE -lWasaRecPS -lWasaRecMDC -lWasaRecCD -lWasaParameter -lWasaAnaRaw -lWasaRecFPC
QT       += core gui
TARGET = Analyser
TEMPLATE = lib
HEADERS += LinkDef.hh analysisjob.hh
SOURCES += main.cc analysisjob.cc
