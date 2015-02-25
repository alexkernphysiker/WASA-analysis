include(../settings.pri)
include(../libs.pri)
QMAKE_CXXFLAGS += -std=c++11
QT             -= gui
TARGET          = analysisjob
TEMPLATE        = lib
SOURCES        += analysisjob.cpp
HEADERS        += analysisjob.h
