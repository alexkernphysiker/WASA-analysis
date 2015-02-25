include(../settings.pri)
include(../libs.pri)
QMAKE_CXXFLAGS += -std=c++11
QT             -= gui
TARGET          = run_analysisjob
TEMPLATE        = app
SOURCES        += main.cpp analysisjob.cpp
HEADERS        += analysisjob.h
