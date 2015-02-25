include(../settings.pri)
include(../libs.pri)
QMAKE_CXXFLAGS += -std=c++11
QT             -= gui
TARGET          = analysis_routines
TEMPLATE        = lib
SOURCES        += routines.cpp
HEADERS        += routines.h
