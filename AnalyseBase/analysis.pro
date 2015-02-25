include(../settings.pri)
include(../libs.pri)
QT             -= gui
TARGET          = analysis_routines
TEMPLATE        = lib
SOURCES        += analysis.cpp wrap.cpp
HEADERS        += analysis.h wrap.h
