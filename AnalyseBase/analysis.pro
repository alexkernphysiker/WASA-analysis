include(../settings.pri)
include(../libs.pri)
QT             -= gui
TARGET          = run_analysisjob
TEMPLATE        = app
SOURCES        += main.cpp analysisjob.cpp
HEADERS        += analysisjob.h
