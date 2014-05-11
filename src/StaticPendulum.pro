#-------------------------------------------------
#
# Project created by QtCreator 2014-04-20T02:18:10
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += xml
TARGET = StaticPendulum
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++1y
QMAKE_CXXFLAGS += -pthread
QMAKE_CXXFLAGS_RELEASE += -ffast-math
QMAKE_CXXFLAGS_RELEASE += -march=native
QMAKE_CXXFLAGS_RELEASE += -funroll-loops
SOURCES += main.cpp

HEADERS += \
    pendulum_system.h \
    pendulum_map.h \
    Integrators/ck45.h \
    Integrators/rk4.h

LIBS += -pthread
