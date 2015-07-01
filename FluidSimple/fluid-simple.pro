TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp

HEADERS += \
    solver.h

LIBS += -L /usr/lib/x86_64-linux-gnu/libGLEW.so
LIBS += -L /usr/lib/x86_64-linux-gnu/libGLU.so
LIBS += -L /usr/lib/x86_64-linux-gnu/libglut.so
LIBS += -L /usr/lib/nvidia-346/libGL.so.1

QMAKE_CXXFLAGS_RELEASE += -O3
