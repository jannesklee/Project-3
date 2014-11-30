TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    singleparticle.cpp \
    manybody.cpp \
    lib.cpp
LIBS += -L/usr/local/Cellar/gcc49/4.9.2/lib/gcc/x86_64-apple-darwin14.0.0/4.9.2
LIBS += -L/usr/local/lib/armadillo-4.400.1 -larmadillo
INCLUDEPATH += /usr/include
DEPENDPATH += /usr/include

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

HEADERS += \
    singleparticle.h \
    manybody.h \
    lib.h
