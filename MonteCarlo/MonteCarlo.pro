TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    singleparticle.cpp \
    manybody.cpp \
    lib.cpp
LIBS += -L/usr/local/lib/armadillo-4.400.1 -larmadillo -openmp
INCLUDEPATH += /usr/include
DEPENDPATH += /usr/include

HEADERS += \
    singleparticle.h \
    manybody.h \
    lib.h
