QMAKE_CXXFLAGS += -Werror
HEADERS       = window.h \
                method.h \
		functions.h

SOURCES       = main.cpp \
                method.cpp \
                window.cpp
QT += opengl
