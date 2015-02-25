INCLUDEPATH += $$ROOT_DIR/include
INCLUDEPATH += $$WASA_DIR/core/include
INCLUDEPATH += $$WASA_DIR/Vt++/include
INCLUDEPATH += $$WASA_DIR/wasa/include
LIBS           += $$WASALIBS -lWasaRecFD -lWasaRecSE -lWasaRecPS -lWasaRecMDC -lWasaRecCD -lWasaParameter -lWasaAnaRaw -lWasaRecFPC -lvt
QMAKE_CXXFLAGS += $$WASACXX -fno-strict-aliasing -std=c++11
QMAKE_CPPFLAGS += $$WASACPP
