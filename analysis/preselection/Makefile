CXX = g++
CPPFLAGS = -g -Wall -O2

ROOTFLAGS = $(shell root-config --libs --cflags --auxcflags) -I./include
ROOTFLAGS += -lRooFit -lRooFitCore -lRooStats -lMinuit -lFoam -lMathCore -lMathMore -lGenVector -lTreePlayer -lTMVA -lXMLIO  -lMLP -lm 

CORRECTIONLIBFLAGS = -I$(shell correction config --incdir --ldflags --rpath) 
BOOSTFLAGS = -I ${BOOSTDIR}/include

MAKEFLAGS := --jobs=8

EXE = bin/runAnalysis
BUILD_DIR = build
SRC_DIR = src
INCL_DIR = include

HEADERS = $(wildcard ${INCL_DIR}/*.h)
SOURCES = $(wildcard ${SRC_DIR}/*.cpp)
OBJECTS = $(SOURCES:${SRC_DIR}/%.cpp=${BUILD_DIR}/%.o)

all: makedir ${HEADERS} ${SOURCES} ${EXE}
	@echo "Built Successfully."

${EXE}: $(OBJECTS)
	@echo "Linking..."
	$(CXX) $(CPPFLAGS) $(ROOTFLAGS) $(CORRECTIONLIBFLAGS) ${BOOSTFLAGS} $(OBJECTS) -o $(EXE)

${BUILD_DIR}/%.o: ${SRC_DIR}/%.cpp  
	@echo "Compiling $<..."
	$(CXX) $(CPPFLAGS) $(ROOTFLAGS) $(CORRECTIONLIBFLAGS) ${BOOSTFLAGS} -c $< -o $@

makedir:
	@mkdir -p $(BUILD_DIR)

# delete the target of a rule if it has changed and its recipe exits with a nonzero exit status
.DELETE_ON_ERROR:

clean:
	@rm -r build/

cleanall:
	@rm -r $(OBJECTS) $(EXE) 
