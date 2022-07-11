CXX = g++
INCLUDE_DIR = -I./include

CXXFLAGS = -g -Wall -std=gnu++14 -O3 $(INCLUDE_DIR)

TARGET = gfastats
TEST_TARGET = gfastats-validate
GENERATE_TARGET = gfastats-generate-tests
RANDOM_FASTA_TARGET = gfastats-generate-random-fasta
BUILD_PATH = build/bin
SOURCE_PATH = src

LIBS += -lz
LDFLAGS= -pthread

$(TARGET): | head validate regenerate random_fasta

head: | main input output functions log struct bed gfa gfa-lines uid-generator
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD_PATH)/$(TARGET) *.o $(LIBS)
	$(RM) -r *.o

main:
	$(CXX) $(CXXFLAGS) -c $(SOURCE_PATH)/$(TARGET).cpp
input:
	$(CXX) $(CXXFLAGS) -c include/gfastats-input.cpp
output:
	$(CXX) $(CXXFLAGS) -c include/gfastats-output.cpp
functions:
	$(CXX) $(CXXFLAGS) -c include/gfastats-functions.cpp
log:
	$(CXX) $(CXXFLAGS) -c include/gfastats-log.cpp
struct:
	$(CXX) $(CXXFLAGS) -c include/gfastats-struct.cpp
bed:
	$(CXX) $(CXXFLAGS) -c include/bed.cpp
gfa:
	$(CXX) $(CXXFLAGS) -c include/gfa.cpp
gfa-lines:
	$(CXX) $(CXXFLAGS) -c include/gfa-lines.cpp
uid-generator:
	$(CXX) $(CXXFLAGS) -c include/uid-generator.cpp

validate: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(GENERATE_TARGET) $(SOURCE_PATH)/$(GENERATE_TARGET).cpp $(LIBS)

random_fasta: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(RANDOM_FASTA_TARGET) $(SOURCE_PATH)/$(RANDOM_FASTA_TARGET).cpp $(LIBS)

$(BUILD_PATH):
	mkdir -p $@
	
debug: CXXFLAGS += -DDEBUG
debug: head

clean:
	$(RM) -r build
