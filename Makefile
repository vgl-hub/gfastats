CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfastats
TEST_TARGET = gfastats-validate
GENERATE_TARGET = gfastats-generate-tests
RANDOM_FASTA_TARGET = gfastats-generate-random-fasta
BUILD_PATH = build/bin
SOURCE_PATH = src

LIBS += -lz
LDFLAGS= -pthread

link = $(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD_PATH)/$(TARGET) $(BUILD_PATH)/o/*.o $(LIBS)

$(TARGET): | $(BUILD_PATH) head validate regenerate random_fasta

head: | $(BUILD_PATH) main input output functions log struct bed gfa gfa-lines uid-generator stream-obj
	$(link)
	
link:
	$(link)

main: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c $(SOURCE_PATH)/$(TARGET).cpp -o $(BUILD_PATH)/o/$(TARGET).o
input: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfastats-input.cpp -o $(BUILD_PATH)/o/gfastats-input.o
output: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfastats-output.cpp -o $(BUILD_PATH)/o/gfastats-output.o
functions: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfastats-functions.cpp -o $(BUILD_PATH)/o/gfastats-functions.o
log: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfastats-log.cpp -o $(BUILD_PATH)/o/gfastats-log.o
struct: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfastats-struct.cpp -o $(BUILD_PATH)/o/gfastats-struct.o
bed: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/bed.cpp -o $(BUILD_PATH)/o/bed.o
gfa: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfa.cpp -o $(BUILD_PATH)/o/gfa.o
gfa-lines: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/gfa-lines.cpp -o $(BUILD_PATH)/o/gfa-lines.o
uid-generator: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/uid-generator.cpp -o $(BUILD_PATH)/o/uid-generator.o
stream-obj: $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -c include/stream-obj.cpp -o $(BUILD_PATH)/o/stream-obj.o

validate: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(GENERATE_TARGET) $(SOURCE_PATH)/$(GENERATE_TARGET).cpp $(LIBS)

random_fasta: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(RANDOM_FASTA_TARGET) $(SOURCE_PATH)/$(RANDOM_FASTA_TARGET).cpp $(LIBS)

$(BUILD_PATH):
	mkdir -p $@ $@/o
	
debug: CXXFLAGS += -DDEBUG -O0
debug: head

clean:
	$(RM) -r build
