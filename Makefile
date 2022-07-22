CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfastats
TEST_TARGET = gfastats-validate
GENERATE_TARGET = gfastats-generate-tests
RANDOM_FASTA_TARGET = gfastats-generate-random-fasta
BUILD_PATH = build/bin
SOURCE_PATH = src
BINDIR := $(BUILD_PATH)/.o

LIBS += -lz
LDFLAGS= -pthread

OBJS := ${TARGET} input input-agp input-filters output functions log struct bed gfa gfa-lines uid-generator stream-obj reads
BINS := $(addprefix $(BINDIR)/, $(OBJS))

head: $(BINS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD_PATH)/$(TARGET) $(wildcard $(BINDIR)/*) $(LIBS)

all: head validate regenerate random_fasta

%: include/%.cpp $(BINDIR)/%
	
$(BINDIR)/%: include/%.cpp include/%.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c include/$(notdir $@).cpp -o $@

$(BINS): | $(BINDIR)

$(BUILD_PATH)/.o/$(TARGET): include/$(TARGET).h include/threadpool.h include/global.h
	$(CXX) $(CXXFLAGS) -c $(SOURCE_PATH)/$(TARGET).cpp -o $@

validate: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(GENERATE_TARGET) $(SOURCE_PATH)/$(GENERATE_TARGET).cpp $(LIBS)

random_fasta: | $(BUILD_PATH)
	$(CXX) $(CXXFLAGS) -o $(BUILD_PATH)/$(RANDOM_FASTA_TARGET) $(SOURCE_PATH)/$(RANDOM_FASTA_TARGET).cpp $(LIBS)

$(BUILD_PATH):
	-mkdir -p $@

$(BINDIR):
	-mkdir -p $@
	
debug: CXXFLAGS += -DDEBUG -O0
debug: head

clean:
	$(RM) -r build
