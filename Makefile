CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfastats
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
RANDOM_FASTA_TARGET = generate-random-fasta
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

LIBS = -lz
LDFLAGS= -pthread

OBJS := main input input-agp input-filters output functions log struct bed gfa gfa-lines uid-generator stream-obj reads
BINS := $(addprefix $(BINDIR)/, $(OBJS))

head: $(INCLUDE)/threadpool.h $(INCLUDE)/global.h $(BINS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(LIBS)

all: head validate regenerate random_fasta

$(OBJS): %: $(SOURCE)/%.cpp $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp $(LIBS)

random_fasta: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(RANDOM_FASTA_TARGET) $(SOURCE)/$(RANDOM_FASTA_TARGET).cpp $(LIBS)

$(BUILD):
	-mkdir -p $@

$(BINDIR):
	-mkdir -p $@
	
debug: CXXFLAGS += -DDEBUG -O0
debug: head

clean:
	$(RM) -r build
