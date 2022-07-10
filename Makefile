CC = g++
INCLUDE_DIR = -I./include

CFLAGS = -g -Wall -std=gnu++14 -O3 $(INCLUDE_DIR)

TARGET = gfastats
TEST_TARGET = gfastats-validate
GENERATE_TARGET = gfastats-generate-tests
RANDOM_FASTA_TARGET = gfastats-generate-random-fasta
BUILD_PATH = build/bin
SOURCE_PATH = src

LIBS += -lz
LDFLAGS= -pthread

$(TARGET): | $(BUILD_PATH) gf validate regenerate random_fasta

head: | $(BUILD_PATH) main input output functions log struct bed gfa gfa-lines uid-generator
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BUILD_PATH)/$(TARGET) *.o $(LIBS)
	$(RM) -r *.o

main:
	$(CC) $(CFLAGS) -c $(SOURCE_PATH)/$(TARGET).cpp
input:
	$(CC) $(CFLAGS) -c include/gfastats-input.cpp
output:
	$(CC) $(CFLAGS) -c include/gfastats-output.cpp
functions:
	$(CC) $(CFLAGS) -c include/gfastats-functions.cpp
log:
	$(CC) $(CFLAGS) -c include/gfastats-log.cpp
struct:
	$(CC) $(CFLAGS) -c include/gfastats-struct.cpp
bed:
	$(CC) $(CFLAGS) -c include/bed.cpp
gfa:
	$(CC) $(CFLAGS) -c include/gfa.cpp
gfa-lines:
	$(CC) $(CFLAGS) -c include/gfa-lines.cpp
uid-generator:
	$(CC) $(CFLAGS) -c include/uid-generator.cpp

validate: | $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(GENERATE_TARGET) $(SOURCE_PATH)/$(GENERATE_TARGET).cpp $(LIBS)

random_fasta: | $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(RANDOM_FASTA_TARGET) $(SOURCE_PATH)/$(RANDOM_FASTA_TARGET).cpp $(LIBS)

$(BUILD_PATH):
	mkdir -p $@

clean:
	$(RM) -r build
