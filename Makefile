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

gf: | $(BUILD_PATH) main input functions log struct
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BUILD_PATH)/$(TARGET) *.o $(LIBS)
	$(RM) -r *.o

main:
	$(CC) $(CFLAGS) -c $(SOURCE_PATH)/$(TARGET).cpp
input:
	$(CC) $(CFLAGS) -c include/gfastats-input.cpp
functions:
	$(CC) $(CFLAGS) -c include/gfastats-functions.cpp
log:
	$(CC) $(CFLAGS) -c include/gfastats-log.cpp
struct:
	$(CC) $(CFLAGS) -c include/gfastats-struct.cpp

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
