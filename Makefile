CC = g++
INCLUDE_DIR = -I./include -I./include/zlib

CFLAGS += -g -Wall -std=gnu++11 -O3 $(INCLUDE_DIR)

TARGET = gfastats
GENERATE_TARGET = gfastats-generate-tests
TEST_TARGET = gfastats-validate
RANDOM_FASTA_TARGET = gfastats-generate-random-fasta
BUILD_PATH = build/bin
SOURCE_PATH = src

LIBS += -lz


$(TARGET): $(SOURCE_PATH)/$(TARGET).cpp
	mkdir -p $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(TARGET) $(SOURCE_PATH)/$(TARGET).cpp $(LIBS)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(GENERATE_TARGET) $(SOURCE_PATH)/$(GENERATE_TARGET).cpp $(LIBS)

validate:
	mkdir -p $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(GENERATE_TARGET) $(SOURCE_PATH)/$(GENERATE_TARGET).cpp $(LIBS)

random_fasta:
	mkdir -p $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(RANDOM_FASTA_TARGET) $(SOURCE_PATH)/$(RANDOM_FASTA_TARGET).cpp $(LIBS)

clean:
	$(RM) -r build
