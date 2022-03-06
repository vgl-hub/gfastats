CC = g++
INCLUDE_DIR = -I./include -I./zlib

CFLAGS += -g -Wall -std=gnu++11 -O3 $(INCLUDE_DIR)

TARGET = gfastats
TEST_TARGET = gfastats-validate
BUILD_PATH = build/bin
SOURCE_PATH = src

LIBS += -lz


$(TARGET): $(SOURCE_PATH)/$(TARGET).cpp
	mkdir -p $(BUILD_PATH)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(TARGET) $(SOURCE_PATH)/$(TARGET).cpp $(LIBS)
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/$(TEST_TARGET) $(SOURCE_PATH)/$(TEST_TARGET).cpp $(LIBS)

clean:
	$(RM) -r build
