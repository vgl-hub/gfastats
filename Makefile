CC = g++

CFLAGS  = -g -Wall -std=c++1y -O3 -I./include -I./zlib

TARGET = gfastats
BUILD_PATH = build
SOURCE_PATH = src

LIBS += -lz


$(TARGET): $(SOURCE_PATH)/$(TARGET).cpp
	mkdir -p $(BUILD_PATH)/bin
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/bin/$(TARGET) $(SOURCE_PATH)/$(TARGET).cpp $(LIBS)

clean:
	$(RM) -r build
