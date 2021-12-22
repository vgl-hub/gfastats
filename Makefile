CC = g++

CFLAGS  = -g -Wall -std=c++11 -O3 -I. -I./include -I./zlib

TARGET = fastats
BUILD_PATH = build
SOURCE_PATH = src

LIBS += -lz


$(TARGET): $(SOURCE_PATH)/$(TARGET).cpp
	mkdir -p $(BUILD_PATH)/bin
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/bin/$(TARGET) $(SOURCE_PATH)/$(TARGET).cpp $(LIBS)

clean:
	$(RM) build/bin