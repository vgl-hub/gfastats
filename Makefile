CC = g++

CFLAGS  = -g -Wall -std=c++11 -O3 -I. -I./include

# The build target 
TARGET = fastats
BUILD_PATH = build
SOURCE_PATH = src

$(TARGET): $(SOURCE_PATH)/$(TARGET).cpp
	$(CC) $(CFLAGS) -o $(BUILD_PATH)/bin/$(TARGET) $(SOURCE_PATH)/$(TARGET).cpp

clean:
	$(RM) $(BUILD_PATH)/bin/$(TARGET)