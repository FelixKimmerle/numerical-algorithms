TARGET_EXEC ?= num

BUILD_DIR ?= ./build
SRC_DIRS ?= src
CXX ?= g++
MKDIR_P ?= mkdir -p

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -std=c++17 -MMD -MP -Wall -Wextra
LDFLAGS ?= -lm

#
# Debug build settings
#
DBGDIR := $(BUILD_DIR)/debug
DBGEXE := $(DBGDIR)/$(TARGET_EXEC)
DBGOBJS := $(SRCS:%=$(DBGDIR)/%.o)
DBGCFLAGS := $(CPPFLAGS) -g -O0 -ftrapv -DDEBUG

#
# Release build settings
# (no -ffast-math: the algorithms compare against exact zeros, e.g. when searching for a pivot)
#
RELDIR := $(BUILD_DIR)/release
RELEXE := $(RELDIR)/$(TARGET_EXEC)
RELOBJS := $(SRCS:%=$(RELDIR)/%.o)
RELCFLAGS := $(CPPFLAGS) -O2 -DRELEASE -DNDEBUG

DEPS := $(DBGOBJS:.o=.d) $(RELOBJS:.o=.d)

.PHONY: all clean release debug run rund remake test

all: release debug
remake: clean all

release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(MKDIR_P) $(dir $@)
	$(CXX) $(RELOBJS) -o $@ $(LDFLAGS)

$(RELDIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(RELCFLAGS) -c $< -o $@

debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(MKDIR_P) $(dir $@)
	$(CXX) $(DBGOBJS) -o $@ $(LDFLAGS)

$(DBGDIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(DBGCFLAGS) -c $< -o $@

clean:
	$(RM) -r $(BUILD_DIR)

run: release
	$(RELEXE)

rund: debug
	$(DBGEXE)

test: debug
	valgrind --leak-check=full $(DBGEXE)

-include $(DEPS)
