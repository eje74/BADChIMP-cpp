# Use this make file by writing:
#
# $ make -f Makefile.terminal
#
# From : Hilton Lipschitz
#

CC := g++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
TARGET := bin/runner

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
#CFLAGS := -std=c++11 -O3 -g -pg # -Wall
CFLAGS := -std=c++11 -O3
LIB := -L lib
INC := -I include

$(TARGET): $(OBJECTS)
	   @echo " Linking..."
	   @echo " $(CC) $(CFLAGS) $^ -o $(TARGET) $(LIB)"; $(CC) $(CFLAGS) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
		 @mkdir -p $(BUILDDIR)
		 @echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester
# Snippets
snippet:
	$(CC) $(CFLAGS) snippets/snippet.cpp $(INC) $(LIB) -o bin/snippet

.PHONY: clean