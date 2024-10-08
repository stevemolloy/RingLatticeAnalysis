CFLAGS = -O0 -Wall -Wpedantic -Wextra -std=c18 -ggdb
CINCLUDES = -I./third_party/
ifeq ($(OS),Windows_NT)
	CC = gcc
	CLIBS = 
else
	CC ?= clang
	CLIBS = -lm
endif

SRC = src
OBJ = obj

SRCS = $(wildcard $(SRC)/*.c)
OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))

TEST_SRCS = $(filter-out $(SRC)/main.c, $(SRCS))
TEST_OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(TEST_SRCS))

MAIN_SRCS = $(filter-out $(SRC)/test.c, $(SRCS))
MAIN_OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(MAIN_SRCS))

BINDIR = bin
BIN = $(BINDIR)/rla
TESTBIN = $(BINDIR)/test

all: $(BIN) $(TESTBIN)

$(BIN): $(MAIN_OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) $^ -o $@ $(CLIBS)

$(TESTBIN): $(TEST_OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) $^ -o $@ $(CLIBS)

$(OBJ)/%.o: $(SRC)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) -c $< -o $@

clean:
	rm -rf $(BINDIR) $(OBJ)

$(OBJ):
	@mkdir -p $@

test: $(TESTBIN)
	$(TESTBIN)

