CC = clang
CFLAGS = -O0 -Wall -Wpedantic -Wextra -std=c18 -ggdb
CINCLUDES = -I./third_party/
CLIBS = -lm

SRC = src
OBJ = obj

SRCS = $(wildcard $(SRC)/*.c)
OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))

BINDIR = bin
BIN = $(BINDIR)/main

all: $(BIN)

$(BIN): $(OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) $(CLIBS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) -c $< -o $@

clean:
	rm -rf $(BINDIR) $(OBJ)

$(OBJ):
	@mkdir -p $@

