CFLAGS = -O0 -Wall -Wpedantic -Wextra -std=c18 -ggdb
CINCLUDES = -I./third_party/
ifeq ($(OS),Windows_NT)
	CC = gcc
	CFLAGS += -D__USE_MINGW_ANSI_STDIO
	CLIBS = -lcblas
else
	CC ?= clang
	CLIBS = -lm -lcblas
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

docs: physics/dispersion_propagation/dispersion.pdf

$(BIN): $(MAIN_OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) $^ -o $@ $(CLIBS)

$(TESTBIN): $(TEST_OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) $^ -o $@ $(CLIBS)

$(OBJ)/%.o: $(SRC)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CINCLUDES) -c $< -o $@

physics/dispersion_propagation/dispersion.pdf: physics/dispersion_propagation/dispersion.tex
	pdflatex -output-directory=physics/dispersion_propagation physics/dispersion_propagation/dispersion.tex
	pdflatex -output-directory=physics/dispersion_propagation physics/dispersion_propagation/dispersion.tex

clean:
	rm -rf $(BINDIR) $(OBJ) tests/*result*

$(OBJ):
	@mkdir -p $@

test: $(TESTBIN)
	$(TESTBIN)

run: $(BIN)
	$(BIN) ./lattices/m4U_f02020101_lattice.mad8 -p 20 -E 3.0

