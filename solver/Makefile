# Courant EM solver core
CC ?= gcc
CFLAGS ?= -O3 -march=native -std=c11 -D_DEFAULT_SOURCE -Wall -Wextra -Wno-unused-parameter
VERSION := $(shell cat ../VERSION 2>/dev/null || echo dev)
CFLAGS += -DCEM_VERSION=\"$(VERSION)\"
INCV := -DCEM_VERSION=\"$(VERSION)\"
OMP ?= -fopenmp
BUILD := build
INC := -Iinclude
SRC := src/sim.c src/mode2d.c src/modeport.c src/dispersive.c src/sibc.c src/export.c
TESTS := test_cavity test_cpml test_sibc test_mode2d test_checkpoint
APPS := gcpw_line condloss2d

.PHONY: all test clean wasm parity
all: $(BUILD)/libcourant.a $(addprefix $(BUILD)/,$(TESTS) $(APPS))

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.c include/courant.h src/sim_internal.h | $(BUILD)
	$(CC) $(CFLAGS) $(OMP) $(INC) -c $< -o $@

$(BUILD)/libcourant.a: $(BUILD)/sim.o $(BUILD)/mode2d.o $(BUILD)/modeport.o $(BUILD)/dispersive.o $(BUILD)/sibc.o $(BUILD)/export.o
	ar rcs $@ $^

$(BUILD)/%: tests/%.c $(BUILD)/libcourant.a
	$(CC) $(CFLAGS) $(OMP) $(INC) $< -L$(BUILD) -lcourant -lm -o $@

$(BUILD)/%: apps/%.c $(BUILD)/libcourant.a
	$(CC) $(CFLAGS) $(OMP) $(INC) $< -L$(BUILD) -lcourant -lm -o $@

test: $(addprefix $(BUILD)/,$(TESTS))
	@set -e; for t in $(TESTS); do echo "== $$t"; ./$(BUILD)/$$t; done

# WebAssembly build with clang + wasi-libc (no Emscripten needed)
WASI_SYSROOT ?= /opt/wasi/wasi-sysroot-25.0
WASI_RT ?= /opt/wasi/libclang_rt.builtins-wasm32-wasi-25.0
wasm: | $(BUILD)
	clang --target=wasm32-wasip1 --sysroot=$(WASI_SYSROOT) -O3 -msimd128 -std=c11 $(INC) $(INCV) -D_DEFAULT_SOURCE \
	  -mexec-model=reactor -Wl,--export-dynamic -Wl,--export=malloc -Wl,--export=free \
	  -Wl,--max-memory=4294967296 -L$(WASI_RT) -lclang_rt.builtins-wasm32 \
	  $(SRC) -o $(BUILD)/courant.wasm
	clang --target=wasm32-wasip1-threads --sysroot=$(WASI_SYSROOT) -O3 -msimd128 -matomics -mbulk-memory -pthread \
	  -std=c11 $(INC) $(INCV) -D_DEFAULT_SOURCE \
	  -mexec-model=reactor -Wl,--export-dynamic -Wl,--export=malloc -Wl,--export=free -Wl,--export=__stack_pointer \
	  -Wl,--import-memory -Wl,--shared-memory -Wl,--max-memory=4294967296 -L$(WASI_RT) -lclang_rt.builtins-wasm32 \
	  $(SRC) -o $(BUILD)/courant-mt.wasm

clean:
	rm -rf $(BUILD)

# Native vs WebAssembly parity on a small scenario covering every kernel feature
parity: | $(BUILD)
	$(CC) -O3 -march=native -std=c11 -D_DEFAULT_SOURCE $(INC) apps/parity.c $(SRC) -lm -o $(BUILD)/parity
	./$(BUILD)/parity > $(BUILD)/parity_native.txt
	clang --target=wasm32-wasip1 --sysroot=$(WASI_SYSROOT) -O3 -msimd128 -std=c11 -D_DEFAULT_SOURCE $(INC) \
	  apps/parity.c $(SRC) -o $(BUILD)/parity.wasm
	node tests/run_wasi.mjs $(BUILD)/parity.wasm > $(BUILD)/parity_wasm.txt
	python3 tests/parity_check.py $(BUILD)/parity_native.txt $(BUILD)/parity_wasm.txt
