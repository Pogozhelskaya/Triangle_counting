GRAPHBLAS=deps/GraphBLAS/build/libgraphblas.a

SOURCEDIR=$(shell pwd -P)
CC_SOURCES = $(wildcard $(SOURCEDIR)/*.c)
CC_SOURCES += $(wildcard $(SOURCEDIR)/timer/*.c)
CC_SOURCES += $(wildcard $(SOURCEDIR)/mytricount/*.c)

run: all

all: $(GRAPHBLAS) $(CC_SOURCES)
	gcc -o main ${CC_SOURCES} -fopenmp $(GRAPHBLAS) -lm

$(GRAPHBLAS):
ifeq (,$(wildcard $(GRAPHBLAS)))
	@$(MAKE) -C deps/GraphBLAS CMAKE_OPTIONS="-DCMAKE_C_COMPILER='gcc' -DCMAKE_CXX_COMPILER='g++'" static_only
endif
.PHONY: $(GRAPHBLAS)