#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "deps/GraphBLAS/Include/GraphBLAS.h"
#include "timer/simple_timer.h"

#define MAX_GRAPH_SIZE 1000
#define MAX_ITEM_NAME_LEN 100

#define DEBUG 0

const char* GRAPH_INPUT_FILE = "input/one_triangle.txt";

GrB_Info info; // Log of GraphBLAS operations

void load_graph(GrB_Matrix* graph, FILE* f) {
    char* line_buf;
    size_t buf_size = 0;

    while (getline(&line_buf, &buf_size, f) != -1) {
        line_buf[strcspn(line_buf, "\n")] = 0;

        char v[MAX_ITEM_NAME_LEN], edge[MAX_ITEM_NAME_LEN], to[MAX_ITEM_NAME_LEN];
        int nitems = sscanf(line_buf, "%s %s", v, to);

        uint64_t v_id = atoll(v) - 1;
        uint64_t to_id = atoll(to) - 1;

        info = GrB_Matrix_setElement(*graph, 1, v_id, to_id);
        assert(info == GrB_SUCCESS && "GraphBlas: failed to set matrix element (v, to)\n");

        info = GrB_Matrix_setElement(*graph, 1, to_id, v_id);  
        assert(info == GrB_SUCCESS && "GraphBlas: failed to set matrix element (to, v)\n");      
    }
}

int main(int argc, char* argv[]) {
    // Load input path from command line
    if (argc == 2) {
        GRAPH_INPUT_FILE = argv[1];
        printf("Input: %s", GRAPH_INPUT_FILE);
    }

    // Initialize GraphBLAS
    GrB_init(GrB_NONBLOCKING);

    // Create graph
    GrB_Matrix graph;
    info = GrB_Matrix_new( // create a new matrix with no entries
        &graph, // handle of matrix to create
        GrB_UINT64, // type of matrix to create
        MAX_GRAPH_SIZE, MAX_GRAPH_SIZE // matrix dimension is nrows-by-ncols
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to construct matrix\n");

    // Open file
    FILE* f = fopen(GRAPH_INPUT_FILE, "r");
    assert(f != NULL);

    // Load graph
    load_graph(&graph, f);
    fclose(f);

    // Create temp matrix for triangle calculation
    GrB_Matrix tri;
    info = GrB_Matrix_dup( // make an exact copy of a matrix
        &tri, // handle of output matrix to create
        graph // input matrix to copy
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to construct temp matrix\n");

    // Create monoid for UINT64 Matrix PLUS
    GrB_Monoid monoid;
    info = GrB_Monoid_new_UINT64( // create a monoid
        &monoid, // handle of monoid to create
        GrB_PLUS_UINT64, // binary operator of the monoid
        0 // identity value of the monoid
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to construct the monoid\n");

    // Create semiring for UINT64 Matrix MULT
    GrB_Semiring semiring;
    info = GrB_Semiring_new( // create a semiring
        &semiring, // handle of semiring to create
        monoid, // add monoid of the semiring
        GrB_TIMES_UINT64 // multiply operator of the semiring
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to construct the semiring\n");

    double timer[2];
    simple_tic(timer);

    for (int i = 0; i < 2; ++i) {        
        info = GrB_mxm( // C<Mask> = accum (C, A*B)
            tri, // input/output matrix for results
            GrB_NULL, // optional mask for C, unused if NULL
            GrB_NULL, // optional accum for Z=accum(C,T)
            semiring, // defines ’+’ and ’*’ for A*B
            tri, // first input: matrix A
            graph, // second input: matrix B
            GrB_NULL // descriptor for C, Mask, A, and B
            );
        assert(info == GrB_SUCCESS && "GraphBlas: failed to multiply matrices\n");
    }

    // Create result
    GrB_Matrix main_diagonal;
    GrB_Matrix_new(&main_diagonal, GrB_UINT64, MAX_GRAPH_SIZE, MAX_GRAPH_SIZE);
    info = GxB_select( // C<Mask> = accum (C, op(A,k)) or op(A’,k)
        main_diagonal, // input/output matrix for results
        GrB_NULL, // optional mask for C, unused if NULL
        GrB_PLUS_UINT64, // optional accum for Z=accum(C,T)
        GxB_DIAG, // operator to apply to the entries
        tri, // first input: matrix A
        GrB_NULL, // optional input for the select operator
        GrB_NULL // descriptor for C, mask, and A
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to select elements from matrix\n");

    uint64_t res = 0;
    
    info = GrB_Matrix_reduce_UINT64( // c = accum (c, reduce_to_scalar (A))
        &res, // result scalar
        GrB_PLUS_UINT64, // optional accum for c=accum(c,t)
        monoid, // monoid to do the reduction
        main_diagonal, // matrix to reduce
        GrB_NULL // descriptor (currently unused)
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to reduce matrix\n");

#ifdef DEBUG
    GxB_print(graph, GxB_COMPLETE);
    GxB_print(tri, GxB_COMPLETE);
    GxB_print(main_diagonal, GxB_COMPLETE);
#endif

    double time_query = simple_toc(timer);

    printf("Number of triangles in graph = %ld\n", res / 6);
    printf("Used time (in seconds): %f\n", time_query);

    return 0;
}