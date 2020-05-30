#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "deps/GraphBLAS/Include/GraphBLAS.h"
#include "mytricount/mytricount.h"
#include "timer/simple_timer.h"

#define MAX_GRAPH_SIZE 2000000
#define MAX_ITEM_NAME_LEN 10

char GRAPH_INPUT_FILE[255];

char METHODS[][255] = {
    "Naive",
    "Burkhardt",
    "Cohen",
    "Sandia",
    "Sandia2",
    "SandiaDot",
    "SandiaDot2"
};

GrB_Info info; // Log of GraphBLAS operations

void load_graph(GrB_Matrix* graph, FILE* f) {
    char* line_buf;
    size_t buf_size = 0;

    while (getline(&line_buf, &buf_size, f) != -1) {
        line_buf[strcspn(line_buf, "\n")] = 0;

        char v[MAX_ITEM_NAME_LEN], to[MAX_ITEM_NAME_LEN];
        int nitems = sscanf(line_buf, "%s %s", v, to);
        assert(nitems == 2);

        uint32_t v_id = atoll(v);
        uint32_t to_id = atoll(to);

        info = GrB_Matrix_setElement(*graph, 1, v_id, to_id);
        assert(info == GrB_SUCCESS && "GraphBlas: failed to set matrix element (v, to)\n");

        info = GrB_Matrix_setElement(*graph, 1, to_id, v_id);  
        assert(info == GrB_SUCCESS && "GraphBlas: failed to set matrix element (to, v)\n");      
    }
}

int main(int argc, char* argv[]) {
    // Load input path from command line
    if (argc == 2) {
        strcpy(GRAPH_INPUT_FILE, argv[1]);
        printf("Input: %s\n", GRAPH_INPUT_FILE);
    }

    // Initialize GraphBLAS
    GrB_init(GrB_NONBLOCKING);

    // Create graph
    GrB_Matrix graph;
    info = GrB_Matrix_new( // create a new matrix with no entries
        &graph, // handle of matrix to create
        GrB_UINT32, // type of matrix to create
        MAX_GRAPH_SIZE, MAX_GRAPH_SIZE // matrix dimension is nrows-by-ncols
        );
    assert(info == GrB_SUCCESS && "GraphBlas: failed to construct matrix\n");

    // Open file
    FILE* f = fopen(GRAPH_INPUT_FILE, "r");
    assert(f != NULL);

    // Load graph
    load_graph(&graph, f);
    fclose(f);

    double timer[2];

    uint64_t res = 0;

    for (int i = 1; i <= 7; ++i) {
        mytricount(&res, i % 7, graph, timer);

        printf("%s number of triangles in graph = %ld\n", METHODS[i % 7], res);
        printf("%s used time (in seconds): %f\n\n", METHODS[i % 7], timer[0] + timer[1]);
        fflush(NULL);
    }

    return 0;
}