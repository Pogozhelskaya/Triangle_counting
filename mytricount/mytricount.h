#pragma once

#include "../deps/GraphBLAS/Demo/Include/demos.h"

GrB_Info mytricount           // count # of triangles
(
    int64_t *p_ntri,        // # of trianagles
    const int method,       // 1 to 6, see above
    const GrB_Matrix A,     // adjacency matrix
    double t [2]            // t [0]: multiply time, t [1]: reduce time
);
