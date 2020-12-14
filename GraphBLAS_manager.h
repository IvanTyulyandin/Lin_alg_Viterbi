#pragma once

extern "C" {
    #include "GraphBLAS.h"
}


void check_for_error(const GrB_Info& info);

// Init GraphBLAS
// Should be done exactly once
void launch_GraphBLAS();


// Finalize GraphBLAS
// Should be done exactly once
void stop_GraphBLAS();
