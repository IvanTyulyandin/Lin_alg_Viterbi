#include "GraphBLAS_manager.h"

#include <iostream>

void check_for_error([[maybe_unused]] const GrB_Info& info) {
#ifndef NDEBUG
    if (!(info == GrB_Info::GrB_SUCCESS || info == GrB_Info::GrB_NO_VALUE)) {
        printf("info: %d error: %s\n", info, GrB_error());
    }
#endif
}

// Init GraphBLAS
// Should be done exactly once
void launch_GraphBLAS() {
    auto info = GrB_init(GrB_Mode::GrB_NONBLOCKING);
    check_for_error(info);
}

// Finalize GraphBLAS
// Should be done exactly once
void stop_GraphBLAS() {
    auto info = GrB_finalize();
    check_for_error(info);
}
