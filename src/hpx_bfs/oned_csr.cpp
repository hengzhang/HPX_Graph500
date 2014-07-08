#include "oned_csr.hpp"
#include "common.hpp"

//!!!!the define is not completed
#define CONV1D_FUNCNAME convert_graph_to_oned_csr_helper

void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g) {
    g->tg = tg; 
    g->nlocaledges = 0;
    convert_graph_to_oned_csr_helper(tg, g);
}
