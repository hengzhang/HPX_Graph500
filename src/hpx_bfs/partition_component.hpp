#if !defined(PARTITION_COMPONENT_HPP)
#define PARTITION_COMPONENT_HPP

#include <vector>
#include <queue>
#include <iostream>
#include <hpx/include/components.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>
#include <hpx/runtime/actions/component_action.hpp>
 #include <hpx/runtime/components/stubs/stub_base.hpp>
#include <hpx/runtime/applier/applier.hpp>
#include <hpx/runtime/components/client_base.hpp>
#include <hpx/include/client.hpp>
#include <hpx/hpx_fwd.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/include/lcos.hpp> 
#include <hpx/include/actions.hpp>  
#include <hpx/runtime/components/component_factory.hpp>
#include <hpx/include/async.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include "common.h"
#include "../generator/make_graph.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

namespace graph500{
  namespace server{

    class HPX_COMPONENT_EXPORT partition: 
      public hpx::components::managed_component_base<partition>
    {
      typedef struct edge_page{
        std::int64_t targets[16];
        struct edge_page* next;
      }edge_page;

      typedef struct adjacency_list {
        std::size_t nvertices; /* User-visible number of vertices */
        std::size_t nvertices_allocated; /* Actual allocated size of data */
        edge_page** data; /* Indexed by vertex */
      } adjacency_list;

      private:
      std::size_t rank;
      std::size_t level;
      //std::size_t parent;
      bool visited;
      int root_ok;
      std::vector<packed_edge> edges;
      std::uint64_t nedges;

      public:
      partition(){}
      
      static inline edge_page* new_edge_page(void){
        edge_page * ep = (edge_page*) malloc(sizeof(edge_page));
        int i;
        ep->next = NULL;
        for(i =0; i< 16; ++i)
          ep->targets[i] = -1;
        return ep;
      }

      static void grow_adj_list(adjacency_list* al, size_t min_nvertices) {
        if (min_nvertices <= al->nvertices) return;
        while (min_nvertices > al->nvertices_allocated) {
          al->nvertices_allocated = (al->nvertices_allocated == 0) ? 16 : (al->nvertices_allocated * 2);
          al->data = (edge_page**)xrealloc(al->data, al->nvertices_allocated * sizeof(edge_page*));
        }
        size_t i;
        for (i = al->nvertices; i < min_nvertices; ++i) {
          al->data[i] = NULL;
        }
        al->nvertices = min_nvertices;
      }
      static void add_adj_list_edge(adjacency_list* al, size_t src, int64_t tgt) {
        grow_adj_list(al, src + 1);
        edge_page** p = al->data + src;
        /* Each page is filled before we allocate another one, so we only need to
         * check the last one in the chain. */
        while (*p && (*p)->next) {p = &((*p)->next);}
        if (*p) {
          assert (!(*p)->next);
          int i;
          for (i = 0; i < 16; ++i) {
            if ((*p)->targets[i] == -1) {
              (*p)->targets[i] = tgt;
              return;
            }
          }
          p = &((*p)->next);
          assert (!*p);
        }
        assert (!*p);
        *p = new_edge_page();
        (*p)->targets[0] = tgt;
      }

      static void clear_adj_list(adjacency_list* al) {
        size_t i;
        for (i = 0; i < al->nvertices; ++i) delete_edge_page(al->data[i]);
        free(al->data);
        al->data = NULL;
        al->nvertices = al->nvertices_allocated = 0;
      }

      //[! Partition_component Action
      //[! Init Action Function
      void init(std::size_t scale, 
          std::int64_t M,  
          std::uint64_t seed1, 
          std::uint64_t seed2, 
          std::size_t rank_in, 
          std::size_t size_in){
        // parent = 9999;
        rank = rank_in;
        root_ok = 0;
        edges = make_graph(SCALE, M, seed1, seed2, rank, size_in); 
        nedges = edges.size();
      }
      //!]

      //[! Convert_graph_to_csr Action Function
      void convert_graph_to_csr(const std::uint64_t nedges_in, 
          const std::vector<packed_edge> edges_in, 
          csr_graph * const g){
        adjacency_list adj_list = {0,0, NULL};

        #define PROCESS_EDGE(src, tgt) \
        /* Do handling for one edge.*/  \
        do{ \
          assert(VERTEX_OWNER((src)) == rank);  \
          add_adj_list_edge(&adj_list, VERTEX_LOCAL((src)), (tgt)); \ 
        } while(0)


      }
      //!]
      //[! Get_root_ok Action Function
       int get_rootok() const {
         return root_ok;
       }
      //!]
      //[! Find_bfs_roots Action Function
      void find_bfs_roots(
          int* num_bfs_roots, const csr_graph* g, 
          const std::uint64_t seed1, 
          const std::uint64_t seed2, 
          const std::int64_t * bfs_roots,
          std::vector<hpx::naming::id_type> const &partitions_components) {
        std::uint64_t counter = 0;
        int bfs_root_idx;
        for(bfs_root_idx = 0; bfs_root_idx < *num_bfs_roots; ++bfs_root_idx){
          std::int64_t root;
          while(1){
            std::vector<double> d;
            d.reserve(2);
            d = make_random_numbers(2, seed1, seed2, counter);
            root = (int64_t)((d[0] + d[1]) * g->nglobalverts) % g->nglobalverts;
            counter += 2;
            if (counter > 2*g->nglobalverts) break;
            int is_duplicate = 0;
            int i;
            for (i = 0; i < bfs_root_idx; ++i) {
              if (root == bfs_roots[i]) {
                is_duplicate = 1;
                break;
              }
            }
            if (is_duplicate) continue;
        //    int root_ok;
            if (VERTEX_OWNER(root) == rank) {
              root_ok = 0;
              size_t ei, ei_end = g->rowstarts[VERTEX_LOCAL(root) + 1];
              for (ei = g->rowstarts[VERTEX_LOCAL(root)]; ei < ei_end; ++ei) {
                if (g->column[ei] != root) {
                  root_ok = 1;
                  break;
                }
              }
            }
            graph500::server::partition::get_rootok_action get_rootok; 
            std::vector<hpx::lcos::future< int > >  lazy_results;
            for (std::size_t i=0;i<partitions_components.size();i++) {
              lazy_results.push_back( hpx::async(get_rootok, partitions_components[i]) );
            }
            std::vector<int > n_ = hpx::util::unwrapped(lazy_results); 
            for(std::vector<int>::iterator it = n_.begin(); it != n_.end(); ++it)
              if(*it) root_ok=*it;
            if(root_ok) break;
          }
          bfs_roots[bfs_root_idx] = root;
        }
        *num_bfs_roots = bfs_root_idx;
      }
      //!]
      //[! Run_hpx_bfs Action Function
      void run_hpx_bfs(const csr_graph* const g, 
          std::int64_t root, 
          std::int64_t* pred, 
          std::int64_t* nvisited) {

        const std::size_t nlocalverts = g->nlocalverts;
        std::uint64_t nvisited_local = 0;

        //set up queues
        std::uint64_t *oldq = 
          (std::uint64_t*)malloc(nlocalverts * sizeof(std::uint64_t));
        std::uint64_t *newq = (std::uint64_t*)malloc(nlocalverts * sizeof(std::uint64_t));
        std::size_t oldq_count = 0;
        std::size_t newq_count = 0;

        //set up visited bitmap
        const int ulong_bits = 
          sizeof(unsigned long) * CHAR_BIT;
        std::uint64_t visited_size = 
          (nlocalverts + ulong_bits -1) /ulong_bits;
        unsigned long *visited = 
          (unsigned long *)calloc(visited_size, sizeof(unsigned long));

#define SET_VISITED(v) do{visited[VERTEX_LOCAL((v))/ulong_bits] |= (1UL<<(VERTEX_LOCAL((v)) % ulong_bits)); }while(0)
#define TEST_VISITED(v) ((visited[VERTEX_LOCAL((v)) / ulong_bits] & (1UL<<(VERTEX_LOCAL((V)) % ulong_bits))) != 0)

        //set up buffers for message coalescing, for localities communication
        /*[!!This buffers need to replace 
          const int coalescing_size = 256;
          */

       /* Termination counter for each level: this variable counts the number of
         * ranks that have said that they are done sending to me in the current
         * level.  This rank can stop listening for new messages when it reaches
         * size. */
        int num_ranks_done;
        /* Set all vertices to "not visited." */
        {std::size_t i; for (i = 0; i < nlocalverts; ++i) pred[i] = -1;}

        /* Mark the root and put it into the queue. */
        if (VERTEX_OWNER(root) == rank) {
          SET_VISITED(root);
          pred[VERTEX_LOCAL(root)] = root;
          oldq[oldq_count++] = root;
        }
        while(1){
          num_ranks_done = 1;//never send to myself, so I'm always done
          //!!!Receive buffer from other localities
          //To Do

          recvreq_active =1;

          //Step through the current level's queue.
          std::size_t i;
          for(i = 0; i<oldq_count; ++i){
            assert(VERTEX_OWNER(oldq[i]) == rank);
            assert(pred[VERTEX_LOCAL(oldq[i])] >= 0 && pred[VERTEX_LOCAL(oldq[i])] < g->nglobalverts);
            ++nvisited_local;
            std::uint64_t src = oldq[i];
            //Iterate through its incident edges.
            std::size_t j, j_end = g->rowstarts[VERTEX_LOCAL(oldq[i]) +1];
            for (j = g->rowstarts[VERTEX_LOCAL(oldq[i])]; j < j_end; ++j) {
              std::uint64_t tgt = g->column[j];
              int owner = VERTEX_OWNER(tgt);
              /* If the other endpoint is mine, update the visited map, predecessor
               * map, and next-level queue locally; otherwise, send the target and
               * the current vertex (its possible predecessor) to the target's owner.
               * */
              if (owner == rank) {
                if (!TEST_VISITED(tgt)) {
                  SET_VISITED(tgt);
                  pred[VERTEX_LOCAL(tgt)] = src;
                  newq[newq_count++] = tgt;
                }
              } else {
                //!!! need to fix, modify with HPX
                std::size_t c = outgoing_counts[owner];
                outgoing[owner * coalescing_size * 2 + c] = tgt;
                outgoing[owner * coalescing_size * 2 + c + 1] = src;
                outgoing_counts[owner] += 2;
                if (outgoing_counts[owner] == coalescing_size * 2) {
                  outgoing_reqs_active[owner] = 1;
                  outgoing_counts[owner] = 0;
                }
              }
            }
          }
          //Flush any coalescing buffers that still have messages.
          int offset;
          for (offset = 1; offset < size; ++offset) {
            int dest = MOD_SIZE(rank + offset);
            if (outgoing_counts[dest] != 0) {
              outgoing_reqs_active[dest] = 1;
              outgoing_counts[dest] = 0;
            }
            /* Wait until all sends to this destination are done. */
            /* Tell the destination that we are done sending to them. */
          }
          //wait until everyone else is done (and thus couldn't send us any more messages)
          /* Test globally if all queues are empty. */
          std::uint64_t global_newq_count;

          /* Quit if they all are empty. */
          if (global_newq_count == 0) break;

          /* Swap old and new queues; clear new queue for next level. */
          {std::uint64_t* temp = oldq; oldq = newq; newq = temp;}
          oldq_count = newq_count;
          newq_count = 0;
        }

        /* Add up total visited-vertex count. */
        *nvisited = nvisited_local;

        free(oldq);
        free(newq);
        free(outgoing_counts);
        free(outgoing_reqs);
        free(outgoing_reqs_active);
        free(visited);
      }
      //!]
      //[!Validate_bfs_result Action Function
      int validate_bfs_result(const csr_graph* g, 
          const std::int64_t root, 
          const std::int64_t* pred, 
          const std::int64_t nvisited) {
        int validation_passed = 1;
        int root_is_mine = (VERTEX_OWNER(root) == rank);

        const size_t nlocalverts = g->nlocalverts;
        const size_t nlocaledges = g->nlocaledges;
        const int64_t nglobalverts = g->nglobalverts;

        /* Check that root is its own parent. */
        if (root_is_mine) {
          if (pred[VERTEX_LOCAL(root)] != root) {
            std::cout<<rank<<":  Validation error: parent of root vertex "<<root<<" is "<<pred[VERTEX_LOCAL(root)]<<std::endl;
            validation_passed = 0;
          }
        }

        /* Check that nothing else is its own parent, and check for in-range
         * values. */
        int any_range_errors = 0;
        size_t i;
        for (i = 0; i < nlocalverts; ++i) {
          int64_t v = VERTEX_TO_GLOBAL(i);
          assert (VERTEX_OWNER(v) == rank);
          assert (VERTEX_LOCAL(v) == i);
          if (v != root && pred[i] == v) {
            std::cout<<rank<<":  Validation error: parent of non-root vertex "<<v<<" is itself."<<std::endl;
            validation_passed = 0;
          }
          if (pred[i] < -1 || pred[i] >= nglobalverts) {
            std::cout<<rank<<":  Validation error: parent of root vertex "<<root<<" is "<<v<<"."<<std::endl;
            validation_passed = 0;
            any_range_errors = 1;
          }
        }

        /* Check that nvisited is correct. */
        int64_t nvisited_actual = 0;
        for (i = 0; i < nlocalverts; ++i) {
          if (pred[i] != -1) ++nvisited_actual;
        }
        if (nvisited_actual != nvisited) {
          std::cout<<rank <<": Validation error: claimed visit count "<<nvisited<< "is different from actual count "<<nvisited_actual<<"."<<std::endl;
          validation_passed = 0;
        }

        if (!any_range_errors) { /* Other parts of validation assume in-range values */

          /* Check that there is an edge from each vertex to its claimed
           * predecessor. */
          size_t i;
          for (i = 0; i < nlocalverts; ++i) {
            int64_t v = VERTEX_TO_GLOBAL(i);
            int64_t p = pred[i];
            if (p == -1) continue;
            int found_pred_edge = 0;
            if (v == p) found_pred_edge = 1; /* Root vertex */
            size_t ei, ei_end = g->rowstarts[i + 1];
            for (ei = g->rowstarts[i]; ei < ei_end; ++ei) {
              int64_t w = g->column[ei];
              if (w == p) {
                found_pred_edge = 1;
                break;
              }
            }
            if (!found_pred_edge) {
              std::cout<<rank<<": Validation error: no graph  edge from vertex "<<v<<" to its parent "<<p<<"."<<std::endl;
              validation_passed = 0;
            }
          }

          /* Create a vertex depth map to use for later validation. */
          int64_t* depth = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t));
          { /* Scope some code that has a lot of temporary variables. */
            int64_t* pred_depth = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Depth of predecessor vertex for each local vertex */
            size_t i;
            for (i = 0; i < nlocalverts; ++i) depth[i] = INT64_MAX;
            if (root_is_mine) depth[VERTEX_LOCAL(root)] = 0;
            /* Send each vertex that appears in the local part of the predecessor map
             * to its owner; record the original locations so we can put the answers
             * into pred_depth. */
            /* Do a histogram sort by owner (this same kind of sort is used other
             * places as well).  First, count the number of vertices going to each
             * destination. */
            int* num_preds_per_owner = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
            for (i = 0; i < nlocalverts; ++i) {
              ++num_preds_per_owner[pred[i] == -1 ? size - 1 : VERTEX_OWNER(pred[i])];
            }
            int64_t* preds_per_owner = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Predecessors sorted by owner */
            int64_t* preds_per_owner_results_offsets = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Indices into pred_depth to write */
            /* Second, do a prefix sum to get the displacements of the different
             * owners in the outgoing array. */
            int* pred_owner_displs = (int*)xmalloc((size + 1) * sizeof(int));
            pred_owner_displs[0] = 0;
            int r;
            for (r = 0; r < size; ++r) {
              pred_owner_displs[r + 1] = pred_owner_displs[r] + num_preds_per_owner[r];
            }
            /* Last, put the vertices into the correct positions in the array, based
             * on their owners and the counts and displacements computed earlier. */
            int* pred_owner_offsets = (int*)xmalloc((size + 1) * sizeof(int));
            memcpy(pred_owner_offsets, pred_owner_displs, (size + 1) * sizeof(int));
            for (i = 0; i < nlocalverts; ++i) {
              int* offset_ptr = &pred_owner_offsets[pred[i] == -1 ? size - 1 : VERTEX_OWNER(pred[i])];
              preds_per_owner[*offset_ptr] = pred[i];
              preds_per_owner_results_offsets[*offset_ptr] = i;
              ++*offset_ptr;
            }
            for (r = 0; r < size; ++r) {
              assert (pred_owner_offsets[r] == pred_owner_displs[r + 1]);
            }
            free(pred_owner_offsets);

            /* Send around the number of vertices that will be sent to each destination. */
            int* num_my_preds_per_sender = (int*)xmalloc(size * sizeof(int));
            int* my_preds_per_sender_displs = (int*)xmalloc((size + 1) * sizeof(int));
            my_preds_per_sender_displs[0] = 0;
            for (r = 0; r < size; ++r) {
              my_preds_per_sender_displs[r + 1] = my_preds_per_sender_displs[r] + num_my_preds_per_sender[r];
            }
            /* Send around the actual vertex data (list of depth requests that will
             * be responded to at each BFS iteration). */
            int64_t* my_depth_requests = (int64_t*)xmalloc(my_preds_per_sender_displs[size] * sizeof(int64_t));
            int64_t* my_depth_replies = (int64_t*)xmalloc(my_preds_per_sender_displs[size] * sizeof(int64_t));

            int64_t* pred_depth_raw = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Depth of predecessor vertex for each local vertex, ordered by source proc */

            /* Do a mini-BFS (naively) over just the predecessor graph (hopefully a
             * tree) produced by the real BFS; fill in the depth map. */
            while (1) {
              int any_changed = 0;
              int i;
              /* Create and send the depth values requested by other nodes.  The list
               * of requests is sent once, and are stored on the receiver so the
               * replies can be sent (possibly with updated depth values) at every
               * iteration. */
              for (i = 0; i < my_preds_per_sender_displs[size]; ++i) {
                my_depth_replies[i] = (my_depth_requests[i] == -1 ? INT64_MAX : depth[VERTEX_LOCAL(my_depth_requests[i])]);
              }
              {
                size_t i;
                /* Put the received depths into the local array. */
                for (i = 0; i < nlocalverts; ++i) {
                  pred_depth[preds_per_owner_results_offsets[i]] = pred_depth_raw[i];
                }
                /* Check those values to determine if they violate any correctness
                 * conditions. */
                for (i = 0; i < nlocalverts; ++i) {
                  int64_t v = VERTEX_TO_GLOBAL(i);
                  if (v == root) {
                    /* The depth and predecessor for this were checked earlier. */
                  } else if (depth[i] == INT64_MAX && pred_depth[i] == INT64_MAX) {
                    /* OK -- depth should be filled in later. */
                  } else if (depth[i] == INT64_MAX && pred_depth[i] != INT64_MAX) {
                    depth[i] = pred_depth[i] + 1;
                    any_changed = 1;
                  } else if (depth[i] != pred_depth[i] + 1) {
                    std::cout<<rank<<": Validation error: BFS predecessors do not form a tree; see vertices "<<v<< " (depth "<<depth[i]<<") and "<<pred[i]<<" (depth "<<pred_depth[i]<<")."<<std::endl;
                    validation_passed = 0;
                  } else {
                    /* Vertex already has its correct depth value. */
                  }
                }
              }
              if (!any_changed) break;
            }

            free(num_preds_per_owner);
            free(num_my_preds_per_sender);
            free(preds_per_owner);
            free(preds_per_owner_results_offsets);
            free(my_preds_per_sender_displs);
            free(my_depth_requests);
            free(my_depth_replies);
            free(pred_owner_displs);
            free(pred_depth);
            free(pred_depth_raw);
          }

          /* Check that all edges connect vertices whose depths differ by at most
           * one. */
          {
            int64_t maxlocaledges = 0;
            const int edge_chunk_size = (1 << 23); /* Reduce memory usage */
            int num_edge_groups = (maxlocaledges + edge_chunk_size - 1) / edge_chunk_size;
            int eg;
            for (eg = 0; eg < num_edge_groups; ++eg) {
              size_t first_edge_index = (size_t)(eg * edge_chunk_size);
              if (first_edge_index > nlocaledges) first_edge_index = nlocaledges;
              size_t last_edge_index = (size_t)((eg + 1) * edge_chunk_size);
              if (last_edge_index > nlocaledges) last_edge_index = nlocaledges;
              /* Sort the edge targets in this chunk by their owners (histogram
               * sort); see the BFS code above for details of the steps of the
               * algorithm. */
              int* num_edge_targets_by_owner = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
              size_t ei;
              for (ei = first_edge_index; ei < last_edge_index; ++ei) {
                ++num_edge_targets_by_owner[VERTEX_OWNER(g->column[ei])];
              }
              int* edge_targets_by_owner_displs = (int*)xmalloc((size + 1) * sizeof(int));
              edge_targets_by_owner_displs[0] = 0;
              int i;
              for (i = 0; i < size; ++i) {
                edge_targets_by_owner_displs[i + 1] = edge_targets_by_owner_displs[i] + num_edge_targets_by_owner[i];
              }
              int64_t* edge_targets_by_owner = (int64_t*)xmalloc(edge_targets_by_owner_displs[size] * sizeof(int64_t));
              int64_t* edge_targets_by_owner_indices = (int64_t*)xmalloc(edge_targets_by_owner_displs[size] * sizeof(int64_t)); /* Source indices for where to write the targets */
              int* edge_targets_by_owner_offsets = (int*)xmalloc((size + 1) * sizeof(int));
              memcpy(edge_targets_by_owner_offsets, edge_targets_by_owner_displs, (size + 1) * sizeof(int));
              for (ei = first_edge_index; ei < last_edge_index; ++ei) {
                edge_targets_by_owner[edge_targets_by_owner_offsets[VERTEX_OWNER(g->column[ei])]] = g->column[ei];
                edge_targets_by_owner_indices[edge_targets_by_owner_offsets[VERTEX_OWNER(g->column[ei])]] = ei;
                ++edge_targets_by_owner_offsets[VERTEX_OWNER(g->column[ei])];
              }
              for (i = 0; i < size; ++i) {
                assert (edge_targets_by_owner_offsets[i] == edge_targets_by_owner_displs[i + 1]);
              }
              free(edge_targets_by_owner_offsets);

              /* Send around the number of data elements that will be sent later. */
              int* num_incoming_targets_by_src = (int*)xmalloc(size * sizeof(int));
              int* incoming_targets_by_src_displs = (int*)xmalloc((size + 1) * sizeof(int));
              incoming_targets_by_src_displs[0] = 0;
              for (i = 0; i < size; ++i) {
                incoming_targets_by_src_displs[i + 1] = incoming_targets_by_src_displs[i] + num_incoming_targets_by_src[i];
              }

              int64_t* target_depth_requests = (int64_t*)xmalloc(incoming_targets_by_src_displs[size] * sizeof(int64_t));
              int64_t* target_depth_replies = (int64_t*)xmalloc(incoming_targets_by_src_displs[size] * sizeof(int64_t));

              /* Send the actual requests for the depths of edge targets. */

              free(edge_targets_by_owner);

              /* Fill in the replies for the requests sent to me. */
              for (i = 0; i < incoming_targets_by_src_displs[size]; ++i) {
                assert (VERTEX_OWNER(target_depth_requests[i]) == rank);
                target_depth_replies[i] = depth[VERTEX_LOCAL(target_depth_requests[i])];
              }

              free(target_depth_requests);

              int64_t* target_depth_raw = (int64_t*)xmalloc((last_edge_index - first_edge_index) * sizeof(int64_t));

              /* Send back the replies. */

              free(target_depth_replies);
              free(num_incoming_targets_by_src);
              free(num_edge_targets_by_owner);
              free(incoming_targets_by_src_displs);
              free(edge_targets_by_owner_displs);

              int64_t* target_depth = (int64_t*)xmalloc((last_edge_index - first_edge_index) * sizeof(int64_t));

              /* Put the replies into the proper order (original order of the edges).
               * */
              for (ei = 0; ei < last_edge_index - first_edge_index; ++ei) {
                target_depth[edge_targets_by_owner_indices[ei] - first_edge_index] = target_depth_raw[ei];
              }

              free(target_depth_raw);
              free(edge_targets_by_owner_indices);

              /* Check the depth relationship of the endpoints of each edge in the
               * current chunk. */
              size_t src_i = 0;
              for (ei = first_edge_index; ei < last_edge_index; ++ei) {
                while (ei >= g->rowstarts[src_i + 1]) {
                  ++src_i;
                }
                int64_t src = VERTEX_TO_GLOBAL(src_i);
                int64_t src_depth = depth[src_i];
                int64_t tgt = g->column[ei];
                int64_t tgt_depth = target_depth[ei - first_edge_index];
                if (src_depth != INT64_MAX && tgt_depth == INT64_MAX) {
                  std::cout<<rank<<": Validation error: edge connects vertex  " <<src<<" in the BFS tree (depth "<<src_depth<<") to vertex "<<tgt<<" outside the tree."<<std::endl;
                  validation_passed = 0;
                } else if (src_depth == INT64_MAX && tgt_depth != INT64_MAX) {
                  /* Skip this for now; this problem will be caught when scanning
                   * reversed copy of this edge.  Set the failure flag, though,
                   * just in case. */
                  validation_passed = 0;
                } else if (src_depth - tgt_depth < -1 ||
                    src_depth - tgt_depth > 1) {
                  std::cout<<rank<<":  Validation error: depths of edge endpoints "<<src<<" (depth "<<src_depth<< ") and "<<tgt<< " (depth "<<tgt_depth<<") are too far apart (abs. val. > 1)."<<std::endl;
                  validation_passed = 0;
                }
              }
              free(target_depth);
            }
          }

          free(depth);
        } /* End of part skipped by range errors */
        /* Collect the global validation result. */
        return validation_passed;
      }
      //!]
      //!]
      HPX_DEFINE_COMPONENT_ACTION(partition, init);
      HPX_DEFINE_COMPONENT_ACTION(partition, run_hpx_bfs);
      HPX_DEFINE_COMPONENT_ACTION(partition, validate_bfs_result);
      HPX_DEFINE_COMPONENT_ACTION(partition, find_bfs_roots);
      HPX_DEFINE_COMPONENT_ACTION(partition, get_rootok);
    };
  }

  namespace stubs{
    struct partition : hpx::components::stub_base<server::partition> {
      static void init(hpx::naming::id_type const & gid,
          std::size_t scale, 
          std::int64_t M,  
          std::uint64_t seed1, 
          std::uint64_t seed2, 
          std::size_t rank_in, 
          std::size_t size_in){
        hpx::apply<server::partition::init_action>(gid, scale, M, seed1, seed2, rank_in, size_in); 
      } 

      static void find_bfs_roots_async(hpx::naming::id_type const& gid_in, 
          int* num_bfs_roots, 
          const csr_graph* g, 
          const std::uint64_t seed1, 
          const std::uint64_t seed2, 
          const std::int64_t * bfs_roots){
        hpx::async<server::partition::find_bfs_roots_action>(gid_in, num_bfs_roots, g, seed1, seed2, bfs_roots);
      }
      static void find_bfs_roots_async(hpx::naming::id_type const& gid_in, 
          int* num_bfs_roots, 
          const csr_graph* g, 
          const std::uint64_t seed1, 
          const std::uint64_t seed2, 
          const std::int64_t * bfs_roots){
        find_bfs_roots_async(gid_in, num_bfs_roots, g, seed1, seed2, bfs_roots);
      }

      static hpx::lcos::future<int> get_rootok_async(hpx::naming::id_type const& gid_in){
        return hpx::async<server::partition::get_rootok>(gid_in);
      }
      static hpx::lcos::future<int> get_rootok(hpx::naming::id_type const& gid_in){
        return get_rootok_async(gid_in);
      }
 
      static void find_bfs_roots_async(hpx::naming::id_type const& gid_in, 
          int* num_bfs_roots, 
          const csr_graph* g, 
          const std::uint64_t seed1, 
          const std::uint64_t seed2, 
          const std::int64_t * bfs_roots){
        find_bfs_roots_async(gid_in, num_bfs_roots, g, seed1, seed2, bfs_roots);
      }


      static void run_hpx_bfs_async(hpx::naming::id_type const& gid_in, 
          const csr_graph* g, 
          const std::int64_t root,
          const std::int64_t* pred,
          const std::int64_t nvisited){
        hpx::async<server::partition::run_hpx_bfs_action>(gid_in, g,root, pred, nvisited);
      }
      static void run_hpx_bfs(hpx::naming::id_type const& gid_in, 
          const csr_graph* g, 
          const std::int64_t root,
          const std::int64_t* pred,
          const std::int64_t nvisited){
        run_hpx_bfs_async(gid_in, g,root, pred, nvisited);
      }

      static hpx::lcos::future<int> validate_bfs_result_async(hpx::naming::id_type const& gid_in, 
          const csr_graph* g, 
          const std::int64_t root,
          const std::int64_t* pred,
          const std::int64_t nvisited){
        return hpx::async<server::partition::validate_bfs_result_action>(gid_in, g,root, pred, nvisited);
      }
      static hpx::lcos::future<int> validate_bfs_result_async(hpx::naming::id_type const& gid_in, 
          const csr_graph* g, 
          const std::int64_t root,
          const std::int64_t* pred,
          const std::int64_t nvisited){
        return validate_bfs_result_async(gid_in, g,root, pred, nvisited);
      }
    };
  }

  class partition : hpx::components::client_base<partition, stubs::partition>{
    typedef hpx::components::client_base<partition, stubs::partition>  base_type;

    public:
    partition(){}

    void init(std::size_t scale, 
        std::int64_t M,  
        std::uint64_t seed1, 
        std::uint64_t seed2, 
        std::size_t rank_in, 
        std::size_t size_in){
      HPX_ASSERT(this->get_gid());
      this->base_type::init(this->get_gid(), scale, M, seed1, seed2, rank_in, size_in);
    }

    int get_rootok_async() {
      HPX_ASSERT(this->get_gid());
      return this->base_type::get_rootok_async(this->get_gid());
    }

    int get_rootok() {
      HPX_ASSERT(this->get_gid());
      return this->base_type::get_rootok(this->get_gid());
    }

    void find_bfs_roots_async(
        int* num_bfs_roots, const csr_graph* g, 
        const std::uint64_t seed1, 
        const std::uint64_t seed2, 
        const std::int64_t * bfs_roots) {
      HPX_ASSERT(this->get_gid());
      this->base_type::find_bfs_roots_async(this->get_gid(), num_bfs_roots, g, seed1, seed2, bfs_roots);
    }

    void find_bfs_roots(
        int* num_bfs_roots, const csr_graph* g, 
        const std::uint64_t seed1, 
        const std::uint64_t seed2, 
        const std::int64_t * bfs_roots) {
      HPX_ASSERT(this->get_gid());
      this->base_type::find_bfs_roots(this->get_gid(), num_bfs_roots, g, seed1, seed2, bfs_roots);
    }

    void run_hpx_bfs_async(const csr_graph* const g, 
        std::int64_t root, 
        std::int64_t* pred, 
        std::int64_t* nvisited){
      HPX_ASSERT(this->get_gid());
      this->base_type::run_hpx_bfs_async(this->get_gid(), g, root, pred, nvisited);
    }
    void run_hpx_bfs(const csr_graph* const g, 
        std::int64_t root, 
        std::int64_t* pred, 
        std::int64_t* nvisited){
      HPX_ASSERT(this->get_gid());
      this->base_type::run_hpx_bfs(this->get_gid(), g, root, pred, nvisited);
    }

    int validate_bfs_result_async(const csr_graph* g, 
        const std::int64_t root, 
        const std::int64_t* pred, 
        const std::int64_t nvisited) {
      HPX_ASSERT(this->get_gid());
      return this->base_type::validate_bfs_result_async(this->get_gid(), g, root, pred, nvisited);
    }

    int validate_bfs_result(const csr_graph* g, 
        const std::int64_t root, 
        const std::int64_t* pred, 
        const std::int64_t nvisited) {
      HPX_ASSERT(this->get_gid());
      return this->base_type::validate_bfs_result(this->get_gid(), g, root, pred, nvisited);
    }

  };
}

HPX_REGISTER_ACTION_DECLARATION(graph500::server::partition::init_action, 
    partition_init_action);

HPX_REGISTER_ACTION_DECLARATION(graph500::server::partition::run_hpx_bfs_action, 
    partition_run_hpx_bfs_action);

HPX_REGISTER_ACTION_DECLARATION(graph500::server::partition::find_bfs_roots_action, 
    partition_find_bfs_roots_action);

HPX_REGISTER_ACTION_DECLARATION(graph500::server::partition::validate_bfs_result_action, 
    partition_validate_bfs_result_action);

#endif
