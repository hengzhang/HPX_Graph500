#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/lcos/future_wait.hpp>
#include <hpx/components/distributing_factory/distributing_factory.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/util/unwrapped.hpp>
#include <boost/ref.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <boost/foreach.hpp>
#include <time.h>

#include <fstream>
#include <iostream>
#include <vector>
//!!from graph500 head
#include "../generator/make_graph.hpp"
#include "partition_component.hpp"
#include "common.hpp"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>



inline void 
init(
    hpx::components::server::distributing_factory::iterator_range_type r,  
    std::vector<hpx::naming::id_type>& p)
{
  BOOST_FOREACH(hpx::naming::id_type const& id, r)
    p.push_back(id);
}

int hpx_main(boost::program_options::variables_map& vm){
  {
    hpx::util::high_resolution_timer t_wall;

    std::size_t const SCALE = vm["SCALE"].as<std::size_t>();
    double const edgefactor = vm["edgefactor"].as<double>();
    std::size_t locality_num = vm["lnum"].as<std::size_t>();
    if(SCALE == 0 || edgefactor == 0){
      std::cout<<"Usage: SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [float, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n";
      return -1;
    }
    hpx::components::distributing_factory factory =  
      hpx::components::distributing_factory::create(hpx::find_here());
    hpx::components::component_type block_type = 
      graph500::server::partition::get_component_type();

    std::vector<hpx::naming::id_type> localities = hpx::find_all_localities(block_type);
    std::cout << " Number of localities: " << localities.size() << std::endl;

    hpx::components::distributing_factory::result_type blocks = 
      factory.create_components(block_type,  locality_num);
    std::vector<hpx::naming::id_type> locality_ids;
    init(hpx::util::locality_results(blocks), locality_ids);

    std::size_t size = localities.size();

    //[! Construct graph
    std::uint64_t seed1 = 2, seed2 = 3;
    const double initiator[4] = {.57, .19, .19, .05};

    std::int64_t nedges;
    std::vector<packed_edge> edges;
    std::int64_t* nedges;
    hpx::util::high_resolution_timer t_make_graph;
 //   make_graph(SCALE, (std::int64_t)(edgefactor * pow(2., SCALE)), seed1, seed2, initiator, &nedges, &edges);//from ../generator/make_graph.hpp
    std::vector<hpx::lcos::future<void> > init_phase;
    graph500::server::partition::init_action  init;
    for(std::size_t i =0; i<locality_num ; i++){
      init_phase.push_back(hpx::async(init, locality_ids[i], SCALE, (std::int64_t)(edgefactor * pow(2., SCALE)), seed1, seed2, i , size));
    }
    hpx::wait_all(init_phase);

    hpx::cout << "Elapsed Execution graph_make time: " << t_make_graph.elapsed() << " [s]" << std::endl; 
    //!]
    //[! CSR graph data structure.
    csr_graph g;
    /* Make CSR data structure, redistributing data using VERTEX_OWNER
     * distribution. */
    hpx::util::high_resolution_timer t_make_csr;
    std::vector<hpx::lcos::future<void> > convert_phase;
    graph500::server::partition::convert_graph_to_csr_action  convert;
    for(std::size_t i =0; i<locality_num ; i++){
      convert_phase.push_back(hpx::async(convert, locality_ids[i], nedges, edges, &g));
    }
    hpx::wait_all(convert_phase);
    hpx::cout << "Elapsed Execution convert_csr time: " << t_make_csr.elapsed() << " [s]" << std::endl; 
    free(edges); edges = NULL;
    //!]
    //[! Get roots for BFS runs. 
    int num_bfs_roots = 64;
    std::int64_t* bfs_roots = (std::int64_t*)malloc(num_bfs_roots * sizeof(std::int64_t));
    std::vector<hpx::lcos::future<void> > find_roots_phase;
    graph500::server::partition::find_bfs_roots  find_roots;
    for(std::size_t i =0; i<locality_num ; i++){
      find_roots_phase.push_back(hpx::async(find_roots, locality_ids[i], &num_bfs_roots, &g, seed1,seed2, bfs_roots,locality_ids));
    }
    hpx::wait_all(find_roots_phase);
    //!]

    /* Number of edges visited in each BFS; a double so get_statistics can be used directly. */
    double* edge_counts = (double*)malloc(num_bfs_roots * sizeof(double));


    //run bfs localities algorithm
    hpx::util::high_resolution_timer t_bfs_all;
    int validation_passed = 1;
    double* bfs_times = (double*)malloc(num_bfs_roots * sizeof(double));
    double* validate_times = (double*)malloc(num_bfs_roots * sizeof(double));
    std::int64_t* pred = (std::int64_t *) malloc(g.nlocalverts *sizeof(std::int64_t));

    int bfs_root_idx;
    for(bfs_root_idx = 0;bfs_root_idx < num_bfs_roots; ++bfs_root_idx){
      std::int64_t root = bfs_roots[bfs_root_idx];
      std::cout<<"Running BFS "<<bfs_root_idx << std::endl;

      memset(pred, 0, g.nlocalverts* sizeof(std::int64_t));
      
      //[! do BFS
      hpx::util::high_resolution_timer bfs_one_round;
      {
        std::int64_t nvisited = 0;
        std::vector<hpx::lcos::future<void> > run_bfs_phase;
        graph500::server::partition::run_hpx_bfs  run_bfs;
        for(std::size_t i =0; i<locality_num ; i++){
          run_bfs_phase.push_back(hpx::async(run_bfs, locality_ids[i], &g, root, &pred[0], &nvisited));
        }
        hpx::wait_all(run_bfs_phase);
      }
      bfs_times[bfs_root_idx] = bfs_one_round.elapsed();
      //!]

      //[! Validate result.
      std::cout<<"Validating BFS "<<bfs_root_idx<<std::endl;
      hpx::util::high_resolution_timer bfs_validate;
      {
        std::vector<hpx::lcos::future<int> > validation_phase;
        std::vector<int> validation_passed_one;
        validation_passed_one.reserve(locality_num);
        graph500::server::partition::validate_bfs_result validation_bfs;
        for(std::size_t i =0; i<locality_num ; i++){
          validation_phase.push_back(hpx::async(validate_bfs, locality_ids[i], &g, root, pred, nvisited));
        }
        hpx::wait_all(run_bfs_phase);
        hpx::lcos::wait(validation_phase, 
            [&](std::size_t, int t){
              validation_passed_one.push_back(t);
            });
      }
      validate_time[bfs_root_idx] = bfs_validate.elapsed();
      for(std::size_t i =0;i<validation_passed_one.size();i++){
        if(!validation_passed_one){
          validation_passed = 0;
          std::cout<< "Validation failed for this BFS root; skipping rest."<<std::endl;
          break;
        }
      }
   }
//[! Output result

//!]
 
    free(bfs_roots);
    free_csr_graph(&g);

    std::cout << "Elapsed BFS all time: " << t_bfs_all.elapsed() << " [s]" << std::endl;
    std::cout << "Elapsed Execution all time: " << t_wall.elapsed() << " [s]" << std::endl;
    free(bfs_times);
    free(validate_times);
  }
  return hpx::finalize(); // Initiate shutdown of the runtime system.
}

int main(int argc, char* argv[]){
  using boost::program_options::value;
  // Configure application-specific options.
  boost::program_options::options_description
    desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

  desc_commandline.add_options()
    ("lnum", value<std::size_t>()->default_value(1),
     "lnum = custom defined locality number[integer,defaults to 1]")
    ("SCALE", value<std::size_t>()->default_value(16),
     "SCALE = log_2(# vertices) [integer, required]")
    ("edgefactor", value<double>()->default_value("16"),
     "edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [float, defaults to 16]");

  return hpx::init(desc_commandline, argc, argv); 
}
