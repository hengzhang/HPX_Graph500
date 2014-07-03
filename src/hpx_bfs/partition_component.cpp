
#include "partition_component.hpp"

HPX_REGISTER_COMPONENT_MODULE();

typedef hpx::components::managed_component<
  graph500::server::partition> partition_type;

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(partition_type, partition);

HPX_REGISTER_ACTION( graph500::server::partition::run_hpx_bfs_action, 
    partition_run_hpx_bfs_action);

HPX_REGISTER_ACTION( graph500::server::partition::validate_bfs_result_action, 
    partition_validate_bfs_result_action);

HPX_REGISTER_ACTION( graph500::server::partition::find_bfs_roots, 
    partition_find_bfs_roots_action);

HPX_REGISTER_ACTION( graph500::server::partition::init_action, 
    partition_init_action);
