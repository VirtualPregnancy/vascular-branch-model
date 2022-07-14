import numpy as np
import os
import placentagen as pg

if __name__ == '__main__':
    
    path = 'data/'
    identifiers = ['F8']
    ####################################################
    
    for sample in identifiers:
        print('Analysing sample',sample)
        import_directory = path + sample + '/' + sample + '-skeleton/'
        export_directory = path + sample + '/' + sample + '-model/'
        if not os.path.exists(export_directory):
            print("Creating export directory")
            os.makedirs(export_directory)
        import_tree = import_directory + sample + '-tree-reordered'
        tree_nodes = pg.import_exnode_tree(import_tree + '.exnode')
        tree_elems = pg.import_exelem_tree(import_tree + '.exelem')
        import_rad = import_directory + sample + '-normal-radius-reordered'
        normal_radii = pg.import_exelem_field(import_rad + '.exelem')
        import_b_id = import_directory + sample + '-branch-id-reordered'
        branch_id = pg.import_exelem_field(import_b_id + '.exelem')
        import_branch = import_directory + sample + '-branch-reordered'   
        branch_elems = pg.import_exelem_tree(import_branch + '.exelem')
        
        
        geom = {}
        geom['nodes']=tree_nodes['nodes']
        geom['elems']=tree_elems['elems']
        geom['radii']=normal_radii
        geom['length']=pg.define_elem_lengths(geom['nodes'], geom['elems'])
        geom['branch id'] = branch_id
        

        connectivity = pg.element_connectivity_1D(geom['nodes'], geom['elems'])
        for ne in range(0,len(geom['elems'])):
            if connectivity['elem_up'][ne,0]==0:
                inlet_element = ne
        
        node_radius = np.zeros(len(geom['nodes']), dtype=float)
        for i in range(0,len(geom['elems'])):
            node1 =  geom['elems'][i,1]
            node2 =  geom['elems'][i,2]
            if i == inlet_element:
                node_radius[node1] = geom['radii'][i]
                print(node1,node_radius[node1])
                node_radius[node2] = geom['radii'][i]
            else:
                node_radius[node2] = geom['radii'][i]
        
        export_tree = export_directory + sample + '-tree'
        pg.export_ip_coords(geom['nodes'][:,1:4],'arteries',export_tree)
        pg.export_ipelem_1d(geom['elems'],'arteries',export_tree)
        pg.export_ipfiel(node_radius,export_tree)

  
  



