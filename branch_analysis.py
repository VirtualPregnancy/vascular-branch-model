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
        export_directory = path + sample + '/' + sample + '-branch/'
        if not os.path.exists(export_directory):
            print("Creating export directory")
            os.makedirs(export_directory)
        import_tree = import_directory + sample + '-tree-reordered'
        tree_nodes = pg.import_exnode_tree(import_tree + '.exnode')
        tree_elems = pg.import_exelem_tree(import_tree + '.exelem')
        import_rad = import_directory + sample + '-euclid-radius-reordered'
        euclid_radii = pg.import_exelem_field(import_rad + '.exelem')
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
        
        branch_geom = {}
        branch_geom['nodes']=tree_nodes['nodes']
        branch_geom['elems']=branch_elems['elems']
        branch_geom['euclidean length'] = pg.define_elem_lengths(geom['nodes'], branch_geom['elems'])
        geom, branch_geom, generation_table,strahler_table,branch_table = pg.analyse_branching(geom,branch_geom,'strahler',1.,1.)
        

        # csv files
        print('Writing files')
        output = export_directory + 'StrahlerTable.csv'
        headerTable = "'Order', 'NumBranches', 'Length(mm)', 'std', 'Diameter(mm)', 'std', 'EuclideanLength(mm)', 'std', 'Len/Diam', 'std', 'Tortuosity', 'std', 'Angles', 'std', 'LenRatio', 'std', 'DiamRatio', 'std'"
        np.savetxt(output, strahler_table, fmt='%.4f', delimiter=',', header=headerTable)
        
        output = export_directory + 'GenerationTable.csv'
        headerTable= "'Gen', 'NumBranches', 'Length(mm)', 'std', 'Diameter(mm)', 'std', 'Euclidean Length(mm)', 'std', 'Len/Diam', 'std', 'Tortuosity', 'std', 'Angles', 'std', 'Minor Angle', 'std', 'Major Angle', 'std', 'LLparent', 'std', 'LminLparent', 'std', 'LmajLparent', 'std', 'LminLmaj', 'std', 'DDparent', 'std', 'DminDparent', 'std', 'DmajDparent', 'std', 'DminDmaj', 'std'"
        np.savetxt(output, generation_table, fmt='%.4f', delimiter=',', header=headerTable)
      
        output = export_directory + 'OverallTable.csv'
        headerTable = "'Num branches', 'Total length','Total vessel volume', 'Total volume', 'vascular span','inlet diameter','num generations', 'num orders', 'ave term gen','std','tortuosity','std','branch length', 'std', 'euc length', 'std', 'diameter','std','L/D','std', 'branch angle', 'std','minor angle','std', 'major angle', 'std', 'D/Dparent', 'std', 'Dmin/Dparent','std', 'Dmaj/Dparent', 'std', 'L/Lparent', 'std','L/Lparent','std', 'Lmin/Lparent','std','Lmaj/Lparent', 'std', 'Lmaj/Lmin','std', 'Rb', 'rsq','Rd','rsq','Rl','rsq'"
        np.savetxt(output, branch_table, fmt='%.4f', delimiter=',', header=headerTable)
  



