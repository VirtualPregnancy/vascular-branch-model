import numpy as np
import os
import SimpleITK as sitk
from skimage import measure, morphology
import placentagen as pg
from skan import skeleton_to_csgraph

if __name__ == '__main__':
    
    path = 'data/'
    identifiers = ['F8']
    write_skeleton=True #Can set to True if looking to see skeletonised image 
    write_distanceimg=True #This is used to create vessel radii and can set to True if you want to see this
    pixel_resolution_um = 6.499919 #set up for isotropic pixel size, here in um. "Scaling dimensions to mm" is the only place this is used, and so could be non-isotrophic
    pixel_resolution_mm = pixel_resolution_um/1000. #models typically use mm units
    ############
    #There is an option here to overide the inlet node number, this code assumes the inlet is the largest 'open ended' vessel which can sometimes not be true. To override you need to set the logical to true and visualise the exported tree (record node number of inlet noew)
    override_inlet = False
    inlet_cmgui = 0 # Needs to be an integer value
    
    #############################################################
    
    for sample in identifiers:#allows you to list many identifiers and loop through images
        print('Analysing sample',sample)
        export_directory = path + sample + '/' + sample + '-skeleton/'
        if not os.path.exists(export_directory):
            print("Creating export directory")
            os.makedirs(export_directory)

        print("Reading segmented image")
        seg= sitk.ReadImage(path + sample + '/'+ sample + '-segment/' + sample + '-segment.nii')
        print("Converting image to array for use in morphology skeletonize")        
        seg_array = sitk.GetArrayFromImage(seg)
        print('Creating Binary skeleton')
        #Ideally we'd use simple ITK where we can, but this binary thinning filter is 2D only, no 3d one released yet
        skel = morphology.skeletonize_3d(seg_array)
        if write_skeleton:
            print('Writing out skeletonised image')
            skel_img = sitk.GetImageFromArray(skel)
            skel_img.CopyInformation(seg)
            sitk.WriteImage(skel_img, export_directory + sample + '-skeleton.nii')
            #Deleting the skeleton image from memory, in case of working on the edge of memory limitations
            del skel_img
        print('Creating graph from skeleton')
        dimensions = 3
        pixel_graph, coordinates = skeleton_to_csgraph(skel)
        #deleting the skeleton from memory, again in case we are dealing with large datasets that take up heaps of memory
        del skel
        

        print('Creating an initial node and element structure from the skeleton')
        output_files = export_directory + sample + '-tree'
        elems,nodes, nodal_degrees = pg.create_graph_structure(pixel_graph, coordinates, dimensions,'arteries',output_files)
        
        print('Creating Euclidean signed distance imaging')
        distance = sitk.SignedMaurerDistanceMapImageFilter()
        distance.InsideIsPositiveOn()
        distance.SquaredDistanceOff()
        distance_img = distance.Execute(seg)
        if write_distanceimg:
            distance_img.CopyInformation(seg)
            sitk.WriteImage(distance_img, export_directory + sample + '-signed-distance.nii')
        print("Creating an initial Euclid distance based radius for each element in the graph")
        euclid_radii = pg.find_radius_euclidean(distance_img, elems,nodes)
        if write_distanceimg:
            output_files = export_directory + sample + '-euclid-radius'
            pg.export_exfield_1d_linear( euclid_radii, 'arteries', 'euradius', output_files)
            
            
        
        #For memory usage we delete the now not needed distance image from the graph
        del distance_img
        
        if not override_inlet:
            print("Finding inlet node")
            length_threshold = 10.
            inlet = pg.find_inlet_auto(elems,nodes,euclid_radii,length_threshold) #finds largest vessel, excluding short stubs < length_threshold 
        else: #cmgui node number is one more than internal python
            inlet = nodes[inlet_cmgui-1,:]

        print("Fixing branch direction to allow more intelligent trimming")
        elems,branch_id,branch_start,branch_end,cycles,seen = pg.fix_elem_direction(inlet[1:4],elems,nodes)
        print('fix direction',len(elems),inlet)
        print("Removing disconnected elements after branch analysis")
        elems, euclid_radii,branch_id = pg.remove_disconnected(elems, euclid_radii, branch_id, seen)
        print('lll',len(elems))
        if cycles.any():
            print("Cutting loops")
            elems, euclid_radii = pg.cut_loops(elems,nodes,branch_id,branch_start,branch_end,cycles,euclid_radii) 
            print("Recalculating element branching")
            elems,branch_id,branch_start,branch_end,cycles,seen = pg.fix_elem_direction(inlet[1:4],elems,nodes)
        print("Reordering branches with inlet as first element")
        elems,elem_map = pg.sort_from_inlet(inlet[1:4],nodes,elems,branch_id,branch_start,branch_end)
        elems,branch_id,branch_start,branch_end,cycles,seen = pg.fix_elem_direction(inlet[1:4],elems,nodes)
        euclid_radii_new = np.zeros(len(elems))
        for ne in range(0,len(elems)):
            euclid_radii_new[ne] = euclid_radii[elem_map[ne]]
            
        print("Removing sub-resolution branches (small radii)")
        elems,euclid_radii =pg.remove_small_radius(elems,euclid_radii_new,branch_id,branch_start,1.5) #Deleting elements with radius less than 1.5 pixel
        print("Recalculating element branching")
        elems,branch_id,branch_start,branch_end,cycles,seen = pg.fix_elem_direction(inlet[1:4],elems,nodes)
        print("Removing disconnected elements after branch analysis")
        elems, euclid_radii,branch_id = pg.remove_disconnected(elems, euclid_radii, branch_id, seen)
        
        
        print("Deleting short order 1 branches connected to the top 2 Strahler orders")
        length_threshold = 10.
        elems, euclid_radii = pg.remove_order1(nodes,elems,branch_id,euclid_radii,length_threshold) 
        print("Recalculating element branching")
        elems,branch_id,branch_start,branch_end,cycles,seen = pg.fix_elem_direction(inlet[1:4],elems,nodes)
        
        print("Deleting unused nodes") #Only do this once happy with branching as it will renumber the whole tree!
        nodes, elems = pg.delete_unused_nodes(nodes,elems)
        print("Calculating radius by normal projection")
        normal_radii = pg.find_radius_normal_projection(seg_array, elems, nodes, euclid_radii)
        
        print("Scaling dimensions to mm")
        normal_radii = normal_radii * pixel_resolution_mm
        euclid_radii = euclid_radii * pixel_resolution_mm
        nodes[:,1:4]=nodes[:,1:4]* pixel_resolution_mm 
        
        branch_elems = np.zeros((len(branch_start),3), dtype=int)
        for nb in range(0,len(branch_start)):
          nnod1 = elems[int(branch_start[nb]),1]
          nnod2 = elems[int(branch_end[nb]),2]
          branch_elems[nb,0] = nb
          branch_elems[nb,1] = nnod1
          branch_elems[nb,2] = nnod2 
          
        
        print("Exporting trees and radii")
        output_files = export_directory + sample + '-tree-reordered'
        pg.export_ex_coords(nodes,'arteries',output_files,'exnode') #Final nodes (all of them)
        pg.export_exelem_1d(elems, 'arteries', output_files) #Final elements (all of them)
        output_files = export_directory + sample + '-branch-id-reordered'
        pg.export_exfield_1d_linear( branch_id, 'arteries', 'branchid', output_files) #BRanch id for each element
        output_files = export_directory + sample + '-euclid-radius-reordered'
        pg.export_exfield_1d_linear( euclid_radii, 'arteries', 'euradius', output_files) #(innacurate) Euclid radius for each element
        output_files = export_directory + sample + '-normal-radius-reordered' #(more accuarate normal projection radius for each element)
        pg.export_exfield_1d_linear( normal_radii, 'arteries', 'normradius', output_files) 
        output_files = export_directory + sample + '-branch-reordered'     
        pg.export_exelem_1d(branch_elems, 'arteries', output_files)        #elements that reflect the lines between branch points only

        
   
        




