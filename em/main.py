
import reader
import processing
import visualizer

#Read molecule map from file
#mapReader = reader.Reader()
#Open file
#mapReader.open("../maps/1010/EMD-1010.map") #level 2.2479
# #mapReader.open("../maps/EMD-2596.map") #level 0.1462
#mapReader.open("../maps/1010/EMD-1010.map") #level 23
#mapReader.open("../maps/1364/EMD-1364.map") #level 39.0086
#mapReader.open("../maps/5017/EMD-5017.map") #level 17.3347


#Get map object
#myMap = mapReader.read()
# Create visualizer with a map surface threshold level
# Otherwise use otsu threshold
#v= Visualizer(myMap, level=9.07)
#v = Visualizer(myMap)
# Watershed 
#v.segmentation(step_sigma=1, steps=3)

# add corresponding atomic structure
#v.add_structure("../maps/1010/pdb1mi6.ent")
#v.map_structure_to_domain("../maps/1010/A-1GQE.aligned.pdb")
#v.map_structure_to_domain("../maps/1010/B-1GQE.aligned.pdb")
#v.map_structure_to_domain("../maps/1010/C-1GQE.aligned.pdb")
#v.map_structure_to_domain("../maps/1010/D-1GQE.aligned.pdb")

#v.add_structure("../maps/1364/pdb1pn6.ent")
#v.map_structure_to_domain("../maps/1364/A-1FNM.aligned.pdb")
#v.map_structure_to_domain("../maps/1364/B-1FNM.aligned.pdb")
#v.map_structure_to_domain("../maps/1364/C-1FNM.aligned.pdb")
#v.map_structure_to_domain("../maps/1364/D-1FNM.aligned.pdb")
#v.map_structure_to_domain("../maps/1364/E-1FNM.aligned.pdb")

#v.add_structure("../maps/5017/pdb3dny.ent")
#v.map_structure_to_domain("../maps/5017/A-1N0U.aligned.pdb")
#v.map_structure_to_domain("../maps/5017/B-1N0U.aligned.pdb")
#v.map_structure_to_domain("../maps/5017/C-1N0U.aligned.pdb")

#v.show_atom_matching()
#v.show()



if __name__ == "__main__":
    main()