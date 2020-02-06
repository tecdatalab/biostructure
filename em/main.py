
import reader
import visualizer
import processing
import model

import configparser
import argparse




def main():
    #config = configparser.ConfigParser()
    #parser = argparse.ArgumentParser(description='Molecule Visualization and segmentation tool')
    
    #Read molecule map from file
    mapReader = reader.Reader()
    #Open file
    mapReader.open("../maps/EMD-9590.map") #level 2.2479
    # #mapReader.open("../maps/EMD-2596.map") #level 0.1462
    #mapReader.open("../maps/1010/EMD-1010.map") #level 23
    #mapReader.open("../maps/1364/EMD-1364.map") #level 39.0086
    #mapReader.open("../maps/5017/EMD-5017.map") #level 17.3347


    #Get map object
    myMap = mapReader.read()
    # Create visualizer with a map surface threshold level
    # Otherwise use otsu threshold
    #v= visualizer.Visualizer(myMap, level=0.05)
    #
    # Watershed 
    contourRatioLvl = [1]
    myModel = model.Model(myMap, 0.05, contourRatioLvl)
    data_labels_list = processing.segment(myModel, step_sigma=1, steps=4)
    
    countourLvl = myModel.getCutoffLevels()
    v1 = visualizer.Visualizer(data_labels_list[0][0], countourLvl[0], data_labels_list[0][1])
    #v2 = visualizer.Visualizer(data_labels_list[1][0], countourLvl[1], data_labels_list[1][1])
    #v3 = visualizer.Visualizer(data_labels_list[2][0], countourLvl[2], data_labels_list[2][1])

    #v1.add_structure("../maps/pdb6acf.ent")
    #v2.add_structure("../maps/pdb6acf.ent")
    #v3.add_structure("../maps/pdb6acf.ent")

    v1.show()
    #v2.show()
    #v3.show()


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