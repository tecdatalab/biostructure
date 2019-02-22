from xml.dom import minidom
from classes.Map_information import Map_information
'''
Created on 22 feb. 2019

@author: luis98
'''

directory = "../db/"

def get_map_information(emd_id):
    doc = minidom.parse(directory +"xml/"+"emd-"+str(emd_id)+".xml")
    
    fileInformationElement = doc.getElementsByTagName("file")[0]
    file_information = [fileInformationElement.getAttribute("format"),
                        fileInformationElement.getAttribute("sizeKb"),
                        fileInformationElement.getAttribute("type"),
                        fileInformationElement.firstChild.data]
    
    dataTypeElement = doc.getElementsByTagName("dataType")[0]
    data_type = dataTypeElement.firstChild.data
    
    numColumnsElement = doc.getElementsByTagName("numColumns")[0]
    num_columns = numColumnsElement.firstChild.data
    
    numRowsElement = doc.getElementsByTagName("numRows")[0]
    num_rows = numRowsElement.firstChild.data
    
    numSectionsElement = doc.getElementsByTagName("numSections")[0]
    num_sections = numSectionsElement.firstChild.data
    
    originColElement = doc.getElementsByTagName("originCol")[0]
    origin_col = originColElement.firstChild.data
    
    originRowElement = doc.getElementsByTagName("originRow")[0]
    origin_row = originRowElement.firstChild.data
    
    originSecElement = doc.getElementsByTagName("originSec")[0]
    origin_sec = originSecElement.firstChild.data
    
    limitColElement = doc.getElementsByTagName("limitCol")[0]
    limit_col = limitColElement.firstChild.data
    
    limitRowElement = doc.getElementsByTagName("limitRow")[0]
    limit_row = limitRowElement.firstChild.data
    
    limitSecElement = doc.getElementsByTagName("limitSec")[0]
    limit_sec = limitSecElement.firstChild.data
    
    spacingColElement = doc.getElementsByTagName("spacingCol")[0]
    spacing_col = spacingColElement.firstChild.data
    
    spacingRowElement = doc.getElementsByTagName("spacingRow")[0]
    spacing_row = spacingRowElement.firstChild.data
    
    spacingSecElement = doc.getElementsByTagName("spacingSec")[0]
    spacing_sec = spacingSecElement.firstChild.data
    
    cellAElement = doc.getElementsByTagName("cellA")[0]
    cell_a = [cellAElement.getAttribute("units"),
              cellAElement.firstChild.data]
    
    cellBElement = doc.getElementsByTagName("cellB")[0]
    cell_b = [cellBElement.getAttribute("units"),
              cellBElement.firstChild.data]
    
    cellCElement = doc.getElementsByTagName("cellC")[0]
    cell_c = [cellCElement.getAttribute("units"),
              cellCElement.firstChild.data]
    
    cellAlphaElement = doc.getElementsByTagName("cellAlpha")[0]
    cell_alpha = [cellAlphaElement.getAttribute("units"),
                  cellAlphaElement.firstChild.data]
    
    cellBetaElement = doc.getElementsByTagName("cellBeta")[0]
    cell_beta = [cellBetaElement.getAttribute("units"),
                 cellBetaElement.firstChild.data]
    
    cellGammaElement = doc.getElementsByTagName("cellGamma")[0]
    cell_gamma = [cellGammaElement.getAttribute("units"),
                  cellGammaElement.firstChild.data]
    
    axisOrderFastElement = doc.getElementsByTagName("axisOrderFast")[0]
    axis_order_fast = axisOrderFastElement.firstChild.data
    
    axisOrderMediumElement = doc.getElementsByTagName("axisOrderMedium")[0]
    axis_order_medium = axisOrderMediumElement.firstChild.data
    
    axisOrderSlowElement = doc.getElementsByTagName("axisOrderSlow")[0]
    axis_order_slow = axisOrderSlowElement.firstChild.data
    
    minimuElement = doc.getElementsByTagName("minimum")[0]
    minimum = minimuElement.firstChild.data
    
    maximumElement = doc.getElementsByTagName("maximum")[0]
    maximum = maximumElement.firstChild.data
    
    averageElement = doc.getElementsByTagName("average")[0]
    average = averageElement.firstChild.data
    
    stdElement = doc.getElementsByTagName("std")[0]
    std = stdElement.firstChild.data
    
    spaceGroupNumberElement = doc.getElementsByTagName("spaceGroupNumber")[0]
    space_group_number = spaceGroupNumberElement.firstChild.data
    
    detailsElement = doc.getElementsByTagName("details")[0]
    details = detailsElement.firstChild.data
    
    pixelXElement = doc.getElementsByTagName("pixelX")[0]
    pixel_x = [pixelXElement.getAttribute("units"),
               pixelXElement.firstChild.data]
    
    pixelYElement = doc.getElementsByTagName("pixelY")[0]
    pixel_y = [pixelYElement.getAttribute("units"),
               pixelYElement.firstChild.data]
    
    pixelZElement = doc.getElementsByTagName("pixelZ")[0]
    pixel_z = [pixelZElement.getAttribute("units"),
               pixelZElement.firstChild.data]
    
    countourLevelElement = doc.getElementsByTagName("contourLevel")[0]
    countour_level = countourLevelElement.firstChild.data
    
    annotationDetailsElement = doc.getElementsByTagName("annotationDetails")[0]
    annotation_details = annotationDetailsElement.firstChild.data
    
    return Map_information(None, 
    file_information,
    data_type,
    num_columns,
    num_rows,
    num_sections,
    origin_col,
    origin_row, 
    origin_sec,
    limit_col,
    limit_row,
    limit_sec,
    spacing_col, 
    spacing_row,
    spacing_sec,
    cell_a,
    cell_b,
    cell_c,
    cell_alpha, 
    cell_beta,
    cell_gamma, 
    axis_order_fast, 
    axis_order_medium, 
    axis_order_slow,
    minimum,
    maximum,
    average,
    std,
    space_group_number, 
    details,
    pixel_x,
    pixel_y,
    pixel_z,
    countour_level, 
    annotation_details) 


if __name__ == '__main__':
    resultado = get_map_information(1364)
    print(resultado)
