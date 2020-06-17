# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:04:52 2020

@author: ctorti
"""

def ParseOutputPoints():
    # Import packages:
    import re

    """ THIS NEEDS TO BE CORRECTED! """
    
    # Define the patterns to search:
    """ 2D example:
        OutputPtsSearch = re.compile("OutputPoint\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s*\]")
    """
    
    Point_pattern = re.compile("Point\s*([\d]+)")
    #InputIndex_pattern = re.compile("InputIndex\s*=\s*\[\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*\]")
    InputIndex_pattern = re.compile("InputIndex\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s*\]")
    
    InputPoint_pattern = re.compile("InputPoint\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]")
    OutputIndexFixed_pattern = re.compile("OutputIndexFixed\s*=\s*\[\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*\]")
    OutputPoint_pattern = re.compile("OutputPoint\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]")
    Deformation_pattern = re.compile("Deformation\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]")
    OutputIndexMoving_pattern = re.compile("OutputIndexMoving\s*=\s*\[\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*\]")
    
    Point = []
    InputIndex = []
    InputPoint = []
    OutputIndexFixed = []
    OutputPoint = []
    Deformation = []
    OutputIndexMoving = []
    
    with open("outputpoints.txt", "r") as file:
        for line in file:
            Parsed = Point_pattern.search(line)
            PointNo = int(Parsed.group(1))
            Point.append(PointNo)
            
            Parsed = InputIndex_pattern.search(line)
            Indices = (int(Parsed.group(1)), int(Parsed.group(2)), int(Parsed.group(3)))
            InputIndex.append(Indices)
        
            Parsed = InputPoint_pattern.search(line)
            Points = (float(Parsed.group(1)), float(Parsed.group(2)), float(Parsed.group(3)))
            InputPoint.append(Points)
            
            Parsed = OutputIndexFixed_pattern.search(line)
            Indices = (int(Parsed.group(1)), int(Parsed.group(2)), int(Parsed.group(3)))
            OutputIndexFixed.append(Indices)
            
            Parsed = OutputPoint_pattern.search(line)
            Points = (float(Parsed.group(1)), float(Parsed.group(2)), float(Parsed.group(3)))
            OutputPoint.append(Points)
            
            Parsed = Deformation_pattern.search(line)
            Points = (float(Parsed.group(1)), float(Parsed.group(2)), float(Parsed.group(3)))
            Deformation.append(Points)
            
            Parsed = OutputIndexMoving_pattern.search(line)
            Indices = (int(Parsed.group(1)), int(Parsed.group(2)), int(Parsed.group(3)))
            OutputIndexMoving.append(Indices)
            
    return {'Point':Point, 'InputIndex':InputIndex, 'InputPoint':InputPoint,
            'OutputIndexFixed':OutputIndexFixed, 'OutputPoint':OutputPoint,
            'Deformation':Deformation, 'OutputIndexMoving':OutputIndexMoving}