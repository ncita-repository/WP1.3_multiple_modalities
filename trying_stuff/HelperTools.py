# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:57:05 2020

@author: ctorti
"""


# Import packages:
import copy



"""
******************************************************************************
******************************************************************************
GENERAL HELPER FUNCTIONS
******************************************************************************
******************************************************************************
"""


def ItemsUniqueToWithin(List, epsilon=1e-5):
    
    UniqueList = []
    
    L = len(List)
    
    if L < 2:
        print(f'The list has {L} items (it must have at least 2 items).')
        
        return UniqueList
    
    else:
        UniqueList.append(List[0])
        
        for i in range(L - 1):
            if abs(List[i] - List[0]) > epsilon:
                UniqueList.append(List[i])
                
        return UniqueList
    



def AddItemToListAndSort(OrigList, Item):
    """
    Append an item to a list and return the sorted list. 
    
    Inputs:
        OrigList - List of items 
        
        Item     - Item to add to List
        
    Returns:
        NewList  - Item appended to OrigList and sorted 
    """
    
    NewList = copy.deepcopy(OrigList)

    NewList.append(Item)
    
    NewList.sort()
    
    return NewList