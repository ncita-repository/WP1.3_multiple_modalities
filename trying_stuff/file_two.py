# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 20:20:28 2021

@author: ctorti
"""

print("File two __name__ is set to: {}" .format(__name__))

def function_three():
   print("Function three is executed")
   
def function_four():
   print("Function four is executed")
   
   
if __name__ == "__main__":
   print("File two executed when ran directly")
else:
   print("File two executed when imported")