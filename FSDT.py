# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 22:34:52 2020

@author: admin
"""

import numpy as np
import math

def Qbar(theta,Q11,Q22,Q12,Q66):
#    Calculates the m and n transformation matrix components.     
    m = math.cos(math.radians(theta))
    n = math.sin(math.radians(theta))
#    Defines the transformation matrix 
    T = np.matrix([[m**2, n**2, 2*m*n],[n**2, m**2, -2*m*n],[-m*n, m*n, (m**2 - n**2)]])
#    Manipulates the transformation matrix
    Ti = np.linalg.inv(T)
    Tt = np.transpose(Ti)

#    Defines the Q matrix
    Q = np.matrix([[Q11 , Q12, 0],[Q12, Q22, 0],[0, 0, Q66]])
#    Computes Qbar matrix and returns it 
    qbar = Ti * Q * Tt
    return(qbar)
    
def Qbarshear(theta,Q44,Q55):
#    Calculates the m and n transformation matrix components.     
    m = math.cos(math.radians(theta))
    n = math.sin(math.radians(theta))
#    Defines the transformation matrix 
    T = np.matrix([[m, -n],[n, m]])
#    Manipulates the transformation matrix
    Ti = np.linalg.inv(T)
    Tt = np.transpose(Ti)
#    Defines the Q matrix
    Q = np.matrix([[Q44 , 0],[0, Q55]])
#    Computes Qbar matrix and returns it 
    qbar = Ti * Q * Tt
    return(qbar)
    
def Function2(E1,E2,V12,G23,G13,G12,pthick,stse,Z):
   "Function 2 computes A, B and D matrix. Function 3 computes Exlam, Eylam, Vxylam and Gxylam. Z refers to the round factor you want applied on final matrix answer. "
#   Setups the variables used throughout the function.
   d = {}
   d1 = {}
   o = 0
   r =len(stse)+ 1    
#   Setups the arrays used throughout the function as arrays of 0 that will be added to.
   ktable = np.zeros((r,5))
   A = np.zeros((3,3))
   B = np.zeros((3,3))
   D = np.zeros((3,3))
   Abar = np.zeros((2,2))
#   Calculates V21
   V21 = (E2 * V12) / E1 
#   Calculates Q11, Q22, Q12 and Q66
   Q11 = E1 / (1 - (V12 * V21))
   Q22 = E2 / (1 - (V12 * V21))
   Q12 = (V12 * E2) / (1 - (V12 * V21))
   Q44 = G23
   Q55 = G13
   Q66 = G12 
#   Generates a dictionary (python data type) of the Qbar matrix for each orientation. Using the Qbar functions. 
   for i in stse:
       o = 1 + o
       d["{0}".format(o)]= Qbar(i,Q11,Q22,Q12,Q66)
   for i in stse:
       o = 1 + o
       d1["{0}".format(o)]= Qbarshear(i,Q44,Q55) 
#   Calculates the positive max of the stacking sequences this is inputted into the first ktable cell.
   ktable[0,0] = sum(pthick)/2
#   Calculates the first row of the ktable
   for i in range(1,r):
       ktable[i,0] = ktable[i - 1,0] - pthick[i-1] 
       for i in range(len(stse)):
#           Calculates the second row of the ktable
           ktable[i,1] = ktable[i+1,0]  
#    Calculates the third, fourth and fifth row of the ktable 
        
   for i in range(r):
       ktable[i,2] = ktable[i,0] - ktable[i,1]
       ktable[i,3] = 0.5 * (((ktable[i,0]) ** 2)- (ktable[i,1]** 2 ))
       ktable[i,4] = (1/3) * ((ktable[i,0] ** 3 )- (ktable[i,1]** 3 ))     
#  Ensures the all the columns but the first one in the last row are set to zero in the ktable.  
   ktable[r-1,1] = 0 
   ktable[r-1,2] = 0 
   ktable[r-1,3] = 0 
   ktable[r-1,4] = 0
#   Using the function input to define which output will be produced. 
#   These could have been joined but this lead to too many outputs to be useful to the user. 
#       Calculate the A, B and D array using the relevant functions
#       A and B are converted to GN(mm^3)/m^2
   A = Aarray(A,d,ktable) #* 10 ** 6
   B = Barray(B,d,ktable) #* 10 ** 3
   D = Darray(D,d,ktable)
   Abar = Abararray(Abar,d1,ktable)
   print(A,B,D)
   print(Abar)
#  This sections join A,B and D together to form the full matrix
#  E = np.concatenate((A,B),axis=1)
#  F = np.concatenate((B,D),axis=1)
#  J = np.concatenate((E, F))
#  A, B, D and the full matrix are returned 
   return(A)
   return(B)
   return(D)
   return(Abar)

       
def Aarray(A,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    a = []
    a1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
#        Indexing   
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values
        a.append(ktable[ty,2] * l[0,0]) 
        a1.append(ktable[ty,2] * l[1,1])
        a2.append(ktable[ty,2] * l[2,2])
        a3.append(ktable[ty,2] * l[1,0])
        a4.append(ktable[ty,2] * l[2,0])
        a5.append(ktable[ty,2] * l[2,1])
#    Sums each list and inputs them into the matrix
    A[0,0] = sum(a)
    A[1,1] = sum(a1)
    A[2,2] = sum(a2)
    A[0,1] = sum(a3)
    A[0,2] = sum(a4)
    A[2,1] = sum(a5)
    A[1,0] = sum(a3)
    A[2,0] = sum(a4)
    A[1,2] = sum(a5)
    return(A)
    
def Barray(B,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    a = []
    a1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
#        Indexing 
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values 
        a.append(ktable[ty,3] * l[0,0]) 
        a1.append(ktable[ty,3] * l[1,1])
        a2.append(ktable[ty,3] * l[2,2])
        a3.append(ktable[ty,3] * l[1,0])
        a4.append(ktable[ty,3] * l[2,0])
        a5.append(ktable[ty,3] * l[2,1])
#    Sums each list and inputs them into the matrix
    B[0,0] = sum(a)
    B[1,1] = sum(a1)
    B[2,2] = sum(a2)
    B[0,1] = sum(a3)
    B[0,2] = sum(a4)
    B[2,1] = sum(a5)
    B[1,0] = sum(a3)
    B[2,0] = sum(a4)
    B[1,2] = sum(a5)
    return(B)
    
def Darray(D,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    k = []
    k1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
#        Indexing 
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values 
        k.append(ktable[ty,4] * l[0,0]) 
        k1.append(ktable[ty,4] * l[1,1])
        a2.append(ktable[ty,4] * l[2,2])
        a3.append(ktable[ty,4] * l[1,0])
        a4.append(ktable[ty,4] * l[2,0])
        a5.append(ktable[ty,4] * l[2,1])
#    Sums each list and inputs them into the matrix
    D[0,0] = sum(k)
    D[1,1] = sum(k1)
    D[2,2] = sum(a2)
    D[0,1] = sum(a3)
    D[0,2] = sum(a4)
    D[2,1] = sum(a5)
    D[1,0] = sum(a3)
    D[2,0] = sum(a4)
    D[1,2] = sum(a5)
    return(D)

def Abararray(Abar,d1,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    a = []
    a1 = []  
    a3 = []
    for i,l in d1.items():
#        Indexing   
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values
        a.append(ktable[ty,2] * l[0,0]) 
        a1.append(ktable[ty,2] * l[1,1])
        a3.append(ktable[ty,2] * l[1,0])
#    Sums each list and inputs them into the matrix
    Abar[0,0] = sum(a)
    Abar[1,1] = sum(a1)
    Abar[0,1] = sum(a3)
    Abar[1,0] = sum(a3)
    return(Abar)    
    
pthick5 = [0.125,0.125,0.125,0.125]
stse5 =  [90, 45, -45, 0] 
#stse5 =  [0, 90, 90, 0] 
print("Question 5 =") 
Function2(54,18,0.28,6,6,6,pthick5,stse5,3)    
    
    
    
    
    