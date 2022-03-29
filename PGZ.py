#!/usr/bin/env python
# coding: utf-8


import numpy as np


class Peterson_Gorenstein_Zierler_decoder:
    
    #We first define the the input variables in the decoder class
    
    def __init__(self,receivedword, omega,t,k):
        
        self.word = receivedword   #This is the received word from the channel        
        self.alpha = omega         #This is the generator of the multiplicative group of non-zero elements
        self.t = t                 #This is the value of max errors that can be corrected
        self.k = k                 #This is the field over which we are working
        self.n = k.cardinality()-1 #This is the number of elements in the multiplicative group of non-zero elements
        
    #We first compute the error syndrome. This is done by evaluating the receiced word at appropriate powers 
    #of the generator of the multiplicative group of non-zero elements
    
    def ErrorSyndrome(self):
        
        #The output of the function is a numpy array comprising the error syndrome
        
        self.syndrome = []
         
        for j in range(1,2*self.t+1):
            
            Sj = self.word(self.alpha**j)
            self.syndrome.append(Sj)
            
        self.syndrome = np.array(self.syndrome)
        
        return self.syndrome
    
    #We now use the error syndrome array to construct the 'syndrome matrix.' This function simply reorganizes the 
    #syndrome array in a syndrome matrix of appropriate size
    
    def SyndromeMatrix(self,mu):
        
        syndrome = self.ErrorSyndrome()
        M = matrix(k, mu, mu, lambda x, y: syndrome[x+y])
                
        return M, syndrome
    
    #Assuming at most t errors have taken place, we now compute the actual number of errors that have taken place.
    #This is done by computing the determinant of the of the syndrome matrix computed above. The syndrome matrix 
    #of appropriate size is constructed until a matrix is found with non-zero determinant. The size of the matrix
    #equals to the number of errors.
    
    def NumberofErrors(self,t):
        
        errors = t
        syndromematrix, syndrome = self.SyndromeMatrix(errors)
        det = syndromematrix.determinant()
        
        while det == 0 and errors > 0:
            
            errors = errors - 1
            syndromematrix, syndrome = self.SyndromeMatrix(errors)
            det = syndromematrix.determinant() 
            
        return errors, syndromematrix, syndrome
    
    #The coefficients of the locator polynomial are now computed. This boils down to solving a linear system 
    #of equations.
    
    def LocatorPolynomial(self,t):
        
        errors, syndromematrix, syndrome = self.NumberofErrors(t)
        
        b = matrix(k, errors, 1, lambda x, y: - syndrome[errors+x])
        error_comp = syndromematrix.inverse()*b
        
        loc_poly = 0
        z = errors
        
        while z >= 0:
            
            if z == 0:
                
                loc_poly = loc_poly[0]+1
                break
            
            loc_poly = loc_poly+error_comp[errors-z]*x^z
            z = z - 1  
        
        return loc_poly
    
    #Using the error locator polynomial, the error polynomial is now computed. This is done by 
    #evaluating the error locator polynomial of appropriate powers of the generator of the 
    #multiplicative group of non-zero elements. The zerors locate the bit flip errors.
    
    def ErrorPolynomial(self,t):
        
        loc_poly = self.LocatorPolynomial(t)
        deg = loc_poly.degree()
        error_poly = 0
        
        for j in range(n):
            
            z = loc_poly(self.alpha**(-j))
            
            if z == 0:
                
                error_poly = error_poly+x^j
                deg = deg - 1
            
            if deg == 0:
                break
        
        return error_poly
    
    #We now decode by returning the the codeword that was sent. If less then t errors indeed took place, the
    #correct codeword will be decoded.
    
    def Codeword(self,t):
        
        error_poly = self.ErrorPolynomial(t)
        return error_poly + self.word


#We now test an example. We work over the finite field GF(2^6), and assume that at most t = 3 can be corrected. In other words, we simulate the decoding 
#algorithm using a distance 7 primitive narrow-sense BCH code. The polynomial received_polynomial is received and is decoded using the codeword
#function in the above class.


prime = 2
dimension = 6
t=3
n = prime**dimension - 1
k = GF((prime,dimension))
omega = k.multiplicative_generator()
R = PolynomialRing(k,'x')
x = R.gen()

received_poly = x^62 + x^60 + x^57 + x^56 + x^54 + x^53 + x^51 + x^49 + x^47 + x^45 + x^43 + x^36 + x^35 + x^34 + x^32 + x^29 +x^12 + x^5


# In[25]:


PGZ = Peterson_Gorenstein_Zierler_decoder(received_poly,omega,3,k)
#PGZ.ErrorSyndrome()
PGZ.Codeword(t)


# In[ ]:




