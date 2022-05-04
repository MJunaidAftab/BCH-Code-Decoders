import numpy as np


class Sugiyama_decoder:
    
    #We first define the the input variables in the decoder class
    
    def __init__(self,receivedword, omega,t,k):
        
        self.word = receivedword   #This is the received word from the channel        
        self.alpha = omega         #This is the generator of the multiplicative group of non-zero elements
        self.t = t                 #This is the value of max errors that can be corrected
        self.k = k                 #This is the field over which we are working
        self.n = k.cardinality()-1 #This is the number of elements in the multiplicative group of non-zero elements
        
    #We first compute the error syndrome. This is done by evaluating the receiced word at appropriate powers 
    #of the generator of the multiplicative group of non-zero elements
    
    def SyndromePolynomial(self):
        
        #The output of the function is the syndrome polynomial
        self.syndrome_poly = 0
         
        for j in range(1,2*self.t+1):
            
            Sj = self.word(self.alpha**j)
            self.syndrome_poly += Sj*x^(j-1)    
        
        return self.syndrome_poly
    
    def EuclideanAlgo(self,t):
        
        a = x^(2*t)
        b = self.SyndromePolynomial()
        
        f_i_minus2 = 1
        f_i_minus1 = 0
        
        g_i_minus2 = 0
        g_i_minus1 = 1
        
        r_i_minus2 = a
        r_i_minus1 = b
        
        while r_i_minus1.degree() >= t:
            
            q_i = r_i_minus2//r_i_minus1
            r_i = r_i_minus2%r_i_minus1
            
            f_i = f_i_minus2 - q_i*f_i_minus1
            g_i = g_i_minus2 - q_i*g_i_minus1
            
            f_i_minus2 = f_i_minus1
            f_i_minus1 = f_i
            
            g_i_minus2 = g_i_minus1
            g_i_minus1 = g_i
            
            r_i_minus2 = r_i_minus1
            r_i_minus1 = r_i
            
        LocatorPolynomial = g_i_minus1
            
        return LocatorPolynomial
    
    #Using the error locator polynomial, the error polynomial is now computed. This is done by 
    #evaluating the error locator polynomial of appropriate powers of the generator of the 
    #multiplicative group of non-zero elements. The zerors locate the bit flip errors.
    
    def ErrorPolynomial(self,t):
        
        loc_poly = self.EuclideanAlgo(t)
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
