#### define CE class
from scipy.optimize import fsolve
from scipy.optimize import fmin
from scipy.optimize import brenth
import numpy as np
import sys
import trans_matrix as tm
import time
import math



def distance_ce(matrixor,matrixfi):
    matrix1,matrix2 = np.zeros((len(matrixor),len(matrixor))),np.zeros((len(matrixor),len(matrixor)))
    for i in range(len(matrix1)):
        matrix1[i] = matrixor[i]/np.sum(matrixor[i])
        matrix2[i] = matrixfi[i]/np.sum(matrixfi[i])
    dis = 0
    for i in range(len(matrix2)):
        for j in range(len(matrix2[i])):
            if matrix1[i,j] == matrix2[i,j]:
                dis += 0
            elif matrix2[i,j]*matrix1[i,j] == 0:
                dis += 9999999999999999999999
            else:
                dis += matrix1[i,j]*math.log(matrix1[i,j]/matrix2[i,j])
    return dis
                
        
class scenario:
    def __init__(self,trans_matrix):
        self.trans = trans_matrix
        self.trans_p = np.zeros((len(trans_matrix),len(trans_matrix)))
        self.area_or = []
        for i in range(len(trans_matrix)):
            self.area_or.append(np.sum(trans_matrix[i]))
            for j in range(len(trans_matrix[i])):
                self.trans_p[i,j] = self.trans[i,j]/float(self.area_or[i])
        self.condition_1 = []
        self.condition_2 = []
    def add_condition(self,range_num,area):
        if type(range_num) == list:
            self.condition_2.append([range_num,area])
        else:
            self.condition_1.append([range_num,area])

### waiting for complete
    def func0(self,x):
        res = []
        n = len(self.area_or)
        m1 = len(self.condition_1)
        for i in range(m1):
            col = self.condition_1[i][0]
            s1 = 0
            for j in range(n):
                row = j
                s1 += self.trans[row,col]/(x[row]+self.area_or[row]*x[n+i])
            res.append(s1-self.condition_1[i][1])            
        for i in range(n):
            s3 = 0
            for j in range(n):
                row,col = i,j
                control = 0
                for k1 in range(m1):
                    if self.condition_1[k1][0] == col:
                        s3 += self.trans[row,col]/(x[row]+self.area_or[row]*x[n+k1])
                        control = 1
                        break                                                   
                if control == 0:
                    s3 += self.trans[row,col]/x[row]
            res.append(s3-self.area_or[i])
        return res
    
    def func1(self,x):
        res = []
        n = len(self.area_or)
        m2 = len(self.condition_2)

        for i in range(m2):
            row,col = self.condition_2[i][0]
            s1 = self.trans[row,col]/(x[row]+self.area_or[row]*x[n+i])
            res.append(s1-self.condition_2[i][1])

            
        for i in range(n):
            s3 = 0
            for j in range(n):
                row,col = i,j
                control = 0
                for k1 in range(m2):
                    if self.condition_2[k1][0] == [row,col]:
                        s3 += self.trans[row,col]/(x[row]+self.area_or[row]*x[n+k1])
                        control = 1
                        break                                                   
                if control == 0:
                    s3 += self.trans[row,col]/x[row]
            res.append(s3-self.area_or[i])
        return res
    
    def func2(self,x):
        res = []
        n = len(self.area_or)
        m1 = len(self.condition_1)
        m2 = len(self.condition_2)

        for i in range(m1):
            col = self.condition_1[i][0]
            s1 = 0
            for j in range(n):
                row = j
                control = 0
                for k in range(m2):
                    if self.condition_2[k][0] == [row,col]:
                        s1 += (self.trans[row,col]/(x[row]+self.area_or[row]*(x[n+m1+k]+x[n+i])))
                        control = 1
                        break
                if control == 0:
                    s1 += (self.trans[row,col]/(x[row]+self.area_or[row]*x[n+i]))
            res.append(s1-self.condition_1[i][1])

        for i in range(m2):
            s2 = 0            
            control = 0
            row,col = self.condition_2[i][0]
            for j in range(m1):
                if col == self.condition_1[j][0]:
                    s2 = self.trans[row,col]/(x[row]+self.area_or[row]*(x[n+m1+i]+x[n+j]))
                    control = 1
                    break
            if control == 0:
                s2 = self.trans[row,col]/(x[row]+self.area_or[row]*x[n+m1+i])
            res.append(s2 - self.condition_2[i][1])

        for i in range(n):
            s3 = 0
            for j in range(n):
                row,col = i,j
                control = 0
                for k1 in range(m1):
                    for k2 in range(m2):
                        if self.condition_2[k2][0] == [row,col] and self.condition_1[k1][0]==col:
                            s3 += (self.trans[row,col]/(x[row]+self.area_or[row]*(x[n+m1+k2]+x[n+k1])))
                            control = 1
                            break
                    if control == 1:
                        break
                if control == 0:
                    for k1 in range(m1):
                        if self.condition_1[k1][0]==col:
                            s3 += (self.trans[row,col]/(x[row]+self.area_or[row]*x[n+k1]))
                            control = 1
                            break
                if control == 0:
                    for k2 in range(m2):
                        if self.condition_2[k2][0] == [row,col]:
                            s3 += (self.trans[row,col]/(x[row]+self.area_or[row]*x[n+m1+k2]))
                            control = 1
                            break                                                                             
                if control == 0:
                    s3 += (self.trans[row,col]/x[row])
            res.append(s3-self.area_or[i])
        return res
    
    def solve(self):
        n = len(self.area_or)
        m1 = len(self.condition_1)
        m2 = len(self.condition_2)
        x = []
        for i in range(len(self.trans_p)):
            x.append(1)
        for i in range(m1+m2):
            x.append(0)
            
        if m1 >0 and m2>0:
            x = fsolve(self.func2,x)
        elif m1 >0 and m2 == 0:
            x = fsolve(self.func0,x)
        elif m1 ==0 and m2 > 0:
            x = fsolve(self.func1,x)
        else:
            x = x
##        print result
        res = np.zeros((len(self.trans),len(self.trans)))
        for i in range(len(res)):
            for j in range(len(res[i])):
                row,col = i,j
                control = 0
                for k1 in range(m1):
                    for k2 in range(m2):
                        if self.condition_2[k2][0] == [row,col] and self.condition_1[k1][0]==col:
                            res[row,col] = (self.trans_p[row,col]/(x[row]+self.area_or[row]*(x[n+m1+k2]+x[n+k1])))
                            control = 1
                            break
                    if control == 1:
                        break
                if control == 0:
                    for k1 in range(m1):
                        if self.condition_1[k1][0]==col:
                            res[row,col] = (self.trans_p[row,col]/(x[row]+self.area_or[row]*x[n+k1]))
                            control = 1
                            break
                if control == 0:
                    for k2 in range(m2):
                        if self.condition_2[k2][0] == [row,col]:
                            res[row,col] = (self.trans_p[row,col]/(x[row]+self.area_or[row]*x[n+m1+k2]))
                            control = 1
                            break                                                                             
                if control == 0:
                    res[row,col] = (self.trans_p[row,col]/x[row])
                    
                    
        return res,x
##    def func(self,x):
##        n = len(self.area_or)
##        m1 = len(self.condition_1)
##        m2 = len(self.condition_2)
##        res = 0
##        for i in range(n):
##            for j in range(n):
##                res += (x[i*n+j])*math.log((x[i*n+j])/(self.trans_p[i,j]))                
##            res += x[n*n+i]*(sum(x[n*i:n*i+n])-1)
##        for k1 in range(m1):
##            col = self.condition_1[k1][0]
##            Area = 0
##            for c in range(n):
##                Area += x[c*n+col]*self.area_or[c]
##            res += x[n*n+n+k1]*(Area-self.condition_1[k1][1])
##        for k2 in range(m2):
##            row,col = self.condition_2[k2][0]
##            res += x[n*n+n+m1+k2]*(x[n*row+col]*self.area_or[row]-self.condition_2[k2][1])
##        print res
##        return res
##    
##    def solve_1(self):
##        n = len(self.area_or)
##        m1 = len(self.condition_1)
##        m2 = len(self.condition_2)
##        x = []
##        self.trans_p = self.trans_p + 0.0000001
##        for i in range(n):
##            for j in range(n):
##                x.append(self.trans_p[i,j])
##        a = [1]*(n+m1+m2)
##        x = x + a
##        print x
        result = fmin(self.func,x)
        res = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                res[i,j] = result[i*n+j]
        return res
