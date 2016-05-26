import sys
sys.path.append("g:/python/arcpy")
import trans_matrix as tm
import numpy as np
import arcpy
import os
import matplotlib.pyplot as plt
import random
##read maps and define the tags
def cf(array,t):
    step_array = np.zeros((len(array),len(array)))
    if t == 0:
        return 1
    if t == 1:
        step_array = array
        return step_array
    elif t > 1:
        out = np.array(array)
        for i in range(t-1):
            step = np.array(out)
            for j in range(len(array)):
                for k in range(len(array)):
                    num = 0.0
                    for m in range(len(array)):
                        num = num + step[j,m]*array[m,k]
                    out[j,k] = num
        return out
    else:
        return "error"

class matrix_divide:
    "create land use transfor probablity for per year"
##    or_map = "g:/CLUE-sII/data/markov/lucc1990"
##    fi_map= "g:/CLUE-sII/data/markov/lucc2010"
##    tag = ['arable land','forest','building','other land']
##    year = 20
##    inter_num = 1000
##    alpha = 0.001
    def __init__(self,matrix,year,tag,inter_num = 1000,alpha= 0.001):
        self.year = year
        self.f = np.float32(matrix)
        self.tag = tag
        self.sum_or,self.sum_fi = [],[]
        self.as_list = []
        self.inter_num = inter_num
        self.alpha = alpha
        for x in self.f:
            self.sum_or.append(np.sum(x))
        for x in self.f.transpose():
            self.sum_fi.append(np.sum(x))
        for i in range(len(self.f)):
            self.f[i] = self.f[i]/np.sum(self.f[i])
        self.p = np.eye(len(self.f))
        for i in range(self.inter_num):
            ascent = (cf(self.p,self.year)-self.f)
            self.p = self.p - self.alpha*ascent
            self.as_list.append(ascent*ascent)
## show processing of iteration
        for i in range(len(self.tag)):
            for j in range(len(self.tag)):
                b = random.uniform(0,1),random.uniform(0,1)
                color = [float(i)/len(self.tag),float(j)/len(self.tag),0]
                ll = []
                for k in range(self.inter_num):
                    ll.append(self.as_list[k][i,j])
                plt.plot(range(self.inter_num),ll,label=self.tag[i]+' to '+self.tag[j])
        plt.legend()

    def area_show(self):
        ##calculte and show area_p for simulation    
        self.area_p = np.zeros((self.year+1,len(self.tag)))
        for i in range(self.year+1):
            self.area_p[i] = np.dot(self.sum_or,cf(self.p,i))/float(sum(self.sum_or))
        self.area_p = self.area_p.transpose()
        bottom = 0
        rec = []
        for i in range(len(tag)):
            r,g,b = random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)
            color = [r,g,b]
            rec.append(plt.Rectangle((0,0),1,1,fc = color))
            top = bottom + self.area_p[i]
#            plt.fill_between(range(self.year+1),bottom,top,facecolor = color)
            plt.fill_between(range(self.year),bottom,top,facecolor = color)
            bottom = top
        plt.legend(rec,self.tag)
        plt.show()
    def area_simulate(self,simulate_year):
        area_p = np.zeros((simulate_year+1,len(self.tag)))
        area_p[0] = self.sum_fi
        for i in range(1,simulate_year+1):
            area_p[i] = np.dot(self.sum_fi,cf(self.p,i))
        area_p = area_p.transpose()
        bottom = 0
        rec = []
        for i in range(len(tag)):
            r,g,b = random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)
            color = [r,g,b]
            rec.append(plt.Rectangle((0,0),1,1,fc = color))
            top = bottom + area_p[i]
#            plt.fill_between(range(self.year+1),bottom,top,facecolor = color)
            plt.fill_between(range(simulate_year+1),bottom,top,facecolor = color)
            bottom = top
        plt.legend(rec,self.tag)
        plt.show()
        return area_p
