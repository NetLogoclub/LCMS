## this script using the following encoding: utf-8
## used for intergrated simulation of CE-markov and spatial logistics
import arcpy
import numpy as np
import random
import os
import sys
import time
import matplotlib.pyplot as plt

def random_list(inlist):
    out_list = []
    for i in range(len(inlist)):
        num = random.randint(0,len(inlist)-1)
        out_list.append(inlist[num])
        inlist.remove(inlist[num])
    return out_list
        
class CLUE:
    def __init__(self,land_map ='G:/CLUE-sII/clue_simulation_s5215/validating_1/validating1.tif',path ='G:/CLUE-sII/clue_simulation_s5215/validating_2',simulate_map = 'simulate_result4311.tif',
                 probabilities=['P_1.tif','P_2.tif','P_3.tif','P_4.tif','P_5.tif','P_6.tif'],matrix=np.eye(5),prt_zone = ''):
        self.land_map = land_map
### check whether the out_raster is existing
        i = 1
        while True:
            if os.path.exists(os.path.join(path,simulate_map)):
                simulate_map = simulate_map + str(i)
                i+= 1
            else:
                break
        self.simulate_map = os.path.join(path,simulate_map)
        self.path = path
        if os.path.exists(self.path) == False:
            os.makedirs(self.path)

        self.probabilities = [os.path.join(path,x) for x in probabilities]
        self.matrix = matrix * 1
        self.prt_zone = prt_zone

        
        if len(self.probabilities) != len(self.matrix):
            print "the number of probabilities({0}) can`t match the number of matrix({1})".format(len(self.probabilities),len(self.matrix))
            sys.exit(0)
            

    def allocate(self,nodataValue = 0,random_effect = 0.5,error_num = 10,max_interation = 100):
        x0 = float(arcpy.GetRasterProperties_management(self.land_map,'LEFT').getOutput(0))
        y0 = float(arcpy.GetRasterProperties_management(self.land_map,'BOTTOM').getOutput(0))
        size_land = float(arcpy.GetRasterProperties_management(self.land_map,'CELLSIZEX').getOutput(0))

        

        land_array = arcpy.RasterToNumPyArray(self.land_map)
        land_index = list(np.unique(land_array))
        try:
            land_index.remove(nodataValue)
        except:
            None
##        print land_index
## add prt_zone to probabilities
        for i in range(len(land_index)):
            try:
                self.probabilities[i] = arcpy.sa.Con(self.prt_zone,arcpy.sa.Con(arcpy.sa.EqualTo(self.land_map,land_index[i]),1,0),self.probabilities[i])
            except:
                None
##            
        or_area = []
        for a in range(len(land_index)):
            area = np.sum(land_array==land_index[a])
            or_area.append(area)
            self.matrix[a] = self.matrix[a]*area

        print self.matrix
        arcpy.AddMessage('allocate matrix:')
        arcpy.AddMessage(self.matrix)
##        print or_area,self.matrix
        
            
        probability_array = [arcpy.RasterToNumPyArray(x) for x in self.probabilities]
### add random factor
        for pp in probability_array:
            for i in range(len(pp)):
                for j in range(len(pp[i])):
                    pp[i,j] = (pp[i,j]+random_effect*random.random())/(1+random_effect)
## start processing            
        out_array = np.int8(np.zeros((len(land_array),len(land_array[0]))))


##        print trans_area,or_area,trans_array
        
## iteration for every categories
        for i in range(len(land_index)):
            land_index1 = [x for x in land_index]
            land_index1.remove(land_index[i])
            land_index1 = [land_index[i]]+land_index1
            for j in range(len(land_index1)-1):
                print i,j,self.matrix[i,j]
                list_error,list_break = [],[]
                break_value = 0.5
                interation_num = 1
                while True:
                    list_break.append(break_value)
                    num = np.sum(((land_array==land_index[i])*(out_array==0)*(probability_array[j]>break_value)))
                    error = (num-self.matrix[i,j])
                    list_error.append(error)
                    if interation_num > max_interation:
                        plt.plot(list_error,label = 'iteration error')
                        plt.plot(list_break,label = 'break_value')
                        plt.legend()
                        plt.title('allocation {0} to {1}, error ={2}'.format(land_index[i],land_index1[j],error))
                        print "unsuccesfully allocate category {0} to {1} with error_rate of {2}, pleast close the figure to continue".format(land_index[i],land_index1[j],error)
                        out_array = out_array+(((land_array==land_index[i])*(out_array==0)*(probability_array[j]>break_value))*land_index[j])
                        plt.savefig(os.path.join(self.path,str(land_index[i])+'to'+str(land_index1[j])))
                        plt.close()
                        break
                    try:
                        if error > error_num:
                            interation_num += 1
                            print break_value,num,self.matrix[i,j]
                            break_value += 1.0/(2**interation_num)                                                                        
                        elif error < (-error_num):
                            interation_num += 1
                            print break_value,num,self.matrix[i,j]
                            break_value -= 1.0/(2**interation_num)                        
                        else:
                            plt.plot(list_error,label = 'iteration error')
                            plt.plot(list_break,label = 'break_value')
                            plt.legend()
                            plt.title('allocation {0} to {1}, error ={2}'.format(land_index[i],land_index1[j],error))
                            print "succesfully allocate category {0} to {1} with error_rate of {2}, pleast close the figure to continue".format(land_index[i],land_index1[j],error)
                            plt.savefig(os.path.join(self.path,str(land_index[i])+'to'+str(land_index1[j])))
                            out_array = out_array+(((land_array==land_index[i])*(out_array==0)*(probability_array[j]>break_value))*land_index[j])
                            plt.close()
                            break
                    except:
                        plt.plot(list_error,label = 'iteration error')
                        plt.plot(list_break,label = 'break_value')
                        plt.legend()
                        plt.title('allocation {0} to {1}, error ={2}'.format(land_index[i],land_index1[j],error))
                        print "unsuccesfully allocate category {0} to {1} with error_rate of {2}, pleast close the figure to continue".format(land_index[i],land_index1[j],error)
                        out_array = out_array+(((land_array==land_index[i])*(out_array==0)*(probability_array[j]>break_value))*land_index[j])
                        plt.savefig(os.path.join(self.path,str(land_index[i])+'to'+str(land_index1[j])))
                        plt.close()
                        break
            j += 1
            num = np.sum(((land_array==land_index[i])*(out_array==0)*(probability_array[j]>0)))
            error = (num-self.matrix[i,j])/np.sum(self.matrix[i])
            print "succesfully allocate category {0} to {1} with error_rate of {2}".format(land_index[i],land_index1[j],error)
            out_array = out_array+(((land_array==land_index[i])*(out_array==0)*(probability_array[j]>0))*land_index[j])

        out_raster = arcpy.NumPyArrayToRaster(out_array,arcpy.Point(x0,y0),size_land,size_land,nodataValue)
        out_raster.save(self.simulate_map)
        arcpy.DefineProjection_management(self.simulate_map,self.land_map)
        del probability_array[:]
        return self.simulate_map
