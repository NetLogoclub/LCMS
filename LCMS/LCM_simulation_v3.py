## this script using the following encoding: utf-8
## used for intergrated simulation of land use change based on logistics and markov

import arcpy
import numpy as np
import random
import os
import sys
import time
import matplotlib.pyplot as plt
import intergrated_ce_logistics_v2 as CLUE
import CEsolve_v3 as CE3
import trans_matrix_divide as TD
import trans_matrix as tm

arcpy.env.overwriteOutput = True

def getMeanOfRaster_noneZero(input_raster):
    array = arcpy.RasterToNumPyArray(input_raster)
    minimum = float(arcpy.GetRasterProperties_management(input_raster,'minimum').getOutput(0))
    maximum = float(arcpy.GetRasterProperties_management(input_raster,'maximum').getOutput(0))
    count = 0.0
    Sum = 0.0
    for i in array:
        for j in i:
            if j <= maximum and j>=minimum and j != 0:
                count = count + 1
                Sum = Sum + j
    return Sum/count
    
def logistics_func(x_array,weight_array):
    weight_array.reshape((len(x_array[0]),1))
    inx = np.dot(x_array,weight_array)
    return 1/(1+np.exp(-inx))


def logistics(factor_rasters,weights,constant,mask):
    arcpy.env.mask = mask
    out_raster = arcpy.sa.Times(mask,0)
    for i in range(len(factor_rasters)):
        min_value = float(arcpy.GetRasterProperties_management(factor_rasters[i],'minimum').getOutput(0))
        max_value = float(arcpy.GetRasterProperties_management(factor_rasters[i],'maximum').getOutput(0))
        raster = arcpy.sa.Divide(arcpy.sa.Minus(factor_rasters[i],min_value),(max_value-min_value))
        out_raster = arcpy.sa.Plus(out_raster,arcpy.sa.Times(raster,weights[i]))
    if constant:
        out_raster = arcpy.sa.Plus(out_raster,weights[i+1])
    out_raster = arcpy.sa.Minus(0,out_raster)
    out_raster = arcpy.sa.Divide(1.0,arcpy.sa.Plus(1,arcpy.sa.Exp(out_raster)))
    return out_raster

def kappa(first_raster,second_raster,nodata_value):
    first_array = arcpy.RasterToNumPyArray(first_raster)
    second_array = arcpy.RasterToNumPyArray(second_raster)
    out_lst = []

    index = list(np.unique(first_array))
    try:
        index.remove(nodata_value)
    except:
        None
    
    num_total = float(np.sum(first_array != nodata_value))
    
    for i in index:
        num1 = np.sum((first_array == second_array)*(first_array != nodata_value))/(num_total)
        num2 = np.sum(first_array == i)/num_total*np.sum(second_array == i)/num_total + (1-np.sum(first_array == i)/num_total)*(1-np.sum(second_array == i)/num_total)
        num = (num1-num2)/(1-num2)
        out_lst.append(num)

    del first_array,second_array
        
    return out_lst

def FRAG(in_raster,nodata_value):
    array = arcpy.RasterToNumPyArray(in_raster)
    index = list(np.unique(array))
    out_lst = []
    try:
        index.remove(nodata_value)
    except:
        None
        
    area,length = [0.0]*len(index),[0.0]*len(index)
    for row in range(len(array)):
        for col in range(len(array[row])):
            for i in range(len(index)):
                if array[row,col] == index[i]:
                    area[i] = area[i] + 1
                    length[i] = length[i] + 4
                    if row > 0 and array[row-1,col] == index[i]:
                        length[i] = length[i] - 1
                    if row < (len(array)-1) and array[row+1,col] == index[i]:
                        length[i] = length[i] - 1
                    if col > 0  and array[row,col-1] == index[i]:
                        length[i] = length[i] - 1
                    if col < (len(array[row])-1) and array[row,col+1] == index[i]:
                        length[i] = length[i] - 1
    for i in range(len(index)):
        out_lst.append(round(float(length[i])/area[i],4))                

    return out_lst
        
        

def ROC_plt(array_x,array_y,weights,num_point = 10):
    R_x = []
    R_y = []

        
    for i in range(num_point+1):
        cut_value = 1-float(i)/num_point
        num_p,num_p_t,num_n,num_n_f = [0.0]*4
        array_res = logistics_func(array_x,weights)>cut_value
        num_p = np.sum(array_y == 1)
        num_n = np.sum(array_y == 0)
        for j in range(len(array_res)):
            if array_y[j][0] == 1 and array_res[j][0] == True:
                num_p_t = num_p_t + 1
            elif array_y[j][0] == 0 and array_res[j][0] == True:
                num_n_f += 1           
        R_x.append(num_n_f/num_n)
        R_y.append(num_p_t/num_p)

    return R_x,R_y

def AUC(array_x,array_y,weights):
    array_score = logistics_func(array_x,weights)
    score_true = []
    score_false = []
    for i in range(len(array_y)):
        if array_y[i] == 1:
            score_true.append(array_score[i])
        else:
            score_false.append(array_score[i])
    AUC_t = 0
    for s_c in score_true:
        AUC_t = AUC_t + np.sum(s_c>np.array(score_false))
    if (len(score_true)*len(score_false)) != 0:
        AUC = float(AUC_t)/(len(score_true)*len(score_false))
    else:
        AUC = 0
    return AUC

def gradAscent(array_x,array_y,min_alpha=0.001):
    m,n = np.shape(array_x)
    weights = np.zeros((n,1))
    array_weights = np.zeros((m,n))
    for k in range(m):
        array_weights[k] = weights.transpose()
        alpha = 1-(1-min_alpha)/m*(k+1)
        error = array_y[k] - 1/(1+np.exp(-np.dot(array_x[k].reshape((1,n)),weights)))
        weights = weights + alpha*array_x[k].reshape((n,1))*error
#        array_weights[k] = weights.transpose()
    return weights,array_weights.transpose()
def arrayToFile(in_array,file_name):
    ff = open(file_name,'w+')
    for row in in_array:
        for col in row:
            ff.write(str(col)+'  ')
        ff.write('\n')
    ff.close()
    print "successfullt translate the input array to "+file_name
            

class simulation:

    def __init__(self,past_map='g:/clue-sii/data/clue_simulation/lucc2000',original_map = 'g:/clue-sii/data/clue_simulation/lucc2010',
                 local_factor = ['g:/clue-sii/data/dis_road_f_f/class_motorway1.tif','g:/clue-sii/data/dis_road_f_f/class_primary1.tif',
                                 'g:/clue-sii/data/dis_road_f_f/class_secondary1.tif','g:/clue-sii/data/factor_rasters/dem',
                                 'g:/clue-sii/data/factor_rasters/dis_platform1.tif','g:/clue-sii/data/factor_rasters/slope'],
                 neighborhood = arcpy.sa.NbrRectangle(10, 10, "CELL"),work_path = 'g:/clue-sii/clue_simulationt',step = 1,error = 1,constant = True,
                 min_alpha = 0.0001, max_cycle = 10000,random_effect = 0.5,prt_zone = '',max_cycle1 = 10000, min_alpha1 = 0.001):
        
        
        self.past_map = past_map
        self.original_map = original_map
        self.local_factor = local_factor
        self.neighborhood = neighborhood
        self.work_path = work_path
        if os.path.exists(work_path) == False:
            os.makedirs(work_path)
        self.error = error
        try:
            self.nodata_value = int(arcpy.Raster(past_map).noDataValue)
        except:
            self.nodata_value = int(past_map.noDataValue)
        self.constant = constant
        self.min_alpha = min_alpha
        self.max_cycle = max_cycle
        self.step = step
        self.weights = []
        self.random_effect = random_effect
        self.prt_zone = prt_zone
        self.min_alpha1 = min_alpha1
        self.max_cycle1 = max_cycle1       

## information writing
        f = open(os.path.join(work_path,'setting.txt'),'w+')
        f.write('past_map: '+ str(past_map)+'\n')
        f.write('original_map: '+ str(original_map)+'\n')
        f.write('local_factor: '+ str(local_factor)+'\n')
        f.write('neighborhood: '+ str(neighborhood)+'\n')
        f.write('work_path: '+ str(work_path)+'\n')
        f.write('error: '+ str(error)+'\n')
        f.write('constant: '+ str(constant)+'\n')
        f.write('min_alpha: '+ str(min_alpha)+'\n')
        f.write('max_cycle: '+ str(max_cycle)+'\n')
        f.write('step: '+ str(step)+'\n')
        f.write('random_effect: '+ str(random_effect)+'\n')
        f.write('protect zone: '+ str(prt_zone)+'\n')
        f.close()
       

    def __neighbor_factor(self,raster):
        input_raster = raster
        neighborhood = self.neighborhood
        out_folder = os.path.join(self.work_path,'neighborhood_'+(os.path.split(raster)[1]).split('.')[0])
        out_list = []
        if os.path.exists(out_folder) == False:
            os.makedirs(out_folder)
        input_array = arcpy.RasterToNumPyArray(input_raster)
    #    des = arcpy.Describe(input_raster+"/Band_1")
    ###define none_value to 0 or others
        none_value = self.nodata_value
        input_class = list(np.unique(input_array))
        input_class.remove(self.nodata_value)

        num_total = np.sum(input_array != none_value)
    ### create dbf file to record average neighborhood enrichment for each class
##        arcpy.CreateTable_management(out_folder,'average.dbf')
##        out_dbf = os.path.join(out_folder,'average.dbf')
##        arcpy.AddField_management(out_dbf,'class_n','LONG')
##        for field in input_class:
##            if field == none_value:
##                continue
##            arcpy.AddField_management(out_dbf,'to_'+str(field),'FLOAT')
##        rows = arcpy.InsertCursor(out_dbf)
    ###calculate total num of pixel

        for c in input_class:
            num_class = 0
            for j in range(len(input_array)):
                num_class = num_class + len([x for x in input_array[j] if x == c])
            arcpy.AddMessage("start processing for class: " + str(c))
            print "start processing for class: " + str(c)
            class_raster = arcpy.sa.EqualTo(input_raster,int(c))
            class_c = arcpy.sa.FocalStatistics(class_raster,neighborhood,"MEAN")
            class_cf = arcpy.sa.Divide(class_c,(float(num_class)/num_total))
            out_file = os.path.join(out_folder,"E_"+str(c))
            out_list.append(out_file)
            class_cf.save(out_file)
    ###calculate and record average neighborhood enrichement of each class for class c
##                arcpy.AddMessage("start calculate neighborhood enrichment of class_" + str(c)+" to each class")
##                row = rows.newRow()
##                row.class_n = int(c) 
##                for cl in input_class:
##                    if cl == none_value:
##                        continue
##                    arcpy.AddMessage("start calculate neighborhood enrichment of class_" + str(c)+" to class: "+str(cl))
##                    c_name = 'to_'+str(cl)
##                    c_value = 0.0
##                    c_count = 0.0
##                    cl_raster = arcpy.sa.EqualTo(arcpy.Raster(input_raster),int(cl))
##                    av_raster = arcpy.sa.Times(cl_raster,class_cf)
##     #               arcpy.AddMessage(111)
##                    no_value_c = arcpy.Describe(av_raster).noDataValue
##                    cl_array = arcpy.RasterToNumPyArray(av_raster)
##    #                arcpy.AddMessage(111)
##                    for av_c in range(len(input_array)):
##                        c_count = c_count + len([x for x in input_array[av_c] if x == cl])
##                    for av in range(len(cl_array)):
##                        c_value = c_value + sum([x for x in cl_array[av] if x != no_value_c])
##                    row.setValue(c_name,c_value/c_count)
##            rows.insertRow(row)
##        del rows
        return out_list

    def probability(self,raster):
        plt.close()
        arcpy.CheckOutExtension('Spatial')
        lucc_raster = raster
        factor_rasters = self.__neighbor_factor(raster)+self.local_factor
        C_value = self.constant
        out_folder = os.path.join(self.work_path,os.path.join('logistics_'+(os.path.split(raster)[1]).split('.')[0]))
        min_alpha = self.min_alpha
        maxCycles = self.max_cycle
        out_list = [[]]

    ## processing the parameters
        if os.path.exists(out_folder) == False:
            os.makedirs(out_folder)
    ### init array_x,array_weight and array_y and set them to zero note that the nodata value must be equal to 0(y) 1(x)
        none_value = self.nodata_value
        array_lucc_raster = arcpy.RasterToNumPyArray(lucc_raster)
        array_lucc = np.delete(np.unique(array_lucc_raster),none_value)
        num_type = len(array_lucc)
        array_y = np.zeros((maxCycles,num_type))
        if C_value:
            num_factors = len(factor_rasters)+1
        else:
            num_factors = len(factor_rasters)
        array_x = np.ones((maxCycles,num_factors))
        array_weights = np.zeros((num_type,num_factors))

        text_fields = []
        for text in factor_rasters:
            text_fields.append(os.path.split(text)[1])
        if C_value:
            text_fields.append('constant value')

    ### randomly get the value for array_y
        array_point = np.zeros((maxCycles,2))
        min_x = float(arcpy.GetRasterProperties_management(lucc_raster,'LEFT').getOutput(0))
        max_x = float(arcpy.GetRasterProperties_management(lucc_raster,'RIGHT').getOutput(0))
        min_y = float(arcpy.GetRasterProperties_management(lucc_raster,'bottom').getOutput(0))
        max_y = float(arcpy.GetRasterProperties_management(lucc_raster,'top').getOutput(0))
        cell_size_x =  float(arcpy.GetRasterProperties_management(lucc_raster,'cellsizex').getOutput(0))
        cell_size_y =  float(arcpy.GetRasterProperties_management(lucc_raster,'cellsizey').getOutput(0))    
        for i in range(maxCycles):
            while True:
                point_x = random.uniform(min_x,max_x)
                point_y = random.uniform(min_y,max_y)
                ncol = int((point_x-min_x)/cell_size_x)
                nrow = int((max_y-point_y)/cell_size_y)
                value = array_lucc_raster[nrow,ncol]
                if value != none_value:
                    array_y[i] = array_lucc == value
                    array_point[i] = np.array([point_x,point_y])                
                    break
          
    ### randomly get the value for array_x
        for j in range(len(factor_rasters)):
            array_f = arcpy.RasterToNumPyArray(factor_rasters[j])
            min_x = float(arcpy.GetRasterProperties_management(factor_rasters[j],'LEFT').getOutput(0))
            max_x = float(arcpy.GetRasterProperties_management(factor_rasters[j],'RIGHT').getOutput(0))
            min_y = float(arcpy.GetRasterProperties_management(factor_rasters[j],'bottom').getOutput(0))
            max_y = float(arcpy.GetRasterProperties_management(factor_rasters[j],'top').getOutput(0))
            cell_size_x =  float(arcpy.GetRasterProperties_management(factor_rasters[j],'cellsizex').getOutput(0))
            cell_size_y =  float(arcpy.GetRasterProperties_management(factor_rasters[j],'cellsizey').getOutput(0))

            max_value =  float(arcpy.GetRasterProperties_management(factor_rasters[j],'maximum').getOutput(0))
            min_value =  float(arcpy.GetRasterProperties_management(factor_rasters[j],'minimum').getOutput(0))
            for k in range(maxCycles):
                point_x = array_point[k,0]
                point_y = array_point[k,1]
                ncol = int((point_x-min_x)/cell_size_x)
                nrow = int((max_y-point_y)/cell_size_y)
                value = array_f[nrow,ncol]

                if value >= min_value and value <= max_value and min_value != max_value:
                    array_x[k,j] = float(value-min_value)/(max_value-min_value)
                else:
                    array_x[k,j] = 0.0

    ### get,write and represent array_weights through gradAscent

        f = open(os.path.join(out_folder,'weights.txt'),'w')
        
        f.write('min_alpha: '+str(min_alpha)+'  maxCycles:  '+str(maxCycles)+'\n')

        f.write('landuse_type  ')

        for text in text_fields:
            f.write(text+'  ')

        f.write('AUC')

        f.write('\n')

        array_y = array_y.transpose()
        arcpy.env.overwriteOutput = True
        color = np.zeros((len(text_fields),3))
        for i in range(len(text_fields)):
            color[i] = np.array([random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)])
        for i in range(num_type):
            result_name = os.path.join(out_folder,'P_'+str(array_lucc[i])+'.tif')
            arcpy.AddMessage('start processing for class: '+str(array_lucc[i]))
            print 'start processing for class: '+str(array_lucc[i])
            grad = gradAscent(array_x,array_y[i],min_alpha)
            array_weights[i] = grad[0].transpose()
            auc = AUC(array_x,array_y[i],array_weights[i])
            f.write(str(array_lucc[i])+'  ')

            for weight in array_weights[i]:
                f.write(str(weight)+'  ')
            f.write(str(auc)+'\n')
    ### represent the iteration of variables

            arcpy.AddMessage('start to record and represent')
            if num_type%2 == 0:
                f_row = num_type/2
            else:
                f_row = num_type/2+1
            
            plt.title("iterations-weights_"+str(array_lucc[i]).zfill(2)+'(AUC = '+str(auc)[:5]+'; n = '+str(int(np.sum(array_y[i])))+')')
            

            for w_x in range(len(text_fields)):
                w_x_r,w_x_g,w_x_b = color[w_x]
                plt.plot(range(1,len(grad[1][w_x])+1),grad[1][w_x],'-',color =(w_x_r,w_x_g,w_x_b))
            plt.savefig(os.path.join(out_folder,"iterations-weights_"+str(array_lucc[i]).zfill(2)+'.jpg'))
            plt.close()
            

            plt.title("weights_"+str(array_lucc[i]).zfill(2)+'(AUC = '+str(auc)[:5]+'; n = '+str(int(np.sum(array_y[i])))+')')
            plt.bar(range(1,num_factors+1),array_weights[i],width= 0.8)
            xticks = range(1,num_factors+1)
            for x in range(num_factors):
                xticks[x] = 'x'+str(xticks[x])
            plt.xticks(np.array(range(1,num_factors+1))+0.4,xticks)
            arcpy.AddMessage(text_fields)
    ### save result rasters
            raster = logistics(factor_rasters,array_weights[i],C_value,lucc_raster)
            raster.save(result_name)
            out_list[0].append(result_name)
            plt.savefig(os.path.join(out_folder,"weights_"+str(array_lucc[i]).zfill(2)+'.jpg'))
            plt.close()
            
        out_list.append(array_weights)               
        f.close()
                
        del array_lucc
        return out_list

    def trans_matrix(self,control=True):
        arcpy.AddMessage("start calculating the trans_matrix between {0} to {1}".format(self.past_map,self.original_map))
        
        if control:
            trans = tm.trans_matrix(self.past_map,self.original_map,self.nodata_value).get_array_p()
            print np.round(trans,4)
        else:
            trans = tm.trans_matrix(self.past_map,self.original_map,self.nodata_value).get_array()
            trans_t = trans.transpose()
            trans_p = tm.trans_matrix(self.past_map,self.original_map,self.nodata_value).get_array_p()
            for i in range(len(trans)):
                trans[i] = np.sum(trans_t[i])*trans_p[i]
            print trans

        
        return trans
        
    def validating(self,raster_name = 'validating',kappa_step = 3,max_radius = 30):
        print "start time of validating: "+ time.ctime()
        trans = self.trans_matrix()
        td = TD.matrix_divide(trans,self.step,['A']*len(trans),self.max_cycle1,self.min_alpha1)
        trans_p = td.p

        print np.round(trans_p,4)
        
        PRB = self.probability(self.past_map)      
        self.weights = PRB[1]
        print self.weights
        
        
        for i in range(self.step):
            arcpy.AddMessage("start processing probability maps for categories with interasion-{0}".format(i+1))
            print "start processing probability maps for categories with interasion-{0}".format(i+1)
            if os.path.exists(os.path.join(self.work_path,raster_name+'_'+str(i+1))) == False:
                os.makedirs(os.path.join(self.work_path,raster_name+'_'+str(i+1)))
            if i == 0:
                probabilities = PRB[0]
                arcpy.AddMessage("finised processing of probability maps")
                arcpy.AddMessage("spatial allocation... ...")
                clue = CLUE.CLUE(land_map = self.past_map,path = os.path.join(self.work_path,raster_name+'_'+str(i+1)),simulate_map = raster_name+str(i+1)+'.tif',
                                 probabilities = probabilities,matrix = trans_p,prt_zone = self.prt_zone)
                out_file = clue.allocate(error_num = self.error,random_effect = self.random_effect)
                arcpy.AddMessage("succefully allocated for all categories")
            else:
                nrb_list = self.__neighbor_factor(out_file)
                probabilities = []
                w = 1
                for weight in self.weights:
                    prb = logistics(nrb_list+self.local_factor,weight,True,out_file)
                    prb_path = os.path.join(os.path.join(self.work_path,raster_name+'_'+str(i+1)),'P_'+str(w)+'.tif')
                    prb.save(prb_path)
                    probabilities.append(prb_path)
                    w += 1
                arcpy.AddMessage("finised processing of probability maps")
                arcpy.AddMessage("spatial allocation... ...")
                clue = CLUE.CLUE(land_map = out_file,path = os.path.join(self.work_path,raster_name+'_'+str(i+1)),simulate_map = raster_name+str(i+1)+'.tif',
                                 probabilities = probabilities,matrix = trans_p,prt_zone = self.prt_zone)

                out_file = clue.allocate(error_num = self.error,random_effect = self.random_effect)
                arcpy.AddMessage("succefully allocated for all categories")
                

## resulting recording
        file_validate = os.path.join(self.work_path,'validating_result')
        os.makedirs(file_validate)        
        array = arcpy.RasterToNumPyArray(self.original_map)
        index = list(np.unique(array))
        try:
            index.remove(self.nodata_value)
        except:
            None
## FRAG calculating
        f0 = open(os.path.join(file_validate,'FRAG.txt'),'w+')
        f0.write('radius'.center(8))
        for x in index:
            f0.write(('Pred_'+str(x)).center(8))
            f0.write(('Null_'+str(x)).center(8))
            f0.write(('Real_'+str(x)).center(8))
        f0.write('\n')
        FRAG_array = np.ones((int(max_radius/kappa_step)+1,len(trans)))
        FRAG_array0 = np.ones((int(max_radius/kappa_step)+1,len(trans)))
        FRAG_array1 = np.ones((int(max_radius/kappa_step)+1,len(trans)))
        FRAG_array[0] = FRAG(out_file,self.nodata_value)
        FRAG_array0[0] = FRAG(self.past_map,self.nodata_value)
        FRAG_array1[0] = FRAG(self.original_map,self.nodata_value)
        f0.write('0'.center(8))
        for j in range(len(trans)):
            f0.write(str(FRAG_array[0,j]).center(8))
            f0.write(str(FRAG_array0[0,j]).center(8))
            f0.write(str(FRAG_array1[0,j]).center(8))
        f0.write('\n')
## kappa statistics with multi fuzzy radious
        f = open(os.path.join(file_validate,'kappa.txt'),'w+')
        f.write('radius'.center(6))
        for x in index:
            f.write(('P_'+str(x)).center(6))
            f.write(('N_'+str(x)).center(6))
        f.write('\n')     
        kappa_array = np.ones((int(max_radius/kappa_step)+1,len(trans)))
        kappa_array0 = np.ones((int(max_radius/kappa_step)+1,len(trans)))
        arcpy.AddMessage("kappa calculating for validating map")
        kappa_array[0] = kappa(self.original_map,out_file,self.nodata_value)
        kappa_array0[0] = kappa(self.original_map,self.past_map,self.nodata_value)
        f.write('0'.center(6))
        for j in range(len(trans)):
            f.write(str(kappa_array[0,j])[:4].center(6))
            f.write(str(kappa_array0[0,j])[:4].center(6))
        f.write('\n')
## multi statistic
        for r in range(kappa_step,max_radius+1,kappa_step):
            k = r/kappa_step
            arcpy.AddMessage('kappa validating radius = '+str(r))
            print 'kappa validating radius = '+str(r)

            first_raster = arcpy.sa.FocalStatistics(self.original_map,arcpy.sa.NbrCircle(r,'CELL'),'MAJORITY','DATA')
            second_raster = arcpy.sa.FocalStatistics(out_file,arcpy.sa.NbrCircle(r,'CELL'),'MAJORITY','DATA')
            third_raster = arcpy.sa.FocalStatistics(self.past_map,arcpy.sa.NbrCircle(r,'CELL'),'MAJORITY','DATA')
            nodata_value = int(first_raster.noDataValue)

            try:
                kappa_array[k] = kappa(first_raster,second_raster,nodata_value)
                kappa_array0[k] = kappa(first_raster,third_raster,nodata_value)
                f.write(str(r).center(6))
                for j in range(len(trans)):
                    f.write(str(kappa_array[k,j])[:4].center(6))
                    f.write(str(kappa_array0[k,j])[:4].center(6))
                f.write('\n')
            except:
                kappa_array = kappa_array[:k]
                kappa_array0 = kappa_array0[:k]
                break
            
            try:
                FRAG_array[k] = FRAG(second_raster,nodata_value)
                FRAG_array0[k] = FRAG(third_raster,nodata_value)
                FRAG_array1[k] = FRAG(first_raster,nodata_value)
                f0.write(str(r).center(8))
                for j in range(len(trans)):
                    f0.write(str(FRAG_array[k,j]).center(8))
                    f0.write(str(FRAG_array0[k,j]).center(8))
                    f0.write(str(FRAG_array1[k,j]).center(8))
                f0.write('\n')
            except:
                FRAG_array = FRAG_array[:k]
                FRAG_array0 = FRAG_array0[:k]
                FRAG_array1 = FRAG_array1[:k]
                break
        f.close()
        f0.close()
        del array
        
        kappa_array = kappa_array.transpose()
        kappa_array0 = kappa_array0.transpose()
        FRAG_array = FRAG_array.transpose()
        FRAG_array0 = FRAG_array0.transpose()
        FRAG_array1 = FRAG_array1.transpose()
        for k in range(len(kappa_array)):
            plt.plot(kappa_array[k],'-',label='LCM model_'+str(index[k]))
            plt.plot(kappa_array0[k],'-',label='NULL model_'+str(index[k]))
            plt.legend(loc=4)
            plt.title('kappa_'+str(index[k])+' radius_step= '+str(kappa_step))
            plt.savefig(os.path.join(file_validate,'kappa_'+str(index[k])+'.png'))
            plt.close()
            
            plt.plot(FRAG_array[k],'-',label='LCM model_'+str(index[k]))
            plt.plot(FRAG_array0[k],'-',label='past Map_'+str(index[k]))
            plt.plot(FRAG_array1[k],'-',label='real MAP_'+str(index[k]))
            plt.legend(loc=0)
            plt.title('FRAG_'+str(index[k])+' radius_step= '+str(kappa_step))
            plt.savefig(os.path.join(file_validate,'FRAG_'+str(index[k])+'.png'))
            plt.close()        
        return out_file

    def predicting(self,condition_file,raster_name = 'predicting',inNum = 1):
        arcpy.AddMessage("start time of predicting: "+ time.ctime())
        if os.path.exists(condition_file):
            con_file = open(condition_file).readlines()
            con = [x.replace('\n','').split() for x in con_file]
            trans = self.trans_matrix(False)
            trans_p = np.array([x/float(np.sum(x)) for x in trans])
            trans_p = TD.cf(trans_p,inNum)
            for i in range(len(trans)):
                trans[i] = trans_p[i]*np.sum(trans[i])            
            print trans
            scean = CE3.scenario(trans)
            for c in con:
                if len(c) == 2:
                    scean.add_condition(int(c[0]),int(c[1]))
                if len(c) == 3:
                    scean.add_condition([int(c[0]),int(c[1])],int(c[2]))
                print 'add condition {0} to trans_matrix'.format(c)
            trans = scean.solve()[0]
            print trans
        else:
            trans = self.trans_matrix()
            trans = TD.cf(trans,inNum)
            arcpy.AddMessage('no condition added to the prediction')
        arrayToFile(trans,os.path.join(self.work_path,'trans_matrix.txt'))
        
        td = TD.matrix_divide(trans,self.step,['A']*len(trans),self.max_cycle1,self.min_alpha1)
        trans_p = td.p
        print trans_p
        
        PRB = self.probability(self.original_map)
        self.weights = PRB[1]
        print self.weights
        
        
        for i in range(self.step):
            arcpy.AddMessage("start processing probability maps for categories with interasion-{0}".format(i+1))
            print "start processing probability maps for categories with interasion-{0}".format(i+1)
            if os.path.exists(os.path.join(self.work_path,raster_name+'_'+str(i+1))) == False:
                os.makedirs(os.path.join(self.work_path,raster_name+'_'+str(i+1)))
            if i == 0:
                probabilities = PRB[0]
                arcpy.AddMessage("finised processing of probability maps")
                arcpy.AddMessage("spatial allocation... ...")
                clue = CLUE.CLUE(land_map = self.original_map,path = os.path.join(self.work_path,raster_name+'_'+str(i+1)),simulate_map = raster_name+str(i+1)+'.tif',
                                 probabilities = probabilities,matrix = trans_p,prt_zone = self.prt_zone)
                out_file = clue.allocate(error_num = self.error,random_effect = self.random_effect)
                arcpy.AddMessage("succefully allocated for all categories")
            else:
                nrb_list = self.__neighbor_factor(out_file)
                probabilities = []
                w = 1
                for weight in self.weights:
                    prb = logistics(nrb_list+self.local_factor,weight,True,out_file)
                    prb_path = os.path.join(os.path.join(self.work_path,raster_name+'_'+str(i+1)),'P_'+str(w)+'.tif')
                    prb.save(prb_path)
                    print prb_path
                    probabilities.append(prb_path)
                    w += 1
                arcpy.AddMessage("finised processing of probability maps")
                arcpy.AddMessage("spatial allocation... ...")
                clue = CLUE.CLUE(land_map = out_file,path = os.path.join(self.work_path,raster_name+'_'+str(i+1)),simulate_map = raster_name+str(i+1)+'.tif',
                                 probabilities = probabilities,matrix = trans_p,prt_zone = self.prt_zone)

                out_file = clue.allocate(error_num = self.error,random_effect = self.random_effect)
                arcpy.AddMessage("succefully allocated for all categories")
                
        return out_file        
    

