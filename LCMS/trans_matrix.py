import arcpy
import numpy as np
import os

class trans_matrix:
    "original mao, final map, out_path,create trans matrix for the original map and final map"

#    or_map = ""
#    fi_map = ""
#    out_path = ""
#    out_name = ""

    def __init__(self,or_map,fi_map,nodata_value):
        self.or_map = or_map
        self.fi_map = fi_map
        self.out_path = "c:/python/temple/trans_matrix"
        self.out_name = "trans_matrix.dbf"
        if os.path.exists(self.out_path) == False:
            os.makedirs(self.out_path)
        self.nodata_value = nodata_value
    def add_field(self,file_name,field_name,type_name):
        for i in range(2):
            try:
                arcpy.AddField_management(file_name,field_name,type_name)
                break
            except:
                arcpy.DeleteField_management(file_name,field_name,type_name)
    def get_array(self):
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            arcpy.AddError("spatial analysi liscense is not available")
            pass
        or_list = arcpy.RasterToNumPyArray(self.or_map)
        or_null = self.nodata_value
        fi_list = arcpy.RasterToNumPyArray(self.fi_map)
        fi_null = self.nodata_value
        or_list = list(np.unique(or_list))
        try:
            or_list.remove(or_null)
        except:
            None
        fi_list = list(np.unique(fi_list))
        try:
            fi_list.remove(fi_null)
        except:
            None
        out_array = np.zeros((len(or_list),len(fi_list)))
        for j in range(len(or_list)):
            for i in range(len(fi_list)):
                map_or = arcpy.sa.EqualTo(self.or_map,int(or_list[j]))
                map_fi = arcpy.sa.EqualTo(self.fi_map,int(fi_list[i]))
                ff = arcpy.sa.Times(map_or,map_fi)
                ff_array = arcpy.RasterToNumPyArray(ff)
                num = 0
                for k in ff_array:
                    num += len([x for x in k if x == 1])
#                print j,i,int(num)
                out_array[j,i] = num
        try:
            arcpy.CheckInExtension("Spatial")
        except:
            pass
        return np.int64(out_array)
    def get_array_p(self):
        array = self.get_array()
        out_array = np.eye(len(array))
        for a in range(len(array)):
            out_array[a] = array[a]/float(np.sum(array[a]))
        return out_array
            
        
    def get_matrix(self,out_path=None,out_name=None):
        if out_path != None:
            self.out_path = out_path
        if out_name != None:
            self.out_name = out_name
        dbf = os.path.join(self.out_path,self.out_name)
        or_list = arcpy.RasterToNumPyArray(self.or_map)
        or_null = 0
        fi_list = arcpy.RasterToNumPyArray(self.fi_map)
        fi_null = 0
        or_list = list(np.unique(or_list))
        try:
            or_list.remove(or_null)
        except:
            None
        fi_list = list(np.unique(fi_list))
        try:
            fi_list.remove(fi_null)
        except:
            None
### create table and calculate
        arcpy.env.overwriteOutput = True
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            arcpy.AddError("spatial analysi liscense is not available")
            pass
        arcpy.CreateTable_management(self.out_path,self.out_name)
        c_field = "ori_map"
        self.add_field(dbf,c_field,"LONG")
        for i in fi_list:
            self.add_field(dbf,'c'+str(i),"LONG")
        Rows = arcpy.InsertCursor(dbf)
        for j in or_list:
            row = Rows.newRow()
            row.setValue(c_field,int(j))
            for i in fi_list:
                n_field = 'c'+str(i)
                map_or = arcpy.sa.EqualTo(self.or_map,int(j))
                map_fi = arcpy.sa.EqualTo(self.fi_map,int(i))
                ff = arcpy.sa.Times(map_or,map_fi)
                ff_array = arcpy.RasterToNumPyArray(ff)
                num = 0
                for k in ff_array:
                    num += len([x for x in k if x == 1])
                row.setValue(n_field,num)
            Rows.insertRow(row)
        del row,Rows
        try:
            arcpy.CheckInExtension("Spatial")
        except:
            pass
        return "trans matrix has been successfully calculated in {0}".format(dbf)
