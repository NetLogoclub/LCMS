import LCMS.LCM_simulation_v3 as LCM3
import arcpy

random_effect = 2.0
work_path = 'g:/clue-sii/clue_simulationtfm_circle5_5f/trans/'+'random_effect='+str(random_effect)+'/natura2030'
yao = LCM3.simulation(past_map='g:/clue-sii/data/clue_simulation/lucc2000_f_m', original_map='g:/clue-sii/data/clue_simulation/lucc2010_f_m',
                      local_factor=['g:/clue-sii/data/dis_road_f_f/class_motorway1.tif', 'g:/clue-sii/data/dis_road_f_f/class_primary1.tif', 'g:/clue-sii/data/dis_road_f_f/class_secondary1.tif',
                                    'g:/clue-sii/data/dis_road_f_f/class_residential1.tif','g:/clue-sii/data/dis_road_f_f/class_services1.tif','g:/clue-sii/data/dis_road_f_f/class_tertiary1.tif',
                                    'g:/clue-sii/data/dis_road_f_f/class_trunk1.tif','g:/clue-sii/data/dis_road_f_f/class_unclassified1.tif',
                                    'g:/clue-sii/data/factor_rasters/dem', 'g:/clue-sii/data/factor_rasters/dis_platform1.tif','g:/clue-sii/data/factor_rasters/slope'],
                      neighborhood=arcpy.sa.NbrCircle(5, 'CELL'), work_path=work_path, step=10, error=1,
                      constant=True, min_alpha=0.0001, max_cycle=10000, random_effect=random_effect)
yao.validating()
yao.predicting('',inNum = 2)
yao.work_path = 'g:/clue-sii/clue_simulationtfm_circle5_5f/trans/'+'random_effect='+str(random_effect)+'/s2030'
yao.predicting('DC_condition.txt',inNum = 2)
