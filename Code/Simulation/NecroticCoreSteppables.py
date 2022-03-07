from PySteppables import *
import CompuCell
import sys
from copy import deepcopy
import random
from math import *
from PlayerPython import *
import os
import numpy as np
from PySteppablesExamples import MitosisSteppableBase

import time
import os.path



V0=16.0
S0=16.0
LBD_V0=15.0
LBD_S0=15.0
V0_nec = 4.0 
S0_nec = 4.0 
LBD_NV0 = 10.0  
LBD_NS0 = 5.0 
incvol=0.2
decvol=0.01
volmaxmit=32.0
Svolmaxmit=32.0
ktgs=4.0
maxdiv=8.0
probstem=0.2
probmut=0.1
PGrThr0=0.032
SGrThr0=0.032
PNeThr0=102.0
QNeThr0=204.0
SNeThr0=408.0
QSNeThr0=916.0
QPThr0=79.5
QSSThr0=79.5
GluD=0.0032
cadhstdev=2.0
def fiber_mapper(x,y,z,map_method='multiline',line_num=12):
    fiber_const = 1
    ECM_no_fiber_const = 0.4
    fiber_conc_boost =3
    fiber_diameter = 1
    #(0,0,1)->(300,300,1)        (-150,-150,1)->(150,150,1)
    if map_method == 'None':
        return fiber_const
    if map_method == 'random':
        return random.random()*ECM_no_fiber_const
    if map_method=='line': 
        if x >= 150 - fiber_diameter and x <= 150 + fiber_diameter:
            return fiber_const
        elif y >= 150 - fiber_diameter and y <= 150 + fiber_diameter:
            return fiber_const
        elif y >= x - fiber_diameter and y <= x + fiber_diameter:
            return fiber_const
        elif y <= 150 - x + fiber_diameter and y >= 150 - x - fiber_diameter:
            return fiber_const
        else:
            return 0
    if map_method == 'multiline':
        for i in xrange(line_num):
            if (y-100) * cos(pi * i/line_num) >= (x-100) * sin(pi * i/line_num) - fiber_diameter and (y-100) * cos(pi * i/line_num) <= (x-100) * sin(pi * i/line_num) + fiber_diameter:
                return fiber_const
        return 0
    if map_method =='lines_aligned':
        for i in xrange(line_num):
            if x >= 300 * i / line_num - fiber_diameter and x<= 300 * i /line_num + fiber_diameter:
                return fiber_const
        return 0

def getCellCOMPoint3D(cell):
    pt=CompuCell.Point3D()
    pt.x=int(round(cell.xCOM))
    pt.y=int(round(cell.yCOM))
    pt.z=int(round(cell.zCOM))
    
    return pt
class FiberAlignment(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        field_fiber      = self.getConcentrationField("fiber")
        field_fiber_cl = self.getConcentrationField('fiber_cl')
        
        for x,y,z in self.everyPixel():   
            field_fiber[x,y,z]  =  fiber_mapper(x,y,z,map_method='none') 
            field_fiber_cl[x,y,z]  = 0
        
    def step(self,mcs):
        fiber_concentration     = 0
        fiber_cl_concentration = 0  
        
#         field_fiber      = self.getConcentrationField("fiber")
#         field_fiber_cl = self.getConcentrationField('fiber_cl')
        
#         for x,y,z in self.everyPixel():
#             fiber_concentration     += field_fiber[x,y,z]
#             fiber_cl_concentration += field_fiber_cl[x,y,z]
                    
#         fileName='FiberAlignment.csv'
#         try:
#             fileHandle,fullFileName=self.openFileInSimulationOutputDirectory(fileName,"a")
#         except IOError:
#             print "Could not open file ", fileName," for writing. "
#             return
            
#         print >>fileHandle,mcs,",", fiber_concentration,",", fiber_cl_concentration
#         fileHandle.close()
            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
class TGrowthSteppable(SteppableBasePy):
    
    def start(self):    
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM):
            cell.targetVolume=V0
            cell.lambdaVolume=LBD_V0
            cell.targetSurface=S0
            cell.lambdaSurface=LBD_S0
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["Counter"] = 0
            cellDict["Health"] = 0
            
        for cell in self.cellListByType(self.NECROTIC):
            cell.targetVolume=V0_nec
            cell.lambdaVolume=LBD_NV0
            cell.targetSurface=S0_nec
            cell.lambdaSurface=LBD_NS0
            
    def step(self,mcs):    
        glucoseField=CompuCell.getConcentrationField(self.simulator,"Glucose")
        
        for cell in self.cellList: 
          pt = getCellCOMPoint3D(cell)
          glucoseField=CompuCell.getConcentrationField(self.simulator,"Glucose")
          conc = glucoseField.get(pt)
          cellDict = CompuCell.getPyAttrib(cell)
          if cell.type == self.NECROTIC:
              cell.targetVolume-=min(decvol, cell.targetVolume)
              cell.targetSurface = ktgs*sqrt(cell.targetVolume)
          if cell.type == self.PCANCER:
              cell.targetVolume+=incvol*max(0, conc-PGrThr0)
              cell.targetSurface =ktgs*sqrt(cell.targetVolume)
          if cell.type == self.PSTEM:
              cell.targetVolume+=incvol*max(0, conc-SGrThr0)
              cell.targetSurface =ktgs*sqrt(cell.targetVolume)   
    

    
class StravHealthCalculator(SteppableBasePy):
     
    def start(self):
            
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM, self.NECROTIC):
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["Starv"]=0
            cellDict["Health"]=0
            
    def MM(self,x,m,k):
        return (m*x/(x+k))
        
    def step(self,mcs):    
        glucoseField=CompuCell.getConcentrationField(self.simulator,"Glucose")
        for cell in self.cellList:
            if cell.type!=self.NECROTIC:
                cellDict = CompuCell.getPyAttrib(cell)
                pt = getCellCOMPoint3D(cell)
                conc = glucoseField.get (pt)
            if cell.type == self.PCANCER: 
               if conc < 0.0032:
                   cellDict["Starv"]+=abs(self.MM(conc,2.25,0.00256)\
                   -self.MM(GluD,2.25,0.00256))
               else:
                   cellDict["Health"]+=self.MM(conc,2.25,0.00256)\
                   -self.MM(GluD,2.25,0.00256)
            if cell.type == self.QCANCER: 
               if conc < 0.0032:
                   cellDict["Starv"]+=abs(self.MM(conc,1.69,0.00256)\
                   -self.MM(GluD,1.69,0.00256))
               else:
                   cellDict["Health"]+=self.MM(conc,1.69,0.00256)\
                   -self.MM(GluD,1.69,0.00256)
            if cell.type == self.PSTEM:
               if conc < 0.0032:
                   cellDict["Starv"]+=abs(self.MM(conc,2.25,0.00256)\
                   -self.MM(GluD, 2.25, 0.00256))
               else:
                   cellDict["Health"]+=self.MM(conc, 2.25, 0.00256)\
                   -self.MM(GluD, 2.25, 0.00256)
            if cell.type == self.QSTEM:
               if conc < 0.0032:
                   cellDict["Starv"]+=abs(self.MM(conc, 1.69, 0.00256)\
                   -self.MM(GluD, 1.69, 0.00256))
               else:
                   cellDict["Health"]+=self.MM(conc, 1.69, 0.00256)\
                   -self.MM(GluD, 1.69, 0.00256)
    def finish(self):
        
        return    
                        
class CellStateTransition(SteppableBasePy):
    
    def start(self):
        
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM):
            cellDict=CompuCell.getPyAttrib(cell) 
            
    def step(self,mcs): 
        for cell in self.cellList: 
          cellDict=CompuCell.getPyAttrib(cell)
          if cell.type == self.PCANCER: 
            if cellDict["Starv"] > PNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
          if cell.type == self. QCANCER: 
            if cellDict["Starv"] > QNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
            if cellDict["Health"] > QPThr0:
              cell.type = self.PCANCER
              cellDict["Health"]=0 
          if cell.type == self. PSTEM: 
            if cellDict["Starv"] > SNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
          if cell.type == self. QSTEM: 
            if cellDict["Starv"] > QSNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
            if cellDict["Health"] > 79.5:
              cell.type = self.PSTEM
              cellDict["Health"]=0   
                     

class MitosisSteppable(MitosisSteppableBase):
    def  start(self):
             
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM):
         
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["Counter"]=0    
        
    
    def step(self,mcs):
        
        cells_to_divide=[]
        
        
        for cell in self.cellList:
            if ((cell.type==self.PCANCER or cell.type==self.QCANCER)\
            and cell.volume>volmaxmit)\
            or ((cell.type==self.PSTEM or cell.type==self.QSTEM)\
            and cell.volume>Svolmaxmit):
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
                self.divideCellRandomOrientation(cell)
          
    def updateAttributes(self):
        
        parentCell = self.mitosisSteppable.parentCell
        childCell = self.mitosisSteppable.childCell
        parentCell.targetVolume=V0
        childCell.targetVolume= V0
        parentCell.targetSurface=ktgs*sqrt(parentCell.targetVolume)
        childCell.targetSurface= ktgs*sqrt(childCell.targetVolume)           
        parentCell.lambdaVolume=LBD_V0
        parentCell.lambdaSurface=LBD_S0
        childCell.lambdaVolume=LBD_V0
        childCell.lambdaSurface=LBD_S0
        parentCellDict=CompuCell.getPyAttrib(parentCell)
        childCellDict=CompuCell.getPyAttrib(childCell)
        
        if parentCell.type==self.PCANCER:
            
          temp = random.gauss(maxdiv,2)
          if parentCellDict["Counter"]<= temp:
             parentCell.type=self.QCANCER #both cells are QC after mitosis 
             childCell.type=self.QCANCER
             parentCellDict["Counter"]+=1
            
          childCellDict["Counter"] =deepcopy(parentCellDict["Counter"])
            
          if parentCellDict["Counter"]> temp:
              parentCell.type=self.NECROTIC
              childCell.type=self.NECROTIC
              
        if parentCell.type == self.PSTEM or parentCell.type == self.QSTEM: 
            parentCellDict["Counter"]+=1
            parentCell.type=self.QSTEM #one is QS
            childCell.type=self.QCANCER #the other is QC
            
            if random.random()<=probstem: #0.2 chance for a 2nd QS
                childCell.type=self.QSTEM
            childCellDict["Counter"]=0
            
        parentCellDict["Starv"]=0
        childCellDict["Starv"]=0
        parentCellDict["Health"]=0
        childCellDict["Health"]=0 
            
        
        
class LogData(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):   
        fileName='CELL_NUM.csv'
        try:
            fileHandle,fullFileName=self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
            print "Could not open file ", fileName," for writing. "
            return
        print >>fileHandle,"mcs",",", "PCANCER_NUM",",","QCANCER_NUM",",","PSTEM_NUM",",","QSTEM_NUM",",","NECROTIC_NUM"
        fileHandle.close()     

        fileName='CELL_VOL.csv'
        try:
            fileHandle,fullFileName=self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
            print "Could not open file ", fileName," for writing. "
            return
        print >>fileHandle,"mcs",",", "PCANCER_VOL",",","QCANCER_VOL",",","PSTEM_VOL",",","QSTEM_VOL",",","NECROTIC_VOL"
        fileHandle.close() 

    def step(self,mcs):
        PCANCER_NUM = 0
        QCANCER_NUM = 0
        PSTEM_NUM = 0
        QSTEM_NUM = 0
        NECROTIC_NUM = 0

        PCANCER_VOL = 0
        QCANCER_VOL = 0
        PSTEM_VOL = 0
        QSTEM_VOL = 0
        NECROTIC_VOL = 0

        for cell in self.cellList: 
            if cell.type == self.PCANCER: 
                PCANCER_NUM+=1
                PCANCER_VOL+=cell.volume
            if cell.type == self. QCANCER:
                QCANCER_NUM+=1
                QCANCER_VOL+=cell.volume
            if cell.type == self.PSTEM: 
                PSTEM_NUM+=1
                PSTEM_VOL+=cell.volume
            if cell.type == self. QSTEM:
                QSTEM_NUM+=1
                QSTEM_VOL+=cell.volume
            if cell.type == self. NECROTIC:
                NECROTIC_NUM+=1
                NECROTIC_VOL+=cell.volume
        

                    
        fileName='CELL_NUM.csv'
        try:
            fileHandle,fullFileName=self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
            print "Could not open file ", fileName," for writing. "
            return
            
        print >>fileHandle,mcs,",", PCANCER_NUM,",",QCANCER_NUM,",",PSTEM_NUM,",",QSTEM_NUM,",",NECROTIC_NUM
        fileHandle.close()    

        fileName='CELL_VOL.csv'
        try:
            fileHandle,fullFileName=self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
            print "Could not open file ", fileName," for writing. "
            return
        print >>fileHandle,mcs,",", PCANCER_VOL,",",QCANCER_VOL,",",PSTEM_VOL,",",QSTEM_VOL,",",NECROTIC_VOL
        fileHandle.close() 


    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
