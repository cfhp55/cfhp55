import sys
from os import environ
from os import getcwd
import math
from math import *
import string
import numpy as np

sys.path.append(environ["PYTHON_MODULE_PATH"])

pi = math.pi

import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
# add extra attributes here
pyAttributeAdder,dictionaryAdder=CompuCellSetup.attachDictionaryToCells(sim)        
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
from PySteppables import SteppableRegistry
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from NecroticCoreSteppables import FiberAlignment
instanceOfFiberAlignment=FiberAlignment(sim,1)
steppableRegistry.registerSteppable(instanceOfFiberAlignment)

from NecroticCoreSteppables import StravHealthCalculator
stravHealthCalculator=StravHealthCalculator(sim,1)
steppableRegistry.registerSteppable(stravHealthCalculator)
        
from NecroticCoreSteppables import CellStateTransition
cellStateTransition=CellStateTransition(sim,1)
steppableRegistry.registerSteppable(cellStateTransition)

from NecroticCoreSteppables import TGrowthSteppable
TGrowthSteppables=TGrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(TGrowthSteppables)
        

from NecroticCoreSteppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(MitosisSteppableInstance)
        

from NecroticCoreSteppables import LogData
LogDataSteppableInstance = LogData(sim,_frequency=1)
steppableRegistry.registerSteppable(LogDataSteppableInstance)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        