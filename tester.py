import numpy as np
import time
import sys
from mapio.shake import ShakeGrid, readShakeFile
from mapio.multiple import MultiGrid

direc = '/Users/sverros/Documents/Modules/MM_CP/input/'
g = open(direc+'grid.xml')
sm_ = readShakeFile(g)
g.close()
