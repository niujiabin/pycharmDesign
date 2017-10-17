#coding=utf-8
import matlab
import numpy as np
import matlab.engine
eng = matlab.engine.start_matlab()
searchAgent = 50
function_name = 'F1'
maxi = float(10);
ret = eng.main(searchAgent, function_name,maxi)
print ret



