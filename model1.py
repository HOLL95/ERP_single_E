#!/usr/bin/env python
import isolver
import math
import numpy as np
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import pints
from single_e_Itot import single_E
capacative_parameters=[0.000133878548046, 0.000653657774506,0.000245772700637,1.10053945995e-06,8.88480830076]
parameters=np.zeros((1,len(capacative_parameters)))
#capacative_parameters=np.array(capacative_parameters)
#x=single_E_capacative(150e-3, -0.1,-0.85,20,27.94e-2,-0.393545799342,10000.0,0.0312279186927,1)
x=single_E_capacative(150e-3, -0.85,-0.1,27.1,24.97e-3,-0.393545799342,10000,0.0312279186927,1)
x.characteristic_val()
parameters=x.non_dim_capacative(capacative_parameters)
time=x.times()
attrs=vars(x)
print len(time)
print attrs
print capacative_parameters
output=x.simulate(capacative_parameters, time)
print len(output)
#output += np.random.normal(size=output.shape)*0.04
#problem=pints.SingleOutputProblem(x,time, output)
#score = pints.SumOfSquaresError(problem)
#initial_guess=np.array([0,0,0,0,math.pi*2.85])
#boundaries=pints.Boundaries([0 ,-1,-0.1,-0.1,0.99*math.pi*2.85],[10,1,0.1,0.1,1.01*math.pi*2.85])
#found_parameters, found_value=pints.optimise(score, initial_guess, boundaries=boundaries, method=pints.CMAES)
#print found_parameters


