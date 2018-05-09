#!/usr/bin/env python
import isolver
import math
import numpy as np
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import pints
from python_wrapper import single_E
capacative_parameters=[0.000133878548046, 0.000653657774506,0.000245772700637,1.10053945995e-06,8.88480830076]
parameters=np.zeros((1,len(capacative_parameters)))
x=single_E(150e-3, -0.85,-0.1,27.1,24.97e-3,-0.393545799342,4000,0.0312279186927,1)
x.characteristic_val()
parameters=x.non_dim_capacative(capacative_parameters)
time=x.times()
attrs=vars(x)
print len(time)
print attrs
print capacative_parameters
output=x.simulate(capacative_parameters, time)
E_0, T_0,I_0 = x.characteristic_val()
exp_data = np.loadtxt('HypD_9Hz_cv_current',skiprows=19)
exp_t = exp_data[:,0]/T_0
exp_I = exp_data[:,1]#/I_0

plt.plot(exp_t,exp_I,label='data')
#plt.plot(time,output,label='sim')
plt.show()



