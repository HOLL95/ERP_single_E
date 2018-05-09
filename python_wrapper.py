#!/usr/bin/env python
import isolver
import math
import numpy as np
from timeit import default_timer as timer
import matplotlib.pyplot as plt
#import pints
start=timer()


end=timer()
print start-end
class single_E:
	def __init__(self,d_E, E_start, E_reverse,Ru,v,E0_mean,k0_mean,E0_sigma,k0_sigma ):
		self.dE=d_E
		self.E_start=E_start
		self.E_reverse=E_reverse
		self.E0_mean=E0_mean
		self.E0_sigma=E0_sigma
		self.k0_mean=k0_mean
		self.k0_sigma=k0_sigma
		self.a=0.03
		self.alpha=0.53
		self.Ru=Ru
		self.v=v
		self.gamma=6.5e-12
		self.E_0, self.T_0,self.I_0=single_E.characteristic_val(self)
	def characteristic_val(self):
		T=273+25
		F=96485.3328959
		R=8.314459848
		E_0=(R*T)/F
		T_0=abs(E_0/self.v)
		I_0=(F*self.a*self.gamma)/T_0
		return E_0, T_0,I_0
	def n_outputs(self):
		return 1
	def n_parameters(self):
		return 5
	def non_dim_capacative(self,parameters):
		parameters[0]=(parameters[0]/self.I_0)/self.T_0*(self.a*self.E_0)
		parameters[4]=parameters[4]*(2*math.pi*self.T_0)
		self.sampling_freq=0.05*((float((2*math.pi)))/float(parameters[4]))
		self.E0_mean=self.E0_mean/self.E_0
		self.k0_mean=self.k0_mean*self.T_0
		self.Ru=self.Ru/self.E_0*self.T_0
		self.E_start=self.E_start/self.E_0
		self.E_reverse=self.E_reverse/self.E_0
		self.dE=self.dE/self.E_0
		return parameters


	def times(self):
		final_time=(self.E_reverse-self.E_start)*2
		return np.linspace(0,final_time,final_time/self.sampling_freq, dtype='double')
	def simulate(self, parameters, time):
		#Parameters 1-4= Cdl, Cdl1, Cdl2, Cdl3, omega
		output=isolver.I_tot_solver(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4],self.v, self.alpha, self.E_start, self.E_reverse, self.dE, self.Ru, self.E0_mean, self.k0_mean, self.E0_sigma, 1)
		output=np.array(output)
		##final_time=(self.E_reverse-self.E_start)*2
		##time=np.linspace(0,final_time,len(output), dtype='double')
		##plt.plot(time,output)
		##plt.xlabel('time')
		##plt.ylabel('Itot')
		##plt.show()
		return output




