#!/usr/bin/env python
import isolver
import math
import numpy as np
from timeit import default_timer as timer
import matplotlib.pyplot as plt
#import pints
start=timer()
time1=np.arange(0,20.05,0.05)





#output=simulate(parameters, time1)
#frequency, filtered, unfiltered=fourier_simulate(parameters, time1)
end=timer()
print start-end
class single_E_capacative:

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
		self.E_0, self.T_0,self.I_0=single_E_capacative.characteristic_val(self)
	def characteristic_val(self):
		T=273+25
		F=96485.3328959
		R=8.314459848
		E_0=(R*T)/F
		T_0=abs(E_0/self.v)
		I_0=(F*self.a*self.gamma)/T_0
		return E_0, T_0,I_0
	#
	def n_outputs(self):
		return 1
	def n_parameters(self):
		return 5
	def non_dim_capacative(self,parameters):
		parameters[0]=(parameters[0]/self.I_0)/self.T_0*(self.a*self.E_0)
		parameters[4]=parameters[4]*(2*math.pi*self.T_0)
		self.sampling_freq=0.005*((float((2*math.pi)))/float(parameters[4]))
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
		output=isolver.I_tot_solver(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4],self.v, self.alpha, self.E_start, self.E_reverse, self.dE, self.Ru, self.E0_mean, self.k0_mean, self.E0_sigma, 1)
		output=np.array(output)
		final_time=(self.E_reverse-self.E_start)*2
		print final_time
		time=np.linspace(0,final_time,len(output), dtype='double')
		Fs=0.005*((float((2*math.pi)))/float(parameters[4]))
		T=1/Fs
		L=len(output)
		Y=(np.fft.fft(output))
		plt.title("input function -capcative and faradic")
		plt.plot(time,output)
		plt.xlabel('time')
		plt.ylabel('Itot')
#		Y=abs(Y)
#		f=(Fs)*np.linspace(0,L,len(Y))
#		#Y=Y[:L/2]
#		plt.subplot(122)
#		plt.plot(f,Y)
#		plt.title("fourier transformed input function")
#		x1=np.arange(4*50,50*13,50)
#		top_hat=np.zeros((len(x1), len(Y)))
#		for i in range(0,len(x1)):
#			top_hat[i,:]=Y;
#			top_hat[i,:][f<(x1[i]-10)] =0
#			top_hat[i,:][f>(x1[i]+10)] =0
#		top_hat=np.sum(top_hat, axis=0)
#		#plt.subplot(223)
		#plt.title("filtered fourier transform")
		#plt.plot(f,top_hat)
		#plt.subplot(224)
		#plt.plot(np.fft.ifft(top_hat))
		#plt.title("filtered inverse fourier transform")
		plt.show()
		return output
	
	

#class single_E_faradic:#
#	def __init__(self,Cdl, CdlE,CdlE2, CdlE3,omega,d_E, E_start, E_reverse,Ru,v):
#		self.Cdl=Cdl
#		self.CdlE=CdlE
#		self.CdlE2=CdlE2
#		self.CdlE3=CdlE3
#		self.omega=omega
#		self.dE=d_E
#		self.E_start=E_start
#		self.E_reverse=E_reverse
#		self.a=0.03
#		self.alpha=0.53
#		self.Ru=Ru
#		self.v=v
#		self.sampling_freq=0.005
#		self.E_0, self.T_0,self.I_0=single_E_faradic.characteristic_val(self)
#	def characteristic_val(self):
#		T=273+25
#		F=96485.3328959
#		R=8.314459848
#		E_0=(R*T)/F
#		T_0=abs(E_0/self.v)
#		gamma=6.5e-12
#		I_0=(F*self.a*gamma)/T_0
#		return E_0, T_0,I_0
#	def non_dim_faradic(self,parameters):
#		self.Cdl=self.Cdl/self.I_0/self.T_0*(a*self.E_0)
#		self.CdlE=self.CdlE/self.I_0/self.T_0*(a*self.E_0)
#		self.CdlE2=self.CdlE2/self.I_0/self.T_0*(a*self.E_0)
#		self.CdlE3=self.CdlE3/self.I_0/self.T_0*(a*self.E_0)
#		self.omega=self.omega*(2*math.pi*self.T_0)
#		self.Ru=self.Ru/self.E_0*self.T_0
#		self.E_start=self.E_start/self.E_0
#		self.E_reverse=self.E_reverse/self.E_0
#		self.dE=self.dE/self.E_0
#		parameters[0]=parameters[0]/self.E_0
#		parameters[1]=parameters[1]*self.T_0
#		return parameters
#	def times(self):
#		final_time=((self.E_reverse-self.E_start)/self.v)*2
#		return np.linspace(0,final_time,(final_time/self.sampling_freq), dtype='double')
#	def simulate(self, parameters, time,flag="optimise"):
#		output=isolver.I_tot_solver(self.Cdl, self.CdlE,self.CdlE2, self.CdlE3,self.omega, self.v,self.alpha, self.E_start, self.E_reverse, self.dE, self.Ru,parameters[0], parameters[1],parameters[2], p##arameters[3])
#		output=np.array(output)
##		if flag=="initial":
#			output += np.random.normal(size=output.shape)*0.04
#		Fs=self.sampling_freq
#		T=1/Fs
#		L=len(time1)
##
#		Y=(np.fft.fft(output))
#		plt.subplot(221)
#		plt.title("input function")
#		plt.plot(time,output)
#		Y=abs(Y)
#		#Y=Y[:L/2]
#		plt.subplot(222)
##		plt.plot((Y))
#		plt.title("fourier transformed input function")
#		f=(Fs)*np.linspace(0,L,len(Y))
#		x1=np.arange(4*50,50*13,50)
#		top_hat=np.zeros((len(x1), len(Y)))
#		for i in range(0,len(x1)):
#			top_hat[i,:]=Y;
#			top_hat[i,:][f<(x1[i]-10)] =0
#			top_hat[i,:][f>(x1[i]+10)] =0
##		top_hat=np.sum(top_hat, axis=0)
#		plt.subplot(223)
#		plt.title("filtered fourier transform")
#		plt.plot(top_hat)
#		plt.subplot(224)
#		plt.plot(np.fft.ifft(top_hat))
#		plt.title("filtered inverse fourier transform")
#		plt.show()
#		return top_hat
#	def n_outputs(self):
#		return 1
#	def n_parameters(self):
#		return 4
###
