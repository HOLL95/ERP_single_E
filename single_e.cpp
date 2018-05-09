#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/distributions/normal.hpp> 
#include <boost/math/distributions/lognormal.hpp>
#include <math.h> 
#include <iostream>
#include <exception>
#include <vector>
#include <time.h>
#include "boost/tuple/tuple.hpp"
#include <boost/math/tools/roots.hpp>
#include </users/henney/Documents/pybind11/include/pybind11/pybind11.h>
#include </users/henney/Documents/pybind11/include/pybind11/stl.h>
namespace py = pybind11;
//using namespace boost::python;



struct param1{
	double delta_E;
	double E_start;
	double E_reverse;
	double Cdl;
	double CdlE;
	double CdlE2;
	double CdlE3;
	double E0;
	double Ru;
	double R;
	double k0;
	double alpha;
	double In_0;
	double dt;
	double gamma;
	double omega;
	double phase;
	double a;
	double v;
	double tr;
	double E0_mean;
	double k0_mean;
	double E0_sigma;
	double k0_sigma;
	int time_end;
	int duration;
	double I_0;
	double I_1;
	double theta_0;
	float t;
};

	

///////////////////////////////////////////////////////////////////////////////SUBROUTINES//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////SOLVER SUBROUTINES///////////////////////////////////////////////////////////////////////////////////////////


double e_t(param1& single_e_param, float t){
	double E_dc;
	double E_t;
	if (t<single_e_param.tr){
		E_dc=single_e_param.E_start+(single_e_param.v*t);
	}else {
		E_dc=single_e_param.E_reverse-(single_e_param.v*(t-single_e_param.tr));
	}
	
	E_t=E_dc+(std::sin((single_e_param.omega*t)+single_e_param.phase));
	
	return E_t;
}

double dEdt(param1& single_e_param,float t){
	double E_dc;
	double dedt;
	//t=t+0.5*single_e_param.dt;
	if (t < single_e_param.tr){
		 E_dc=single_e_param.v;
	}else {
		 E_dc=-single_e_param.v;
	}
   	dedt=E_dc+(single_e_param.delta_E*single_e_param.omega*std::cos(single_e_param.omega*t+single_e_param.phase));
	
	return dedt;
}


double theta_1(param1& single_e_param,float t, double In1, float theta_0,double E){
		const double Ereduced = E- (single_e_param.Ru*In1);
		const double expval1 = Ereduced - single_e_param.E0_mean;
		double exp11 = std::exp((1.0-single_e_param.alpha)*expval1);
		double exp12 = std::exp(-single_e_param.alpha*expval1);
		const double u1n1_top = single_e_param.dt*single_e_param.k0_mean*exp11 + theta_0;
		const double denom = (single_e_param.dt*single_e_param.k0_mean*exp11 + single_e_param.dt*single_e_param.k0_mean*exp12 + 1);
		const double tmp = 1.0/denom;
		double u1n1 = u1n1_top*tmp;
		return u1n1;
}
double dtheta_1(param1& single_e_param,float t, double In1, double theta_0, double E){
		const double Ereduced = E - (single_e_param.Ru*In1);
		const double expval1 = Ereduced - single_e_param.E0_mean;
		double exp11 = std::exp((1.0-single_e_param.alpha)*expval1);
		double exp12 = std::exp(-single_e_param.alpha*expval1);

		double dexp11 = -single_e_param.Ru*(1.0-single_e_param.alpha)*exp11;
		double dexp12 = single_e_param.Ru*single_e_param.alpha*exp11;

		const double u1n1_top = single_e_param.dt*single_e_param.k0_mean*exp11 + theta_0;
		const double du1n1_top = single_e_param.dt*single_e_param.k0_mean*dexp11;
		const double denom = (single_e_param.dt*single_e_param.k0_mean*exp11 + single_e_param.dt*single_e_param.k0_mean*exp12 + 1);
		const double ddenom = single_e_param.dt*single_e_param.k0_mean*(dexp11 + dexp12);
		const double tmp = 1.0/denom;
		const double tmp2 = pow(tmp,2);
		double du1n1 = -(u1n1_top*ddenom + du1n1_top*denom)*tmp2;
		return du1n1;
}
double poly_er(param1& single_e_param, double E, double I_1){
	double Er=E-(single_e_param.Ru*I_1);
	double Er2=Er*Er;
	double Er3=Er2*Er;
	double Cdlp=single_e_param.Cdl*(1+(single_e_param.CdlE*Er)+(single_e_param.CdlE2*Er2)+(single_e_param.CdlE3*Er3));
	return Cdlp;
}
double residual(param1& single_e_param, double In1, double In0,float t, double theta_0, double E, double dE) {
        	return single_e_param.gamma*(theta_1(single_e_param, t, In1,theta_0,E)-theta_0)+single_e_param.dt*single_e_param.R*(E-single_e_param.Ru*In1) - single_e_param.dt*In1+poly_er(single_e_param, E,In1)*((single_e_param.dt*dE)-single_e_param.Ru*(In1-In0));

}
double residual_gradient(param1& single_e_param, const double In1,float t, float theta_0, double E) {
	        return single_e_param.gamma*dtheta_1(single_e_param, t, In1,theta_0,E)- single_e_param.dt*single_e_param.R*single_e_param.Ru - single_e_param.dt-poly_er(single_e_param, E,In1)*single_e_param.Ru;
     
}




double newton_raphson(param1& single_e_param,double I_0,double theta_0, float t,double E, double dE){	
	double newton_estimate=I_0-(residual(single_e_param, I_0, I_0, t, theta_0,E,dE)/residual_gradient(single_e_param, I_0, t, theta_0,E));
	int i=0;
	while(  i<100 ){ //abs((newton_estimate[i]-newton_estimate[i-1]))>10e-15 &&
		i++;	
		newton_estimate=I_0-(residual(single_e_param, newton_estimate,I_0, t, theta_0, E,dE)/residual_gradient(single_e_param, newton_estimate, t,theta_0, E));
		
		
	}
	
	return newton_estimate;
}

std::vector<double> non_linear_I_solver(param1& single_e_param){
	float t=0;
	std::vector<double> I_matrix(single_e_param.duration, 0);
	std::vector<double> CDL_matrix(single_e_param.duration, 0);
	I_matrix[0]=0;
	double E=0;
	double dE=0;
	double theta_0=0;
	CDL_matrix[0]=0;
	for(int i=1; i<single_e_param.duration; i++){
		t=t+single_e_param.dt;
		E=e_t(single_e_param,t);
		dE=dEdt(single_e_param,t);
		I_matrix[i]=newton_raphson(single_e_param, I_matrix[i-1], theta_0,t,E,dE);
		//CDL_matrix[i]=poly_er(single_e_param, E, I_matrix[i])*  ((dE)-single_e_param.Ru*(I_matrix[i]-I_matrix[i-1]));//
		//std::cout<<residual(single_e_param, I_matrix[i], I_matrix[i-1], t, theta_0, E,dE)<<"\n";
		theta_0=theta_1(single_e_param,t,I_matrix[i], theta_0,E);
	}
	return I_matrix;
}

///////////////////////////////////////////////////////////////////////////////DISPERSION SUBROUTINES//////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>>weightmatrix(param1& single_e_param,int bins, std::vector<double> E0_disp,std::vector<double> k0_disp ){
	double weight_E;
	double weight_k;
	double x1;
	boost::math::normal_distribution<double> E0(single_e_param.E0_mean, single_e_param.E0_sigma);
	boost::math::lognormal_distribution<double> k0(single_e_param.k0_mean, single_e_param.k0_sigma);
	std::vector< std::vector<double>> weight_matrix((bins), std::vector< double >((bins)));
	for(int i=0; i<bins; i++){
		if(i==0){
			weight_E=cdf(E0,E0_disp[0]);
		}else{
			weight_E=cdf(E0,E0_disp[i])-cdf(E0,E0_disp[i-1]);
		}
		
		for(int j=0; j<bins; j++){
			
			
			if(j==0){
				weight_k=cdf(k0, k0_disp[0]);
			}else{
				weight_k=cdf(k0,k0_disp[j])-cdf(k0, k0_disp[j-1]);
			}	
			x1=x1+weight_E*weight_k;
			weight_matrix[i][j]=weight_E*weight_k;
		}
	
	}
	
	return weight_matrix;
	
}		

std::vector<double> dispersion_solver(param1& single_e_param, std::vector<std::vector<double>>weight_matrix,std::vector<double> E0_disp,std::vector<double> k0_disp ,int bins){
	std::vector<double> I_matrix(single_e_param.duration);
	std::vector<double> I_disp(single_e_param.duration);
	for(int i=0; i<(bins); i++){
		single_e_param.E0=E0_disp[i];
		for(int j=0; j<bins; j++){
		 	single_e_param.k0=k0_disp[j];
			I_matrix=non_linear_I_solver(single_e_param);
	 		for(int k=0; k<single_e_param.duration; k++){
				I_disp[k]=I_disp[k]+(I_matrix[k]*weight_matrix[i][j]);
			}
		}
	}
	return I_disp;
}

///////////////////////////////////////////////////////////////////////////////PYTHON WRAPPER SUBROUTINES//////////////////////////////////////////////////////////////////////////////////



py::object I_tot_solver(double Cdl, double CdlE1, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, double delta_E, double Ru, double dt, double time_end, int time_length, double E0_mean=0, double k0_mean=0, double E0_sigma=0.1, double k0_sigma=1) {	
	param1 single_e_param;
	single_e_param.v=v;
	single_e_param.Ru=Ru;
	single_e_param.delta_E=delta_E;
	single_e_param.E_start=E_start;//E_start;
	single_e_param.E_reverse=E_reverse;//E_reverse;
	single_e_param.Cdl=Cdl;//0.000133878548046;//
	single_e_param.CdlE=CdlE1;//0.000653657774506;//
	single_e_param.CdlE2=CdlE2;//0.000245772700637;//
	single_e_param.CdlE3=CdlE3;//1.10053945995e-06;//
	single_e_param.omega=omega;//boost::math::constants::pi<double>()*2;//
	single_e_param.alpha=0.53;
	single_e_param.gamma=6.5e-12;//6.5e-12/0.03;
	single_e_param.R=0;
	single_e_param.phase=0;
	single_e_param.E0_mean=E0_mean;  
	single_e_param.k0_mean=k0_mean;
	single_e_param.E0_sigma=E0_sigma;//
	single_e_param.k0_sigma=k0_sigma;//
	single_e_param.time_end=time_end;
	single_e_param.dt=dt;
	single_e_param.duration=time_length;	
	single_e_param.tr=((single_e_param.E_reverse-single_e_param.E_start)/single_e_param.v);
	//WEIGHT CALCULATIONS
	int bins=16;
	std::vector<double>I_disp(single_e_param.duration,0);
	std::vector<double>E0_disp(bins,0);
	std::vector<double>k0_disp(bins,0);
	float E0_interval=(abs(-0.4-0.4))/bins;
	float k0_interval=(0+50)/bins;	
	E0_disp[0]=-0.3;
	k0_disp[0]=0;
	for(int i=1; i<bins; i++){
		E0_disp[i]=E0_disp[i-1]+E0_interval;
		k0_disp[i]=k0_disp[i-1]+k0_interval;
	}
	std::vector< std::vector< double > > weight_matrix( bins, std::vector< double >( bins ) );
	//weight_matrix=weightmatrix(single_e_param, bins, E0_disp,k0_disp);
	//I_disp=dispersion_solver(single_e_param, weight_matrix, E0_disp, k0_disp, bins);
	I_disp=non_linear_I_solver(single_e_param);
	return py::cast(I_disp);

}

PYBIND11_MODULE(isolver, m) {
	m.def("I_tot_solver", &I_tot_solver, "solve for I_tot with dispersion");
	
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){ 
	//std::vector<double>I_disp(800,0);
	//I_disp=I_tot_solver(0.000133878548046, 0.000653657774506,0.000245772700637,1.10053945995e-06,boost::math::constants::pi<double>()*2,0,0,0.1,1,40,800);
	
	//clock_t t;
	//t = clock();
	//I_matrix=non_linear_I_solver(single_e_pa);	 
	//for(int i=0; i<800; i++){
	//	std::cout<<I_disp[i] << "\n";
	//}
	//t= clock() - t;
	//std::cout<<(t*1.0/CLOCKS_PER_SEC)<<"\n";	//
	
	
}
	

	

	

	
	

	
		


	
