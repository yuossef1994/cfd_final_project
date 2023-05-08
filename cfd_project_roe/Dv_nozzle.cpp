//
//  Dv_nozzle.cpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#include "Dv_nozzle.hpp"
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;



void Dv_nozzle::set_geometry (int imax, double x_max, double x_min)
{
    double i_n_kasi, i_n_eta;
    
    
    //imax is number of cells
    //_imax is number of faces
    
     im_f = _imax_f-1;
     im_c = _imax_c-1;
     i_n_kasi=im_c+1;
     i_n_eta=im_c+1;
   
    
    for(int j=_imin;j<=i_n_eta;j++)
    {
        for( int i =_imin; i <= i_n_kasi; i++)
        {
            x_n[i][j]=   x_min + (float(i)/float(im_f))*(x_max- x_min) ;
            y_n[i][j]=   x_min + (float(j)/float(im_f))*(x_max- x_min) ;
            x_c[i]=   x_min + (float(i)/float(im_f))*(x_max- x_min) ;
            //std::cout<<"point is "<<" ( "<<x_n[i][j]<<","<<y_n[i][j]<<" ) "<<std::endl;
            
        }
        
    }
    
   
}
    void Dv_nozzle::initialization (double P_0, double T_0)
{
        
        
        
        _P_0= P_0;
        _T_0= T_0;
        
        
        
        
        
        
        
        
        
        // 0 is rho, 1 is velocity,  2 is pressure
        //values calculated at cell as cell average
        for (int j= _imin; j<=im_c; j++)
        {
            for (int i= _imin; i<= im_c; i++)
            {
                
                _M_0 = 0.85*x_c[i] + 1;
                epsi_1= 1.0 + ((gamma-1.0)*_M_0*_M_0)/2.0;
                
                T[i][j] = _T_0/epsi_1;
                
                
                _V[1][i][j] = _M_0 * sqrt(gamma*R_air*abs(T[i][j]));
                _V[2][i][j] = _M_0 * sqrt(gamma*R_air*abs(T[i][j]));
                
                _V[3][i][j]= _P_0/pow(epsi_1,gamma/(gamma-1));
                
                _V[0][i][j]= _V[3][i][j]/(R_air*T[i][j]);
                
                et[i][j] = (R_air/(gamma-1))*T[i][j] + 0.5* (_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j]);
                
                _U[0][i][j]= _V[0][i][j];
                
                _U[1][i][j]= _V[0][i][j]*_V[1][i][j];
                
                _U[2][i][j]= _V[0][i][j]*_V[2][i][j];
                
                _U[3][i][j]=_V[0][i][j]*et[i][j];
                
                
                
                _mach[i][j]=_M_0;
                
            }
            
        }
        
        //calculate epsi- kasi direction kasi corresponds to x kinda
        
        for (int j=2;j<im_c;j++)
        {
            for (int i=2;i<im_c;i++)
            {
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i-1][j]),1e-6),_V[0][i][j]-_V[0][i-1][j]);
                r_plus= (_V[0][i+1][j]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i-1][j]-_V[0][i-2][j])/(denm);
                epsi_plus_x[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i-1][j]),1e-6),_V[1][i][j]-_V[1][i-1][j]);
                r_plus= (_V[1][i+1][j]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i-1][j]-_V[1][i-2][j])/(denm);
                epsi_plus_x[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i-1][j]),1e-6),_V[2][i][j]-_V[2][i-1][j]);
                r_plus= (_V[2][i+1][j]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i-1][j]-_V[2][i-2][j])/(denm);
                epsi_plus_x[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i-1][j]),1e-6),_V[3][i][j]-_V[3][i-1][j]);
                r_plus= (_V[3][i+1][j]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i-1][j]-_V[3][i-2][j])/(denm);
                epsi_plus_x[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                
            }
            
        }
        
        
        //calculate epsi-eta direction  corresponds to y kinda
        
        for (int j=2;j<im_c;j++)
        {
            for (int i=2;i<im_c;i++)
            {
                
              
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i][j-1]),1e-6),_V[0][i][j]-_V[0][i][j-1]);
                r_plus= (_V[0][i][j+1]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i][j-1]-_V[0][i][j-2])/(denm);
                epsi_plus_y[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
               
                
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i][j-1]),1e-6),_V[1][i][j]-_V[1][i][j-1]);
                r_plus= (_V[1][i][j+1]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i][j-1]-_V[1][i][j-2])/(denm);
                epsi_plus_y[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i][j-1]),1e-6),_V[2][i][j]-_V[2][i][j-1]);
                r_plus= (_V[2][i][j+1]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i][j-1]-_V[2][i][j-2])/(denm);
                epsi_plus_y[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i][j-1]),1e-6),_V[3][i][j]-_V[3][i][j-1]);
                r_plus= (_V[3][i][j+1]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i][j-1]-_V[3][i][j-2])/(denm);
                epsi_plus_y[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                
            }
            
        }
        
        roe_flux();
        
        
        
        
        
        
    }



void Dv_nozzle::set_boundary_cond()
{
   
    for (int j=0;j<=im_c;j++)
        
    {
        //inflow boundary conditions
        mach_in = 0.5*(3.0*_mach[0][j] - _mach[1][j] );
        epsi_1= 1 + ((gamma-1)*mach_in*mach_in)/2;
        T_in = _T_0/epsi_1;
        
        
        u_in =  mach_in * sqrt(gamma*R_air*abs(T_in));
        v_in=0;
        P_in= _P_0/pow(epsi_1,gamma/(gamma-1));
        
        rho_in= P_in/(R_air*T_in);
        
        // et_in = (R_air/(gamma-1))*T_in + 0.5*u_in*u_in;
        
        //extrapolate U and V for 1st ghost cell
        /*
         _U_ghost_inflow[0]=2*_U[0][_imin]-_U[0][_imin-1];
         _U_ghost_inflow[1]=2*_U[1][_imin]-_U[1][_imin-1];
         _U_ghost_inflow[2]=2*_U[2][_imin]-_U[2][_imin-1];
         */
        
        //first ghoast cell
        
        _V_ghost_inflow[0][j]=2*rho_in-_V[0][_imin][j];
        _V_ghost_inflow[1][j]=2*u_in-_V[1][_imin][j];
        _V_ghost_inflow[2][j]=2*v_in-_V[2][_imin][j];
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_inflow[3][j]=2*P_in-_V[3][_imin][j];
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow[5][j]=2*T_in-T[_imin][j];    //T
        
        //second ghost cell
        _V_ghost_inflow_2[0][j]=2*_V_ghost_inflow[0][j]-_V[0][_imin][j];
        _V_ghost_inflow_2[1][j]=2*_V_ghost_inflow[1][j]-_V[1][_imin][j];
        
        _V_ghost_inflow_2[2][j]=2*_V_ghost_inflow[2][j]-_V[2][_imin][j];
        _V_ghost_inflow_2[3][j]=2*_V_ghost_inflow[3][j]-_V[3][_imin][j];
        //  _V_ghost_inflow_2[2][j]=0;
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow_2[5][j]=2*_V_ghost_inflow[4][j]-T[_imin][j];
        
        //third ghost cell
        _V_ghost_inflow_3[0][j]=2*_V_ghost_inflow_2[0][j]-_V_ghost_inflow[0][j];
        _V_ghost_inflow_3[1][j]=2*_V_ghost_inflow_2[1][j]-_V_ghost_inflow[1][j];
        
        _V_ghost_inflow_3[2][j]=2*_V_ghost_inflow_2[2][j]-_V_ghost_inflow[2][j];
        _V_ghost_inflow_3[3][j]=2*_V_ghost_inflow_2[3][j]-_V_ghost_inflow[3][j];
        // _V_ghost_inflow_3[2][j]=0;
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow_3[5][j]=2*_V_ghost_inflow_2[4][j]-_V_ghost_inflow[4][j];
        
        
        
        
        
        
        
        
        
        //calculate epsi at 1
        
        
        
        //for density
        double denm= copysign(max(abs(_V[0][1][j]-_V[0][0][j]),1e-6),_V[0][1][j]-_V[0][0][j]);
        r_plus= (_V[0][2][j]-_V[0][1][j])/(denm);
        r_minus= (_V[0][0][j]-_V_ghost_inflow[0][j])/(denm);
        epsi_plus_x[0][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[0][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u- velocity
        denm= copysign(max(abs(_V[1][1][j]-_V[1][0][j]),1e-6),_V[1][1][j]-_V[1][0][j]);
        r_plus= (_V[1][2][j]-_V[1][1][j])/(denm);
        r_minus= (_V[1][0][j]-_V_ghost_inflow[1][j])/(denm);
        epsi_plus_x[1][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[1][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v velocity
        denm= copysign(max(abs(_V[2][1][j]-_V[2][0][j]),1e-6),_V[2][1][j]-_V[2][0][j]);
        r_plus= (_V[2][2][j]-_V[2][1][j])/(denm);
        r_minus= (_V[2][0][j]-_V_ghost_inflow[2][j])/(denm);
        epsi_plus_x[2][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[2][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        denm= copysign(max(abs(_V[3][1][j]-_V[3][0][j]),1e-6),_V[3][1][j]-_V[3][0][j]);
        r_plus= (_V[3][2][j]-_V[3][1][j])/(denm);
        r_minus= (_V[3][0][j]-_V_ghost_inflow[3][j])/(denm);
        epsi_plus_x[3][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[3][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
        //calculate epsi at 0
        
        //for density
        denm= copysign(max(abs(_V[0][0][j]-_V_ghost_inflow[0][j]),1e-6),_V[0][0][j]-_V_ghost_inflow[0][j]);
        r_plus= (_V[0][1][j]-_V[0][0][j])/(denm);
        r_minus= (_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j])/(denm);
        epsi_plus_x[0][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[0][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u velocity
        denm= copysign(max(abs(_V[1][0][j]-_V_ghost_inflow[1][j]),1e-6),_V[1][0][j]-_V_ghost_inflow[1][j]);
        r_plus= (_V[1][1][j]-_V[1][0][j])/(denm);
        r_minus= (_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j])/(denm);
        epsi_plus_x[1][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[1][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v- velocity
        denm= copysign(max(abs(_V[2][0][j]-_V_ghost_inflow[2][j]),1e-6),_V[2][0][j]-_V_ghost_inflow[2][j]);
        r_plus= (_V[2][1][j]-_V[2][0][j])/(denm);
        r_minus= (_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j])/(denm);
        epsi_plus_x[2][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[2][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        denm= copysign(max(abs(_V[3][0][j]-_V_ghost_inflow[3][j]),1e-6),_V[3][0][j]-_V_ghost_inflow[3][j]);
        r_plus= (_V[3][1][j]-_V[3][0][j])/(denm);
        r_minus= (_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j])/(denm);
        epsi_plus_x[3][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[3][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        //calculate epsi at imin-1
        
        //for density
        denm= copysign(max(abs(_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j]),1e-6),_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j]);
        r_plus= (_V[0][0][j]-_V_ghost_inflow[0][j])/(denm);
        r_minus= (_V_ghost_inflow_2[0][j]-_V_ghost_inflow_3[0][j])/(denm);
        epsi_plus_x_ghost[0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u-velocity
        denm= copysign(max(abs(_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j]),1e-6),_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j]);
        r_plus= (_V[1][0][j]-_V_ghost_inflow[1][j])/(denm);
        r_minus= (_V_ghost_inflow_2[1][j]-_V_ghost_inflow_3[1][j])/(denm);
        epsi_plus_x_ghost[1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v -velocity
        
        denm= copysign(max(abs(_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j]),1e-6),_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j]);
        r_plus= (_V[2][0][j]-_V_ghost_inflow[2][j])/(denm);
        r_minus= (_V_ghost_inflow_2[2][j]-_V_ghost_inflow_3[2][j])/(denm);
        epsi_plus_x_ghost[2][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[2][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        //for pressure
        denm= copysign(max(abs(_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j]),1e-6),_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j]);
        r_plus= (_V[3][0][j]-_V_ghost_inflow[3][j])/(denm);
        r_minus= (_V_ghost_inflow_2[3][j]-_V_ghost_inflow_3[3][j])/(denm);
        epsi_plus_x_ghost[3][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[3][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
        //flux at face 1
        
        
        
        //density
        V_L[0]= _V[0][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][0][j]*(_V[0][0][j]-_V_ghost_inflow[0][j]) + (1+kappa)*epsi_minus_x[0][1][j]*(_V[0][1][j]-_V[0][0][j])        );
        
        V_R[0]= _V[0][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][2][j]*(_V[0][2][j]-_V[0][1][j]) + (1+kappa)*epsi_plus_x[0][1][j]*(_V[0][1][j]-_V[0][0][j])        );
        //velocity
        V_L[1]= _V[1][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][0][j]*(_V[1][0][j]-_V_ghost_inflow[1][j]) + (1+kappa)*epsi_minus_x[1][1][j]*(_V[1][1][j]-_V[1][0][j])        );
        
        V_R[1]= _V[1][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][2][j]*(_V[1][2][j]-_V[1][1][j]) + (1+kappa)*epsi_plus_x[1][1][j]*(_V[1][1][j]-_V[1][0][j])        );
        //pressure
        V_L[2]= _V[2][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][0][j]*(_V[2][0][j]-_V_ghost_inflow[2][j]) + (1+kappa)*epsi_minus_x[2][1][j]*(_V[2][1][j]-_V[2][0][j])        );
        
        V_R[2]= _V[2][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][2][j]*(_V[2][2][j]-_V[2][1][j]) + (1+kappa)*epsi_plus_x[2][1][j]*(_V[2][1][j]-_V[2][0][j])        );
        //pressure
        V_L[3]= _V[3][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][0][j]*(_V[3][0][j]-_V_ghost_inflow[3][j]) + (1+kappa)*epsi_minus_x[3][1][j]*(_V[3][1][j]-_V[3][0][j])        );
        
        V_R[3]= _V[3][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][2][j]*(_V[3][2][j]-_V[3][1][j]) + (1+kappa)*epsi_plus_x[3][1][j]*(_V[3][1][j]-_V[3][0][j])        );
        
        ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + V_L[1]*V_L[1]*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + V_R[1]*V_R[1]*0.5;
        
        roe_boundary_in_out(1,j);
        if((j>1)&&(j<im_f-1))
        {
            
            roe_boundary_F_eta(1, j);
            
        }
        
        
        
        
        
        // flux at 0
        
        
        //density
        V_L[0]= _V_ghost_inflow[0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[0][j]*(_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j]) + (1+kappa)*epsi_minus_x[0][0][j]*(_V[0][0][j]-_V_ghost_inflow[0][j]));
        
        V_R[0]= _V[0][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][1][j]*(_V[0][1][j]-_V[0][0][j]) + (1+kappa)*epsi_plus_x[0][0][j]*(_V[0][0][j]-_V_ghost_inflow[0][j])        );
        
        
        
        //u-velocity
        V_L[1]= _V_ghost_inflow[1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[1][j]*(_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j]) + (1+kappa)*epsi_minus_x[1][0][j]*(_V[1][0][j]-_V_ghost_inflow[1][j]));
        
        V_R[1]= _V[1][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][1][j]*(_V[1][1][j]-_V[1][0][j]) + (1+kappa)*epsi_plus_x[1][0][j]*(_V[1][0][j]-_V_ghost_inflow[1][j])        );
        //v-velocity
        V_L[2]= _V_ghost_inflow[2][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[2][j]*(_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j]) + (1+kappa)*epsi_minus_x[2][0][j]*(_V[2][0][j]-_V_ghost_inflow[2][j]));
        
        V_R[2]= _V[2][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][1][j]*(_V[2][1][j]-_V[2][0][j]) + (1+kappa)*epsi_plus_x[2][0][j]*(_V[2][0][j]-_V_ghost_inflow[2][j])        );
        //pressure
        V_L[3]= _V_ghost_inflow[3][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[3][j]*(_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j]) + (1+kappa)*epsi_minus_x[3][0][j]*(_V[3][0][j]-_V_ghost_inflow[3][j]));
        
        V_R[3]= _V[3][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][1][j]*(_V[3][1][j]-_V[3][0][j]) + (1+kappa)*epsi_plus_x[3][0][j]*(_V[3][0][j]-_V_ghost_inflow[3][j])        );
        
        //total enthalpy
        ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + V_L[1]*V_L[1]*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + V_R[1]*V_R[1]*0.5;
        
        
        
        
        roe_boundary_in_out(0,j);
        
        if((j>1)&&(j<im_f-1))
        {
            
            roe_boundary_F_eta(0, j);
            
        }
        
        
        

        
        
        
       
        
        //******************************************
        // outflow boundary conditions
        
        
        
            // P_out = 125000;
            P_out = 0.5*(3*_V[3][im_c][j]-_V[3][im_c-1][j]);
            // std::cout<<"back pressure is "<<P_out<<std::endl;
            
            
            rho_out = 0.5*(3*_U[0][im_c][j]-_U[0][im_c-1][j]);
            
            
            
            u_out = (0.5*(3*_U[1][im_c][j]-_U[1][im_c-1][j]))/rho_out;
            v_out = (0.5*(3*_U[2][im_c][j]-_U[2][im_c-1][j]))/rho_out;
            T_out=P_out/(rho_out*R_air);
          
            
            
           // ///here
            
            
            
            
            //ghost cells
            //1st ghost
            //
            //
            _V_ghost_outflow[0][j]=2*_V[0][im_c][j]-_V[0][im_c-1][j];
            _V_ghost_outflow[1][j]=2*_V[1][im_c][j]-_V[1][im_c-1][j];
            
            _V_ghost_outflow[2][j]=2*_V[2][im_c][j]-_V[2][im_c-1][j];
            
            _V_ghost_outflow[3][j]=2*P_out-_V[3][im_c][j];
            
            //_V_ghost_outflow[3]=2*et[_imin]-et[_imin];  //et
            _V_ghost_outflow[5][j]=2*T[im_c][j]-T[im_c-1][j];   //T
            
            ///second ghost cell
            ///
            ///
            _V_ghost_outflow_2[0][j]=2*_V_ghost_outflow[0][j]- _V[0][im_c][j];
            _V_ghost_outflow_2[1][j]=2*_V_ghost_outflow[1][j]- _V[1][im_c][j];
            _V_ghost_outflow_2[2][j]=2*_V_ghost_outflow[2][j]- _V[2][im_c][j];
            _V_ghost_outflow_2[3][j]=2*_V_ghost_outflow[3][j]-_V[3][im_c][j];
            
            
            //  _V_ghost_outflow_2[3]=2*_V_ghost_outflow[3]-et[im_c];
            _V_ghost_outflow_2[5][j]=2*_V_ghost_outflow[5][j]-T[im_c][j];
            
            //third ghost cell
            ///
            ///
            _V_ghost_outflow_3[0][j]=2*_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j];
            _V_ghost_outflow_3[1][j]=2*_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j];
            _V_ghost_outflow_3[2][j]=2*_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j];
            _V_ghost_outflow_3[3][j]=2*_V_ghost_outflow_2[3][j]-_V_ghost_outflow[3][j];
            //  _V_ghost_outflow_2[3]=2*_V_ghost_outflow[3]-et[im_c];
            _V_ghost_outflow_3[5][j]=2*_V_ghost_outflow_2[5][j]-_V_ghost_outflow[5][j];
            
            
            
            
            //calculate epsi at imf-1 or imc
            
            
            //for density
            denm= copysign(max(abs(_V[0][im_c][j]-_V[0][im_c-1][j]),1e-6),_V[0][im_c][j]-_V[0][im_c-1][j]);
            r_plus= (_V_ghost_outflow[0][j]-_V[0][im_c][j])/(denm);
            r_minus= (_V[0][im_c-1][j]-_V[0][im_c-2][j])/(denm);
            epsi_plus_x[0][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[0][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for velocity
            denm= copysign(max(abs(_V[1][im_c][j]-_V[1][im_c-1][j]),1e-6),_V[1][im_c][j]-_V[1][im_c-1][j]);
            r_plus= (_V_ghost_outflow[1][j]-_V[1][im_c][j])/(denm);
            r_minus= (_V[1][im_c-1][j]-_V[1][im_c-2][j])/(denm);
            epsi_plus_x[1][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[1][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            denm= copysign(max(abs(_V[2][im_c][j]-_V[2][im_c-1][j]),1e-6),_V[2][im_c][j]-_V[2][im_c-1][j]);
            r_plus= (_V_ghost_outflow[2][j]-_V[2][im_c][j])/(denm);
            r_minus= (_V[2][im_c-1][j]-_V[2][im_c-2][j])/(denm);
            epsi_plus_x[2][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[2][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            
            
            
            //calculate epsi for imf or imc+1
            
            
            //for density
            denm= copysign(max(abs(_V_ghost_outflow[0][j]-_V[0][im_c][j]),1e-6),_V_ghost_outflow[0][j]-_V[0][im_c][j]);
            r_plus= (_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j])/(denm);
            r_minus= (_V[0][im_c][j]-_V[0][im_c-1][j])/(denm);
            epsi_plus_x[0][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[0][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for velocity
            denm= copysign(max(abs(_V_ghost_outflow[1][j]-_V[1][im_c][j]),1e-6),_V_ghost_outflow[1][j]-_V[1][im_c][j]);
            r_plus= (_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j])/(denm);
            r_minus= (_V[1][im_c]-_V[1][im_c-1])/(denm);
            epsi_plus_x[1][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[1][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            denm= copysign(max(abs(_V_ghost_outflow[2][j]-_V[2][im_c][j]),1e-6),_V_ghost_outflow[2][j]-_V[2][im_c][j]);
            r_plus= (_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j])/(denm);
            r_minus= (_V[2][im_c][j]-_V[2][im_c-1][j])/(denm);
            epsi_plus_x[2][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[2][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            
            
            //calculate epsi for imf+1 or imc+2
            
            
            //for density
            denm=copysign(max(abs(_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j]),1e-6),_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j]);
            
            r_plus= (_V_ghost_outflow_3[0][j]-_V_ghost_outflow_2[0][j])/(denm);
            r_minus= (_V_ghost_outflow[0][j]-_V[0][im_c][j])/(denm);
            
            epsi_plus_x[0][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[0][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for velocity
            denm
            =copysign(max(abs(_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j]),1e-6),_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j]);
            
            
            r_plus= (_V_ghost_outflow_3[1][j]-_V_ghost_outflow_2[1][j])/(denm);
            r_minus= (_V_ghost_outflow[1][j]-_V[1][im_c][j])/(denm);
            
            epsi_plus_x[1][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[1][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            denm
            =copysign(max(abs(_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j]),1e-6),_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j]);
            
            
            r_plus= (_V_ghost_outflow_3[2][j]-_V_ghost_outflow_2[2][j])/(denm);
            r_minus= (_V_ghost_outflow[2][j]-_V[2][im_c][j])/(denm);
            
            epsi_plus_x[2][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[2][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            
            //calculate fluxes at imf-1 or imc
            
            
            //mistake here
            
            //density
            V_L[0]= _V[0][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][im_c-1][j]*(_V[0][im_c-1][j]-_V[0][im_c-2][j]) + (1+kappa)*epsi_minus_x[0][im_c][j]*(_V[0][im_c][j]-_V[0][im_c-1][j])        );
            
            V_R[0]= _V[0][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][im_c+1][j]*(_V_ghost_outflow[0][j]-_V[0][im_c][j]) + (1+kappa)*epsi_plus_x[0][im_c][j]*(_V[0][im_c][j]-_V[0][im_c-1][j])        );
            //velocity
            V_L[1]= _V[1][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][im_c-1][j]*(_V[1][im_c-1][j]-_V[1][im_c-2][j]) + (1+kappa)*epsi_minus_x[1][im_c][j]*(_V[1][im_c][j]-_V[1][im_c-1][j])        );
            
            V_R[1]=_V[1][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][im_c+1][j]*(_V_ghost_outflow[1][j]-_V[1][im_c][j]) + (1+kappa)*epsi_plus_x[1][im_c][j]*(_V[1][im_c][j]-_V[1][im_c-1][j])        );
            //pressure
            V_L[2]= _V[2][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][im_c-1][j]*(_V[2][im_c-1][j]-_V[2][im_c-2][j]) + (1+kappa)*epsi_minus_x[2][im_c][j]*(_V[2][im_c][j]-_V[2][im_c-1][j])        );
            
            V_R[2]=_V[2][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][im_c+1][j]*(_V_ghost_outflow[2][j]-_V[2][im_c][j]) + (1+kappa)*epsi_plus_x[2][im_c][j]*(_V[2][im_c][j]-_V[2][im_c-1][j])        );
            
            ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + V_L[1]*V_L[1]*0.5;
            ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + V_R[1]*V_R[1]*0.5;
            
            roe_boundary_in_out(im_f-1,j);
            
        if((j>1)&&(j<im_f-1))
        {
            
            roe_boundary_F_eta(im_c-1, j);
            
        }
            
            
            
            
            
            // flux at im_f
            
            
            //density
            V_L[0]= _V[0][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][im_c][j]*(_V[0][im_c][j]-_V[0][im_c-1][j]) + (1+kappa)*epsi_minus_x[0][im_c+1][j]*(_V_ghost_outflow[0][j]-_V[0][im_c][j])        );
            
            V_R[0]= _V_ghost_outflow[0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][im_c+2][j]*(_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j]) + (1+kappa)*epsi_plus_x[0][im_c+1][j]*(_V_ghost_outflow[0][j]-_V[0][im_c][j])        );
            
            
            //velocity
            V_L[1]= _V[1][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][im_c][j]*(_V[1][im_c][j]-_V[1][im_c-1][j]) + (1+kappa)*epsi_minus_x[1][im_c+1][j]*(_V_ghost_outflow[1][j]-_V[1][im_c][j])        );
            
            V_R[1]= _V_ghost_outflow[1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][im_c+2][j]*(_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j]) + (1+kappa)*epsi_plus_x[1][im_c+1][j]*(_V_ghost_outflow[1]-_V[1][im_c])        );
            //pressure
            V_L[2]= _V[2][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][im_c][j]*(_V[2][im_c][j]-_V[2][im_c-1][j]) + (1+kappa)*epsi_minus_x[2][im_c+1][j]*(_V_ghost_outflow[2][j]-_V[2][im_c][j])        );
            
            V_R[2]= _V_ghost_outflow[2][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][im_c+2][j]*(_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j]) + (1+kappa)*epsi_plus_x[2][im_c+1][j]*(_V_ghost_outflow[2][j]-_V[2][im_c][j])        );
            
            ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + V_L[1]*V_L[1]*0.5;
            ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + V_R[1]*V_R[1]*0.5;
            
            
            roe_boundary_in_out(im_f,j);
            
        if((j>1)&&(j<im_f-1))
        {
            
            roe_boundary_F_eta(im_c, j);
            
        }
        
    }

    
    
    roe_boundary_walls(0);
   
    roe_boundary_walls(1);
    
    roe_boundary_walls(im_f);
    
    roe_boundary_walls(im_f-1);
    
    
    
    
    
    
}

void Dv_nozzle::euler_explicit()
{
    double a1, a2,a3,a4;
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            
            
            
            a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
            
            a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
            a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
           
            
            //continuity eqn
            
             _U[0][i][j]=((-_F[0][i+1][j])*a1+a2*_F[0][i][j]-a3*_F_eta[0][i][j+1]+(_F_eta[0][i][j])*a4)*(delta_t_gl/(cell_vol(a1, a2))) + _U[0][i][j];
            
            
            //x-mtm
               
            _U[1][i][j]=((-_F[1][i+1][j])*a1+a2*_F[1][i][j]-a3*_F_eta[1][i][j+1]+(_F_eta[1][i][j])*a4)*(delta_t_gl/(cell_vol(a1, a2))) + _U[1][i][j];
                  
            //y-mtm
            _U[2][i][j]=((-_F[2][i+1][j])*a1+a2*_F[2][i][j]-a3*_F_eta[2][i][j+1]+(_F_eta[2][i][j])*a4)*(delta_t_gl/(cell_vol(a1, a2))) + _U[2][i][j];
            
            
            //energy eqn
                     
            _U[3][i][j]=((-_F[3][i+1][j])*a1+a2*_F[3][i][j]-a3*_F_eta[3][i][j+1]+(_F_eta[3][i][j])*a4)*(delta_t_gl/(cell_vol(a1, a2))) + _U[3][i][j];
            
            
            
        }
    }
    
    
    
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            
            _V[0][i][j]=_U[0][i][j];
            _V[1][i][j]= _U[1][i][j]/_V[0][i][j];
            _V[2][i][j]= _U[2][i][j]/_V[0][i][j];
            //  if (_V[1][i]<0)_V[1][i]=1;
            et[i][j] = _U[3][i][j]/ _V[0][i][j];
            //  if (et[i]<200)et[i]=200;
            T[i][j] = (et[i][j] - 0.5* (_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j] ))*((gamma-1)/R_air);
            //  if (T[i]<50)T[i]=50;
            
            
            _V[3][i][j]=_V[0][i][j]*R_air*T[i][j];
            //   if (_V[2][i]<3000)_V[2][i]=3000;
            ht[i][j]=(gamma/(gamma-1))*(_V[3][i][j]/_V[0][i][j]) + (_V[2][i][j]*_V[2][i][j]+_V[1][i][j]*_V[1][i][j])*0.5;
            _mach[i][j]=sqrt(_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j])/sqrt(gamma*R_air*abs(T[i][j]));
            
            
            // cout<< "pressure is    "<<_V[2][i]<<endl;
        }
    }
        //calculate epsi- kasi direction kasi corresponds to x kinda
        
        for (int j=2;j<im_c;j++)
        {
            for (int i=2;i<im_c;i++)
            {
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i-1][j]),1e-6),_V[0][i][j]-_V[0][i-1][j]);
                r_plus= (_V[0][i+1][j]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i-1][j]-_V[0][i-2][j])/(denm);
                epsi_plus_x[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i-1][j]),1e-6),_V[1][i][j]-_V[1][i-1][j]);
                r_plus= (_V[1][i+1][j]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i-1][j]-_V[1][i-2][j])/(denm);
                epsi_plus_x[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i-1][j]),1e-6),_V[2][i][j]-_V[2][i-1][j]);
                r_plus= (_V[2][i+1][j]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i-1][j]-_V[2][i-2][j])/(denm);
                epsi_plus_x[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i-1][j]),1e-6),_V[3][i][j]-_V[3][i-1][j]);
                r_plus= (_V[3][i+1][j]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i-1][j]-_V[3][i-2][j])/(denm);
                epsi_plus_x[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                
            }
            
        }
        
        
        //calculate epsi- eta direction kasi corresponds to y kinda
        
        for (int j=2;j<im_c;j++)
        {
            for (int i=2;i<im_c;i++)
            {
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i][j-1]),1e-6),_V[0][i][j]-_V[0][i][j-1]);
                r_plus= (_V[0][i][j+1]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i][j-1]-_V[0][i][j-2])/(denm);
                epsi_plus_y[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i][j-1]),1e-6),_V[1][i][j]-_V[1][i][j-1]);
                r_plus= (_V[1][i][j+1]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i][j-1]-_V[1][i][j-2])/(denm);
                epsi_plus_y[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i][j-1]),1e-6),_V[2][i][j]-_V[2][i][j-1]);
                r_plus= (_V[2][i][j+1]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i][j-1]-_V[2][i][j-2])/(denm);
                epsi_plus_y[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i][j-1]),1e-6),_V[3][i][j]-_V[3][i][j-1]);
                r_plus= (_V[3][i][j+1]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i][j-1]-_V[3][i][j-2])/(denm);
                epsi_plus_y[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                
            }
            
        }
   
    
    
    
    
    roe_flux();
    
    
    
    
    
}
double Dv_nozzle::area(double x)
{
    _area = 0.2 + 0.4*(1 + sin(3.14159265359*(x-0.5)));
    
    
    return _area;
    
}
double Dv_nozzle::d_area(double x)
{
    _d_area = 0.4*3.14159265359*cos(3.14159265359*(x-0.5));
    
    
    return _d_area;
    
}
void Dv_nozzle::time_step()
{
     double a1, a2, a3, a4;
                
     for (int j=0;j<=im_c;j++)
         for (int i =0;i<=im_c;i++)
         {
             
             
             
              a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
              a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
              a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
              a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
             
             
             
             
             
         //for two x faces
             
             //for kase
             nx=n_x(y_n[i+1][j+1], y_n[i+1][j],a1);
             
             nx= 0.5*(nx + n_x(y_n[i][j+1],y_n[i][j],a2));
             
             
             ny=n_y(x_n[i+1][j+1], x_n[i+1][j], a1);
             
             ny= 0.5*(ny + n_y(x_n[i][j+1], x_n[i][j], a2));
             a[i][j]=sqrt(abs((gamma-1)*(ht[i][j]-0.5*(_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j]))));
             
             
            
             eigen_max[0]=abs(_V[1][i][j]*nx+_V[2][i][j]*ny) + a[i][j];
             
             // for eta
             nx=n_x(y_n[i+1][j+1], y_n[i][j+1], a3);
             
             nx= 0.5*(nx + n_x(y_n[i+1][j],y_n[i][j],a4));
             
             
             ny=n_y(x_n[i+1][j+1], x_n[i][j+1], a3);
             
             ny= 0.5*(ny + n_y(x_n[i+1][j], x_n[i][j], a4));
             
             
             eigen_max[1]=abs(_V[1][i][j]*nx+_V[2][i][j]*ny) + a[i][j];
             
             
             
             delta_t[i][j] = CFL*(cell_vol(a1, a2)/(eigen_max[0]*cell_vol(a1, a2)  + eigen_max[1]*cell_vol(a3, a4)    )     );
             
             
             
             delta_t_gl=min(delta_t[i][j],delta_t_gl);
             
             
         }
     
     
 }


void Dv_nozzle::print_res(){
    rho_vs_x<<"######### x ######   "<<"#########rho######   "<<std::endl;
    u_vs_x<<"######### x ######   "<<"#########ux######   "<<std::endl;
    p_vs_x<<"######### x ######   "<<"#########P######   "<<std::endl;
    for (int i= 0;i<=im_c;i++)
    {
        
        rho_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[0][i]<<std::endl;
        
        
        
        u_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(14)<< _V[1][i]<<std::endl;
        
        
        
        p_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[2][i]<<std::endl;
        
        
        
    }
    
    
    
}
double Dv_nozzle::L2norm(int iter){
    
    
    rL2norm[0]=0;
    rL2norm[1]=0;
    rL2norm[2]=0;
    rL2norm[3]=0;
    
    double a1, a2,a3,a4;
    
    
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            
            
            a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
            
            a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
            a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
        
            _R_iter[0]=-1*(((-_F[0][i+1][j])*a1+a2*_F[0][i][j]-a3*_F_eta[0][i][j+1]+(_F_eta[0][i][j])*a4));
            
            rL2norm[0] = rL2norm[0] +_R_iter[0]*_R_iter[0];
            
            
            _R_iter[1]=-1*(((-_F[1][i+1][j])*a1+a2*_F[1][i][j]-a3*_F_eta[1][i][j+1]+(_F_eta[1][i][j])*a4)   )  ;
            
            rL2norm[1] = rL2norm[1] +  _R_iter[1]* _R_iter[1];
            
            _R_iter[2]=-1*(((-_F[2][i+1][j])*a1+a2*_F[2][i][j]-a3*_F_eta[2][i][j+1]+(_F_eta[2][i][j])*a4) )  ;
            
            rL2norm[2] = rL2norm[2] +  _R_iter[2]* _R_iter[2];
            
            _R_iter[3]=-1*(((-_F[3][i+1][j])*a1+a2*_F[3][i][j]-a3*_F_eta[3][i][j+1]+(_F_eta[3][i][j])*a4))  ;
            
            rL2norm[3] = rL2norm[3] +  _R_iter[3]* _R_iter[3];
            
            
        }
    }
    double imax_c = _imax_c*1.0;
    
    rL2norm[0]=sqrt(rL2norm[0]/(imax_c*imax_c))/_rL2initial[0];
    rL2norm[1]=sqrt(rL2norm[1]/(imax_c*imax_c))/_rL2initial[1];
    rL2norm[2]=sqrt(rL2norm[2]/(imax_c*imax_c))/_rL2initial[2];
    rL2norm[3]=sqrt(rL2norm[3]/(imax_c*imax_c))/_rL2initial[3];
    
    /*
     std::cout<< " R eqn1 "<< _R_iter[0]<<std::endl;
     std::cout<< " R eqn2 "<< _R_iter[1]<<std::endl;
     std::cout<< " R eqn3 "<< _R_iter[2]<<std::endl;
     std::cout<< " L eqn1 "<< rL2norm[0]<<std::endl;
     std::cout<< " L eqn2 "<< rL2norm[1]<<std::endl;
     std::cout<< " L eqn3 "<< rL2norm[2]<<std::endl;
     std::cout<< " L inti eqn1 "<< _rL2initial[0]<<std::endl;
     std::cout<< " L inti eqn2 "<< _rL2initial[1]<<std::endl;
     std::cout<< " L inti eqn3 "<< _rL2initial[2]<<std::endl;
     std::cout<< " u is  "<< _V[1][9]<<std::endl;
     */
    if(iter%10==0)L_vs_Iter<<std::setprecision(5)<<iter<<","<<std::setprecision(5)<<rL2norm[0]<<","<<std::setprecision(5)<<rL2norm[1]<<","<<std::setprecision(5)<<rL2norm[2]<<std::endl;
    
 return std::max(rL2norm[3],max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])));
    
    
}
void Dv_nozzle::rL2initial()


{
    double a1, a2,a3,a4;
    _rL2initial[0]=0;
    _rL2initial[1]=0;
    _rL2initial[2]=0;
    _rL2initial[3]=0;
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
            
            a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
            a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
        
            _R[0]= -1*(((-_F[0][i+1][j])*a1+a2*_F[0][i][j]-a3*_F_eta[0][i][j+1]+(_F_eta[0][i][j])*a4)  )  ;
            
            
            _rL2initial[0] = _rL2initial[0] + _R[0]*_R[0];
            
            
            _R[1]= -1*(((-_F[1][i+1][j])*a1+a2*_F[1][i][j]-a3*_F_eta[1][i][j+1]+(_F_eta[1][i][j])*a4)  )  ;
            
            _rL2initial[1] = _rL2initial[1] + _R[1]*_R[1];
            
            _R[2]= -1*(((-_F[2][i+1][j])*a1+a2*_F[2][i][j]-a3*_F_eta[2][i][j+1]+(_F_eta[2][i][j])*a4) )  ;
            
            _rL2initial[2] = _rL2initial[2] + _R[2]*_R[2];
            _R[3]= -1*(((-_F[3][i+1][j])*a1+a2*_F[3][i][j]-a3*_F_eta[3][i][j+1]+(_F_eta[3][i][j])*a4))  ;
            
            _rL2initial[3] = _rL2initial[3] + _R[3]*_R[3];
        }
    }
    double imax_c = _imax_c*1.0;
    _rL2initial[0] = sqrt(_rL2initial[0]/(imax_c*imax_c));
    _rL2initial[1] = sqrt(_rL2initial[1]/(imax_c*imax_c));
    _rL2initial[2] = sqrt(_rL2initial[2]/(imax_c*imax_c));
    _rL2initial[3] = sqrt(_rL2initial[3]/(imax_c*imax_c));
    
    std::cout<<" r initial is "<<_rL2initial[0]<<std::endl;
    
    
}
void Dv_nozzle::roe_flux()
{
    
    //calculate the flux
for (int j=2;j<=im_c-1;j++){
       for (int i=2;i<=im_c-1;i++)
       {
           
           
           //density
           V_L[0]= _V[0][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][i-1][j]*(_V[0][i-1][j]-_V[0][i-2][j]) + (1+kappa)*epsi_minus_x[0][i][j]*(_V[0][i][j]-_V[0][i-1][j])        );
           
           V_R[0]= _V[0][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][i+1][j]*(_V[0][i+1][j]-_V[0][i][j]) + (1+kappa)*epsi_plus_x[0][i][j]*(_V[0][i][j]-_V[0][i-1][j])        );
           //x-velocity
           V_L[1]= _V[1][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][i-1][j]*(_V[1][i-1][j]-_V[1][i-2][j]) + (1+kappa)*epsi_minus_x[1][i][j]*(_V[1][i][j]-_V[1][i-1][j])        );
           
           V_R[1]= _V[1][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][i+1][j]*(_V[1][i+1][j]-_V[1][i][j]) + (1+kappa)*epsi_plus_x[1][i][j]*(_V[1][i][j]-_V[1][i-1][j])        );
           //y-velocity
           V_L[2]= _V[2][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][i-1][j]*(_V[2][i-1][j]-_V[2][i-2][j]) + (1+kappa)*epsi_minus_x[2][i][j]*(_V[2][i][j]-_V[2][i-1][j])        );
           
           V_R[2]= _V[2][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][i+1][j]*(_V[2][i+1][j]-_V[2][i][j]) + (1+kappa)*epsi_plus_x[2][i][j]*(_V[2][i][j]-_V[2][i-1][j])        );
           //pressure
           V_L[3]= _V[3][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][i-1][j]*(_V[3][i-1][j]-_V[3][i-2][j]) + (1+kappa)*epsi_minus_x[3][i][j]*(_V[3][i][j]-_V[3][i-1][j])        );
           
           V_R[3]= _V[3][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][i+1][j]*(_V[3][i+1][j]-_V[3][i][j]) + (1+kappa)*epsi_plus_x[3][i][j]*(_V[3][i][j]-_V[3][i-1][j])        );
           //total enthalpy
           ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + ((V_L[1]*V_L[1])+(V_L[2]*V_L[2]))*0.5;
           ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + ((V_R[1]*V_R[1])+(V_R[2]*V_R[2]))*0.5;
           
           R= sqrt(V_R[0]/V_L[0]);
           rho_roe=R*V_L[0];
           u_roe=(R*V_R[1]+V_L[1])/(R+1);
           v_roe=(R*V_R[2]+V_L[2])/(R+1);
           ht_roe= (R*ht_R+ht_L)/(R+1);
           a_roe=sqrt(abs((gamma-1)*(ht_roe-0.5*(u_roe*u_roe+v_roe*v_roe))));
           
           double a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
           nx=n_x(y_n[i][j+1], y_n[i][j], a2);
           ny=n_y(x_n[i][j+1], x_n[i][j], a2);
           T_vel_roe=nx*u_roe+ny*v_roe;
           
           
           r_1[0]=1;
           r_1[1]=u_roe;
           r_1[2]=v_roe;
           r_1[3]=0.5*(u_roe*u_roe+v_roe*v_roe);
           
           
           
           r_2[0]=0;
           r_2[1]=ny*rho_roe;
           r_2[2]=-nx*rho_roe;
           r_2[3]=rho_roe*(ny*u_roe-nx*v_roe);
           
           
           
           
           
           r_3[0]=rho_roe/(2*a_roe);
           r_3[1]=r_3[0]*(u_roe+nx*a_roe);
           r_3[2]=r_3[0]*(v_roe+ny*a_roe);
           r_3[3]=r_3[0]*(ht_roe+T_vel_roe*a_roe);
           
           
           r_4[0]=-rho_roe/(2*a_roe);
           r_4[1]=r_4[0]*(u_roe-nx*a_roe);
           r_4[2]=r_4[0]*(v_roe-ny*a_roe);
           r_4[3]=r_4[0]*(ht_roe-T_vel_roe*a_roe);
           
           lambda[0]=T_vel_roe;
           lambda[1]=T_vel_roe;
           lambda[2]=T_vel_roe+a_roe;
           lambda[3]=T_vel_roe-a_roe;
           
           
           if(abs(lambda[0])<0.1*2*a_roe)lambda[0]=(lambda[0]*lambda[0])/(0.1*4*a_roe)+0.1*a_roe;
           if(abs(lambda[1])<0.1*2*a_roe)lambda[1]=(lambda[1]*lambda[1])/(0.1*4*a_roe)+0.1*a_roe;
           if(abs(lambda[2])<0.1*2*a_roe)lambda[2]=(lambda[2]*lambda[2])/(0.1*4*a_roe)+0.1*a_roe;
           if(abs(lambda[3])<0.1*2*a_roe)lambda[3]=(lambda[3]*lambda[3])/(0.1*4*a_roe)+0.1*a_roe;
           
           delta_rho=V_R[0]-V_L[0];
           delta_u=V_R[1]-V_L[1];
           delta_v=V_R[2]-V_L[2];
           delta_P=V_R[3]-V_L[3];
           
           
           delta_w[0]= delta_rho-delta_P/(a_roe*a_roe);
           delta_w[1]=ny*delta_u-nx*delta_v;
           delta_w[2]=nx*delta_u+ny*delta_v+(delta_P)/(rho_roe*a_roe);
           delta_w[3]=nx*delta_u+ny*delta_v-(delta_P)/(rho_roe*a_roe);
           
           T_vel_L=V_L[1]*nx+V_L[2]*ny;
           T_vel_R=V_R[1]*nx+V_R[2]*ny;
           
           _F_L[0][i][j] =V_L[0]* T_vel_L ;
           _F_L[1][i][j]= V_L[0]*V_L[1]*T_vel_L + V_L[3]*nx;
           _F_L[2][i][j]= V_L[0]*V_L[2]*T_vel_L + V_L[3]*ny;
           _F_L[3][i][j]= V_L[0]*ht_L*T_vel_L;
           
           
           _F_R[0][i][j] =V_R[0]* T_vel_R ;
           _F_R[1][i][j]= V_R[0]*V_R[1]*T_vel_R + V_R[3]*nx;
           _F_R[2][i][j]= V_R[0]*V_R[2]*T_vel_R + V_R[3]*ny;
           _F_R[3][i][j]= V_R[0]*ht_R*T_vel_R;
           
       
           
           
          //flux for kasi faces
           
           
           _F[0][i][j]=0.5*(_F_L[0][i][j]+_F_R[0][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[0] + abs(lambda[1])*delta_w[1]*r_2[0] + abs(lambda[2])*delta_w[2]*r_3[0] + abs(lambda[3])*delta_w[3]*r_4[0]     );
           
           
           
           _F[1][i][j]=0.5*(_F_L[1][i][j]+_F_R[1][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[1] + abs(lambda[1])*delta_w[1]*r_2[1] + abs(lambda[2])*delta_w[2]*r_3[1] + abs(lambda[3])*delta_w[3]*r_4[1]     );
           
           
           
           _F[2][i][j]=0.5*(_F_L[2][i][j]+_F_R[2][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[2] + abs(lambda[1])*delta_w[1]*r_2[2] + abs(lambda[2])*delta_w[2]*r_3[2] + abs(lambda[3])*delta_w[3]*r_4[2]     );
           
           _F[3][i][j]=0.5*(_F_L[3][i][j]+_F_R[3][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[3] + abs(lambda[1])*delta_w[1]*r_2[3] + abs(lambda[2])*delta_w[2]*r_3[3]  + abs(lambda[3])*delta_w[3]*r_4[3]    );
           
           
           
      // flux for eta faces
           
           
           //density
           V_D[0]= _V[0][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[0][i][j-1]*(_V[0][i][j-1]-_V[0][i][j-2]) + (1+kappa)*epsi_minus_y[0][i][j]*(_V[0][i][j]-_V[0][i][j-1])        );
           
           V_U[0]= _V[0][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[0][i][j+1]*(_V[0][i][j+1]-_V[0][i][j]) + (1+kappa)*epsi_plus_y[0][i][j]*(_V[0][i][j]-_V[0][i][j-1])        );
           //x-velocity
           V_D[1]= _V[1][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[1][i][j-1]*(_V[1][i][j-1]-_V[1][i][j-2]) + (1+kappa)*epsi_minus_y[1][i][j]*(_V[1][i][j]-_V[1][i][j-1])        );
           
           V_U[1]= _V[1][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[1][i][j+1]*(_V[1][i][j+1]-_V[1][i][j]) + (1+kappa)*epsi_plus_y[1][i][j]*(_V[1][i][j]-_V[1][i][j-1])        );
           //y-velocity
           V_D[2]= _V[2][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[2][i][j-1]*(_V[2][i][j-1]-_V[2][i][j-2]) + (1+kappa)*epsi_minus_y[2][i][j]*(_V[2][i][j]-_V[2][i][j-1])        );
           
           V_U[2]= _V[2][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[2][i][j+1]*(_V[2][i][j+1]-_V[2][i][j]) + (1+kappa)*epsi_plus_y[2][i][j]*(_V[2][i][j]-_V[2][i][j-1])        );
           //pressure
           V_D[3]= _V[3][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[3][i][j-1]*(_V[3][i][j-1]-_V[3][i][j-2]) + (1+kappa)*epsi_minus_y[3][i][j]*(_V[3][i][j]-_V[3][i][j-1])        );
           
           V_U[3]= _V[3][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[3][i][j+1]*(_V[3][i][j+1]-_V[3][i][j]) + (1+kappa)*epsi_plus_y[3][i][j]*(_V[3][i][j]-_V[3][i][j-1])        );
           //total enthalpy
           ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
           ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
           
           
           
           R= sqrt(V_U[0]/V_D[0]);
           rho_roe=R*V_D[0];
           u_roe=(R*V_U[1]+V_D[1])/(R+1);
           v_roe=(R*V_U[2]+V_D[2])/(R+1);
           ht_roe= (R*ht_U+ht_D)/(R+1);
           a_roe=sqrt(abs((gamma-1)*(ht_roe-0.5*(u_roe*u_roe+v_roe*v_roe))));
           
           double a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
           
           
           nx=n_x(y_n[i+1][j], y_n[i][j], a4);
           ny=n_y(x_n[i+1][j], x_n[i][j], a4);
           
           T_vel_roe=nx*u_roe+ny*v_roe;
           
           
           r_1[0]=1;
           r_1[1]=u_roe;
           r_1[2]=v_roe;
           r_1[3]=0.5*(u_roe*u_roe+v_roe*v_roe);
           
           
           
           r_2[0]=0;
           r_2[1]=ny*rho_roe;
           r_2[2]=-nx*rho_roe;
           r_2[3]=rho_roe*(ny*u_roe-nx*v_roe);
           
           
           
           
           
           r_3[0]=rho_roe/(2*a_roe);
           r_3[1]=r_3[0]*(u_roe+nx*a_roe);
           r_3[2]=r_3[0]*(v_roe+ny*a_roe);
           r_3[3]=r_3[0]*(ht_roe+T_vel_roe*a_roe);
           
           
           r_4[0]=-rho_roe/(2*a_roe);
           r_4[1]=r_4[0]*(u_roe-nx*a_roe);
           r_4[2]=r_4[0]*(v_roe-ny*a_roe);
           r_4[3]=r_4[0]*(ht_roe-T_vel_roe*a_roe);
           
           lambda[0]=T_vel_roe;
           lambda[1]=T_vel_roe;
           lambda[2]=T_vel_roe+a_roe;
           lambda[3]=T_vel_roe-a_roe;
           
           
           if(abs(lambda[0])<0.1*2*a_roe)lambda[0]=(lambda[0]*lambda[0])/(0.1*4*a_roe)+0.1*a_roe;
           if(abs(lambda[1])<0.1*2*a_roe)lambda[1]=(lambda[1]*lambda[1])/(0.1*4*a_roe)+0.1*a_roe;
           if(abs(lambda[2])<0.1*2*a_roe)lambda[2]=(lambda[2]*lambda[2])/(0.1*4*a_roe)+0.1*a_roe;
           if(abs(lambda[3])<0.1*2*a_roe)lambda[3]=(lambda[3]*lambda[3])/(0.1*4*a_roe)+0.1*a_roe;
    
           delta_rho=V_U[0]-V_D[0];
           delta_u=V_U[1]-V_D[1];
           delta_v=V_U[2]-V_D[2];
           delta_P=V_U[3]-V_D[3];
           
           
           delta_w[0]= delta_rho-delta_P/(a_roe*a_roe);
           delta_w[1]=ny*delta_u-nx*delta_v;
           delta_w[2]=nx*delta_u+ny*delta_v+(delta_P)/(rho_roe*a_roe);
           delta_w[3]=nx*delta_u+ny*delta_v-(delta_P)/(rho_roe*a_roe);
           
           T_vel_D=V_D[1]*nx+V_D[2]*ny;
           T_vel_R=V_U[1]*nx+V_U[2]*ny;
           
           _F_D[0][i][j] =V_D[0]* T_vel_D ;
           _F_D[1][i][j]= V_D[0]*V_D[1]*T_vel_D + V_D[3]*nx;
           _F_D[2][i][j]= V_D[0]*V_D[2]*T_vel_D + V_D[3]*ny;
           _F_D[3][i][j]= V_D[0]*ht_D*T_vel_D;
           
           
           _F_U[0][i][j] =V_U[0]*T_vel_U ;
           _F_U[1][i][j]= V_U[0]*V_U[1]*T_vel_U + V_U[3]*nx;
           _F_U[2][i][j]= V_U[0]*V_U[2]*T_vel_U + V_U[3]*ny;
           _F_U[3][i][j]= V_U[0]*ht_U*T_vel_U;
           
       
           
           _F_eta[0][i][j]=0.5*(_F_D[0][i][j]+_F_U[0][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[0] + abs(lambda[1])*delta_w[1]*r_2[0] + abs(lambda[2])*delta_w[2]*r_3[0] + abs(lambda[3])*delta_w[3]*r_4[0]     );
           
           
           
           _F_eta[1][i][j]=0.5*(_F_D[1][i][j]+_F_U[1][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[1] + abs(lambda[1])*delta_w[1]*r_2[1] + abs(lambda[2])*delta_w[2]*r_3[1] + abs(lambda[3])*delta_w[3]*r_4[1]     );
           
           
           
           _F_eta[2][i][j]=0.5*(_F_D[2][i][j]+_F_U[2][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[2] + abs(lambda[1])*delta_w[1]*r_2[2] + abs(lambda[2])*delta_w[2]*r_3[2] + abs(lambda[3])*delta_w[3]*r_4[2]     );
           
           _F_eta[3][i][j]=0.5*(_F_D[3][i][j]+_F_U[3][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[3] + abs(lambda[1])*delta_w[1]*r_2[3] + abs(lambda[2])*delta_w[2]*r_3[3]  + abs(lambda[3])*delta_w[3]*r_4[3]    );
           
           
           
           
           
           
           
           
       }
   }
    
}
void Dv_nozzle::roe_boundary_in_out (int i,int j)

{
    
    
    //flux in kasi direction
    R= sqrt(V_R[0]/V_L[0]);
    rho_roe=R*V_L[0];
    u_roe=(R*V_R[1]+V_L[1])/(R+1);
    v_roe=(R*V_R[2]+V_L[2])/(R+1);
    ht_roe= (R*ht_R+ht_L)/(R+1);
    a_roe=sqrt(abs((gamma-1)*(ht_roe-0.5*(u_roe*u_roe+v_roe*v_roe))));
    
    double a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
    nx=n_x(y_n[i][j+1], y_n[i][j], a2);
    ny=n_y(x_n[i][j+1], x_n[i][j], a2);
    T_vel_roe=nx*u_roe+ny*v_roe;
    
    
    r_1[0]=1;
    r_1[1]=u_roe;
    r_1[2]=v_roe;
    r_1[3]=0.5*(u_roe*u_roe+v_roe*v_roe);
    
    
    
    r_2[0]=0;
    r_2[1]=ny*rho_roe;
    r_2[2]=-nx*rho_roe;
    r_2[3]=rho_roe*(ny*u_roe-nx*v_roe);
    
    
    
    
    
    r_3[0]=rho_roe/(2*a_roe);
    r_3[1]=r_3[0]*(u_roe+nx*a_roe);
    r_3[2]=r_3[0]*(v_roe+ny*a_roe);
    r_3[3]=r_3[0]*(ht_roe+T_vel_roe*a_roe);
    
    
    r_4[0]=-rho_roe/(2*a_roe);
    r_4[1]=r_4[0]*(u_roe-nx*a_roe);
    r_4[2]=r_4[0]*(v_roe-ny*a_roe);
    r_4[3]=r_4[0]*(ht_roe-T_vel_roe*a_roe);
    
    lambda[0]=T_vel_roe;
    lambda[1]=T_vel_roe;
    lambda[2]=T_vel_roe+a_roe;
    lambda[3]=T_vel_roe-a_roe;
    
    
    if(abs(lambda[0])<0.1*2*a_roe)lambda[0]=(lambda[0]*lambda[0])/(0.1*4*a_roe)+0.1*a_roe;
    if(abs(lambda[1])<0.1*2*a_roe)lambda[1]=(lambda[1]*lambda[1])/(0.1*4*a_roe)+0.1*a_roe;
    if(abs(lambda[2])<0.1*2*a_roe)lambda[2]=(lambda[2]*lambda[2])/(0.1*4*a_roe)+0.1*a_roe;
    if(abs(lambda[3])<0.1*2*a_roe)lambda[3]=(lambda[3]*lambda[3])/(0.1*4*a_roe)+0.1*a_roe;
    
    delta_rho=V_R[0]-V_L[0];
    delta_u=V_R[1]-V_L[1];
    delta_v=V_R[2]-V_L[2];
    delta_P=V_R[3]-V_L[3];
    
    
    delta_w[0]= delta_rho-delta_P/(a_roe*a_roe);
    delta_w[1]=ny*delta_u-nx*delta_v;
    delta_w[2]=nx*delta_u+ny*delta_v+(delta_P)/(rho_roe*a_roe);
    delta_w[3]=nx*delta_u+ny*delta_v-(delta_P)/(rho_roe*a_roe);
    
    T_vel_L=V_L[1]*nx+V_L[2]*ny;
    T_vel_R=V_R[1]*nx+V_R[2]*ny;
    
    _F_L[0][i][j] =V_L[0]* T_vel_L ;
    _F_L[1][i][j]= V_L[0]*V_L[1]*T_vel_L + V_L[3]*nx;
    _F_L[2][i][j]= V_L[0]*V_L[2]*T_vel_L + V_L[3]*ny;
    _F_L[3][i][j]= V_L[0]*ht_L*T_vel_L;
    
    
    _F_R[0][i][j] =V_R[0]* T_vel_R ;
    _F_R[1][i][j]= V_R[0]*V_R[1]*T_vel_R + V_R[3]*nx;
    _F_R[2][i][j]= V_R[0]*V_R[2]*T_vel_R + V_R[3]*ny;
    _F_R[3][i][j]= V_R[0]*ht_R*T_vel_R;
    
    //flux for kasi faces
     
     
     _F[0][i][j]=0.5*(_F_L[0][i][j]+_F_R[0][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[0] + abs(lambda[1])*delta_w[1]*r_2[0] + abs(lambda[2])*delta_w[2]*r_3[0] + abs(lambda[3])*delta_w[3]*r_4[0]     );
     
     
     
     _F[1][i][j]=0.5*(_F_L[1][i][j]+_F_R[1][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[1] + abs(lambda[1])*delta_w[1]*r_2[1] + abs(lambda[2])*delta_w[2]*r_3[1] + abs(lambda[3])*delta_w[3]*r_4[1]     );
     
     
     
     _F[2][i][j]=0.5*(_F_L[2][i][j]+_F_R[2][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[2] + abs(lambda[1])*delta_w[1]*r_2[2] + abs(lambda[2])*delta_w[2]*r_3[2] + abs(lambda[3])*delta_w[3]*r_4[2]     );
     
     _F[3][i][j]=0.5*(_F_L[3][i][j]+_F_R[3][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[3] + abs(lambda[1])*delta_w[1]*r_2[3] + abs(lambda[2])*delta_w[2]*r_3[3]  + abs(lambda[3])*delta_w[3]*r_4[3]    );
    
   
    
}
                        
void Dv_nozzle::roe_boundary_walls(int j)
{
                
    for (int i=0;i<=im_c;i++)
                
  {
      
    
      
      double P_wall=0;
      
      
      if(j==0)
      {
          
          P_wall =_V[3][i][0] - 0.5*(_V[3][i][1]-_V[3][i][0]);
          
        
          
          
          double a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
          nx= n_x(y_n[i+1][j], y_n[i][j], a4);
          ny= n_y(x_n[i+1][j], x_n[i][j], a4);
          
          _F_eta[0][i][j]=0;
          
          _F_eta[1][i][j]=nx*P_wall;
          
          _F_eta[2][i][j]=ny*P_wall;
          
          _F_eta[3][i][j]=0;
          
          if((i>1)&&(i<im_f-1))
          {
              
              roe_boundary_F_kasi(i, 0);
              
          }
      }
      
         
          if (j==1){
              
              
              
              _F_eta[0][i][1]= 0.5*(_F_eta[0][i][0]+_F_eta[0][i][2]);
              _F_eta[1][i][1]= 0.5*(_F_eta[1][i][0]+_F_eta[1][i][2]);
              _F_eta[2][i][1]= 0.5*(_F_eta[2][i][0]+_F_eta[2][i][2]);
              _F_eta[3][i][1]= 0.5*(_F_eta[3][i][0]+_F_eta[3][i][2]);
              if((i>1)&&(i<im_f-1))
              {
                  
                  roe_boundary_F_kasi(i, 1);
                  
              }
              
          }
          if(j==im_f)
          
          {
              
              
              P_wall =_V[3][i][im_c] - 0.5*(_V[3][i][im_c-1]-_V[3][i][im_c]);
              
              
              double a4=face_area(x_n[i+1][im_f], x_n[i][im_f],y_n[i+1][im_f], y_n[i][im_f]);
              nx= n_x(y_n[i+1][im_f], y_n[i][im_f], a4);
              ny= n_y(x_n[i+1][im_f], x_n[i][im_f], a4);
              
              _F_eta[0][i][im_f]=0;
              
              _F_eta[1][i][im_f]=nx*P_wall;
              
              _F_eta[2][i][im_f]=ny*P_wall;
              
              _F_eta[3][i][im_f]=0;
              
              if((i>1)&&(i<im_f-1))
              {
                  
                  roe_boundary_F_kasi(i, im_c);
                  
              }
          }
          if (j==im_f-1)
          {
              
              
              
              
              _F_eta[0][i][im_f-1]= 0.5*(_F_eta[0][i][im_f]+_F_eta[0][i][im_f-2]);
              _F_eta[1][i][im_f-1]= 0.5*(_F_eta[1][i][im_f]+_F_eta[1][i][im_f-2]);
              _F_eta[2][i][im_f-1]= 0.5*(_F_eta[2][i][im_f]+_F_eta[2][i][im_f-2]);
              _F_eta[3][i][im_f-1]= 0.5*(_F_eta[3][i][im_f]+_F_eta[3][i][im_f-2]);
              if((i>1)&&(i<im_f-1))
              {
                  
                  roe_boundary_F_kasi(i, im_c-1);
                  
              }
              
          }
          
      
          
          
     
          
          
          
          
      
      
      
      
  }
      

                
                
}

void Dv_nozzle::mesh_nodes()
{
   
    
    
    
    
    
    
}


double Dv_nozzle:: face_area(double x_1,double x_2,double y_1,double y_2)


{
    
    double f_area;
    
    f_area= sqrt((x_1-x_2)*(x_1-x_2)+(y_1-y_2)*(y_1-y_2));
    
    return f_area;
    
    
}
double Dv_nozzle:: n_x(double y_1, double y_2,double f_area)

{
    double n_x;
    
    n_x = (y_1-y_2)/f_area;
    
    return n_x;
    
}
double Dv_nozzle:: n_y(double x_1, double x_2,double f_area)

{
    double n_y;
    
    n_y = (x_1-x_2)/f_area;
    
    return n_y;
    
}
double Dv_nozzle:: cell_vol(double area_1,double area_2)


{
    double vol;
    
    vol = 0.5*(area_1+area_2);
    
    return vol;
    
}
void Dv_nozzle::roe_boundary_F_kasi (int i,int j)

{
    
    //density
    V_L[0]= _V[0][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][i-1][j]*(_V[0][i-1][j]-_V[0][i-2][j]) + (1+kappa)*epsi_minus_x[0][i][j]*(_V[0][i][j]-_V[0][i-1][j])        );
    
    V_R[0]= _V[0][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][i+1][j]*(_V[0][i+1][j]-_V[0][i][j]) + (1+kappa)*epsi_plus_x[0][i][j]*(_V[0][i][j]-_V[0][i-1][j])        );
    //x-velocity
    V_L[1]= _V[1][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][i-1][j]*(_V[1][i-1][j]-_V[1][i-2][j]) + (1+kappa)*epsi_minus_x[1][i][j]*(_V[1][i][j]-_V[1][i-1][j])        );
    
    V_R[1]= _V[1][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][i+1][j]*(_V[1][i+1][j]-_V[1][i][j]) + (1+kappa)*epsi_plus_x[1][i][j]*(_V[1][i][j]-_V[1][i-1][j])        );
    //y-velocity
    V_L[2]= _V[2][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][i-1][j]*(_V[2][i-1][j]-_V[2][i-2][j]) + (1+kappa)*epsi_minus_x[2][i][j]*(_V[2][i][j]-_V[2][i-1][j])        );
    
    V_R[2]= _V[2][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][i+1][j]*(_V[2][i+1][j]-_V[2][i][j]) + (1+kappa)*epsi_plus_x[2][i][j]*(_V[2][i][j]-_V[2][i-1][j])        );
    //pressure
    V_L[3]= _V[3][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][i-1][j]*(_V[3][i-1][j]-_V[3][i-2][j]) + (1+kappa)*epsi_minus_x[3][i][j]*(_V[3][i][j]-_V[3][i-1][j])        );
    
    V_R[3]= _V[3][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][i+1][j]*(_V[3][i+1][j]-_V[3][i][j]) + (1+kappa)*epsi_plus_x[3][i][j]*(_V[3][i][j]-_V[3][i-1][j])        );
    //total enthalpy
    ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + ((V_L[1]*V_L[1])+(V_L[2]*V_L[2]))*0.5;
    ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + ((V_R[1]*V_R[1])+(V_R[2]*V_R[2]))*0.5;
    
    
    //flux in kasi direction
    R= sqrt(V_R[0]/V_L[0]);
    rho_roe=R*V_L[0];
    u_roe=(R*V_R[1]+V_L[1])/(R+1);
    v_roe=(R*V_R[2]+V_L[2])/(R+1);
    ht_roe= (R*ht_R+ht_L)/(R+1);
    a_roe=sqrt(abs((gamma-1)*(ht_roe-0.5*(u_roe*u_roe+v_roe*v_roe))));
    
    double a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
    nx=n_x(y_n[i][j+1], y_n[i][j], a2);
    ny=n_y(x_n[i][j+1], x_n[i][j], a2);
    T_vel_roe=nx*u_roe+ny*v_roe;
    
    
    r_1[0]=1;
    r_1[1]=u_roe;
    r_1[2]=v_roe;
    r_1[3]=0.5*(u_roe*u_roe+v_roe*v_roe);
    
    
    
    r_2[0]=0;
    r_2[1]=ny*rho_roe;
    r_2[2]=-nx*rho_roe;
    r_2[3]=rho_roe*(ny*u_roe-nx*v_roe);
    
    
    
    
    
    r_3[0]=rho_roe/(2*a_roe);
    r_3[1]=r_3[0]*(u_roe+nx*a_roe);
    r_3[2]=r_3[0]*(v_roe+ny*a_roe);
    r_3[3]=r_3[0]*(ht_roe+T_vel_roe*a_roe);
    
    
    r_4[0]=-rho_roe/(2*a_roe);
    r_4[1]=r_4[0]*(u_roe-nx*a_roe);
    r_4[2]=r_4[0]*(v_roe-ny*a_roe);
    r_4[3]=r_4[0]*(ht_roe-T_vel_roe*a_roe);
    
    lambda[0]=T_vel_roe;
    lambda[1]=T_vel_roe;
    lambda[2]=T_vel_roe+a_roe;
    lambda[3]=T_vel_roe-a_roe;
    
    
    if(abs(lambda[0])<0.1*2*a_roe)lambda[0]=(lambda[0]*lambda[0])/(0.1*4*a_roe)+0.1*a_roe;
    if(abs(lambda[1])<0.1*2*a_roe)lambda[1]=(lambda[1]*lambda[1])/(0.1*4*a_roe)+0.1*a_roe;
    if(abs(lambda[2])<0.1*2*a_roe)lambda[2]=(lambda[2]*lambda[2])/(0.1*4*a_roe)+0.1*a_roe;
    if(abs(lambda[3])<0.1*2*a_roe)lambda[3]=(lambda[3]*lambda[3])/(0.1*4*a_roe)+0.1*a_roe;
    
    delta_rho=V_R[0]-V_L[0];
    delta_u=V_R[1]-V_L[1];
    delta_v=V_R[2]-V_L[2];
    delta_P=V_R[3]-V_L[3];
    
    
    delta_w[0]= delta_rho-delta_P/(a_roe*a_roe);
    delta_w[1]=ny*delta_u-nx*delta_v;
    delta_w[2]=nx*delta_u+ny*delta_v+(delta_P)/(rho_roe*a_roe);
    delta_w[3]=nx*delta_u+ny*delta_v-(delta_P)/(rho_roe*a_roe);
    
    T_vel_L=V_L[1]*nx+V_L[2]*ny;
    T_vel_R=V_R[1]*nx+V_R[2]*ny;
    
    _F_L[0][i][j] =V_L[0]* T_vel_L ;
    _F_L[1][i][j]= V_L[0]*V_L[1]*T_vel_L + V_L[3]*nx;
    _F_L[2][i][j]= V_L[0]*V_L[2]*T_vel_L + V_L[3]*ny;
    _F_L[3][i][j]= V_L[0]*ht_L*T_vel_L;
    
    
    _F_R[0][i][j] =V_R[0]* T_vel_R ;
    _F_R[1][i][j]= V_R[0]*V_R[1]*T_vel_R + V_R[3]*nx;
    _F_R[2][i][j]= V_R[0]*V_R[2]*T_vel_R + V_R[3]*ny;
    _F_R[3][i][j]= V_R[0]*ht_R*T_vel_R;
    
    //flux for kasi faces
    
    
    _F[0][i][j]=0.5*(_F_L[0][i][j]+_F_R[0][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[0] + abs(lambda[1])*delta_w[1]*r_2[0] + abs(lambda[2])*delta_w[2]*r_3[0] + abs(lambda[3])*delta_w[3]*r_4[0]     );
    
    
    
    _F[1][i][j]=0.5*(_F_L[1][i][j]+_F_R[1][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[1] + abs(lambda[1])*delta_w[1]*r_2[1] + abs(lambda[2])*delta_w[2]*r_3[1] + abs(lambda[3])*delta_w[3]*r_4[1]     );
    
    
    
    _F[2][i][j]=0.5*(_F_L[2][i][j]+_F_R[2][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[2] + abs(lambda[1])*delta_w[1]*r_2[2] + abs(lambda[2])*delta_w[2]*r_3[2] + abs(lambda[3])*delta_w[3]*r_4[2]     );
    
    _F[3][i][j]=0.5*(_F_L[3][i][j]+_F_R[3][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[3] + abs(lambda[1])*delta_w[1]*r_2[3] + abs(lambda[2])*delta_w[2]*r_3[3]  + abs(lambda[3])*delta_w[3]*r_4[3]    );
    
    
    
}

void Dv_nozzle::roe_boundary_F_eta (int i,int j)
{
    
    
    
    // flux for eta faces
         
         
         //density
         V_D[0]= _V[0][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[0][i][j-1]*(_V[0][i][j-1]-_V[0][i][j-2]) + (1+kappa)*epsi_minus_y[0][i][j]*(_V[0][i][j]-_V[0][i][j-1])        );
         
         V_U[0]= _V[0][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[0][i][j+1]*(_V[0][i][j+1]-_V[0][i][j]) + (1+kappa)*epsi_plus_y[0][i][j]*(_V[0][i][j]-_V[0][i][j-1])        );
         //x-velocity
         V_D[1]= _V[1][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[1][i][j-1]*(_V[1][i][j-1]-_V[1][i][j-2]) + (1+kappa)*epsi_minus_y[1][i][j]*(_V[1][i][j]-_V[1][i][j-1])        );
         
         V_U[1]= _V[1][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[1][i][j+1]*(_V[1][i][j+1]-_V[1][i][j]) + (1+kappa)*epsi_plus_y[1][i][j]*(_V[1][i][j]-_V[1][i][j-1])        );
         //y-velocity
         V_D[2]= _V[2][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[2][i][j-1]*(_V[2][i][j-1]-_V[2][i][j-2]) + (1+kappa)*epsi_minus_y[2][i][j]*(_V[2][i][j]-_V[2][i][j-1])        );
         
         V_U[2]= _V[2][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[2][i][j+1]*(_V[2][i][j+1]-_V[2][i][j]) + (1+kappa)*epsi_plus_y[2][i][j]*(_V[2][i][j]-_V[2][i][j-1])        );
         //pressure
         V_D[3]= _V[3][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[3][i][j-1]*(_V[3][i][j-1]-_V[3][i][j-2]) + (1+kappa)*epsi_minus_y[3][i][j]*(_V[3][i][j]-_V[3][i][j-1])        );
         
         V_U[3]= _V[3][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[3][i][j+1]*(_V[3][i][j+1]-_V[3][i][j]) + (1+kappa)*epsi_plus_y[3][i][j]*(_V[3][i][j]-_V[3][i][j-1])        );
         //total enthalpy
         ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
         ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
         
         
         
         R= sqrt(V_U[0]/V_D[0]);
         rho_roe=R*V_D[0];
         u_roe=(R*V_U[1]+V_D[1])/(R+1);
         v_roe=(R*V_U[2]+V_D[2])/(R+1);
         ht_roe= (R*ht_U+ht_D)/(R+1);
         a_roe=sqrt(abs((gamma-1)*(ht_roe-0.5*(u_roe*u_roe+v_roe*v_roe))));
         
         double a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
         
         
         nx=n_x(y_n[i+1][j], y_n[i][j], a4);
         ny=n_y(x_n[i+1][j], x_n[i][j], a4);
         
         T_vel_roe=nx*u_roe+ny*v_roe;
         
         
         r_1[0]=1;
         r_1[1]=u_roe;
         r_1[2]=v_roe;
         r_1[3]=0.5*(u_roe*u_roe+v_roe*v_roe);
         
         
         
         r_2[0]=0;
         r_2[1]=ny*rho_roe;
         r_2[2]=-nx*rho_roe;
         r_2[3]=rho_roe*(ny*u_roe-nx*v_roe);
         
         
         
         
         
         r_3[0]=rho_roe/(2*a_roe);
         r_3[1]=r_3[0]*(u_roe+nx*a_roe);
         r_3[2]=r_3[0]*(v_roe+ny*a_roe);
         r_3[3]=r_3[0]*(ht_roe+T_vel_roe*a_roe);
         
         
         r_4[0]=-rho_roe/(2*a_roe);
         r_4[1]=r_4[0]*(u_roe-nx*a_roe);
         r_4[2]=r_4[0]*(v_roe-ny*a_roe);
         r_4[3]=r_4[0]*(ht_roe-T_vel_roe*a_roe);
         
         lambda[0]=T_vel_roe;
         lambda[1]=T_vel_roe;
         lambda[2]=T_vel_roe+a_roe;
         lambda[3]=T_vel_roe-a_roe;
         
         
         if(abs(lambda[0])<0.1*2*a_roe)lambda[0]=(lambda[0]*lambda[0])/(0.1*4*a_roe)+0.1*a_roe;
         if(abs(lambda[1])<0.1*2*a_roe)lambda[1]=(lambda[1]*lambda[1])/(0.1*4*a_roe)+0.1*a_roe;
         if(abs(lambda[2])<0.1*2*a_roe)lambda[2]=(lambda[2]*lambda[2])/(0.1*4*a_roe)+0.1*a_roe;
         if(abs(lambda[3])<0.1*2*a_roe)lambda[3]=(lambda[3]*lambda[3])/(0.1*4*a_roe)+0.1*a_roe;
  
         delta_rho=V_U[0]-V_D[0];
         delta_u=V_U[1]-V_D[1];
         delta_v=V_U[2]-V_D[2];
         delta_P=V_U[3]-V_D[3];
         
         
         delta_w[0]= delta_rho-delta_P/(a_roe*a_roe);
         delta_w[1]=ny*delta_u-nx*delta_v;
         delta_w[2]=nx*delta_u+ny*delta_v+(delta_P)/(rho_roe*a_roe);
         delta_w[3]=nx*delta_u+ny*delta_v-(delta_P)/(rho_roe*a_roe);
         
         T_vel_D=V_D[1]*nx+V_D[2]*ny;
         T_vel_R=V_U[1]*nx+V_U[2]*ny;
         
         _F_D[0][i][j] =V_D[0]* T_vel_D ;
         _F_D[1][i][j]= V_D[0]*V_D[1]*T_vel_D + V_D[3]*nx;
         _F_D[2][i][j]= V_D[0]*V_D[2]*T_vel_D + V_D[3]*ny;
         _F_D[3][i][j]= V_D[0]*ht_D*T_vel_D;
         
         
         _F_U[0][i][j] =V_U[0]*T_vel_U ;
         _F_U[1][i][j]= V_U[0]*V_U[1]*T_vel_U + V_U[3]*nx;
         _F_U[2][i][j]= V_U[0]*V_U[2]*T_vel_U + V_U[3]*ny;
         _F_U[3][i][j]= V_U[0]*ht_U*T_vel_U;
         
     
         
         _F_eta[0][i][j]=0.5*(_F_D[0][i][j]+_F_U[0][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[0] + abs(lambda[1])*delta_w[1]*r_2[0] + abs(lambda[2])*delta_w[2]*r_3[0] + abs(lambda[3])*delta_w[3]*r_4[0]     );
         
         
         
         _F_eta[1][i][j]=0.5*(_F_D[1][i][j]+_F_U[1][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[1] + abs(lambda[1])*delta_w[1]*r_2[1] + abs(lambda[2])*delta_w[2]*r_3[1] + abs(lambda[3])*delta_w[3]*r_4[1]     );
         
         
         
         _F_eta[2][i][j]=0.5*(_F_D[2][i][j]+_F_U[2][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[2] + abs(lambda[1])*delta_w[1]*r_2[2] + abs(lambda[2])*delta_w[2]*r_3[2] + abs(lambda[3])*delta_w[3]*r_4[2]     );
         
         _F_eta[3][i][j]=0.5*(_F_D[3][i][j]+_F_U[3][i][j]) -0.5*(abs(lambda[0])*delta_w[0]*r_1[3] + abs(lambda[1])*delta_w[1]*r_2[3] + abs(lambda[2])*delta_w[2]*r_3[3]  + abs(lambda[3])*delta_w[3]*r_4[3]    );
         
         
         
         
    
    
}
