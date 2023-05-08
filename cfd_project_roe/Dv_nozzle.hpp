//
//  Dv_nozzle.hpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#ifndef Dv_nozzle_hpp
#define Dv_nozzle_hpp

#include <stdio.h>
#include <fstream>


class Dv_nozzle {
    
    
private:
    
    static int const _imax_f=21, _imax_c=20, _imin=0;
    int im_f, im_c;
    double A, B, C, D, E, F;  // just variables for calculations
    double k_2=0.2, k_4=0.02;
    double x_f[_imax_f];
    double y_f[_imax_f];
    double r_1[4],r_2[4],r_3[4],r_4[4];
    double lambda[4];
    double eigen_max[2];
    double R;
    double rho_roe, u_roe,v_roe, ht_roe, a_roe;
    double T_vel_L,T_vel_R,T_vel_D,T_vel_U,T_vel_roe;
    double nx,ny;
    
    double delta_u,delta_v, delta_P, delta_rho;
    //variables for grid nodes and geometry terms
    static int const  i_n_max=21, j_n_max=21, k_n_max=2;
    double x_n[i_n_max][j_n_max];
    double y_n[i_n_max][j_n_max];
    
    
    
    
    
    
    
    double delta_w[4];
    double _F_L[4][_imax_f][_imax_f] , _F_R[4][_imax_f][_imax_f];
    double _F_D[4][_imax_f][_imax_f] , _F_U[4][_imax_f][_imax_f];
    double avg_area;
    double avg_d_area;
    double x_c[_imax_c];
    double y_c[_imax_c];
    double delta_t[_imax_c][_imax_c];
    double delta_t_gl=999;
    double CFL=0.01;
    double area_f[_imax_f][_imax_f];
    double _area, _d_area;
    double mu[_imax_c];
    double epsilon_2;
    double _R[4]={0,0,0,0};
    double _R_iter[4]={0,0,0,0};
    double eigen_v_avg;
    double _delta_x, _delta_y;
    double _P_0,  _T_0, mach_in, mach_out, T_in, u_in, P_in, rho_in, et_in, T_out, u_out, P_out, rho_out, et_out, v_in,v_out;
    double _M_0;
    double _mach[_imax_c][_imax_c];
    double T[_imax_c][_imax_c];
    double et[_imax_c][_imax_c];
    double ht[_imax_c][_imax_c];
    double _V[4][_imax_c][_imax_c];
    double _U[4][_imax_c][_imax_c];
    double _F[4][_imax_f][_imax_f];
    double _F_eta[4][_imax_f][_imax_f];
    double _U_ghost_outflow[4][_imax_c];
    double _V_ghost_outflow[6][_imax_c];
    double _V_ghost_outflow_2[6][_imax_c];
    double _V_ghost_outflow_3[6][_imax_c];
    double _U_ghost_inflow[4][_imax_c];
    double _V_ghost_inflow[6][_imax_c];
    double _V_ghost_inflow_2[6][_imax_c];
    double _V_ghost_inflow_3[6][_imax_c];
    //upwind schemes variables
    
    double kappa=-1;
    double epsilon_upwind=0;
    double _F_C[4][_imax_f][_imax_f];
    double _F_P[4][_imax_f][_imax_f];
    
    double c_plus;
    double c_minus;
    double alpha_plus;
    double alpha_minus;
    double beta_R;
    double beta_L;
    double a[_imax_c][_imax_c];
    double a_L;
    double a_R;
    
    double D_plus;
    double D_minus;
    double P_plus;
    double P_minus;
    double V_L[4];
    double V_R[4];
    double V_D[4];
    double V_U[4];
    double T_L;
    double T_R;
    double T_D;
    double T_U;
    double r_plus;
    double r_minus;
    double epsi_plus_x[4][_imax_f+2][_imax_f+2];
    double epsi_minus_x[4][_imax_f+2][_imax_f+2];
    double epsi_plus_y[4][_imax_f+2][_imax_f+2];
    double epsi_minus_y[4][_imax_f+2][_imax_f+2];
    double epsi_minus_x_ghost[4][_imax_f+2];
    double epsi_plus_x_ghost[4][_imax_f+2];
    double epsi_minus_y_ghost[4][_imax_f+2];
    double epsi_plus_y_ghost[4][_imax_f+2];
    double ht_L,ht_R,ht_D,ht_U;
    
    
    double mach_knight;
    double M_plus;
    double M_minus;
    double M_L;
    double M_R;
    
    
    
    double _D[3][_imax_f][_imax_f]; //JST dissipation scheme
    double gamma = 1.4, R_u=  8314, M_air= 28.96; //gamma air, universal gas const, molecular weight of air
    double R_air = R_u/M_air;
    double epsi_1;
    double rL2norm[4];
    double _rL2initial[4]={0,0,0,0};
    std::fstream L_vs_Iter; //creation of file to write the results
    
    std::fstream rho_vs_x; //creation of file to write the results
    std::fstream u_vs_x; //creation of file to write the results
    
    std::fstream p_vs_x; //creation of file to write the results
    
public:
    //constructor to open the file for write
    Dv_nozzle(){
        
     
        
        rho_vs_x.open("/Users/cringedaddy/CFD class/project_roe/hw_3_roe/hw_3_roe/cfd_project_roe/cfd_project_roe/results/rho.csv", std::ios::trunc | std::ios::out);
        
        u_vs_x.open("/Users/cringedaddy/CFD class/project_roe/hw_3_roe/hw_3_roe/cfd_project_roe/cfd_project_roe/results/u.csv", std::ios::trunc | std::ios::out);
        p_vs_x.open("/Users/cringedaddy/CFD class/project_roe/hw_3_roe/hw_3_roe/cfd_project_roe/cfd_project_roe/results/p.csv", std::ios::trunc | std::ios::out);
        L_vs_Iter.open("/Users/cringedaddy/CFD class/project_roe/hw_3_roe/hw_3_roe/cfd_project_roe/cfd_project_roe/results/L2.csv", std::ios::trunc | std::ios::out);
        L_vs_Iter<<"#iter#"<<"##L2_1#"<<"##L2_2#"<<"##L2_3#"<<std::endl;
        
    };
    // destructor to close the file
    ~Dv_nozzle()
    {
        
        rho_vs_x.close();
        u_vs_x.close();
        p_vs_x.close();
        L_vs_Iter.close();
        
    };

   
    void set_geometry(int,double,double);
    void initialization(double,double);
    void set_boundary_cond();
    void euler_explicit();
    double area(double);
    double d_area(double);
    void time_step();
    void print_res();
    void rL2initial();
    double L2norm(int);
    void roe_flux();
    void roe_boundary_in_out(int,int);
    void roe_boundary_in_out_F_kasi(int ,int);
    void roe_boundary_walls(int);
    void mesh_nodes();
    // to calculate face area, cell volume (for the eigen values) and outpointing normals
    double face_area(double,double,double,double);
    double n_x(double, double,double);
    double n_y(double,double,double);
    double cell_vol(double,double);
    
    
    
};




#endif /* Dv_nozzle_hpp */
