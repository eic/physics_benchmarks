#include <iostream>
//#include <math>
#define SUBCELLS 12
struct hexplit_opts {
  // longitudinal space between layers
  double layer_spacing;
  // side length of a cell
  double side_length;
  // MIP value for the detector
  double MIP;
  // minimum energy cut
  double Emin=0;
  // maximum time cut
  double tmax=1000;
};
struct position_recon_opts{
  //use the sampling fraction to get the recon energy,
  // which is used in the parameterization for the w0.
  double sampling_fraction;
  
  // parameterization of values to determine which
  // value of w0 to use
  double w0_a;
  double w0_b;
  double w0_c;
  double E0=50;
  
};
const double pi=M_PI;
/*
phi=np.linspace(0, np.pi*5/3, 6)
        cph=np.cos(phi)
        sph=np.sin(phi)
    
        offsetsx=np.concatenate((1.5*cph, -sqrt3/2*sph))*sl
        offsetsy=np.concatenate((1.5*sph, sqrt3/2*cph))*sl
    
        #determine the positions of the vertices
        subcell_offsetsx=[]
        subcell_offsetsy=[]
        subcell_offsetsx+=[[sl/2*cph[k], sl/2*(cph[k]+cph[(k+1)%6]), 
                            sl*cph[k], sl/2*(cph[k]+cph[(k+5)%6])] for k in range(6)]
        subcell_offsetsx+=[[0, sl/4*cph[k]-sqrt3*sl*sph[k]/4, 
                            -sqrt3*sl*sph[k]/2, -sl/4*cph[k]-sqrt3*sl*sph[k]/4] for k in range(6)]
        
        subcell_offsetsy+=[[sl/2*sph[k], sl/2*(sph[k]+sph[(k+1)%6]), 
                            sl*sph[k], sl/2*(sph[k]+sph[(k+5)%6])] for k in range(6)]
        subcell_offsetsy+=[[0, sqrt3*sl*cph[k]/4+sl/4*sph[k], 
                            sqrt3*sl*cph[k]/2, sqrt3*sl*cph[k]/4-sl/4*sph[k]] for k in range(6)]*/

//positions where the overlapping cells are relative to a given cell (in units of hexagon side length)
double neighbor_offsets_x[]={1.5*cos(0), 1.5*cos(pi/3), 1.5*cos(2*pi/3),1.5*cos(3*pi/3), 1.5*cos(4*pi/3), 1.5*cos(5*pi/3),
			     -sqrt(3)/2.*sin(0),-sqrt(3)/2.*sin(pi/3),-sqrt(3)/2.*sin(2*pi/3),-sqrt(3)/2.*sin(3*pi/3),-sqrt(3)/2.*sin(4*pi/3),-sqrt(3)/2.*sin(5*pi/3)};
double neighbor_offsets_y[]={1.5*sin(0), 1.5*sin(pi/3), 1.5*sin(2*pi/3),1.5*sin(3*pi/3), 1.5*sin(4*pi/3), 1.5*sin(5*pi/3),
                              sqrt(3)/2.*cos(0), sqrt(3)/2.*cos(pi/3), sqrt(3)/2.*cos(2*pi/3), sqrt(3)/2.*cos(3*pi/3), sqrt(3)/2.*cos(4*pi/3), sqrt(3)/2.*cos(5*pi/3)};

//indices of the neighboring cells which overlap to produce a given subcell
int neighbor_indices[12][3]={{0, 11,10}, {1, 6, 11},{2, 7, 6}, {3,8,7}, {4,9,8}, {5,10,9},
			     {6, 11, 7}, {7, 6, 8}, {8, 7, 9}, {9,8,10},{10,9,11},{11,10,6}};

//positions of the centers of subcells
double subcell_offsets_x[]={0.75*cos(0), 0.75*cos(pi/3), 0.75*cos(2*pi/3), 0.75*cos(3*pi/3), 0.75*cos(4*pi/3), 0.75*cos(5*pi/3),
			    -sqrt(3)/4*cos(0),-sqrt(3)/4*cos(pi/3),-sqrt(3)/4*cos(2*pi/3),-sqrt(3)/4*cos(3*pi/3),-sqrt(3)/4*cos(4*pi/3),-sqrt(3)/4*cos(5*pi/3)};
double subcell_offsets_y[]={0.75*sin(0), 0.75*sin(pi/3), 0.75*sin(2*pi/3), 0.75*sin(3*pi/3), 0.75*sin(4*pi/3), 0.75*sin(5*pi/3),
                             sqrt(3)/4*sin(0), sqrt(3)/4*sin(pi/3), sqrt(3)/4*sin(2*pi/3), sqrt(3)/4*sin(3*pi/3), sqrt(3)/4*sin(4*pi/3), sqrt(3)/4*sin(5*pi/3)};
			    
void subcell_reweight(TTreeReaderArray<float>&E,TTreeReaderArray<float>&t,TTreeReaderArray<float>&x,TTreeReaderArray<float>&y,TTreeReaderArray<float>&z, //number of hits, energies, time, and local x,y,z positions of the hits
		      vector<double> & subcellE, vector<double>& subcellx, vector<double> & subcelly, vector<double> & subcellz,  //returned position of the cluster
		  hexplit_opts opts) {
  int nhits=E.GetSize();
  double sl=opts.side_length;
  double layer_spacing=opts.layer_spacing;
  double Emin=opts.Emin;
  double tmax=opts.tmax;
  double MIP=opts.MIP;

  double Esum=0;
  for(int i=0; i<nhits; i++){
    //skip hits that do not pass E and t cuts
    if (E[i]<Emin || t[i]>tmax)
      continue;
    
    //keep track of the energy in each neighboring cell
    double Eneighbors[SUBCELLS];
    for (int j=0; j<SUBCELLS; j++)
      Eneighbors[j]=0;
      
    for (int j=0; j<nhits; j++){
      //only look at hits nearby within two layers of the current layer
      if (abs(z[i]-z[j])>2.5*layer_spacing || z[i]==z[j])
	continue;
      //cout << "z cut" << endl;
      if (E[j]<Emin || t[j]>tmax)
	continue;
      //difference in transverse position (in units of side lengths)
      double dx=(x[j]-x[i])/sl;
      double dy=(y[j]-y[i])/sl;
      if (abs(dx)>2 || abs(dy)>sqrt(3))
	continue;
      //cout << "dx, dy cut passed: " << dx << " " << dy <<  endl;
      
      //loop over locations of the neighboring cells
      //and check if the jth hit matches this location
      double tol=0.01; //tolerance for rounding errors
      for(int k=0;k<SUBCELLS;k++){
	//cout << "neighbor pos"<< neighbor_offsets_x[k] << " " << neighbor_offsets_y[k] << endl;
	if(abs(dx-neighbor_offsets_x[k])<tol && abs(dy-neighbor_offsets_y[k])<tol){
	  //cout << "found neighbor " << k;
	  Eneighbors[k]+=E[j];
	  break;
	}
      }
    }
    double weights[SUBCELLS];
    for(int k=0; k<SUBCELLS; k++){
      Eneighbors[k]=max(Eneighbors[k],MIP);
    }
    double sum_weights=0;
    for(int k=0; k<SUBCELLS; k++){
      weights[k]=Eneighbors[neighbor_indices[k][0]]*Eneighbors[neighbor_indices[k][1]]*Eneighbors[neighbor_indices[k][2]];
      sum_weights+=weights[k];
    }
    for(int k=0; k<SUBCELLS;k++){
      subcellE.push_back(E[i]*weights[k]/sum_weights);
      subcellx.push_back(x[i]+subcell_offsets_x[k]*sl);
      subcelly.push_back(y[i]+subcell_offsets_y[k]*sl);
      subcellz.push_back(z[i]);
    }
  }
}
// Optional subroutine to merge all hits or subcells at a given transverse position together after doing
// subcell reweighting.  This method was found to not improve the angle resolution, so it is not
// used.  
void flatten(vector<double>&E, vector<double>& x, vector<double>& y, vector<double>& z,
	     vector<double>&Enew,vector<double>& xnew, vector<double>& ynew, vector<double>& znew){
  Enew.clear();
  xnew.clear();
  ynew.clear();
  znew.clear();
  for(int i = 0; i<E.size(); i++){
    int found=0;
    for(int j=0; j<Enew.size(); j++){
      double tol=0.01;
      if(abs(xnew[j]-x[i])<tol && abs(ynew[j]-y[i])<tol){
	Enew[j]+=E[i];
	znew[i]+=z[i]*E[i];
	found=1;
	break;
      }
    }
    if (!found){
      Enew.push_back(E[i]);
      xnew.push_back(x[i]);
      ynew.push_back(y[i]);
      znew.push_back(z[i]*E[i]);
    }
  }
  for(int j = 0; j<Enew.size();j++){
    znew[j]/=Enew[j];
  }
}


void recon_position(vector<double>& E, vector<double>& x, vector<double>& y, vector<double>& z, double & E_recon, double & x_recon,
		    double & y_recon, double & z_recon, position_recon_opts& opts){

  //first determine the total energy of the particle and of the shower
  E_recon=0;
  double E_shower=0;
  for (int i = 0; i<E.size();i++){
    E_recon += E[i]/opts.sampling_fraction;
    E_shower += E[i];
  }
  //now determine the value of w0 to use
  double w0=opts.w0_a+opts.w0_b*log(E_recon/opts.E0)+opts.w0_c*pow(log(E_recon/opts.E0),2);  
  //cout <<"w0" << w0 << endl;
  double sum_weights=0;
  //position of the cluster
  x_recon=0;y_recon=0;z_recon=0;
  for (int i=0; i<E.size(); i++){    
    double weight=max(log(E[i]/E_shower)+w0, 0.);
    sum_weights+=weight;
    x_recon+=x[i]*weight;
    y_recon+=y[i]*weight;
    z_recon+=z[i]*weight;
  }
  x_recon/=sum_weights;
  y_recon/=sum_weights;
  z_recon/=sum_weights;  
}
