#include "hexplit.h"
TLorentzVector ZDC_neutron_recon(TTreeReaderArray<float>&E,TTreeReaderArray<float>&t,TTreeReaderArray<float>&x,TTreeReaderArray<float>&y,
				 TTreeReaderArray<float>&z){ 
  //cout << "hit count" << E.GetSize() << endl;
  vector<double> subcellE, subcellx, subcelly, subcellz;
  hexplit_opts opts1;
  opts1.MIP=0.000470;
  opts1.side_length=31.3;
  opts1.layer_spacing=25.1;
  opts1.Emin=opts1.MIP*0.1;
  opts1.tmax=320;
  
  subcell_reweight(E, t, x, y, z, subcellE, subcellx, subcelly, subcellz, opts1);

  //now determine the position from the subcells
  position_recon_opts opts2;
  opts2.sampling_fraction=0.0203;


  /*vector<double> flatE;flatE.reserve(2048*10);
  vector<double> flatx;flatx.reserve(2048*10);
  vector<double> flaty;flaty.reserve(2048*10);
  vector<double> flatz;flatz.reserve(2048*10);
  
  flatten(subcellE, subcellx, subcelly, subcellz, flatE, flatx,	flaty, flatz);
  subcellE=flatE;
  subcellx=flatx;
  subcelly=flaty;
  subcellz=flatz;*/

  
  // parameterization of values to determine which
  // value of w0 to use
  opts2.w0_a=5.0;
  opts2.w0_b=0.65;
  opts2.w0_c=0.31;
  opts2.E0=50;

  /*opts2.w0_a=3.0;
  opts2.w0_b=0;
  opts2.w0_c=0;*/
  
  //obtains the local position
  double E_recon, x_recon, y_recon, z_recon;  
  
  recon_position(subcellE, subcellx, subcelly, subcellz, E_recon, x_recon,
		 y_recon, z_recon, opts2);

  //translate and then rotate
  double z_offset=36601;  //z' position of the center of the ZDC SiPM-on-tile detector
  double cross=-0.025;
  z_recon=z_recon+z_offset;
  double z_recon_new=z_recon*cos(cross)-x_recon*sin(cross);
  double x_recon_new=x_recon*cos(cross)+z_recon*sin(cross);

  x_recon=x_recon_new;
  z_recon=z_recon_new;
  double r_recon=sqrt(x_recon*x_recon+y_recon*y_recon+z_recon*z_recon);
  double m_n=0.939565;
  double p_recon=sqrt(E_recon*E_recon-m_n*m_n);
  //cout << "x, y, z recon "  <<x_recon << ", " << y_recon << ", " << z_recon << endl;
  //cout << "p_recon " << p_recon << endl;
  TLorentzVector a(p_recon*x_recon/r_recon, p_recon*y_recon/r_recon, p_recon*z_recon/r_recon, E_recon);
  return a;
}
