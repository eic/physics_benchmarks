#pragma once
#include "Indicing.h"
#include "Constants.h"

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

void BasicKinematics(){}

//namespace rad{
//RDF interpreter can only take 1 level of namespace...
namespace rad{
    ///\brief Helper functions and functors for RDF processing
    ///       combining momentum components into 4-vectors
    using ROOT::Math::PxPyPzMVector ;
    using ROOT::Math::XYZVector ;
    using ROOT::RVecF;
    using ROOT::RVecI;
    using ROOT::RVec;
    using ROOT::Math::VectorUtil::boost;

  ///\brief functor returning 4-vector of fixed components
  // class FixedP4 {
    
  // public:
  //   // FixedP4()=default;
  //   FixedP4(double x, double y, double z, double m):
  //     _p4{x,y,z,m}{};
    
  //   const PxPyPzMVector& operator()() const{
  //     return _p4;
  //   }
    
  // private:
  //   const PxPyPzMVector _p4;
    
  // };
  
 
   ///\brief return 4-vector of particle idx
    template<typename Tp, typename Tm>
    PxPyPzMVector FourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      return PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]);
    }
  
    ///\brief add 4-vectors of particles ip to p4
    template<typename Tp, typename Tm>
    void SumFourVector(PxPyPzMVector& p4, const RVecI &ip,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
    
      auto np = ip.size();
      for (size_t i = 0; i < np; ++i) {
	p4 += PxPyPzMVector(px[ip[i]], py[ip[i]], pz[ip[i]], m[ip[i]]);
      }
    }
  
    ///\brief subtract 4-vectors of particles ip from p4
    template<typename Tp, typename Tm>
    void SubtractFourVector(PxPyPzMVector& p4, const RVecI &ip,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto np = ip.size();
      for (size_t i = 0; i <np ; ++i) {
	p4 -= PxPyPzMVector(px[ip[i]], py[ip[i]], pz[ip[i]], m[ip[i]]);
      }
    }
  
    ///\brief return 4-vector of summed particles ipart
    template<typename Tp, typename Tm>
    PxPyPzMVector FourVector(const RVecI &ipart,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
  
      PxPyPzMVector psum(0,0,0,0);
      SumFourVector(psum,ipart,px,py,pz,m);
      return psum;
    
    }
    ///\brief return mass of combined 4-vectors, adding particles with indices ipos and subtracting ineg
    template<typename Tp, typename Tm>
    double FourVectorMassCalc(const RVecI &ipos, const RVecI &ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      /*
	//Optional we could add in check like this for indices
      if(indice::InvalidIndices(ipos)) return constant::InvalidEntry();
      if(indice::InvalidIndices(ineg)) return constant::InvalidEntry();
      */
      
      PxPyPzMVector psum(0,0,0,0);
      SumFourVector(psum,ipos,px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M();
    }
    ///\brief return magnitude of momentum
    template<typename T>
      RVec<double> ThreeVectorMag(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      return sqrt(x * x + y * y + z * z);
    }
    ///\brief return eta of momentum
    template<typename T>
      RVec<double> ThreeVectorTheta(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      // RVec<T> ThreeVectorTheta(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      auto mag = ThreeVectorMag(x,y,z);
      auto costh = z/mag;
      return acos(costh);
    }
    ///\brief return eta of momentum
    template<typename T>
      RVec<double> ThreeVectorPhi(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      //RVec<T> ThreeVectorPhi(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      //std::cout<<" ThreeVectorPhi "<<x.size()<<std::endl;
      return atan2(y,x); //will use vectorised version
    }
   ///\brief return eta of momentum
    template<typename T>
      RVec<double> ThreeVectorEta(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      //RVec<T> ThreeVectorEta(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      auto theta = ThreeVectorTheta(x,y,z);
      return -log(tan(0.5 * theta));//will use vectorised version
    }
    ///\brief return x-component
    template<typename T>
      RVec<double> ThreeVectorX(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      //RVec<T> ThreeVectorX(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      return p*sin(theta)*cos(phi);
    }
    ///\brief return y-component
    template<typename T>
      RVec<double> ThreeVectorY(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      //RVec<T> ThreeVectorY(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      return p*sin(theta)*sin(phi);
    }
    ///\brief return z-component
    template<typename T>
      RVec<double> ThreeVectorZ(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      //RVec<T> ThreeVectorZ(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      return p*cos(theta);
    }

  ///\brief print all particles for event
  template<typename Tpid, typename Tp, typename Tm>
  bool PrintParticles(const std::string& type, ULong64_t entry,const RVec<Tpid> &pid,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
    
    std::cout<< type <<"  PrintParticles Event = "<< entry <<std::endl;
    for(size_t  idx=0; idx<px.size();++idx){
      std::cout<< " "<<pid[idx]<<"\t"<<PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx])<<" pmag "<< ThreeVectorMag(px, py, pz)[idx] << " theta "<<ThreeVectorTheta(px, py, pz)[idx]<<"\n";
    }
    return true;
  }
  
    /* NOTE : we might want to change to using edm4hep functions. Then VecMag would change to 
       template <typename T>
       auto VecMag = [](ROOT::VecOps::RVec<T> momenta) {
       return ROOT::VecOps::Map(momenta, [](const T& p) { return edm4hep::utils::magnitude(p.momentum); });
       };
       //and we would call like Define("pmag", VecMag<edm4hep::MCParticleData> ,{"MCParticles"})
       */


  }//compute
//}//rad
