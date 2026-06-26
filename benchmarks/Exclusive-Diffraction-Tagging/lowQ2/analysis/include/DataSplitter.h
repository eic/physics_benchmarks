
#pragma once


/*!
 *  Class to split data in multiple dimensions by
 *  defined bins in DataFrame columns
 */

#include <TMath.h>
#include <vector>
#include <deque>

using AxisData = std::vector<double>; //should contain all edges including upper, N+1 entries where N is number of bins

namespace rad{
  
  namespace histo{

    /*
     * Helper function to create AxisData (aka bin edges)
     * For case of equal bin widths in a given range
     * Given number of bins (N), and range , (lower, upper)
     */
    AxisData RegularSplits(short N,double lower, double upper){
      
      AxisData edges(N+1);
      double width = (upper-lower)/N;
      
      edges[0]=lower;
      for(short i=1;i<=N;++i){
	edges[i]=edges[i-1] + width;
      }

      return edges;
    }
  
  
    class DataSplitter {

      using SplitPtrs = std::vector<std::shared_ptr<DataSplitter> >;

    public:
      DataSplitter(){
	/////////
      }

      
    DataSplitter(std::string name,double low=0, double width=0):_nameDimension{name},_low{low},_width{width},_centre{low+width/2}{
	/////////
      }

      /*
       *  Case all elements of the current split level
       *  have same binning scheme
       */
      void AddRegularDimension(std::string name,const AxisData& bins){
     	_dimensionNames.push_back(name);
	//	std::cout<<"AddRegularDimension "<<name<<" "<<_NDimension<<" "<<_N<<" "<<_dimensionNames.size()<<" elements "<<_elements.size()<<std::endl;
	if(_NDimension==1){ //at deepest dimension
	  std::vector<AxisData> repeat_bins(_N,bins);
	  AddDimension(name,repeat_bins);
	}
	else{
	  for(ushort i = 0 ; i <_N; ++i){
	    _elements[i]->AddRegularDimension(name,bins);
	  }
	}
     }
      
      /*
       *  Case all elements of the current split level
       *  may not have same binning scheme
       *  Define splits for each element in all_bins vector
       */
      void AddDimension(std::string name,const std::vector<AxisData>& all_bins){
	//	std::cout<<"AddDimension "<<name<<" "<<_NDimension<<" "<<_N<<" "<<_dimensionNames.size()<<std::endl;
 	if(all_bins.size()!=_N){
	  std::cerr<< "DataSplitter::AddDimension need an axis for each element, N elements = "<<_N<<" N given = "<< all_bins.size() <<"...exiting..."<<std::endl;
	  exit(0);
	}
	
	if(_elements.empty()){//first dimension
	  AddElements(name,all_bins[0]);
	  return;
	}
	
	for(ushort i = 0 ; i <_N; ++i){
	  _elements[i]->AddElements(name,all_bins[i]);
	  //	_iN[i] = _elements[i]->N();//number of elements in each element!
	}
  
      }
      /*
       * Add bins/splits to a particular element
       */
      void AddElements(std::string name,const AxisData& edges){
	//	std::cout<<"AddElements " <<name<<" "<<_dimensionNames.size()<<" "<<_NDimension<<"  N "<<_N<<std::endl;
	_elements.resize(edges.size() - 1);
	_edges=edges; //keep last high bin
      
	_N = _elements.size();
      	_iN.resize(_N);

	for(short i=0;i<_N;++i){
	  _elements[i]=std::move(std::make_shared<DataSplitter>(name,edges[i],(edges[i+1]-edges[i])));
	}
	_NDimension++;

      }

      /*
       * Create partial bins sums
       * and unique name for each bin in all dimensions
       */
      void Init(){

	SumElements();
	_binNames = ConstructBinNames();
      }

      /*
       * create partial sums of bins
       * so can work out unique bin ID
       */
      uint SumElements(){
	if(_elements.empty()) return 1;
	if(_NDimension==1) return _N;
 
	auto sum = 0;
	for(short i=0; i<_N ; ++i){	  
	  _iN[i]=_elements[i]->SumElements();
	  sum+=_iN[i];
	}
	return sum;
      }
      /*
       * Get deepest level bin name
       */
      std::string ElementName(ushort ib) const {
	return _elements[ib]->DimensionName()+Form("_%1.4f",_elements[ib]->Centre());
      }

      /*
       *  Find bin/split for a given event
       *  xs is the values of each dimension
       *  in order added
       */
      short GetBin(std::deque<double> xs) const {
       	//std::cout<<" GetBin "<<xs.size()<<" "<<_nameDimension<<" "<<_NDimension<<std::endl;
	if(xs.empty()==true) return 1;
	//copy xs as going to strip it
	auto x = xs[0];
	xs.pop_front();
	auto bin0 = GetBin(x); //bin in this axis
	if(bin0==-1) return bin0;
      
	if(xs.empty()==true) {
	  //furthest dimension, can return
	  return bin0;
	}
      
	//now get number bin in next dimension
	auto bin1 = _elements[bin0]->GetBin(xs);
	if(bin1==-1) return bin1;

	//get Number of bins to here
	//summing over elements of previous bins
	ushort sumBins = 0;
	for(auto i =0;i<bin0;++i){
	  sumBins+=_iN[i];
	}
	return sumBins+bin1;
      }
      /*
       * Get bin in particular dimension
       * with current variable value x
       */
      short GetBin(double x) const {
	//std::cout<<NTotal()<<" dim "<<" "<<_NDimension<<" "<<_nameDimension<<" x "<<x<<" from "<<_N<<_edges.front()<<" "<<_edges.back()<<endl;
	//first bin returns 0, second 1,...last _N-1
	if(NTotal()==1) return 0;
	if(x<_edges.front()) return -1;
	if(x>_edges.back()) return -1;
	return TMath::BinarySearch(_N,_edges.data(),x);
      }
    
      void Print(std::string opt="",std::deque<double> xs={}){
	if(_N==1) return;
	opt+="\t";
	std::cout<< "\n"<<opt<<"Print dimension "<<_NDimension<<", "<< _N <<" bin edges :"<<std::endl;
	std::cout<<opt;
	for(const auto& el:_elements){
	  xs.push_back(el->Centre());
	  std::cout<<" "<<el->LowEdge();
	  el->Print(opt,xs);
	  std::cout<<opt;
	  xs.pop_back();
	}
	std::cout<<endl;
      }
      /*
       *  Create a unique bin name containing all dimensions
       *  and limits of each bin
       */
      std::vector<std::string> ConstructBinNames(){
	std::vector<std::string> bin_names;
	auto temp = _NTotal;
	_NTotal=1000;//make sure we get a bin, NTotal is set later in this function
	for(const auto& el:_elements){
	  auto base=ElementName(GetBin(el->Centre()));
	  auto higher_parts = el->ConstructBinNames();
	  if( higher_parts.empty()) bin_names.push_back(base);
	  for(const auto& name:higher_parts){
	    bin_names.push_back(base+name);
	  }
	}

	//case there are no elements, => just single bin
	_NTotal=bin_names.size();
	if(_NTotal==0){
	  bin_names.push_back("");
	  _NTotal=1;
	}
	return bin_names;
      }

      
      double LowEdge() const {return _low;}
      double Centre() const {return _centre;}
      double Width() const {return _width;}
      ushort N()const {return _N;}
      uint NTotal()const {return _NTotal;}
    
      const string& DimensionName(){return  _nameDimension;}
      const std::vector<std::string>&  GetDimensionNames(){return _dimensionNames;}
      const std::vector<std::string>&  GetBinNames(){return _binNames;}
  
    private:

      double _low=0.;
      double _centre=0.;
      double _width=0.;

      //store element data for fast search
      AxisData _edges;
 
      //store element data for configuring
      SplitPtrs _elements;//multi-dimensions

      std::vector<std::string> _dimensionNames;
      std::vector<std::string> _binNames;
      std::string _nameDimension;
      std::vector<ushort> _iN;
      uint _NTotal=1;
      ushort _N=1;
      ushort _NDimension=1;
    
    };

  }
}
