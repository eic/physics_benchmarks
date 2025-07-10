#pragma once

#include "RVecHelpers.h"
#include "Constants.h"

namespace rad{
  namespace indice{

    /**
     * returns the indice of the (n_occurance)th time value is in values
     */ 
    auto useNthOccurance(const int n_occurance,const int value) {
      return  [n_occurance,value]
	(const ROOT::RVecI& values){
 	return helpers::findNthIndex(values,n_occurance,value);
      };
     }

    /**
     * return nth occurance of type value ordered by values of vector sorter
     * unforunately need to give explicit type for sorter as compiler cannot infer it
     * instead just have I and F versions....
     */
    
   template <typename T>
   auto  useNthOccuranceSortedBy(const int n_occurance,const int value) {
      return [n_occurance,value]
	(const ROOT::RVecI& values,const ROOT::RVec<T>& sorter){
	//get sorted indices of sorter
	auto sortIndices = ROOT::VecOps::StableArgsort(sorter);
	//reorder values
	auto reorderValues = ROOT::VecOps::Take(values, sortIndices);
	//find position of nth occurance of value
	return helpers::findNthIndex(reorderValues,n_occurance,value);
      };
    }
    auto  useNthOccuranceSortedByF(const int n_occurance,const int value) {
      return [n_occurance,value]
	(const ROOT::RVecI& values,const ROOT::RVecF& sorter){
	//get sorted indices of sorter
	auto sortIndices = ROOT::VecOps::StableArgsort(sorter);
	//reorder values
	auto reorderValues = ROOT::VecOps::Take(values, sortIndices);
	//find position of nth occurance of value
	return helpers::findNthIndex(reorderValues,n_occurance,value);
      };
    }
    auto  useNthOccuranceSortedByI(const int n_occurance,const int value) {
      return [n_occurance,value]
	(const ROOT::RVecI& values,const ROOT::RVecI& sorter){
	//get sorted indices of sorter
	auto sortIndices = ROOT::VecOps::StableArgsort(sorter);
	//reorder values
	auto reorderValues = ROOT::VecOps::Take(values, sortIndices);
	//find position of nth occurance of value
	return helpers::findNthIndex(reorderValues,n_occurance,value);
      };
    }
    /**
     * return nth occurance of type value ordered by values of vector sorter
     */
    // auto  useNthOccuranceWithCondition(const int n_occurance,const int value) {
    //   return [n_occurance,value]
    // 	(const ROOT::RVecI& values,const ROOT::RVecI& sorter){
    // 	//get sorted indices of sorter
    // 	auto sortIndices = ROOT::VecOps::StableArgsort(sorter);
    // 	//reorder values
    // 	auto reorderValues = ROOT::VecOps::Take(values, sortIndices);
    // 	//find position of nth occurance of value
    // 	return helpers::findNthIndex(reorderValues,n_occurance,value);
    //   };
    // }

    /**
     *simple general function to return zeroth index of branch
     *subtracts 2 to accounts for beam indices being removed
     *from list
     *offset exists incase submodule (i.e. epic-rad) remove/add
     *particles in the list ordering.
     */
    auto UseAsID(int entry, int offset=0) {
      return [entry, offset](const ROOT::RVec<int>& id) -> int {
    	return id[entry]-offset;
      };
    }
    /**
     * Function to check if all indices in a vec are valid
     * i.e. != -1
     */
    template<typename T>
    bool InvalidIndices(const ROOT::RVec<T>& indices){
      return ROOT::VecOps::Any( ( indices==-1 ) );
    }

   /**
    * Invalidate index if it is not in lost of valid indices
    */
    template<typename T,typename Ti>
    void InvalidateIndices(const ROOT::RVec<T>& valid,Ti& index){
      if( ROOT::VecOps::Any((valid==index)))
	return;
      index=constant::InvalidIndex();
    }
  }//indice
}//rad
