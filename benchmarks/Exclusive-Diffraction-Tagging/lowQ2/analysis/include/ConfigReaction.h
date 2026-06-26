#pragma once

//!  Base configuration class for electron scattering reactions

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
*/
#include "DefineNames.h"
#include "RVecHelpers.h"
#include "ReactionUtilities.h"
#include "StringUtilities.h"
#include "Random.h"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>


namespace rad{
  namespace config{
    using rad::names::data_type::Rec;
    class ParticleCreator;
    
    
    //! Code simplifications
  
    using ROOT::RVecI;
  
    /*! RDFstep is used to update the dataframe handle after each filter or define step */
    // using RDFstep = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> ;
    using RDFstep = ROOT::RDF::RNode ;

    /*! RVecIndexMap is gives access to specific particle indices in action functions.
      Use reaction part index getters e.g. ScatEleIdx() to get the index from the map.
      e.g. rdf::RVecIndexMap& react; ... ; auto index = react[ScatEleIdx()]; 
    */
    using RVecIndexMap = ROOT::RVec<ROOT::RVecI>;
    
    void PrintDefinedColumnNames(RDFstep  df){
      std::cout<<"Print Column Names : ";
      auto cols =  df.GetDefinedColumnNames();
      for(auto& col:cols){
	std::cout<<col<<", ";
      }
      cout<<"\n";
    }
   void PrintAllColumnNames(RDFstep  df){
      std::cout<<"Print Column Names : ";
      auto cols =  df.GetColumnNames();
      for(auto& col:cols){
	std::cout<<col<<", ";
      }
      cout<<"\n";
    }

    bool ColumnExists(const string& col,RDFstep  df){
      auto cols =  df.GetDefinedColumnNames();
      if(std::find(cols.begin(),cols.end(),col)==cols.end()) return false;
      return true;
    }
    
    
    std::string as_string(std::string_view v) { 
      return {v.data(), v.size()}; 
    }
    
    //! Class definition

    class ConfigReaction {
    
    public:

      //     ConfigReaction(const std::string& treeName, const std::string& fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ) : _orig_df{treeName,{fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data()},columns},_curr_df{_orig_df}{
      ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ) : _orig_df{treeName,{fileNameGlob.data()},columns},_curr_df{_orig_df},_base_df{_orig_df},_treeName{treeName},_fileName{fileNameGlob}{
	if (fileNameGlob.empty()) {
	  throw std::invalid_argument("ConfigReaction: fileNameGlob cannot be empty.");
	}
	_orig_col_names = _orig_df.GetColumnNames();

    }
      ConfigReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ) : _orig_df{treeName,filenames,columns},_curr_df{_orig_df},_base_df{_orig_df},_treeName{treeName},_fileNames{filenames}{
	if (filenames.empty()) {
	  throw std::invalid_argument("ConfigReaction: fileNameGlob cannot be empty.");
	}
	_orig_col_names = _orig_df.GetColumnNames();
 	
      }

      //if creating from alternative data source
      ConfigReaction(ROOT::RDataFrame rdf ) : _orig_df{rdf},_curr_df{rdf},_base_df{rdf}{
	_orig_col_names = _orig_df.GetColumnNames();

   }
      
      ~ConfigReaction(){ 
	//	std::cout << "ConfigReaction destructor: " <<_triggerSnapshots.size() << std::endl;
	for (auto& trigger : _triggerSnapshots) {
	  if (trigger) trigger();
	}
	
	//std::cout << "ConfigReaction destructor: " << CurrFrame().Count().GetValue() << std::endl;
      }
      
      /** 
       * Get the current dataframe to add further actions
       */
      RDFstep CurrFrame(){return _curr_df;}
      /** 
       * Set the current dataframe after adding further actions
       */
      void setCurrFrame(RDFstep  df){ _curr_df = df ;}
      
      /**
       * Interface to RDataFrame Define
       */
      void Define(const string_view name,const string& expression){
	setCurrFrame(CurrFrame().Define(name,expression));
      }
      //  void Define(const TString& name,const string& expression){
      // 	setCurrFrame(CurrFrame().Define(name,expression));
      // }
      template<typename Lambda>
      void Define(const string_view name,Lambda&& func,const ROOT::RDF::ColumnNames_t&  columns ){
	setCurrFrame(CurrFrame().Define(name,func,columns));
      }

      /**
       * Call Define for predefined data types, pprepend type string
       */
        void DefineForAllTypes(const string& name,const string& expression){
	for(auto &atype:_type_comps){
	  if (atype.second.find("components_p4") == atype.second.end() ||
	      atype.second.find("components_p3") == atype.second.end()) {
            throw std::runtime_error("DefineForAllTypes: Missing 'components_p4' or 'components_p3' for type: " + atype.first);
	  }
	  TString type_expr = expression.data();
	  type_expr.ReplaceAll("components_p4",atype.second["components_p4"]);
	  type_expr.ReplaceAll("components_p3",atype.second["components_p3"]);
	  Define(atype.first + name.data(),type_expr.Data());
	}
      }
      template<typename Lambda>
      void DefineForAllTypes(const string& name,Lambda&& func,const ROOT::RDF::ColumnNames_t&  columns ){
	for(auto &atype:_type_comps){
	  if (atype.second.find("components_p4") == atype.second.end() ||
	      atype.second.find("components_p3") == atype.second.end()) {
            throw std::runtime_error("DefineForAllTypes: Missing 'components_p4' or 'components_p3' for type: " + atype.first);
	  }
	  ROOT::RDF::ColumnNames_t type_cols;
	  for(auto& acol:columns){
	    type_cols.push_back(atype.first + acol);
	  }
	  Define(atype.first + name.data(),func,type_cols);
	}
      }
      /**
       * @brief make a new column for each particle based on applying func_name
       *        the new columns will called name_particle
       * @param name base name if new variable
       * @param particles list of particles to define this variable from
       * @param func_name the function to apply, should be defined in a .h file
       * @param values variables to be used by func to define new column
       */
      void DefineForParticles(const string& name,const ROOT::RDF::ColumnNames_t &particles,const ROOT::RDF::ColumnNames_t &values, const std::string& func_name){
	//loop over particles
	for(const auto& p : particles){
	  auto selected_entries = values;
	  //create string selecting particle entry from the values arrays
	  std::for_each(selected_entries.begin(), selected_entries.end(),
			[&p](std::string& s) {
			  s += '[';
			  s += p;
			  s += ']';
			});
	  //create the function for the Define call
	  auto selected_func = rad::utils::createFunctionCallStringFromVec(func_name,selected_entries);
	  cout<<"DefineForParticles "<<name+"_"+p<<" "<<selected_func<<endl;
	  //Define this variable for this particle
	  Define(name+"_"+p,selected_func);
	}
      }

      /**
       * Interface to RDataFrame Redefine
       */
      void RedefineExpr(const string& name,const string& expression){
	setCurrFrame(CurrFrame().Redefine(name,expression));
      }
      void Redefine(const string& name,const string& expression){
	setCurrFrame(CurrFrame().Redefine(name,expression));
      }
      template<typename Lambda>
      void Redefine(const string& name,Lambda&& func,const ROOT::RDF::ColumnNames_t& columns = {}){
	setCurrFrame(CurrFrame().Redefine(name,func,columns));
      }
      /**
       * Interface to RDataFrame Redefine via any aliases that may be used
       */
      //The redefine with alias does not work as the branchname
      //is not a mathematical expression. This is not checked for
      //the Lambda version below and works.
      // void RedefineViaAlias(const string& alias,const string& expression){
      // 	RedefineExpr(_aliasMap[alias],expression);
      // }
      template<typename Lambda>
      void RedefineViaAlias(const string& alias,Lambda&& func,const ROOT::RDF::ColumnNames_t& columns ){
	//	Redefine(_aliasMap[alias],func,columns);
	
	auto it = _aliasMap.find(alias);
	if (it == _aliasMap.end()) {
	  throw std::invalid_argument("RedefineViaAlias: alias '" + alias + "' does not exist in _aliasMap.");
	}
	Redefine(it->second, std::forward<Lambda>(func), columns);
      }
 
      /** 
       * Add an alias for a branch and update the current frame to the aliased one
       */
      void setBranchAlias(const string& old_name,const string& new_name){
	//Check if the old_name column exists before creating the alias
	if (!OriginalColumnExists(old_name)) {
	  throw std::invalid_argument("setBranchAlias: Source column '" + old_name + "' does not exist in the DataFrame.");
	}
	_aliasMap[new_name] = old_name;
	setCurrFrame(CurrFrame().Alias(new_name,old_name));
      }
      /**
       * Interface to RDataFrame Filter
       */
      template<typename Lambda>
      void Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {},std::string name = "" ){
	setCurrFrame(CurrFrame().Filter(func,columns,name));
      }
      void Filter(const std::string& expression,const std::string& 	name ){
	setCurrFrame(CurrFrame().Filter(expression,name));
      }
     
      /**
       * Make a snapshot of newly defined columns
       */
      void BookLazySnapshot(const string& filename){
	try {
	  ROOT::RDF::RSnapshotOptions opts;
	  opts.fLazy = true;
	  auto cols = CurrFrame().GetDefinedColumnNames();
	  RemoveSnapshotColumns(cols);
	  auto snapshot_result = CurrFrame().Snapshot("rad_tree", filename, cols, opts);
	  _triggerSnapshots.emplace_back([snapshot = std::move(snapshot_result)]() mutable{});
	} catch (const std::exception& ex) {
	  std::cerr << "BookLazySnapshot failed: " << ex.what() << std::endl;
	  throw; // or handle gracefully
	}
      
      }
      
      
      void Snapshot(const string& filename){
	try {
	  auto cols = CurrFrame().GetDefinedColumnNames();
	  RemoveSnapshotColumns(cols);
	  CurrFrame().Snapshot("rad_tree", filename, cols);
	} catch (const std::exception& ex) {
	  std::cerr << "Snapshot failed: " << ex.what() << std::endl;
	  throw;
	}
      }

      virtual void RemoveSnapshotColumns(std::vector<string>& cols){
	cols.erase(std::remove(cols.begin(), cols.end(), names::ReactionMap() ), cols.end());

	//remove any columns with the DoNotWriteTag
	auto tag  = DoNotWriteTag();
	cols.erase( std::remove_if( cols.begin(), cols.end(),
				    [&tag]( const string& col ) -> bool
				    { return col.find(tag) != std::string::npos; } ),
		    cols.end() );

      }
      /** 
       * Set constant index in collection for particle
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setParticleIndex(const string& particle, const int idx, int pdg=0 ){
	// Check if particle index column already exists to avoid accidental overwrite
	if (ColumnExists(particle, CurrFrame())) {
	  throw std::invalid_argument("setParticleIndex: Index column for particle '" + particle + "' already exists!");
	}
	Define(particle,[idx](){return idx;},{});
	if(pdg!=0){
	  if(ColumnExists(Rec()+"pid",CurrFrame())||CheckAlias(Rec()+"pid")){
	    Define(particle+"_OK",Form("%spid[%s]==%d",Rec().data(),particle.data(),pdg));
	  }
	}
 	AddParticleName(particle);
 	AddFinalParticleName(particle);
     }
      
      /** 
       * Set function, func,  which defines variable index in collection for particle
       * The given function func is responsible for giving the correct position of particle in the collection
       * Names of required branches or identifiers given as vector of column names : columns
       * and update the current frame to the aliased one
       */
      template<typename Lambda>
      void setParticleIndex(const string& particle, Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {}, int pdg=0 ){
	// Check if particle index column already exists to avoid accidental overwrite
	if (ColumnExists(particle, CurrFrame())) {
	  throw std::invalid_argument("setParticleIndex: Index column for particle '" + particle + "' already exists!");
	}
	Define(particle,func,columns);

	if(pdg!=0){
	  if(ColumnExists(Rec()+"pid",CurrFrame())||CheckAlias(Rec()+"pid")){
	    Define(string(particle)+"_OK",Form("%spid[%s]==%d",Rec().data(),particle.data(),pdg));
	  }
	}
	AddParticleName(particle);
	AddFinalParticleName(particle);
      }

      /**
       * Make map that links particle names to indices in user functions
       * in C++ functions you can use the RVecIndexMap object indexed by 
       * name of the reaction component you need
       */
      virtual void makeParticleMap() {
	std::string particle_func("1E6+");
	for(auto& part : _particleNames){
	  particle_func+=part+"+";
	}
	particle_func.pop_back(); //remove last +

	Filter(particle_func.data(),"particle_list");

	PostParticles();
      }
      /**
       *Any additional stuff to be done after all particles have been indiced
       */
      virtual void PostParticles(){

      }
     // /**
     //   * Make map that links beam names to indices in user functions
     //   * in C++ functions you can use the RVecIndexMap object indexed by 
     //   * name of the reaction component you need
     //   *
     //   * This function must be implmented by a derived class
     //   */
     //  virtual void makeBeamMap() = 0;


       /**
       * Collect constant indices for final state mesons and baryons
       */
      void setMesonIndices(const RVecI& indices){
	Define(as_string(names::Mesons()),[indices](){return indices;},{});
      } 
      void setBaryonIndices(const RVecI& indices){
	Define(as_string(names::Baryons()),[indices](){return indices;},{});
      }

       /**
       * Collect variable indices for final state mesons and baryons
       */
      void setMesonParticles(const  ROOT::RDF::ColumnNames_t& particles){
	if(particles.empty()==true){
	  std::cout<<"setBaryonParticles "<<as_string(names::Mesons())<<std::endl;
	  Define(as_string(names::Mesons()),[](){return RVecI{-1};},{});
	  return;
	}
	setGroupParticles(as_string(names::Mesons()),particles);
      }
      void setBaryonParticles(const  ROOT::RDF::ColumnNames_t& particles){
	if(particles.empty()==true){
	  std::cout<<"setBaryonParticles "<<as_string(names::Baryons())<<std::endl;
	  Define(as_string(names::Baryons()),[](){return RVecI{-1};},{});
	  return;
	}
	setGroupParticles(as_string(names::Baryons()),particles);
      }
     
      /**
       * Generic function to collect variable indices for groups of particle
       * name = identifier of the indice collection, 
       * particles = list of particle names known to ConfigureReaction
       */

      /*
      //Would like the code below to work. Compiles OK, but the dataframe infers 0 arguments, rather than particles.size()
      template <typename... Args>
      void setGroupParticles(const string& name,const  ROOT::RDF::ColumnNames_t& particles){
	auto ncols=particles.size();
	//create lambda to return current indices of particles
	//args will be the data in the columns defined by particles
	auto func = [ncols](Args&&... args) {
	  RVecI vec;
	  vec.reserve(ncols);
	  (vec.emplace_back(std::forward<Args>(args)), ...);
	  return vec;};

	//need to test which of these works fastest
	// auto func = [ncols](Args&&... args) {
	//   RVecI vec;
	//   vec.reserve(ncols);
	//   for (auto i = 0UL; i < size; ++i) {
	//     ret.emplace_back(args[i]...);
	//   }
	  
	// auto func = [ncols](Args&&... args) {
	//   RVecI vec(ncols);
	//   for (auto i = 0UL; i < size; ++i) {
	//     ret[i]=(args[i]...);
	//   }
	//  return vec;};
	
   	setCurrFrame(CurrFrame().Define(names::Baryons(),func,particles));
   
      }
      */
       void setGroupParticles(const string& name,const ROOT::RDF::ColumnNames_t &particles){

	 auto pstring =reaction::util::ColumnsToString(particles); //"{p1,p2,p3,p4,...}"
	 pstring = pstring.substr(1,pstring.size() - 2); //remove {}
	 Define( name,Form("rad::helpers::Group<int>(%s)",pstring.data()) );
	 return;

	/*
	switch( particles.size() ){
	case 1 : 
	  Define(name,[](const int p0){return RVecI{p0};},particles);
	  break;
	case 2 : 
	  Define(name,[](const int p0,const int p1){return RVecI{p0,p1};},particles);
	  break;
	case 3 : 
	  Define(name,[](const int p0,const int p1,const int p2){return RVecI{p0,p1,p2};},particles);
	  break;
	case 4 : 
	  Define(name,[](const int p0,const int p1,const int p2,const int p3){return RVecI{p0,p1,p2,p3};},particles);
	  break;
	case 5 : 
	  Define(name,[](const int p0,const int p1,const int p2,const int p3,const int p4){return RVecI{p0,01,p2,p3,p4};},particles);
	  break;
	case 6 : 
	  Define(name,[](const int p0,const int p1,const int p2,const int p3,const int p4,const int p5){return RVecI{p0,01,p2,p3,p4,p5};},particles);
	  break;
	default:
	  std::cerr<<"setGroupParticles only defined up to 6 particles! "<<std::endl;
	  exit(0);
	  break;
	}
	*/   
      }
      /**
       * Set all Pid (aka PDG) values to -1 so particles ignored
       */
      template<typename Lambda>
      void  ApplyPidMask(const string& mask_name,const string& type,Lambda&& func,const ROOT::RDF::ColumnNames_t  columns = {}){
	//Define a column for the mask variable
	Define(mask_name,func,columns);
	//Apply mask to PID column, so entries indexed by -1 will not be used
	RedefineViaAlias(type+"_"+"pid",[](const RVecI& pid,const RVecI& mask){
	  return ROOT::VecOps::Where(mask,pid,-1);},
	  {type+"_"+"pid",mask_name.data()});
      }
      /**
       * create shortcut string for 3 and 4 momentum components
       */
      void AddType(const string& atype){
	if(_primary_type.empty()==true) _primary_type=atype;
	_type_comps[atype]["components_p4"] = Form("%spx,%spy,%spz,%sm",atype.data(),atype.data(),atype.data(),atype.data());
	_type_comps[atype]["components_p3"] = Form("%spx,%spy,%spz",atype.data(),atype.data(),atype.data());
     }
      /**
       *  change name of this first type column to same name without type prefix
       */
      void AliasToPrimaryType(const string& name){
	if(_primary_type.empty()==true) return;
	std::string fullName = _primary_type + name;
	if (!OriginalColumnExists(fullName)) {
	  throw std::invalid_argument("AliasToPrimaryType: Column '" + fullName + "' does not exist.");
	}
	setBranchAlias(_primary_type+name,name);
      }
     /**
      *  copy name of this first type column to same name without type prefix
       */
      void CopyToPrimaryType(const string& name){
	if(_primary_type.empty()==true) return;
	std::cout<<"CopyToPrimaryType "<<_primary_type<<" "<<name<<std::endl;
	Define(name,Form("return %s;",(_primary_type+name).data() ) );
      }


      /**
      * check if alias is used
      */
      bool CheckAlias(const string& alias){
	if(_aliasMap.find(alias) != _aliasMap.end()) return true;
	else return false;
      }
     
  

      std::map<string, std::map<string,string>> GetTypes() const {return _type_comps;}

      /**
       * Check if currently using type
       */
      bool CheckForType(const string& type){
	if(_type_comps.find(type) != _type_comps.end()) return true;
	else return false;
 
      }
      /**
       * get or set the base df
       * this should include the configuration of the input
       * data to rad format (arrays of px,px,pz,m)
       * sharing this allows parallisation between reactions
       */
      void setBaseFrame(RDFstep step){
	_base_df = step;
      }
      void setMyBaseFrame(){
	_base_df = CurrFrame();
      }
      RDFstep getBaseFrame() const{
	return _base_df;
      }
      RDFstep getOrigFrame() const{
	return _orig_df;
      }


      std::string GetTreeName() const {
	return _treeName;
      }
      std::string GetFileName() const {
	return _fileName;
      }
      std::vector<std::string> GetFileNames() const {
	return _fileNames;
      }
      
      bool OriginalColumnExists(const string& col){
	if(std::find(_orig_col_names.begin(),_orig_col_names.end(),col)==_orig_col_names.end()){
	  return false;
	}
	return true;
      }
 
      void AddParticleName(const std::string& particle){_particleNames.push_back(particle);}
      void AddFinalParticleName(const std::string& particle){_finalNames.push_back(particle);}
      const ROOT::RDF::ColumnNames_t& ParticleNames() const {return _particleNames;}
      const ROOT::RDF::ColumnNames_t& FinalParticleNames() const {return _finalNames;}
      
      const std::map<string,string>& AliasMap() const {return _aliasMap;}

      const std::string DoNotWriteTag(){return "__dnwtag";};

        void InitRandom(size_t seed){
	 if(_initRandom==true){
	  std::cout<<"Warning, ConfigReaction::InitRandom() random generator already initialised"<<std::endl;
	  return;
	}
	auto frame = CurrFrame();
	frame = frame.DefineSlot("RDF_Internal_RNG_Init", [seed](unsigned int slot) {
	  rad::random::initializeAllThreadRNGs(slot,seed);
	  return true; // Return dummy value
	}).Filter("RDF_Internal_RNG_Init");
	
	setCurrFrame(frame);
      }
      
    protected:

      bool _useBeamsFromMC=false; 
      
    
    private:
    
      /**
       * _bare_df
       * Base dataframe constructed prior to adding actions
       */
      ROOT::RDataFrame _orig_df;
      /** 
       * _curr_df 
       * Handle for the current dataframe state. 
       * This includes all added defines and filters so far.
       * Additional actions will be applied to this.
       */
      RDFstep _curr_df;
      /** 
       * _base_df 
       * Handle for the basic dataframe state. 
       * This will have columns organised for processing.
       * Typically after a SetAlias call.
       */
      RDFstep _base_df;
      
      /**
       * Store name of any aliased branches
       */
      std::map<string,string> _aliasMap;
      /**
       * type of data linked to component names e.g._type_comps["tru"]["components_p4"] = "tru_px,tru_py,tru_pz,tru_m"
       */
      std::map<string, std::map<string,string>> _type_comps;
      /*
       * First/Primary type name
       */
      std::string _primary_type;

      //keep dataset info 
      std::vector<std::string> _fileNames;//if given list of files
      std::string _fileName;//if single file (or wildcards)
      std::string _treeName;
      ROOT::RDF::ColumnNames_t  _particleNames; //list of all particles, so index calculation can be enforced at start of operations
      ROOT::RDF::ColumnNames_t  _finalNames; //list of detectable particles
      ROOT::RDF::ColumnNames_t _orig_col_names; //columns in origin tree
      
      //snapshot
      std::vector<std::function<void()>> _triggerSnapshots;

      bool _initRandom = false; //flag to ensure random generator only init once
      
    };//ConfigReaction

  }//config
}//rad
