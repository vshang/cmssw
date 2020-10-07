//This class holds functional blocks of a sector

 
#ifndef FPGASECTOR_H
#define FPGASECTOR_H

#include "FPGAInputLink.hh"
#include "FPGAAllStubs.hh"
#include "FPGAVMStubsTE.hh"
#include "FPGAVMStubsME.hh"
#include "FPGAStubPairs.hh"
#include "FPGAStubTriplets.hh"
#include "FPGATrackletParameters.hh"
#include "FPGATrackletProjections.hh"
#include "FPGAAllProjections.hh"
#include "FPGAVMProjections.hh"
#include "FPGACandidateMatch.hh"
#include "FPGAFullMatch.hh"
#include "FPGATrackFit.hh"
#include "FPGACleanTrack.hh"

#include "FPGAVMRouter.hh"
#include "FPGATrackletEngine.hh"
#include "FPGATrackletEngineDisplaced.hh"
#include "FPGATripletEngine.hh"
#include "FPGATrackletCalculator.hh"
#include "FPGATrackletCalculatorDisplaced.hh"
#include "FPGAProjectionRouter.hh"
#include "FPGAMatchEngine.hh"
#include "FPGAMatchCalculator.hh"
#include "FPGAMatchProcessor.hh"
#include "FPGAFitTrack.hh"
#include "FPGAPurgeDuplicate.hh"
#include "FPGAUtil.hh"

using namespace std;

class FPGASector{

public:

  FPGASector(unsigned int i){
    isector_=i;
    double dphi=2*M_PI/NSector;
    double dphiHG=0.5*dphisectorHG-M_PI/NSector;
    phimin_=isector_*dphi-dphiHG;
    phimax_=phimin_+dphi+2*dphiHG;
    phimin_-=M_PI/NSector;
    phimax_-=M_PI/NSector;
    phimin_=FPGAUtil::phiRange(phimin_);
    phimax_=FPGAUtil::phiRange(phimax_);
    if (phimin_>phimax_)  phimin_-=2*M_PI;
  }

  ~FPGASector() {
    for (unsigned i = 0; i < TRE_.size(); i++)
      delete TRE_[i];
  }

  bool addStub(L1TStub stub, string dtc) {

    bool add=false;
    
    double phi=stub.phi();
    //cout << "FPGASector::addStub layer phi phimin_ phimax_ : "<<stub.layer()+1<<" "<<" "<<dtc<<" "<<phi<<" "<<phimin_<<" "<<phimax_<<endl;
    double dphi=0.5*dphisectorHG-M_PI/NSector;
    
    static std::map<string,std::vector<int> > ILindex;
    std::vector<int>& tmp=ILindex[dtc];
    if (tmp.size()==0){
      //cout << "Adding entries for dtc : "<<dtc;
      for (unsigned int i=0;i<IL_.size();i++){
	if (IL_[i]->getName().find("_"+dtc)!=string::npos){
	  //cout << " Adding link "<<IL_[i]->getName()<<endl;
	  tmp.push_back(i);
	}
      }
	//cout << endl;
    }
    
    if (((phi>phimin_-dphi)&&(phi<phimax_+dphi))||
	((phi>2*M_PI+phimin_-dphi)&&(phi<2*M_PI+phimax_+dphi))) {
      FPGAStub fpgastub(stub,phimin_,phimax_);
      std::vector<int>& tmp=ILindex[dtc];
      assert(tmp.size()!=0);
      for (unsigned int i=0;i<tmp.size();i++){
	//cout << "Add stub to link"<<IL_[tmp[i]]->getName()<<endl;
	if (IL_[tmp[i]]->addStub(stub,fpgastub,dtc)) add=true;
      }
    }
    
    return add;
    
  }
    
  void addMem(string memType,string memName){
    if (memType=="InputLink:") {
      IL_.push_back(new FPGAInputLink(memName,isector_,phimin_,phimax_));
      Memories_[memName]=IL_.back();
      MemoriesV_.push_back(IL_.back());
    } else if (memType=="AllStubs:") {
      AS_.push_back(new FPGAAllStubs(memName,isector_,phimin_,phimax_));
      Memories_[memName]=AS_.back();
      MemoriesV_.push_back(AS_.back());
    } else if (memType=="VMStubsTE:") {
      VMSTE_.push_back(new FPGAVMStubsTE(memName,isector_,phimin_,phimax_));
      Memories_[memName]=VMSTE_.back();
      MemoriesV_.push_back(VMSTE_.back());
    } else if (memType=="VMStubsME:") {
      VMSME_.push_back(new FPGAVMStubsME(memName,isector_,phimin_,phimax_));
      Memories_[memName]=VMSME_.back();
      MemoriesV_.push_back(VMSME_.back());
    } else if (memType=="StubPairs:"||
               memType=="StubPairsDisplaced:") {
      SP_.push_back(new FPGAStubPairs(memName,isector_,phimin_,phimax_));
      Memories_[memName]=SP_.back();
      MemoriesV_.push_back(SP_.back());
    } else if (memType=="StubTriplets:") {
      ST_.push_back(new FPGAStubTriplets(memName,isector_,phimin_,phimax_));
      Memories_[memName]=ST_.back();
      MemoriesV_.push_back(ST_.back());
    } else if (memType=="TrackletParameters:") {
      TPAR_.push_back(new FPGATrackletParameters(memName,isector_,phimin_,phimax_));
      Memories_[memName]=TPAR_.back();
      MemoriesV_.push_back(TPAR_.back());
    } else if (memType=="TrackletProjections:") {
      TPROJ_.push_back(new FPGATrackletProjections(memName,isector_,phimin_,phimax_));
      Memories_[memName]=TPROJ_.back();
      MemoriesV_.push_back(TPROJ_.back());
    } else if (memType=="AllProj:") {
      AP_.push_back(new FPGAAllProjections(memName,isector_,phimin_,phimax_));
      Memories_[memName]=AP_.back();
      MemoriesV_.push_back(AP_.back());
    } else if (memType=="VMProjections:") {
      VMPROJ_.push_back(new FPGAVMProjections(memName,isector_,phimin_,phimax_));
      Memories_[memName]=VMPROJ_.back();
      MemoriesV_.push_back(VMPROJ_.back());
    } else if (memType=="CandidateMatch:") {
      CM_.push_back(new FPGACandidateMatch(memName,isector_,phimin_,phimax_));
      Memories_[memName]=CM_.back();
      MemoriesV_.push_back(CM_.back());
    } else if (memType=="FullMatch:") {
      FM_.push_back(new FPGAFullMatch(memName,isector_,phimin_,phimax_));
      Memories_[memName]=FM_.back();
      MemoriesV_.push_back(FM_.back());
    } else if (memType=="TrackFit:") {
      TF_.push_back(new FPGATrackFit(memName,isector_,phimin_,phimax_));
      Memories_[memName]=TF_.back();
      MemoriesV_.push_back(TF_.back());
    } else if (memType=="CleanTrack:") {
      CT_.push_back(new FPGACleanTrack(memName,isector_,phimin_,phimax_));
      Memories_[memName]=CT_.back();
      MemoriesV_.push_back(CT_.back());
    } else {
      cout << "Don't know of memory type: "<<memType<<endl;
      exit(0);
    }
    
  }

  void addProc(string procType,string procName){
    if (procType=="VMRouter:") {
      VMR_.push_back(new FPGAVMRouter(procName,isector_));
      Processes_[procName]=VMR_.back();
    } else if (procType=="TrackletEngine:") {
      TE_.push_back(new FPGATrackletEngine(procName,isector_));
      Processes_[procName]=TE_.back();
    } else if (procType=="TrackletEngineDisplaced:") {
      TED_.push_back(new FPGATrackletEngineDisplaced(procName,isector_));
      Processes_[procName]=TED_.back();
    } else if (procType=="TripletEngine:") {
      TRE_.push_back(new FPGATripletEngine(procName,isector_));
      Processes_[procName]=TRE_.back();
    } else if (procType=="TrackletCalculator:"||
	       procType=="TrackletDiskCalculator:") {
      TC_.push_back(new FPGATrackletCalculator(procName,isector_));
      Processes_[procName]=TC_.back();
    } else if (procType=="TrackletCalculatorDisplaced:") {
      TCD_.push_back(new FPGATrackletCalculatorDisplaced(procName,isector_));
      Processes_[procName]=TCD_.back();
    } else if (procType=="ProjectionRouter:") {
      PR_.push_back(new FPGAProjectionRouter(procName,isector_));
      Processes_[procName]=PR_.back();
    } else if (procType=="MatchEngine:") {
      ME_.push_back(new FPGAMatchEngine(procName,isector_));
      Processes_[procName]=ME_.back();
    } else if (procType=="MatchCalculator:"||
               procType=="DiskMatchCalculator:") {
      MC_.push_back(new FPGAMatchCalculator(procName,isector_));
      Processes_[procName]=MC_.back();
    } else if (procType=="MatchProcessor:") {
      MP_.push_back(new FPGAMatchProcessor(procName,isector_));
      Processes_[procName]=MP_.back();
    } else if (procType=="FitTrack:") {
      FT_.push_back(new FPGAFitTrack(procName,isector_));
      Processes_[procName]=FT_.back();
    } else if (procType=="PurgeDuplicate:") {
      PD_.push_back(new FPGAPurgeDuplicate(procName,isector_));
      Processes_[procName]=PD_.back();
    } else {
      cout << "Don't know of processing type: "<<procType<<endl;
      exit(0);      
    }
  }

  void addWire(string mem,string procinfull,string procoutfull){

    //cout << "Mem : "<<mem<<" input from "<<procinfull
    //	 << " output to "<<procoutfull<<endl;

    stringstream ss1(procinfull);
    string procin, output;
    getline(ss1,procin,'.');
    getline(ss1,output);

    stringstream ss2(procoutfull);
    string procout, input;
    getline(ss2,procout,'.');
    getline(ss2,input);

    //cout << "Procin  : "<<procin<<" "<<output<<endl;
    //cout << "Procout : "<<procout<<" "<<input<<endl;

    FPGAMemoryBase* memory=getMem(mem);

    if (procin!="") {
      FPGAProcessBase* inProc=getProc(procin);
      inProc->addOutput(memory,output);
      }

    if (procout!="") {
      FPGAProcessBase* outProc=getProc(procout);
      outProc->addInput(memory,input);
    }



  }

  FPGAProcessBase* getProc(string procName){

    map<string, FPGAProcessBase*>::iterator it=Processes_.find(procName);

    if (it!=Processes_.end()) {
      return it->second;
    }
    cout << "Could not find process with name : "<<procName<<endl;
    assert(0);
    return 0;
  }

  FPGAMemoryBase* getMem(string memName){

    map<string, FPGAMemoryBase*>::iterator it=Memories_.find(memName);

    if (it!=Memories_.end()) {
      return it->second;
    }
    cout << "Could not find memory with name : "<<memName<<endl;
    assert(0);
    return 0;
  }

  void writeInputStubs(bool first) {
    for (unsigned int i=0;i<IL_.size();i++){
      IL_[i]->writeStubs(first);
    }
  }

  void writeInputStubs_in2(bool first) {
    for (unsigned int i=0;i<IL_.size();i++){
      IL_[i]->writeInputStubs(first, writestubs_in2,padding);
    }
  }

  void writeVMSTE(bool first) {
    for (unsigned int i=0;i<VMSTE_.size();i++){
      VMSTE_[i]->writeStubs(first);
    }
  }
  
  void writeVMSME(bool first) {
    for (unsigned int i=0;i<VMSME_.size();i++){
      VMSME_[i]->writeStubs(first);
    }
  }

  void writeAS(bool first) {
    for (unsigned int i=0;i<AS_.size();i++){
      AS_[i]->writeStubs(first);
    }
  }

  void writeSP(bool first) {
    for (unsigned int i=0;i<SP_.size();i++){
      SP_[i]->writeSP(first);
    }
  }

  void writeST(bool first) {
    for (unsigned int i=0;i<ST_.size();i++){
      ST_[i]->writeST(first);
    }
  }

  void writeTPAR(bool first) {
    for (unsigned int i=0;i<TPAR_.size();i++){
      TPAR_[i]->writeTPAR(first);
    }
  }

  void writeTPROJ(bool first) {
    for (unsigned int i=0;i<TPROJ_.size();i++){
      TPROJ_[i]->writeTPROJ(first);
    }
  }

  void writeAP(bool first) {
    for (unsigned int i=0;i<AP_.size();i++){
      AP_[i]->writeAP(first);
    }
  }

  void writeVMPROJ(bool first) {
    for (unsigned int i=0;i<VMPROJ_.size();i++){
      VMPROJ_[i]->writeVMPROJ(first);
    }
  }

  void writeCM(bool first) {
    for (unsigned int i=0;i<CM_.size();i++){
     CM_[i]->writeCM(first);
    }
  }

  void writeMC(bool first) {
    for (unsigned int i=0;i<FM_.size();i++){
      FM_[i]->writeMC(first);
    }
  }

  void writeTF(bool first){
    for(unsigned int i=0; i<TF_.size(); ++i){
      TF_[i]->writeTF(first);
    }
  }

  void writeCT(bool first) {
    for(unsigned int i=0; i<CT_.size(); ++i){
      CT_[i]->writeCT(first);
    }
  }

  void clean() {
    if (writeNMatches) {
      int matchesL1=0;
      int matchesL3=0;
      int matchesL5=0;
      for(unsigned int i=0;i<TPAR_.size();i++) {
	TPAR_[i]->writeMatches(matchesL1,matchesL3,matchesL5);
      }
      static ofstream out("nmatchessector.txt");
      out <<matchesL1<<" "<<matchesL3<<" "<<matchesL5<<endl;
    }
    
    
    for(unsigned int i=0;i<MemoriesV_.size();i++) {
      MemoriesV_[i]->clean();
    }
  }

  void executeVMR(){

    if (writeIL) {
      static ofstream out("inputlink.txt");
      for (unsigned int i=0;i<IL_.size();i++){
	out<<IL_[i]->getName()<<" "<<IL_[i]->nStubs()<<endl;
      } 
    }
    
    for (unsigned int i=0;i<VMR_.size();i++){
      VMR_[i]->execute();
    }
  }

  void executeTE(){
    for (unsigned int i=0;i<TE_.size();i++){
      TE_[i]->execute();
    }
  }

  void executeTED(){
    for (unsigned int i=0;i<TED_.size();i++){
      TED_[i]->execute();
    }
  }

  void executeTRE(){
    for (unsigned int i=0;i<TRE_.size();i++){
      TRE_[i]->execute();
    }
  }

  void executeTC(){
    for (unsigned int i=0;i<TC_.size();i++){
      TC_[i]->execute();
    }
  }

  void executeTCD(){
    for (unsigned int i=0;i<TCD_.size();i++){
      TCD_[i]->execute();
    }
  }

  void executePR(){
    for (unsigned int i=0;i<PR_.size();i++){
      PR_[i]->execute();
    }
  }

  void executeME(){
    for (unsigned int i=0;i<ME_.size();i++){
      ME_[i]->execute();
    }
  }

  void executeMC(){
    for (unsigned int i=0;i<MC_.size();i++){
      MC_[i]->execute();
    }
  }

  void executeMP(){
    for (unsigned int i=0;i<MP_.size();i++){
      MP_[i]->execute();
    }
  }

  void executeFT(){
    for (unsigned int i=0;i<FT_.size();i++){
      FT_[i]->execute();
    }
  }

  void executePD(std::vector<FPGATrack*>& tracks){
    for (unsigned int i=0;i<PD_.size();i++){
      PD_[i]->execute(tracks);
    }
  }


  bool foundTrack(ofstream& outres, L1SimTrack simtrk){
    bool match=false;
    for (unsigned int i=0;i<TF_.size();i++){
      if (TF_[i]->foundTrack(outres,simtrk)) match=true;
    }
    return match;
  }

  std::vector<FPGATracklet*> getAllTracklets() {
    std::vector<FPGATracklet*> tmp;
    for(unsigned int i=0;i<TPAR_.size();i++){
      for(unsigned int j=0;j<TPAR_[i]->nTracklets();j++){
	tmp.push_back(TPAR_[i]->getFPGATracklet(j));
      }
    }
    return tmp;
  }

  std::vector<std::pair<FPGAStub*,L1TStub*> > getStubs() const {
    std::vector<std::pair<FPGAStub*,L1TStub*> > tmp;

    for(unsigned int imem=0;imem<IL_.size();imem++){
      for(unsigned int istub=0;istub<IL_[imem]->nStubs();istub++){
	tmp.push_back(IL_[imem]->getStub(istub));
      }
    }
    
    return tmp;
  }

  std::set<int> seedMatch(int itp) {
    std::set<int> tmpSeeds;
    for(unsigned int i=0;i<TPAR_.size();i++) {
      unsigned int nTracklet=TPAR_[i]->nTracklets();
      for(unsigned int j=0;j<nTracklet;j++){
	if (TPAR_[i]->getFPGATracklet(j)->tpseed()==itp) {
	  tmpSeeds.insert(TPAR_[i]->getFPGATracklet(j)->getISeed());
	}
      }
    }
    return tmpSeeds;
  }

  double phimin() const {return phimin_;}
  double phimax() const {return phimax_;}

private:

  int isector_;
  double phimin_;
  double phimax_;


  std::map<string, FPGAMemoryBase*> Memories_;
  std::vector<FPGAMemoryBase*> MemoriesV_;
  std::vector<FPGAInputLink*> IL_;
  std::vector<FPGAAllStubs*> AS_;
  std::vector<FPGAVMStubsTE*> VMSTE_;
  std::vector<FPGAVMStubsME*> VMSME_;
  std::vector<FPGAStubPairs*> SP_;
  std::vector<FPGAStubTriplets*> ST_;
  std::vector<FPGATrackletParameters*> TPAR_;
  std::vector<FPGATrackletProjections*> TPROJ_;
  std::vector<FPGAAllProjections*> AP_;
  std::vector<FPGAVMProjections*> VMPROJ_;
  std::vector<FPGACandidateMatch*> CM_;
  std::vector<FPGAFullMatch*> FM_;
  std::vector<FPGATrackFit*> TF_;
  std::vector<FPGACleanTrack*> CT_;
  
  std::map<string, FPGAProcessBase*> Processes_;
  std::vector<FPGAVMRouter*> VMR_;
  std::vector<FPGATrackletEngine*> TE_;
  std::vector<FPGATrackletEngineDisplaced*> TED_;
  std::vector<FPGATripletEngine*> TRE_;
  std::vector<FPGATrackletCalculator*> TC_;
  std::vector<FPGATrackletCalculatorDisplaced*> TCD_;
  std::vector<FPGAProjectionRouter*> PR_;
  std::vector<FPGAMatchEngine*> ME_;
  std::vector<FPGAMatchCalculator*> MC_;
  std::vector<FPGAMatchProcessor*> MP_;
  std::vector<FPGAFitTrack*> FT_;
  std::vector<FPGAPurgeDuplicate*> PD_;


  
};

#endif
