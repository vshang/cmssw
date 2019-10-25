#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseII.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

L1Analysis::L1AnalysisPhaseII::L1AnalysisPhaseII()
{
}

L1Analysis::L1AnalysisPhaseII::~L1AnalysisPhaseII()
{

}

void L1Analysis::L1AnalysisPhaseII::SetVertices(float z0Puppi, float z0VertexTDR, const edm::Handle< std::vector<l1t::Vertex> > l1vertices, const edm::Handle< std::vector<l1t::L1TkPrimaryVertex> > l1TkPrimaryVertex){
      edm::Handle<std::vector<l1t::Vertex> > l1vertextdr;
      l1extra_.z0Puppi=z0Puppi;
      l1extra_.z0VertexTDR=z0VertexTDR;

      for (unsigned int i=0; i<l1vertices->size() ; i++){
                  l1extra_.z0Vertices.push_back(l1vertices->at(i).z0());
                  l1extra_.nVertices++;
      }
      for (unsigned int i=0; i<l1TkPrimaryVertex->size() ; i++){
                  l1extra_.z0L1TkPV.push_back(l1TkPrimaryVertex->at(i).getZvertex());
                  l1extra_.sumL1TkPV.push_back(l1TkPrimaryVertex->at(i).getSum());
                  l1extra_.nL1TkPVs++;
      }


}


void L1Analysis::L1AnalysisPhaseII::SetCaloTau(const edm::Handle<l1t::TauBxCollection> calotau, unsigned maxL1Extra)
{
  for (int ibx = calotau->getFirstBX(); ibx <= calotau->getLastBX(); ++ibx) {
    for (l1t::TauBxCollection::const_iterator it=calotau->begin(ibx); it!=calotau->end(ibx) && l1extra_.nCaloTaus<maxL1Extra; it++){
      if (it->pt() > 0){
	l1extra_.caloTauEt .push_back(it->et());
	l1extra_.caloTauEta.push_back(it->eta());
	l1extra_.caloTauPhi.push_back(it->phi());
	l1extra_.caloTauIEt .push_back(it->hwPt());
	l1extra_.caloTauIEta.push_back(it->hwEta());
	l1extra_.caloTauIPhi.push_back(it->hwPhi());
	l1extra_.caloTauIso.push_back(it->hwIso());
	l1extra_.caloTauBx .push_back(ibx);
	l1extra_.caloTauTowerIPhi.push_back(it->towerIPhi());
	l1extra_.caloTauTowerIEta.push_back(it->towerIEta());
	l1extra_.caloTauRawEt.push_back(it->rawEt());
	l1extra_.caloTauIsoEt.push_back(it->isoEt());
	l1extra_.caloTauNTT.push_back(it->nTT());
	l1extra_.caloTauHasEM.push_back(it->hasEM());
	l1extra_.caloTauIsMerged.push_back(it->isMerged());
	l1extra_.caloTauHwQual.push_back(it->hwQual());
	l1extra_.nCaloTaus++;
      }
    }
  }
}


void L1Analysis::L1AnalysisPhaseII::SetCaloJet(const edm::Handle<l1t::JetBxCollection> jet, unsigned maxL1Extra, float caloJetHTT)
{

  double mHT30_px=0, mHT30_py=0, HT30=0;
  double mHT30_3p5_px=0, mHT30_3p5_py=0, HT30_3p5=0;

  for (int ibx = jet->getFirstBX(); ibx <= jet->getLastBX(); ++ibx) {
    for (l1t::JetBxCollection::const_iterator it=jet->begin(ibx); it!=jet->end(ibx) && l1extra_.nCaloJets<maxL1Extra; it++){
      if (it->pt() > 0){
      l1extra_.caloJetEt .push_back(it->et());
      l1extra_.caloJetEta.push_back(it->eta());
      l1extra_.caloJetPhi.push_back(it->phi());
      l1extra_.caloJetBx .push_back(ibx);
      l1extra_.nCaloJets++;

      if(it->et()>30 && fabs(it->eta())<2.4) {
                  HT30+=it->et();
                  mHT30_px+=it->px();
                  mHT30_py+=it->py();
      }
      if(it->et()>30 && fabs(it->eta())<3.5) {
                  HT30_3p5+=it->et();
                  mHT30_3p5_px+=it->px();
                  mHT30_3p5_py+=it->py();
      }
      }
    }
  }

  l1extra_.caloJetHTDefault= caloJetHTT;

  l1extra_.caloJetMHTEt.push_back( sqrt(mHT30_px*mHT30_px+mHT30_py*mHT30_py) );
  l1extra_.caloJetMHTPhi.push_back( atan(mHT30_py/mHT30_px) );
  l1extra_.caloJetHT.push_back( HT30 );

  l1extra_.caloJetMHTEt.push_back( sqrt(mHT30_3p5_px*mHT30_3p5_px+mHT30_3p5_py*mHT30_3p5_py) );
  l1extra_.caloJetMHTPhi.push_back( atan(mHT30_3p5_py/mHT30_3p5_px) );
  l1extra_.caloJetHT.push_back( HT30_3p5 );

  l1extra_.nCaloJetMHT=2;




}


void L1Analysis::L1AnalysisPhaseII::SetMuon(const edm::Handle<l1t::MuonBxCollection> muon, unsigned maxL1Extra)
{
  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
    for (l1t::MuonBxCollection::const_iterator it=muon->begin(ibx); it!=muon->end(ibx) && l1extra_.nGlobalMuons<maxL1Extra; it++){
      if (it->pt() > 0){
	l1extra_.globalMuonPt .push_back(it->et());
	l1extra_.globalMuonEta.push_back(it->eta());
	l1extra_.globalMuonPhi.push_back(it->phi());
	l1extra_.globalMuonEtaAtVtx.push_back(it->etaAtVtx());
	l1extra_.globalMuonPhiAtVtx.push_back(it->phiAtVtx());
	l1extra_.globalMuonIEt .push_back(it->hwPt());
	l1extra_.globalMuonIEta.push_back(it->hwEta());
	l1extra_.globalMuonIPhi.push_back(it->hwPhi());
	l1extra_.globalMuonIEtaAtVtx.push_back(it->hwEtaAtVtx());
	l1extra_.globalMuonIPhiAtVtx.push_back(it->hwPhiAtVtx());
	l1extra_.globalMuonIDEta.push_back(it->hwDEtaExtra());
	l1extra_.globalMuonIDPhi.push_back(it->hwDPhiExtra());
	l1extra_.globalMuonChg.push_back(it->charge());
	l1extra_.globalMuonIso.push_back(it->hwIso());
	l1extra_.globalMuonQual.push_back(it->hwQual());
	l1extra_.globalMuonTfMuonIdx.push_back(it->tfMuonIndex());
	l1extra_.globalMuonBx .push_back(ibx);
	l1extra_.nGlobalMuons++;
      }
    }
  }
}

void L1Analysis::L1AnalysisPhaseII::SetMuonKF(const edm::Handle<l1t::RegionalMuonCandBxCollection> standaloneMuon, unsigned maxL1Extra, unsigned int muonDetector)
{
  for (int ibx = standaloneMuon->getFirstBX(); ibx <= standaloneMuon->getLastBX(); ++ibx) {
    for (l1t::RegionalMuonCandBxCollection::const_iterator it=standaloneMuon->begin(ibx); it!=standaloneMuon->end(ibx) && l1extra_.nStandaloneMuons<maxL1Extra; it++){
      if (it->hwPt() > 0){
//      std::cout<<"hwPt vs hwPt2?"<<it->hwPt()*0.5<<" "<<it->hwPt2()<<"   "<<it->hwSign()<<"   "<<muonDetector<<std::endl;
      l1extra_.standaloneMuonPt .push_back(it->hwPt()*0.5);
      l1extra_.standaloneMuonPt2 .push_back(it->hwPt2());
      l1extra_.standaloneMuonDXY .push_back(it->hwDXY());
      l1extra_.standaloneMuonEta.push_back(it->hwEta()*0.010875);
      l1extra_.standaloneMuonPhi.push_back(l1t::MicroGMTConfiguration::calcGlobalPhi( it->hwPhi(), it->trackFinderType(), it->processor() )*2*M_PI/576);
      l1extra_.standaloneMuonChg.push_back( pow(-1,it->hwSign() ) );
      l1extra_.standaloneMuonQual.push_back(it->hwQual());
      l1extra_.standaloneMuonRegion.push_back(muonDetector);
      l1extra_.standaloneMuonBx .push_back(ibx);
      l1extra_.nStandaloneMuons++;
      }
    }
  }
}
  // RegionalMuons are a bit ugly... why not global muons?? 
  /// Get compressed pT (returned int * 0.5 = pT (GeV))
  //    const int hwPt() const { return m_hwPt; };
  //        /// Get compressed local phi (returned int * 2*pi/576 = local phi in rad)
  //            const int hwPhi() const { return m_hwPhi; };
  //                /// Get compressed eta (returned int * 0.010875 = eta)
  //                    const int hwEta() const { return m_hwEta; };
  //                        /// Get charge sign bit (charge = (-1)^(sign))
  //                        const int hwSign() const { return m_hwSign; };


//EG (seeded by Phase 2 Objects )
void L1Analysis::L1AnalysisPhaseII::SetEG(const edm::Handle<l1t::EGammaBxCollection> EG, const edm::Handle<l1t::EGammaBxCollection> EGHGC,unsigned maxL1Extra)
{
  for(l1t::EGammaBxCollection::const_iterator it=EG->begin(); it!=EG->end() && l1extra_.nEG<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.EGEt .push_back(it->et());
    l1extra_.EGEta.push_back(it->eta());
    l1extra_.EGPhi.push_back(it->phi());
    l1extra_.EGIso.push_back(it->isoEt());
    l1extra_.EGHwQual.push_back(it->hwQual());
    l1extra_.EGBx.push_back(0);//it->bx());
    l1extra_.EGHGC.push_back(0);
    bool quality= ( ( it->hwQual() >> 1 ) & 1   ) > 0 ;
    l1extra_.EGPassesLooseTrackID.push_back(quality); 
    quality= ( ( it->hwQual() >> 2 ) & 1 ) > 0 ;
    l1extra_.EGPassesPhotonID.push_back(quality);
    l1extra_.nEG++;
  }
  }

  for(l1t::EGammaBxCollection::const_iterator it=EGHGC->begin(); it!=EGHGC->end() && l1extra_.nEG<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.EGEt .push_back(it->et());
    l1extra_.EGEta.push_back(it->eta());
    l1extra_.EGPhi.push_back(it->phi());
    l1extra_.EGIso.push_back(it->isoEt());
    l1extra_.EGHwQual.push_back(it->hwQual());
    l1extra_.EGBx.push_back(0);//it->bx());
    l1extra_.EGHGC.push_back(1);
    bool quality= (it->hwQual() ==5 ) ;
    l1extra_.EGPassesLooseTrackID.push_back(quality);  
    l1extra_.EGPassesPhotonID.push_back(quality);
    l1extra_.nEG++;
  }
  }
}

// TrkEG (seeded by Phase 2 Objects)
void L1Analysis::L1AnalysisPhaseII::SetTkEG(const edm::Handle<l1t::L1TkElectronParticleCollection> tkElectron, const edm::Handle<l1t::L1TkElectronParticleCollection> tkElectronHGC,  unsigned maxL1Extra)
{
  for(l1t::L1TkElectronParticleCollection::const_iterator it=tkElectron->begin(); it!=tkElectron->end() && l1extra_.nTkElectrons<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.tkElectronEt .push_back(it->et());
    l1extra_.tkElectronEta.push_back(it->eta());
    l1extra_.tkElectronPhi.push_back(it->phi());
        int chargeFromCurvature = (it->trackCurvature() > 0)? 1 : -1 ; // ThisIsACheck
    l1extra_.tkElectronChg.push_back(chargeFromCurvature);
    l1extra_.tkElectronzVtx.push_back(it->getTrkzVtx());
    l1extra_.tkElectronTrkIso.push_back(it->getTrkIsol());
    l1extra_.tkElectronHwQual.push_back(it->getEGRef()->hwQual());
    l1extra_.tkElectronEGRefPt.push_back(it->getEGRef()->et());
    l1extra_.tkElectronEGRefEta.push_back(it->getEGRef()->eta());
    l1extra_.tkElectronEGRefPhi.push_back(it->getEGRef()->phi());
    l1extra_.tkElectronBx.push_back(0);//it->bx());
    l1extra_.tkElectronHGC.push_back(0);
    bool quality=  ( ( it->getEGRef()->hwQual() >> 1 ) & 1  )> 0; // LooseTrackID should be the second bit 
    l1extra_.tkElectronPassesLooseTrackID.push_back(quality);
    quality=  ( ( it->getEGRef()->hwQual() >> 2 ) & 1  )> 0; // LooseTrackID should be the second bit 
    l1extra_.tkElectronPassesPhotonID.push_back(quality);
    l1extra_.nTkElectrons++;
  }}

  for(l1t::L1TkElectronParticleCollection::const_iterator it=tkElectronHGC->begin(); it!=tkElectronHGC->end() && l1extra_.nTkElectrons<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.tkElectronEt .push_back(it->et());
    l1extra_.tkElectronEta.push_back(it->eta());
    l1extra_.tkElectronPhi.push_back(it->phi());
        int chargeFromCurvature = (it->trackCurvature() > 0)? 1 : -1 ; // ThisIsACheck
    l1extra_.tkElectronChg.push_back(chargeFromCurvature);
    l1extra_.tkElectronzVtx.push_back(it->getTrkzVtx());
    l1extra_.tkElectronTrkIso.push_back(it->getTrkIsol());
    l1extra_.tkElectronHwQual.push_back(it->getEGRef()->hwQual());
    l1extra_.tkElectronEGRefPt.push_back(it->getEGRef()->et());
    l1extra_.tkElectronEGRefEta.push_back(it->getEGRef()->eta());
    l1extra_.tkElectronEGRefPhi.push_back(it->getEGRef()->phi());
    l1extra_.tkElectronBx.push_back(0);//it->bx());
    l1extra_.tkElectronHGC.push_back(1);
    bool quality= (it->getEGRef()->hwQual() ==5 ) ;
    l1extra_.tkElectronPassesLooseTrackID.push_back(quality);
    l1extra_.tkElectronPassesPhotonID.push_back(quality);
    l1extra_.nTkElectrons++;
  }}

}

void L1Analysis::L1AnalysisPhaseII::SetTkEGV2(const edm::Handle<l1t::L1TkElectronParticleCollection> tkElectronV2, const edm::Handle<l1t::L1TkElectronParticleCollection> tkElectronV2HGC,unsigned maxL1Extra)
{
  for(l1t::L1TkElectronParticleCollection::const_iterator it=tkElectronV2->begin(); it!=tkElectronV2->end() && l1extra_.nTkElectronsV2<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.tkElectronV2Et .push_back(it->et());
    l1extra_.tkElectronV2Eta.push_back(it->eta());
    l1extra_.tkElectronV2Phi.push_back(it->phi());
        int chargeFromCurvature = (it->trackCurvature() > 0)? 1 : -1 ; // ThisIsACheck
    l1extra_.tkElectronV2Chg.push_back(chargeFromCurvature);
    l1extra_.tkElectronV2zVtx.push_back(it->getTrkzVtx());
    l1extra_.tkElectronV2TrkIso.push_back(it->getTrkIsol());
    l1extra_.tkElectronV2HwQual.push_back(it->getEGRef()->hwQual());
    l1extra_.tkElectronV2EGRefPt.push_back(it->getEGRef()->et());
    l1extra_.tkElectronV2EGRefEta.push_back(it->getEGRef()->eta());
    l1extra_.tkElectronV2EGRefPhi.push_back(it->getEGRef()->phi());
    l1extra_.tkElectronV2Bx.push_back(0);//it->bx());
    l1extra_.tkElectronV2HGC.push_back(0);
    bool quality=( ( it->getEGRef()->hwQual() >> 1 ) & 1   )> 0;  
    l1extra_.tkElectronV2PassesLooseTrackID.push_back(quality);
    quality=( ( it->getEGRef()->hwQual() >> 2 ) & 1   )> 0;  
    l1extra_.tkElectronV2PassesPhotonID.push_back(quality);
    l1extra_.nTkElectronsV2++;
  }}

  for(l1t::L1TkElectronParticleCollection::const_iterator it=tkElectronV2HGC->begin(); it!=tkElectronV2HGC->end() && l1extra_.nTkElectronsV2<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.tkElectronV2Et .push_back(it->et());
    l1extra_.tkElectronV2Eta.push_back(it->eta());
    l1extra_.tkElectronV2Phi.push_back(it->phi());
        int chargeFromCurvature = (it->trackCurvature() > 0)? 1 : -1 ; // ThisIsACheck
    l1extra_.tkElectronV2Chg.push_back(chargeFromCurvature);
    l1extra_.tkElectronV2zVtx.push_back(it->getTrkzVtx());
    l1extra_.tkElectronV2TrkIso.push_back(it->getTrkIsol());
    l1extra_.tkElectronV2HwQual.push_back(it->getEGRef()->hwQual());
    l1extra_.tkElectronV2EGRefPt.push_back(it->getEGRef()->et());
    l1extra_.tkElectronV2EGRefEta.push_back(it->getEGRef()->eta());
    l1extra_.tkElectronV2EGRefPhi.push_back(it->getEGRef()->phi());
    l1extra_.tkElectronV2Bx.push_back(0);//it->bx());
    l1extra_.tkElectronV2HGC.push_back(1);
    bool quality= (it->getEGRef()->hwQual() ==5 ) ;
    l1extra_.tkElectronV2PassesLooseTrackID.push_back(quality);
    l1extra_.tkElectronV2PassesPhotonID.push_back(quality);
    l1extra_.nTkElectronsV2++;
  }}


}

void L1Analysis::L1AnalysisPhaseII::SetTkEM(const edm::Handle<l1t::L1TkEmParticleCollection> tkPhoton, const edm::Handle<l1t::L1TkEmParticleCollection> tkPhotonHGC, unsigned maxL1Extra)
{
  for(l1t::L1TkEmParticleCollection::const_iterator it=tkPhoton->begin(); it!=tkPhoton->end() && l1extra_.nTkPhotons<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.tkPhotonEt .push_back(it->et());
    l1extra_.tkPhotonEta.push_back(it->eta());
    l1extra_.tkPhotonPhi.push_back(it->phi());
    l1extra_.tkPhotonTrkIso.push_back(it->getTrkIsol());
    l1extra_.tkPhotonTrkIsoPV.push_back(it->getTrkIsolPV());
    l1extra_.tkPhotonBx.push_back(0);//it->bx());
    l1extra_.tkPhotonHwQual.push_back(it->getEGRef()->hwQual());
    l1extra_.tkPhotonEGRefPt.push_back(it->getEGRef()->et());
    l1extra_.tkPhotonEGRefEta.push_back(it->getEGRef()->eta());
    l1extra_.tkPhotonEGRefPhi.push_back(it->getEGRef()->phi());
    l1extra_.tkPhotonHGC.push_back( 0 );
    bool quality= ( ( it->getEGRef()->hwQual() >> 1 ) & 1   ) > 0 ; 
    l1extra_.tkPhotonPassesLooseTrackID.push_back(quality);
    quality= ( ( it->getEGRef()->hwQual() >> 2 ) & 1   ) > 0 ; // Photon Id should be the third bit 
    l1extra_.tkPhotonPassesPhotonID.push_back(quality);
    l1extra_.nTkPhotons++;
  }}
  for(l1t::L1TkEmParticleCollection::const_iterator it=tkPhotonHGC->begin(); it!=tkPhotonHGC->end() && l1extra_.nTkPhotons<maxL1Extra; it++){
    if (it->et() > 5){
    l1extra_.tkPhotonEt .push_back(it->et());
    l1extra_.tkPhotonEta.push_back(it->eta());
    l1extra_.tkPhotonPhi.push_back(it->phi());
    l1extra_.tkPhotonTrkIso.push_back(it->getTrkIsol());
    l1extra_.tkPhotonTrkIsoPV.push_back(it->getTrkIsolPV());
    l1extra_.tkPhotonBx.push_back(0);//it->bx());
    l1extra_.tkPhotonHwQual.push_back(it->getEGRef()->hwQual());
    l1extra_.tkPhotonEGRefPt.push_back(it->getEGRef()->et());
    l1extra_.tkPhotonEGRefEta.push_back(it->getEGRef()->eta());
    l1extra_.tkPhotonEGRefPhi.push_back(it->getEGRef()->phi());
    l1extra_.tkPhotonHGC.push_back( 1 );
    bool quality= (it->getEGRef()->hwQual() ==5 ) ;
    l1extra_.tkPhotonPassesLooseTrackID.push_back(quality);
    l1extra_.tkPhotonPassesPhotonID.push_back(quality);
    l1extra_.nTkPhotons++;
  }}
}

void L1Analysis::L1AnalysisPhaseII::SetTrkTau(const edm::Handle<l1t::L1TrkTauParticleCollection> tkTau, unsigned maxL1Extra)
{

  for(l1t::L1TrkTauParticleCollection::const_iterator it=tkTau->begin(); it!=tkTau->end() && l1extra_.nTkTau<maxL1Extra; it++){

    l1extra_.tkTauEt.push_back(it->et());
    l1extra_.tkTauEta.push_back(it->eta());
    l1extra_.tkTauPhi.push_back(it->phi());
    l1extra_.tkTauTrkIso.push_back(it->getIso());
    l1extra_.tkTauBx.push_back(0);//it->bx());
    l1extra_.nTkTau++;
  }
}

void L1Analysis::L1AnalysisPhaseII::SetCaloTkTau(const edm::Handle<l1t::L1CaloTkTauParticleCollection> caloTkTau, unsigned maxL1Extra)
{

  for(l1t::L1CaloTkTauParticleCollection::const_iterator it=caloTkTau->begin(); it!=caloTkTau->end() && l1extra_.nCaloTkTau<maxL1Extra; it++){

    l1extra_.caloTkTauEt.push_back(it->et());
    l1extra_.caloTkTauEta.push_back(it->eta());
    l1extra_.caloTkTauPhi.push_back(it->phi());
    l1extra_.caloTkTauTrkIso.push_back(it->getVtxIso());
    l1extra_.caloTkTauBx.push_back(0);//it->bx());
    l1extra_.nCaloTkTau++;
  }
}

void L1Analysis::L1AnalysisPhaseII::SetTkEGTau(const edm::Handle<l1t::L1TkEGTauParticleCollection> tkEGTau, unsigned maxL1Extra)
{

  for(l1t::L1TkEGTauParticleCollection::const_iterator it=tkEGTau->begin(); it!=tkEGTau->end() && l1extra_.nTkEGTau<maxL1Extra; it++){

    l1extra_.tkEGTauEt.push_back(it->et());
    l1extra_.tkEGTauEta.push_back(it->eta());
    l1extra_.tkEGTauPhi.push_back(it->phi());
    l1extra_.tkEGTauTrkIso.push_back(it->getIso());
    l1extra_.tkEGTauBx.push_back(0);//it->bx());
    l1extra_.nTkEGTau++;
  }
}


// TkJet
void L1Analysis::L1AnalysisPhaseII::SetTkJet(const edm::Handle<l1t::L1TkJetParticleCollection> trackerJet, unsigned maxL1Extra)
{

  for(l1t::L1TkJetParticleCollection::const_iterator it=trackerJet->begin(); it!=trackerJet->end() && l1extra_.nTrackerJets<maxL1Extra; it++){
    l1extra_.trackerJetEt .push_back(it->et());
    l1extra_.trackerJetEta.push_back(it->eta());
    l1extra_.trackerJetPhi.push_back(it->phi());
    l1extra_.trackerJetzVtx.push_back(it->getJetVtx());
    l1extra_.trackerJetBx .push_back(0);//it->bx());
    l1extra_.nTrackerJets++;
  }
}

void L1Analysis::L1AnalysisPhaseII::SetTkCaloJet(const edm::Handle<l1t::L1TkJetParticleCollection> tkCaloJet, unsigned maxL1Extra)
{

  for(l1t::L1TkJetParticleCollection::const_iterator it=tkCaloJet->begin(); it!=tkCaloJet->end() && l1extra_.nTkCaloJets<maxL1Extra; it++){
    l1extra_.tkCaloJetEt .push_back(it->et());
    l1extra_.tkCaloJetEta.push_back(it->eta());
    l1extra_.tkCaloJetPhi.push_back(it->phi());
    l1extra_.tkCaloJetzVtx.push_back(it->getJetVtx());
    l1extra_.tkCaloJetBx .push_back(0);//it->bx());
    l1extra_.nTkCaloJets++;
  }
}




void L1Analysis::L1AnalysisPhaseII::SetTkMuon(const edm::Handle<l1t::L1TkMuonParticleCollection> muon, unsigned maxL1Extra)
{

  for(l1t::L1TkMuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && l1extra_.nTkMuons<maxL1Extra; it++){

    l1extra_.tkMuonPt .push_back(it->pt());
    l1extra_.tkMuonEta.push_back(it->eta());
    l1extra_.tkMuonPhi.push_back(it->phi());
    int chargeFromCurvature = (it->trackCurvature() > 0)? 1 : -1 ; // ThisIsACheck
    l1extra_.tkMuonChg.push_back( chargeFromCurvature);
    l1extra_.tkMuonTrkIso.push_back(it->getTrkIsol());
    if(it->muonDetector()!=3){
    l1extra_.tkMuonMuRefPt.push_back(it->getMuRef()->hwPt()*0.5);
    l1extra_.tkMuonMuRefEta.push_back(it->getMuRef()->hwEta()*0.010875);
    l1extra_.tkMuonMuRefPhi.push_back(l1t::MicroGMTConfiguration::calcGlobalPhi( it->getMuRef()->hwPhi(), it->getMuRef()->trackFinderType(), it->getMuRef()->processor() )*2*M_PI/576);
    l1extra_.tkMuonDRMuTrack.push_back(it->dR());
    l1extra_.tkMuonNMatchedTracks.push_back(it->nTracksMatched());
    l1extra_.tkMuonQual .push_back(it->getMuRef()->hwQual());
    l1extra_.tkMuonMuRefChg.push_back(pow(-1,it->getMuRef()->hwSign() ) );
    }else {
    l1extra_.tkMuonMuRefPt.push_back(-777);
    l1extra_.tkMuonMuRefEta.push_back(-777);
    l1extra_.tkMuonMuRefPhi.push_back(-777);
    l1extra_.tkMuonDRMuTrack.push_back(-777);
    l1extra_.tkMuonNMatchedTracks.push_back(0);
    l1extra_.tkMuonQual .push_back(999);
    l1extra_.tkMuonMuRefChg.push_back(0);
    }
    l1extra_.tkMuonRegion.push_back(it->muonDetector());
    l1extra_.tkMuonzVtx.push_back(it->getTrkzVtx());
    l1extra_.tkMuonBx .push_back(0); //it->bx());
    l1extra_.nTkMuons++;
  }
}


void L1Analysis::L1AnalysisPhaseII::SetTkMuonStubs(const edm::Handle<l1t::L1TkMuonParticleCollection> muon, unsigned maxL1Extra, unsigned int muonDetector)
{

  for(l1t::L1TkMuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && l1extra_.nTkMuonStubs<maxL1Extra; it++){

    if (muonDetector==1 && fabs(it->eta())>=0.9) continue;
    if (muonDetector==3 && fabs(it->eta())<1.2) continue;

    l1extra_.tkMuonStubsPt .push_back( it->pt());
    l1extra_.tkMuonStubsEta.push_back(it->eta());
    l1extra_.tkMuonStubsPhi.push_back(it->phi());
    l1extra_.tkMuonStubsChg.push_back(it->charge());
    l1extra_.tkMuonStubsTrkIso.push_back(it->getTrkIsol());
    l1extra_.tkMuonStubszVtx.push_back(it->getTrkzVtx());
    l1extra_.tkMuonStubsBx .push_back(0); //it->bx());
    l1extra_.tkMuonStubsQual .push_back(1);
    l1extra_.tkMuonStubsBarrelStubs.push_back(it->getMatchedStubs().size());
    l1extra_.tkMuonStubsRegion.push_back(muonDetector);
    l1extra_.nTkMuonStubs++;
  }
}

void L1Analysis::L1AnalysisPhaseII::SetTkMuonStubsOMTF(const edm::Handle<l1t::BayesMuCorrTrackBxCollection> muon, unsigned maxL1Extra, unsigned int muonDetector)
{

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
   for(l1t::BayesMuCorrTrackBxCollection::const_iterator it=muon->begin(ibx); it!=muon->end(ibx) && l1extra_.nTkMuonStubs<maxL1Extra; it++){

    // filtering to avoid overlaps
    if (fabs(it->eta())<0.9 || fabs(it->eta())>=1.2) continue;

    l1extra_.tkMuonStubsPt .push_back( it->pt());
    l1extra_.tkMuonStubsEta.push_back(it->eta());
    l1extra_.tkMuonStubsPhi.push_back(it->phi());
    l1extra_.tkMuonStubsChg.push_back(it->hwSign());
    l1extra_.tkMuonStubsTrkIso.push_back(0);
    l1extra_.tkMuonStubszVtx.push_back(it->getTtTrackPtr()->getPOCA(4).z());
    l1extra_.tkMuonStubsBx .push_back(ibx);
    l1extra_.tkMuonStubsQual .push_back(1);
    l1extra_.tkMuonStubsBarrelStubs.push_back(0);
    l1extra_.tkMuonStubsRegion.push_back(muonDetector);
    l1extra_.nTkMuonStubs++;
   }
  }

}

void L1Analysis::L1AnalysisPhaseII::SetHSCP(const edm::Handle<l1t::BayesMuCorrTrackBxCollection> muon, unsigned maxL1Extra)
{

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
   for(l1t::BayesMuCorrTrackBxCollection::const_iterator it=muon->begin(ibx); it!=muon->end(ibx) && l1extra_.nTkHSCPs<maxL1Extra; it++){

    l1extra_.tkHSCPsPt .push_back( it->pt());
    l1extra_.tkHSCPsEta.push_back(it->eta());
    l1extra_.tkHSCPsPhi.push_back(it->phi());
    l1extra_.tkHSCPsChg.push_back(it->hwSign());
    l1extra_.tkHSCPszVtx.push_back(it->getTtTrackPtr()->getPOCA(4).z());
    l1extra_.tkHSCPsBx .push_back(ibx);
    l1extra_.nTkHSCPs++;
   }
  }

}




void L1Analysis::L1AnalysisPhaseII::SetTkGlbMuon(const edm::Handle<l1t::L1TkGlbMuonParticleCollection> muon, unsigned maxL1Extra)
{
  for(l1t::L1TkGlbMuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && l1extra_.nTkGlbMuons<maxL1Extra; it++){

    l1extra_.tkGlbMuonPt .push_back( it->pt());
    l1extra_.tkGlbMuonEta.push_back(it->eta());
    l1extra_.tkGlbMuonPhi.push_back(it->phi());
    l1extra_.tkGlbMuonChg.push_back(it->charge());
    l1extra_.tkGlbMuonTrkIso.push_back(it->getTrkIsol());
    l1extra_.tkGlbMuonMuRefPt.push_back(it->getMuRef()->pt());
    l1extra_.tkGlbMuonMuRefEta.push_back(it->getMuRef()->eta());
    l1extra_.tkGlbMuonMuRefPhi.push_back(it->getMuRef()->phi());
    l1extra_.tkGlbMuonDRMuTrack.push_back(it->dR());
    l1extra_.tkGlbMuonNMatchedTracks.push_back(it->nTracksMatched());
    l1extra_.tkGlbMuonzVtx.push_back(it->getTrkzVtx());
    l1extra_.tkGlbMuonBx .push_back(0); //it->bx());
    l1extra_.tkGlbMuonQual .push_back(it->getMuRef()->hwQual());
    l1extra_.nTkGlbMuons++;
  }
}

// trackerMet
void L1Analysis::L1AnalysisPhaseII::SetTkMET(const edm::Handle<l1t::L1TkEtMissParticleCollection> trackerMets)
{
  for(l1t::L1TkEtMissParticleCollection::const_iterator it=trackerMets->begin(); it!=trackerMets->end(); it++) {
    l1extra_.trackerMetSumEt.    push_back( it->etTotal() );
    l1extra_.trackerMetEt.   push_back( it->etMiss() );
    l1extra_.trackerMetPhi.push_back( it->phi() );
    l1extra_.trackerMetBx. push_back( it->bx() );
    l1extra_.nTrackerMet++;
  }
}


void L1Analysis::L1AnalysisPhaseII::SetTkMHT(const edm::Handle<l1t::L1TkHTMissParticleCollection > trackerMHTs){
   // Hardcoding it like this, but this needs to be changed to a vector

  for(l1t::L1TkHTMissParticleCollection::const_iterator it=trackerMHTs->begin(); it!=trackerMHTs->end(); it++) {
    l1extra_.trackerHT.    push_back( it->EtTotal() );
    l1extra_.trackerMHT.   push_back( it->EtMiss() );
    l1extra_.trackerMHTPhi.push_back( it->phi() );
    l1extra_.nTrackerMHT++;
  }

}

void L1Analysis::L1AnalysisPhaseII::SetL1PfPhase1L1TJet(const      edm::Handle< std::vector<reco::CaloJet> >  l1L1PFPhase1L1Jet,    unsigned maxL1Extra){

   double mHT30_px=0, mHT30_py=0, HT30=0;
  double mHT30_3p5_px=0, mHT30_3p5_py=0, HT30_3p5=0; 


    for (reco::CaloJetCollection::const_iterator it=l1L1PFPhase1L1Jet->begin(); it!=l1L1PFPhase1L1Jet->end() && l1extra_.nPfPhase1L1Jets<maxL1Extra; it++){
      if (it->pt() > 0){
      l1extra_.pfPhase1L1JetEt .push_back(it->et());
      l1extra_.pfPhase1L1JetEta.push_back(it->eta());
      l1extra_.pfPhase1L1JetPhi.push_back(it->phi());
//      l1extra_.pfPhase1L1JetBx .push_back(0);
      l1extra_.nPfPhase1L1Jets++;
 
    if(it->et()>30 && fabs(it->eta())<2.4) {
                  HT30+=it->et();
                  mHT30_px+=it->px();
                  mHT30_py+=it->py();
      }
    if(it->et()>30 && fabs(it->eta())<3.5) {
                  HT30_3p5+=it->et();
                  mHT30_3p5_px+=it->px();
                  mHT30_3p5_py+=it->py();
      }



   }  
  }  

  l1extra_.nPfPhase1L1MHT=2;

  l1extra_.pfPhase1L1MHTEt.push_back( sqrt(mHT30_px*mHT30_px+mHT30_py*mHT30_py) );
  l1extra_.pfPhase1L1MHTPhi.push_back( atan(mHT30_py/mHT30_px) );
  l1extra_.pfPhase1L1HT.push_back( HT30 );

  l1extra_.pfPhase1L1MHTEt.push_back( sqrt(mHT30_3p5_px*mHT30_3p5_px+mHT30_3p5_py*mHT30_3p5_py) );
  l1extra_.pfPhase1L1MHTPhi.push_back( atan(mHT30_3p5_py/mHT30_3p5_px) );
  l1extra_.pfPhase1L1HT.push_back( HT30_3p5 );


}


/*
void L1Analysis::L1AnalysisPhaseII::SetPFJetForMET(const edm::Handle<l1t::PFJetCollection> PFJet, unsigned maxL1Extra)
{

  for(l1t::PFJetCollection::const_iterator it=PFJet->begin(); it!=PFJet->end() && l1extra_.nPuppiJetForMETs<maxL1Extra; it++){
    l1extra_.puppiJetForMETEt .push_back(it->pt());
    l1extra_.puppiJetForMETEtUnCorr .push_back(it->rawPt());
    l1extra_.puppiJetForMETEta.push_back(it->eta());
    l1extra_.puppiJetForMETPhi.push_back(it->phi());
//    l1extra_.puppiJetForMETzVtx.push_back(it->getJetVtx());
    l1extra_.puppiJetForMETBx .push_back(0);//it->bx());
    l1extra_.nPuppiJetForMETs++;
  }
}
*/

void L1Analysis::L1AnalysisPhaseII::SetPFJet(const edm::Handle<l1t::PFJetCollection> PFJet, unsigned maxL1Extra)
{
  double mHT30_px=0, mHT30_py=0, HT30=0;
  double mHT30_3p5_px=0, mHT30_3p5_py=0, HT30_3p5=0;

  for(l1t::PFJetCollection::const_iterator it=PFJet->begin(); it!=PFJet->end() && l1extra_.nPuppiJets<maxL1Extra; it++){
    l1extra_.puppiJetEt .push_back(it->pt());
    l1extra_.puppiJetEtUnCorr .push_back(it->rawPt());
    l1extra_.puppiJetEta.push_back(it->eta());
    l1extra_.puppiJetPhi.push_back(it->phi());
//    l1extra_.puppiJetzVtx.push_back(it->getJetVtx());
    l1extra_.puppiJetBx .push_back(0);//it->bx());
    l1extra_.nPuppiJets++;

    if(it->pt()>30 && fabs(it->eta())<2.4) { 
                  HT30+=it->pt();
                  mHT30_px+=it->px();
                  mHT30_py+=it->py();
      }
    if(it->pt()>30 && fabs(it->eta())<3.5) {
                  HT30_3p5+=it->pt();
                  mHT30_3p5_px+=it->px();
                  mHT30_3p5_py+=it->py();
      }


  }
  l1extra_.puppiMHTEt.push_back( sqrt(mHT30_px*mHT30_px+mHT30_py*mHT30_py) );
  l1extra_.puppiMHTPhi.push_back( atan(mHT30_py/mHT30_px) );
  l1extra_.puppiHT.push_back( HT30 );

  l1extra_.puppiMHTEt.push_back( sqrt(mHT30_3p5_px*mHT30_3p5_px+mHT30_3p5_py*mHT30_3p5_py) );
  l1extra_.puppiMHTPhi.push_back( atan(mHT30_3p5_py/mHT30_3p5_px) );
  l1extra_.puppiHT.push_back( HT30_3p5 );

  l1extra_.nPuppiMHT=2;

}

void L1Analysis::L1AnalysisPhaseII::SetL1METPF(const edm::Handle< std::vector<reco::PFMET> > l1MetPF)
{
  reco::PFMET met=l1MetPF->at(0);
  l1extra_.puppiMETEt = met.et();
  l1extra_.puppiMETPhi = met.phi();
}

void L1Analysis::L1AnalysisPhaseII::SetPFObjects(const edm::Handle< vector<l1t::PFCandidate> >  l1pfCandidates,  unsigned maxL1Extra)
{
      for (unsigned int i=0; i<l1pfCandidates->size() && l1extra_.nPFMuons<maxL1Extra; i++){
           //enum Kind { ChargedHadron=0, Electron=1, NeutralHadron=2, Photon=3, Muon=4 };
            if(abs(l1pfCandidates->at(i).id())==4){
                  l1extra_.pfMuonPt.push_back(l1pfCandidates->at(i).pt()); 
                  l1extra_.pfMuonChg.push_back(l1pfCandidates->at(i).charge());
                  l1extra_.pfMuonEta.push_back(l1pfCandidates->at(i).eta());
                  l1extra_.pfMuonPhi.push_back(l1pfCandidates->at(i).phi());
                  l1extra_.pfMuonzVtx.push_back(l1pfCandidates->at(i).pfTrack()->track()->getPOCA(4).z()); // check with Giovanni, there has to be a cleaner way to do this. nParam_=4 should not be hardcoded
                  l1extra_.nPFMuons++;
            }
      }

      for (unsigned int i=0; i<l1pfCandidates->size(); i++){
           //enum Kind { ChargedHadron=0, Electron=1, NeutralHadron=2, Photon=3, Muon=4 };
            if(abs(l1pfCandidates->at(i).id())!=4){
              //  std::cout<<"pf cand id: "<<l1pfCandidates->at(i).id()<<std::endl;
                  l1extra_.pfCandId.push_back(l1pfCandidates->at(i).id()); 
                  l1extra_.pfCandEt.push_back(l1pfCandidates->at(i).pt()); 
                  l1extra_.pfCandChg.push_back(l1pfCandidates->at(i).charge());
                  l1extra_.pfCandEta.push_back(l1pfCandidates->at(i).eta());
                  l1extra_.pfCandPhi.push_back(l1pfCandidates->at(i).phi());
                  if (l1pfCandidates->at(i).id()==0) {
                      l1extra_.pfCandzVtx.push_back(l1pfCandidates->at(i).pfTrack()->track()->getPOCA(4).z()); 
                  } else {
                      l1extra_.pfCandzVtx.push_back(9999.0);
                  };

                  l1extra_.nPFCands++;
            }
      }


}

void L1Analysis::L1AnalysisPhaseII::SetPFTaus(const edm::Handle< vector<l1t::L1PFTau> >  l1pfTaus,  unsigned maxL1Extra)
{

      for (unsigned int i=0; i<l1pfTaus->size() && l1extra_.nPFTaus<maxL1Extra; i++){
                   if(l1pfTaus->at(i).pt()<10) continue;
                   l1extra_.pfTauEt.push_back(l1pfTaus->at(i).pt());
                   l1extra_.pfTauEta.push_back(l1pfTaus->at(i).eta());
                   l1extra_.pfTauPhi.push_back(l1pfTaus->at(i).phi());
                   l1extra_.pfTauChg.push_back(l1pfTaus->at(i).charge());
                   l1extra_.pfTauType.push_back(l1pfTaus->at(i).tauType());
                   l1extra_.pfTauChargedIso.push_back(l1pfTaus->at(i).chargedIso());
                   unsigned int isoflag=l1pfTaus->at(i).tauIsoQuality();
                   l1extra_.pfTauIsoFlag.push_back(isoflag);
                   //std::cout<<l1pfTaus->at(i).pt()<<"   "<<l1pfTaus->at(i).chargedIso()<<"  "<<l1pfTaus->at(i).passTightIso()<<"  "<<l1extra_.pfTauIsoFlag[l1extra_.nPFTaus]<<"  "<<isoflag<<std::endl;
                   // VeryLoose: <50; Loose < 20; Medium<10; Tight<5 
                   isoflag=l1pfTaus->at(i).tauRelIsoQuality();
                   l1extra_.pfTauRelIsoFlag.push_back(isoflag);
                   l1extra_.pfTauPassesMediumIso.push_back(l1pfTaus->at(i).passMediumIso());

                   l1extra_.pfTauZ0.push_back(l1pfTaus->at(i).z0());
                   l1extra_.nPFTaus++;
      }

}


void L1Analysis::L1AnalysisPhaseII::SetHPSPFTaus(const edm::Handle<l1t::L1HPSPFTauCollection >  l1HPSPFTaus,  unsigned maxL1Extra)
{     

      for (unsigned int i=0; i<l1HPSPFTaus->size() && l1extra_.nHPSTaus<maxL1Extra; i++){
                   if(l1HPSPFTaus->at(i).pt()<10) continue;
                   l1extra_.hpsTauEt.push_back(l1HPSPFTaus->at(i).pt());
                   l1extra_.hpsTauEta.push_back(l1HPSPFTaus->at(i).eta());
                   l1extra_.hpsTauPhi.push_back(l1HPSPFTaus->at(i).phi());
                   l1extra_.hpsTauChg.push_back(l1HPSPFTaus->at(i).charge());
                   l1extra_.hpsTauType.push_back(l1HPSPFTaus->at(i).tauType());
                   l1extra_.hpsTauPassTightRelIso.push_back(l1HPSPFTaus->at(i).passTightRelIso());
                   l1extra_.hpsTauZ0.push_back(l1HPSPFTaus->at(i).primaryVertex()->z0()); 
                   l1extra_.nHPSTaus++;
      }

}



void L1Analysis::L1AnalysisPhaseII::SetNNTaus(const edm::Handle< vector<l1t::PFTau> >  l1nnTaus,  unsigned maxL1Extra)
{
 
      for (unsigned int i=0; i<l1nnTaus->size() && l1extra_.nNNTaus<maxL1Extra; i++){
                   if(l1nnTaus->at(i).pt()<10) continue;
                   l1extra_.nnTauEt.push_back(l1nnTaus->at(i).pt());
                   l1extra_.nnTauEta.push_back(l1nnTaus->at(i).eta());
                   l1extra_.nnTauPhi.push_back(l1nnTaus->at(i).phi());
                   l1extra_.nnTauChg.push_back(l1nnTaus->at(i).charge());
                   l1extra_.nnTauChargedIso.push_back(l1nnTaus->at(i).chargedIso());
                   l1extra_.nnTauFullIso.push_back(l1nnTaus->at(i).fullIso());
                   l1extra_.nnTauID.push_back(l1nnTaus->at(i).id());
                   l1extra_.nnTauPassLooseNN.push_back(l1nnTaus->at(i).passLooseNN());
                   l1extra_.nnTauPassLoosePF.push_back(l1nnTaus->at(i).passLoosePF());
                   l1extra_.nnTauPassTightPF.push_back(l1nnTaus->at(i).passTightPF());
                   l1extra_.nnTauPassTightNN.push_back(l1nnTaus->at(i).passTightNN());
                   l1extra_.nNNTaus++;
      }
 
}



void L1Analysis::L1AnalysisPhaseII::SetBsCands (const      edm::Handle< std::vector<l1t::L1TkBsCandidate>>  l1TkBs,    unsigned maxL1Extra, int kind)
{
      for (unsigned int i=0; i<l1TkBs->size() && l1extra_.nTkBsCands<maxL1Extra; i++){
           l1extra_.tkBsCandPt.push_back(l1TkBs->at(i).pt());
           l1extra_.tkBsCandMass.push_back(l1TkBs->at(i).p4().M());
           l1extra_.tkBsCandEta.push_back(l1TkBs->at(i).eta());
           l1extra_.tkBsCandPhi.push_back(l1TkBs->at(i).phi());
           l1extra_.tkBsCandPhi1Pt.push_back(l1TkBs->at(i).getPhiCandidate(0).pt());
           l1extra_.tkBsCandPhi2Pt.push_back(l1TkBs->at(i).getPhiCandidate(1).pt());
           l1extra_.tkBsCandPhi1Mass.push_back(l1TkBs->at(i).getPhiCandidate(0).p4().M());
           l1extra_.tkBsCandPhi2Mass.push_back(l1TkBs->at(i).getPhiCandidate(1).p4().M());
           l1extra_.tkBsCandPhi1Phi.push_back(l1TkBs->at(i).getPhiCandidate(0).phi());
           l1extra_.tkBsCandPhi2Phi.push_back(l1TkBs->at(i).getPhiCandidate(1).phi());
           l1extra_.tkBsCandPhi1Eta.push_back(l1TkBs->at(i).getPhiCandidate(0).eta());
           l1extra_.tkBsCandPhi2Eta.push_back(l1TkBs->at(i).getPhiCandidate(1).eta());
           l1extra_.tkBsCandDRPhiPair.push_back(l1TkBs->at(i).dRPhiPair());
           l1extra_.tkBsCandDxyPhiPair.push_back(l1TkBs->at(i).dxyPhiPair());
           l1extra_.tkBsCandDzPhiPair.push_back(l1TkBs->at(i).dRPhiPair());
           l1extra_.tkBsCandKind.push_back(kind);
           l1extra_.nTkBsCands++; 
      }


}

/*
void L1Analysis::L1AnalysisPhaseII::SetDiMuonTk(const edm::Handle<l1t::L1TkMuonParticleCollection> muon, unsigned maxL1Extra)
{

  for(l1t::L1TkMuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end(); it++){
       for(l1t::L1TkMuonParticleCollection::const_iterator it2=it+1; it2!=muon->end(); it2++){
            if (it==it2) continue;
            double dz=it->getTrkzVtx()-it2->getTrkzVtx();

            // Compute the centers of the circles
            double ra=1./it->getTrkPtr()->getRInv(); // I am assuming this is correct, not sure 
            // centres and radii of track circles
              // double phi = momentum.phi();
              //   // xc = signed_r*sin(phi) + position.x();
              //   // yc = -signed_r*cos(phi) + position.y();
            double xca=ra*sin(it->getTrkPtr()->getMomentum(4).phi())+it->getTrkPtr()->getPOCA(4).x(); // that 4 should not be harcoded 
            double yca=-ra*cos(it->getTrkPtr()->getMomentum(4).phi())+it->getTrkPtr()->getPOCA(4).y();      

            double rb=1./it2->getTrkPtr()->getRInv(); // I am assuming this is correct, not sure 
            double xcb=rb*sin(it2->getTrkPtr()->getMomentum(4).phi())+it2->getTrkPtr()->getPOCA(4).x();
            double ycb=-rb*cos(it2->getTrkPtr()->getMomentum(4).phi())+it2->getTrkPtr()->getPOCA(4).y();   

            // points of closest approach in transverse plane
            // First check the circles cross
               double xg1, yg1, xg2, yg2;

             // transverseCoord  start function - move to separate cc/h
            int flag = transverseCoord(xca, yca, ra, xcb, ycb, rb, xg1, yg1, xg2, yg2); 

            if (flag==0) continue;

            GlobalPoint positionA= it->getTrkPtr()->getPOCA(4) ;//GlobalPoint(xga, yga, zga);
            GlobalPoint positionB= it2->getTrkPtr()->getPOCA(4) ; //GlobalPoint(xga, yga, zga);
            GlobalVector momentumA=it->getTrkPtr()->getMomentum(4);
            GlobalVector momentumB=it2->getTrkPtr()->getMomentum(4);
            double xga, yga, zga, xgb, ygb, zgb;

            if (flag == 1) {
             // two crossing points on each track in transverse plane
             // select point for which z-coordinates on the 2 tracks are the closest
             double za1 = zCoord(momentumA, positionA, ra, xca, yca, xg1, yg1);
             double zb1 = zCoord(momentumB, positionB, rb, xcb, ycb, xg1, yg1);
             double za2 = zCoord(momentumA, positionA, ra, xca, yca, xg2, yg2);
             double zb2 = zCoord(momentumB, positionB, rb, xcb, ycb, xg2, yg2);

             if (abs(zb1 - za1) < abs(zb2 - za2)) {
               xga = xg1; yga = yg1; zga = za1; zgb = zb1;
             }
             else {
               xga = xg2; yga = yg2; zga = za2; zgb = zb2;
             }
             xgb = xga; ygb = yga;
           }
           else {
             // one point of closest approach on each track in transverse plane
                xga = xg1; yga = yg1;
                zga = zCoord(momentumA, positionA, ra, xca, yca, xga, yga);
                xgb = xg2; ygb = yg2;
                zgb = zCoord(momentumB, positionB, rb, xcb, ycb, xgb, ygb);
              }

            GlobalPoint  posA = GlobalPoint(xga, yga, zga);
            GlobalPoint  posB = GlobalPoint(xgb, ygb, zgb);

//            std::cout<<dz<<"  "<<it->pt()<<"  "<<ra<<"  "<<xca<<"   "<<yca<<" ; "<<it2->pt()<<"  "<<rb<<"  "<<xcb<<"   "<<ycb<<"  ->"<< flag<<"  --->" <<(posB - posA).mag()<<std::endl;  // this is not working ? 
  }
 }
}

int 
L1Analysis::L1AnalysisPhaseII::transverseCoord(double cxa, double cya, double ra, 
				       double cxb, double cyb, double rb, 
				       double & xg1, double & yg1, 
				       double & xg2, double & yg2)
{
  int flag = 0;
  double x1, y1, x2, y2;

  // new reference frame with origin in (cxa, cya) and x-axis 
  // directed from (cxa, cya) to (cxb, cyb)

  double d_ab = sqrt((cxb - cxa)*(cxb - cxa) + (cyb - cya)*(cyb - cya));
  if (d_ab == 0) { // concentric circles
    return 0;
  }
  // elements of rotation matrix
  double u = (cxb - cxa) / d_ab;
  double v = (cyb - cya) / d_ab;

  // conditions for circle intersection
  if (d_ab <= ra + rb && d_ab >= abs(rb - ra)) {

    // circles cross each other
    flag = 1;

    // triangle (ra, rb, d_ab)
    double cosphi = (ra*ra - rb*rb + d_ab*d_ab) / (2*ra*d_ab);
    double sinphi2 = 1. - cosphi*cosphi;
    if (sinphi2 < 0.) { sinphi2 = 0.; cosphi = 1.; }

    // intersection points in new frame
    double sinphi = sqrt(sinphi2);
    x1 = ra*cosphi; y1 = ra*sinphi; x2 = x1; y2 = -y1;
  } 
  else if (d_ab > ra + rb) {

    // circles are external to each other
    flag = 2;

    // points of closest approach in new frame 
    // are on line between 2 centers
    x1 = ra; y1 = 0; x2 = d_ab - rb; y2 = 0;
  }
  else if (d_ab < abs(rb - ra)) {

    // circles are inside each other
    flag = 2;

    // points of closest approach in new frame are on line between 2 centers
    // choose 2 closest points
    double sign = 1.;
    if (ra <= rb) sign = -1.;
    x1 = sign*ra; y1 = 0; x2 = d_ab + sign*rb; y2 = 0;
  }
  else {
    return 0;
  }

  // intersection points in global frame, transverse plane
  xg1 = u*x1 - v*y1 + cxa; yg1 = v*x1 + u*y1 + cya;
  xg2 = u*x2 - v*y2 + cxa; yg2 = v*x2 + u*y2 + cya;

  return flag;
}


double 
L1Analysis::L1AnalysisPhaseII::zCoord(const GlobalVector& mom, 
			      const GlobalPoint& pos, 
			      double r, double xc, double yc, 
			      double xg, double yg)
{

  // starting point
  double x = pos.x(); double y = pos.y(); double z = pos.z();

  double px = mom.x(); double py = mom.y(); double pz = mom.z();

  // rotation angle phi from starting point to crossing point (absolute value)
  // -- compute sin(phi/2) if phi smaller than pi/4, 
  // -- cos(phi) if phi larger than pi/4
  double phi = 0.;
  double sinHalfPhi = sqrt((x-xg)*(x-xg) + (y-yg)*(y-yg))/(2*r);
  if (sinHalfPhi < 0.383) { // sin(pi/8)
    phi = 2*asin(sinHalfPhi);
  }
  else {
    double cosPhi = ((x-xc)*(xg-xc) + (y-yc)*(yg-yc))/(r*r);
    if (std::abs(cosPhi) > 1) cosPhi = (cosPhi > 0 ? 1 : -1);
    phi = abs(acos(cosPhi));
  }
  // -- sign of phi
  double signPhi = ((x - xc)*(yg - yc) - (xg - xc)*(y - yc) > 0) ? 1. : -1.;

  // sign of track angular momentum
  // if rotation is along angular momentum, delta z is along pz
  double signOmega = ((x - xc)*py - (y - yc)*px > 0) ? 1. : -1.;

  // delta z
  // -- |dz| = |cos(theta) * path along helix|
  //         = |cos(theta) * arc length along circle / sin(theta)|
  double dz = signPhi*signOmega*(pz/mom.transverse())*phi*r;

  return z + dz;
}
*/


