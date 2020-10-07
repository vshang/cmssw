#include "L1Trigger/VertexFinder/interface/RecoVertexWithTP.h"

namespace l1tVertexFinder {

RecoVertexWithTP::RecoVertexWithTP(const float z0) :
  z0_(z0),
  pT_(0.0),
  met_(-999.)
{
}


RecoVertexWithTP::RecoVertexWithTP(RecoVertex& vertex, std::map<const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>>, const L1TrackTruthMatched*> trackAssociationMap)
{
  z0_ = -999.;
  pT_ = -9999.;
  met_ = -999.;


  // loop over base fitted tracks in reco vertex and find the corresponding TP
  // track using the TTTrack - L1TrackTruthMatched map from above
  for (const auto& trackIt : vertex.tracks()) {
    // using insert ensures that true tracks are also stored in vertex object
    insert(trackAssociationMap[trackIt->getTTTrackPtr()]);
  }
}

void RecoVertexWithTP::insert(const L1TrackTruthMatched* fitTrack)
{
  tracks_.push_back(fitTrack);
  if (fitTrack->getMatchedTP() != nullptr and fitTrack->getMatchedTP()->physicsCollision())
    trueTracks_.insert(fitTrack->getMatchedTP());
}

void RecoVertexWithTP::computeParameters(bool weightedmean)
{
  pT_ = 0.;
  z0_ = 0.;
  met_ = 0.;
  metX_ = 0.;
  metY_ = 0.;
  float z0square = 0.;
  highPt_ = false;
  highestPt_ = 0.;
  numHighPtTracks_ = 0;
  // unsigned int overflows = 0;
  float SumZ_pT = 0.;
  float SumZ = 0.;
  for (const L1TrackTruthMatched* track : tracks_) {
    // if(track->pt() < 100.){
    pT_ += track->pt();
    SumZ += track->z0();
    SumZ_pT += track->z0() * track->pt();
    z0square += track->z0() * track->z0();
    // if(track->getNumStubs() > 4){

    if (track->pt() > 999. and track->getNumStubs() < 6) {
      metX_ += 100. * cos(track->phi0());
      metY_ += 100. * sin(track->phi0());
    }
    else {
      metX_ += track->pt() * cos(track->phi0());
      metY_ += track->pt() * sin(track->phi0());
    }

    // }
    if (track->pt() > 15.) {
      highPt_ = true;
      highestPt_ = track->pt();
      // 	numHighPtTracks_++;
    }
    // } else{
    // 	overflows++;
    // }
  }
  // unsigned int divider = tracks_.size() - overflows;
  // if(divider > 0)	{

  if (weightedmean) {
    z0_ = SumZ_pT / pT_;
  }
  else {
    z0_ = SumZ / tracks_.size();
  }

  met_ = sqrt(metX_ * metX_ + metY_ * metY_);
  z0square /= tracks_.size();
  z0width_ = sqrt(fabs(z0_ * z0_ - z0square));
}

} // end ns l1tVertexFinder
