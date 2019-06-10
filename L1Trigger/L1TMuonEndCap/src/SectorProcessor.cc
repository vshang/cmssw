#include "L1Trigger/L1TMuonEndCap/interface/SectorProcessor.h"

SectorProcessor::SectorProcessor() {}

SectorProcessor::~SectorProcessor() {}

void SectorProcessor::configure(const EMTFSetup* setup, int verbose, int endcap, int sector) {
  emtf_assert(setup != nullptr);
  emtf_assert(emtf::MIN_ENDCAP <= endcap && endcap <= emtf::MAX_ENDCAP);
  emtf_assert(emtf::MIN_TRIGSECTOR <= sector && sector <= emtf::MAX_TRIGSECTOR);

  setup_ = setup;
  verbose_ = verbose;
  endcap_ = endcap;
  sector_ = sector;
}

void SectorProcessor::process(const edm::EventID& event_id,
                              const TriggerPrimitiveCollection& muon_primitives,
                              EMTFHitCollection& out_hits,
                              EMTFTrackCollection& out_tracks) const {
  auto cfg = setup_->getVersionControl();

  // ___________________________________________________________________________
  // List of converted hits, extended from previous BXs
  // deque (double-ended queue) is similar to a vector, but allows insertion or deletion of elements at both beginning and end
  std::deque<EMTFHitCollection> extended_conv_hits;

  // List of best track candidates, extended from previous BXs
  std::deque<EMTFTrackCollection> extended_best_track_cands;

  // Map of pattern detector --> lifetime, tracked across BXs
  std::map<pattern_ref_t, int> patt_lifetime_map;

  // ___________________________________________________________________________
  // Run each sector processor for every BX, taking into account the BX window

  int delayBX = cfg.bxWindow_ - 1;

  for (int bx = cfg.minBX_; bx <= cfg.maxBX_ + delayBX; ++bx) {
    if (verbose_ > 0) {  // debug
      std::cout << "Event: " << event_id << " Endcap: " << endcap_ << " Sector: " << sector_ << " BX: " << bx
                << std::endl;
    }

    process_single_bx(
        bx, muon_primitives, out_hits, out_tracks, extended_conv_hits, extended_best_track_cands, patt_lifetime_map);

    // Drop earliest BX outside of BX window
    if (bx >= cfg.minBX_ + delayBX) {
      extended_conv_hits.pop_front();

      int n = emtf::zone_array<int>().size();
      extended_best_track_cands.erase(extended_best_track_cands.end() - n,
                                      extended_best_track_cands.end());  // pop_back
    }
  }  // end loop over bx


  // ___________________________________________________________________________
  // Debug
  bool dump_emtf_hits = false;

  if (dump_emtf_hits && endcap_ == 2 && sector_ == 6) {
    int i = 0;
    const char subsystem_names[][4] = {"DT","CSC","RPC","GEM","ME0"};
    for (const auto& h : out_hits) {
      std::cout << "conv_hit #" << i++ << " [" << subsystem_names[h.Subsystem()] << "]: BX " << h.BX()
        << ", endcap " << h.Endcap() << ", sector " << h.PC_sector()
        << ", station " << h.Station() << ", ring " << h.Ring()
        << ", chamber " << h.Chamber() << ", roll " << h.Roll()
        << ", subsector " << h.Subsector() << ", CSC ID " << h.CSC_ID()
        << ", strip " << h.Strip() << ", wire " << h.Wire() << ", pattern " << h.Pattern() << ", bend " << h.Bend() << ", quality " << h.Quality()
        << ", phi_fp " << h.Phi_fp() << ", theta_fp " << h.Theta_fp() << ", phi_sim " << h.Phi_sim() << ", theta_sim " << h.Theta_sim()
        << std::endl;
    }
  }

  return;
}

void SectorProcessor::process_single_bx(int bx,
                                        const TriggerPrimitiveCollection& muon_primitives,
                                        EMTFHitCollection& out_hits,
                                        EMTFTrackCollection& out_tracks,
                                        std::deque<EMTFHitCollection>& extended_conv_hits,
                                        std::deque<EMTFTrackCollection>& extended_best_track_cands,
                                        std::map<pattern_ref_t, int>& patt_lifetime_map) const {
  auto cfg = setup_->getVersionControl();

  auto tp_geom_ = &(setup_->getGeometryTranslator());
  auto pc_lut_ = &(setup_->getSectorProcessorLUT());
  auto pt_assign_engine_ = setup_->getPtAssignmentEngine();

  // ___________________________________________________________________________
  // Configure

  PrimitiveSelection prim_sel;
  prim_sel.configure(verbose_,
                     endcap_,
                     sector_,
                     bx,
                     cfg.bxShiftCSC_,
                     cfg.bxShiftRPC_,
                     cfg.bxShiftGEM_,
                     cfg.bxShiftME0_,
                     cfg.includeNeighbor_,
                     cfg.duplicateTheta_,
                     cfg.bugME11Dupes_);

  PrimitiveConversion prim_conv;
  prim_conv.configure(tp_geom_,
                      pc_lut_,
                      verbose_,
                      endcap_,
                      sector_,
                      bx,
                      cfg.bxShiftCSC_,
                      cfg.bxShiftRPC_,
                      cfg.bxShiftGEM_,
                      cfg.bxShiftME0_,
                      cfg.zoneBoundaries_,
                      cfg.zoneOverlap_,
                      cfg.duplicateTheta_,
                      cfg.fixZonePhi_,
                      cfg.useNewZones_,
                      cfg.fixME11Edges_,
                      cfg.bugME11Dupes_);

  PatternRecognition patt_recog;
  patt_recog.configure(verbose_,
                       endcap_,
                       sector_,
                       bx,
                       cfg.bxWindow_,
                       cfg.pattDefinitions_,
                       cfg.symPattDefinitions_,
                       cfg.useSymPatterns_,
                       cfg.maxRoadsPerZone_,
                       cfg.useSecondEarliest_);

  PrimitiveMatching prim_match;
  prim_match.configure(
      verbose_, endcap_, sector_, bx, cfg.fixZonePhi_, cfg.useNewZones_, cfg.bugSt2PhDiff_, cfg.bugME11Dupes_);

  AngleCalculation angle_calc;
  angle_calc.configure(verbose_,
                       endcap_,
                       sector_,
                       bx,
                       cfg.bxWindow_,
                       cfg.thetaWindow_,
                       cfg.thetaWindowZone0_,
                       cfg.bugME11Dupes_,
                       cfg.bugAmbigThetaWin_,
                       cfg.twoStationSameBX_);

  BestTrackSelection btrack_sel;
  btrack_sel.configure(verbose_,
                       endcap_,
                       sector_,
                       bx,
                       cfg.bxWindow_,
                       cfg.maxRoadsPerZone_,
                       cfg.maxTracks_,
                       cfg.useSecondEarliest_,
                       cfg.bugSameSectorPt0_);

  SingleHitTrack single_hit;
  single_hit.configure(verbose_, endcap_, sector_, bx, cfg.maxTracks_, cfg.useSingleHits_);

  PtAssignment pt_assign;
  pt_assign.configure(pt_assign_engine_,
                      verbose_,
                      endcap_,
                      sector_,
                      bx,
                      cfg.readPtLUTFile_,
                      cfg.fixMode15HighPt_,
                      cfg.bug9BitDPhi_,
                      cfg.bugMode7CLCT_,
                      cfg.bugNegPt_,
                      cfg.bugGMTPhi_,
                      cfg.promoteMode7_,
                      cfg.modeQualVer_);

  std::map<int, TriggerPrimitiveCollection> selected_dt_map;
  std::map<int, TriggerPrimitiveCollection> selected_csc_map;
  std::map<int, TriggerPrimitiveCollection> selected_rpc_map;
  std::map<int, TriggerPrimitiveCollection> selected_gem_map;
  std::map<int, TriggerPrimitiveCollection> selected_me0_map;
  std::map<int, TriggerPrimitiveCollection> selected_prim_map;
  std::map<int, TriggerPrimitiveCollection> inclusive_selected_prim_map;

  EMTFHitCollection conv_hits;  // "converted" hits converted by primitive converter
  EMTFHitCollection inclusive_conv_hits;

  emtf::zone_array<EMTFRoadCollection> zone_roads;  // each zone has its road collection

  emtf::zone_array<EMTFTrackCollection> zone_tracks;  // each zone has its track collection

  EMTFTrackCollection best_tracks;  // "best" tracks selected from all the zones

  // ___________________________________________________________________________
  // Process

  // Select muon primitives that belong to this sector and this BX.
  // Put them into maps with an index that roughly corresponds to
  // each input link.
  // From src/PrimitiveSelection.cc
  prim_sel.process(emtf::DTTag(), muon_primitives, selected_dt_map);
  prim_sel.process(emtf::CSCTag(), muon_primitives, selected_csc_map);
  prim_sel.process(emtf::RPCTag(), muon_primitives, selected_rpc_map);
  prim_sel.process(emtf::GEMTag(), muon_primitives, selected_gem_map);
  prim_sel.process(emtf::ME0Tag(), muon_primitives, selected_me0_map);
  prim_sel.merge(
      selected_dt_map, selected_csc_map, selected_rpc_map, selected_gem_map, selected_me0_map, selected_prim_map);

  // Convert trigger primitives into "converted" hits
  // A converted hit consists of integer representations of phi, theta, and zones
  // From src/PrimitiveConversion.cc
#ifdef PHASE_TWO_TRIGGER
  // Exclude Phase 2 trigger primitives before running the rest of EMTF
  prim_conv.process(selected_prim_map, conv_hits);
  EMTFHitCollection tmp_conv_hits;
  for (const auto& conv_hit : conv_hits) {
    if (prim_conv.is_valid_for_run2(conv_hit)) {
      tmp_conv_hits.push_back(conv_hit);
    }
  }
  extended_conv_hits.push_back(tmp_conv_hits);
#else
  prim_conv.process(selected_prim_map, conv_hits);
  extended_conv_hits.push_back(conv_hits);
#endif

  {
    // Keep all the converted hits for the use of data-emulator comparisons.
    // They include the extra ones that are not used in track building and the subsequent steps.
    prim_sel.merge_no_truncate(selected_dt_map,
                               selected_csc_map,
                               selected_rpc_map,
                               selected_gem_map,
                               selected_me0_map,
                               inclusive_selected_prim_map);
    prim_conv.process(inclusive_selected_prim_map, inclusive_conv_hits);

    // Clear the input maps to save memory
    selected_dt_map.clear();
    selected_csc_map.clear();
    selected_rpc_map.clear();
    selected_gem_map.clear();
    selected_me0_map.clear();
  }

  // Detect patterns in all zones, find 3 best roads in each zone
  // From src/PatternRecognition.cc
  patt_recog.process(extended_conv_hits, patt_lifetime_map, zone_roads);

  // Match the trigger primitives to the roads, create tracks
  // From src/PrimitiveMatching.cc
  prim_match.process(extended_conv_hits, zone_roads, zone_tracks);

  // Calculate deflection angles for each track and fill track variables
  // From src/AngleCalculation.cc
  angle_calc.process(zone_tracks);
  extended_best_track_cands.insert(
      extended_best_track_cands.begin(), zone_tracks.begin(), zone_tracks.end());  // push_front

  // Select 3 "best" tracks from all the zones
  // From src/BestTrackSelection.cc
  btrack_sel.process(extended_best_track_cands, best_tracks);

  // Insert single LCTs from station 1 as tracks
  // From src/SingleHitTracks.cc
#ifdef PHASE_TWO_TRIGGER
  // Do not make single-hit tracks
#else
  single_hit.process(conv_hits, best_tracks);
#endif

  // Construct pT address, assign pT, calculate other GMT quantities
  // From src/PtAssignment.cc
  pt_assign.process(best_tracks);

  // ___________________________________________________________________________
  // Output

  out_hits.insert(out_hits.end(), inclusive_conv_hits.begin(), inclusive_conv_hits.end());
  out_tracks.insert(out_tracks.end(), best_tracks.begin(), best_tracks.end());

  return;
}
