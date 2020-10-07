#include "L1Trigger/L1TMuonEndCap/interface/experimental/Phase2SectorProcessor.h"

#include "L1Trigger/L1TMuonEndCap/interface/TrackTools.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"


// _____________________________________________________________________________
// This implements a TEMPORARY version of the Phase 2 EMTF sector processor.
// It is supposed to be replaced in the future. It is intentionally written
// in a monolithic fashion to allow easy replacement.
//


namespace experimental {

void Phase2SectorProcessor::configure(
    // Object pointers
    const GeometryTranslator* geom,
    const ConditionHelper* cond,
    const SectorProcessorLUT* lut,
    PtAssignmentEngine* pt_assign_engine,
    // Sector processor config
    int verbose, int endcap, int sector, int bx,
    int bxShiftCSC, int bxShiftRPC, int bxShiftGEM, int bxShiftME0,
    std::string era
) {
  assert(emtf::MIN_ENDCAP <= endcap && endcap <= emtf::MAX_ENDCAP);
  assert(emtf::MIN_TRIGSECTOR <= sector && sector <= emtf::MAX_TRIGSECTOR);

  assert(geom != nullptr);
  assert(cond != nullptr);
  assert(lut  != nullptr);
  assert(pt_assign_engine != nullptr);

  geom_             = geom;
  cond_             = cond;
  lut_              = lut;
  pt_assign_engine_ = pt_assign_engine;

  verbose_    = verbose;
  endcap_     = endcap;
  sector_     = sector;
  bx_         = bx;

  bxShiftCSC_ = bxShiftCSC;
  bxShiftRPC_ = bxShiftRPC;
  bxShiftGEM_ = bxShiftGEM;
  bxShiftME0_ = bxShiftME0;

  era_        = era;
}

void Phase2SectorProcessor::process(
    // Input
    const edm::Event& iEvent, const edm::EventSetup& iSetup,
    const TriggerPrimitiveCollection& muon_primitives,
    // Output
    EMTFHitCollection& out_hits,
    EMTFTrackCollection& out_tracks
) const {

  // ___________________________________________________________________________
  // Primitive selection & primitive conversion
  // (shared with current EMTF)

  bool includeNeighbor  = true;
  bool duplicateTheta   = true;
  bool bugME11Dupes     = false;

  std::vector<int> zoneBoundaries = {0, 41, 49, 87, 127};
  int zoneOverlap       = 2;
  bool fixZonePhi       = true;
  bool useNewZones      = false;
  bool fixME11Edges     = true;

  PrimitiveSelection prim_sel;
  prim_sel.configure(
      verbose_, endcap_, sector_, bx_,
      bxShiftCSC_, bxShiftRPC_, bxShiftGEM_, bxShiftME0_,
      includeNeighbor, duplicateTheta,
      bugME11Dupes
  );

  PrimitiveConversion prim_conv;
  prim_conv.configure(
      geom_, lut_,
      verbose_, endcap_, sector_, bx_,
      bxShiftCSC_, bxShiftRPC_, bxShiftGEM_, bxShiftME0_,
      zoneBoundaries, zoneOverlap,
      duplicateTheta, fixZonePhi, useNewZones, fixME11Edges,
      bugME11Dupes
  );

  // ___________________________________________________________________________
  // Input

  EMTFHitCollection conv_hits;     // "converted" hits converted by primitive converter
  std::vector<Track> best_tracks;  // "best" tracks selected from all the zones. 'Track' is an internal class

  std::map<int, TriggerPrimitiveCollection> selected_dt_map;
  std::map<int, TriggerPrimitiveCollection> selected_csc_map;
  std::map<int, TriggerPrimitiveCollection> selected_rpc_map;
  std::map<int, TriggerPrimitiveCollection> selected_gem_map;
  std::map<int, TriggerPrimitiveCollection> selected_me0_map;
  std::map<int, TriggerPrimitiveCollection> selected_prim_map;
  std::map<int, TriggerPrimitiveCollection> inclusive_selected_prim_map;

  // Select muon primitives that belong to this sector and this BX.
  // Put them into maps with an index that roughly corresponds to
  // each input link.
  prim_sel.process(DTTag(), muon_primitives, selected_dt_map);
  prim_sel.process(CSCTag(), muon_primitives, selected_csc_map);
  prim_sel.process(RPCTag(), muon_primitives, selected_rpc_map);
  prim_sel.process(GEMTag(), muon_primitives, selected_gem_map);
  prim_sel.process(ME0Tag(), muon_primitives, selected_me0_map);
  prim_sel.merge_no_truncate(selected_dt_map, selected_csc_map, selected_rpc_map, selected_gem_map, selected_me0_map, selected_prim_map);

  // Convert trigger primitives into "converted" hits
  // A converted hit consists of integer representations of phi, theta, and zones
  prim_conv.process(selected_prim_map, conv_hits);

  {
    // Clear the input maps to save memory
    selected_dt_map.clear();
    selected_csc_map.clear();
    selected_rpc_map.clear();
    selected_gem_map.clear();
    selected_me0_map.clear();
  }

  // ___________________________________________________________________________
  // Build

  build_tracks(conv_hits, best_tracks);

  // ___________________________________________________________________________
  // Output

  EMTFTrackCollection best_emtf_tracks;
  convert_tracks(conv_hits, best_tracks, best_emtf_tracks);

  out_hits.insert(out_hits.end(), conv_hits.begin(), conv_hits.end());
  out_tracks.insert(out_tracks.end(), best_emtf_tracks.begin(), best_emtf_tracks.end());
  return;
}


// _____________________________________________________________________________
// Specific data formats
// (adapted from rootpy_trackbuilding9.py)

constexpr int NLAYERS = 16;      // 5 (CSC) + 4 (RPC) + 3 (GEM) + 4 (DT)
constexpr int NFEATURES = 36;    // NN features
constexpr int NPREDICTIONS = 9;  // NN outputs

constexpr int PATTERN_BANK_NPT = 18;   // straightness
constexpr int PATTERN_BANK_NETA = 7;   // zone
constexpr int PATTERN_BANK_NLAYERS = NLAYERS;
constexpr int PATTERN_BANK_NVARS = 3;  // min, med, max

constexpr int PATTERN_X_CENTRAL = 31;  // pattern bin number 31 is the central
constexpr int PATTERN_X_SEARCH_MIN = 33;
constexpr int PATTERN_X_SEARCH_MAX = 154-10;
//constexpr int PATTERN_X_SEARCH_MAX = 154-10+12;  // account for DT


class Hit {
public:
  // id = (type, station, ring, endsec, fr, bx)
  using hit_id_t = std::array<int32_t, 6>;
  constexpr hit_id_t id() const {
    hit_id_t ret {{type, station, ring, endsec, fr, bx}};
    return ret;
  }

  explicit Hit(int16_t vh_type, int16_t vh_station, int16_t vh_ring,
               int16_t vh_endsec, int16_t vh_fr, int16_t vh_bx,
               int32_t vh_emtf_layer, int32_t vh_emtf_phi, int32_t vh_emtf_theta,
               int32_t vh_emtf_bend, int32_t vh_emtf_qual, int32_t vh_emtf_time,
               int32_t vh_old_emtf_phi, int32_t vh_old_emtf_bend,
               int32_t vh_sim_tp, int32_t vh_ref)
  {
    type             = vh_type;
    station          = vh_station;
    ring             = vh_ring;
    endsec           = vh_endsec;
    fr               = vh_fr;
    bx               = vh_bx;
    emtf_layer       = vh_emtf_layer;
    emtf_phi         = vh_emtf_phi;
    emtf_theta       = vh_emtf_theta;
    emtf_bend        = vh_emtf_bend;
    emtf_qual        = vh_emtf_qual;
    emtf_time        = vh_emtf_time;
    old_emtf_phi     = vh_old_emtf_phi;
    old_emtf_bend    = vh_old_emtf_bend;
    sim_tp           = vh_sim_tp;
    ref              = vh_ref;
  }

  // Properties
  int16_t type;           // DT=0,CSC=1,RPC=2,GEM=3,ME0=4
  int16_t station;        // 1 to 4
  int16_t ring;           // 1 to 4
  int16_t endsec;         // 0 to 5: endcap 1 sector 1-6; 6 to 11: endcap 2 sector 1-6
  int16_t fr;             // 0: rear CSC chamber; 1: front CSC chamber
  int16_t bx;             // -3 to 3
  int32_t emtf_layer;     // 0 to 4: CSC stations; 5 to 8: RPC stations; 9 to 11: GEM stations; 12 to 15: DT stations
  int32_t emtf_phi;       // 13-bit integer (0 to 8191)
  int32_t emtf_theta;     // 7-bit integer (0 to 127)
  int32_t emtf_bend;      // 6-bit integer (0 to 63) if no DT; 10-bit integer (-512 to 511) with DT
  int32_t emtf_qual;      // 1 to 6: number of layers in CSC or ME0; includes sign: +/- for F/R
  int32_t emtf_time;      // not currently used
  int32_t old_emtf_phi;   // used only for sotfware
  int32_t old_emtf_bend;  // used only for software
  int32_t sim_tp;         // used only for software
  int32_t ref;            // used only for software
};

class Road {
public:
  using road_hits_t = std::vector<Hit>;

  // id = (endcap, sector, ipt, ieta, iphi)
  using road_id_t = std::array<int32_t, 5>;
  constexpr road_id_t id() const {
    road_id_t ret {{endcap, sector, ipt, ieta, iphi}};
    return ret;
  }

  // Provide hash function for road_id_t
  struct Hasher {
    inline std::size_t operator()(const road_id_t& road_id) const noexcept {
      int32_t endcap = road_id[0];
      int32_t sector = road_id[1];
      int32_t endsec = (endcap == 1) ? (sector - 1) : (sector - 1 + 6);

      std::size_t seed = 0;
      seed |= (road_id[4] << (0));        // allocates 256 (1<<8) entries for iphi (needs ~160)
      seed |= (road_id[3] << (0+8));      // allocates 8 (1<<3) entries for ieta (needs 7)
      seed |= (road_id[2] << (0+8+3));    // allocates 32 (1<<5) entries for ipt (needs 18)
      seed |= (endsec     << (0+8+3+5));
      return seed;
    }
  };

  explicit Road(int16_t vr_endcap, int16_t vr_sector, int16_t vr_ipt, int16_t vr_ieta, int16_t vr_iphi,
                const road_hits_t& vr_hits, int16_t vr_mode, int16_t vr_quality,
                int16_t vr_sort_code, int32_t vr_phi_median, int32_t vr_theta_median)
  {
    endcap       = vr_endcap;
    sector       = vr_sector;
    ipt          = vr_ipt;
    ieta         = vr_ieta;
    iphi         = vr_iphi;
    hits         = vr_hits;
    mode         = vr_mode;
    quality      = vr_quality;
    sort_code    = vr_sort_code;
    phi_median   = vr_phi_median;
    theta_median = vr_theta_median;
  }

  // Properties
  int16_t endcap;         // +1: positive; -1: negative
  int16_t sector;         // 1 to 6
  int16_t ipt;            // 0 to 8: prompt; 9 to 17: displaced
  int16_t ieta;           // 0 to 6: zone 0-6
  int16_t iphi;           // 0 to 159: quadstrip number
  road_hits_t hits;       // hits that belong to this road
  int16_t mode;           // 4-bit word: see create_road()
  int16_t quality;        // 0 to 9: see find_emtf_road_quality()
  int16_t sort_code;      // 10-bit word: see find_emtf_road_sort_code()
  int32_t phi_median;     // 13-bit integer (0 to 8191): median phi of the hits
  int32_t theta_median;   // 7-bit integer (0 to 127): median theta of the hits
};

class Track {
public:
  using road_hits_t = std::vector<Hit>;

  // id = (endcap, sector, ipt, ieta, iphi)
  using road_id_t = std::array<int32_t, 5>;
  constexpr road_id_t id() const {
    road_id_t ret {{endcap, sector, ipt, ieta, iphi}};
    return ret;
  }

  explicit Track(int16_t vt_endcap, int16_t vt_sector, int16_t vt_ipt, int16_t vt_ieta, int16_t vt_iphi,
                 const road_hits_t& vt_hits, int16_t vt_mode, int16_t vt_quality, int16_t vt_sort_code,
                 float vt_xml_pt, float vt_pt, int16_t vt_q, float vt_y_pred, float vt_y_discr,
                 float vt_y_displ, float vt_d0_displ, float vt_pt_displ, int32_t vt_emtf_phi, int32_t vt_emtf_theta)
  {
    endcap     = vt_endcap;
    sector     = vt_sector;
    ipt        = vt_ipt;
    ieta       = vt_ieta;
    iphi       = vt_iphi;
    hits       = vt_hits;
    mode       = vt_mode;
    quality    = vt_quality;
    sort_code  = vt_sort_code;
    xml_pt     = vt_xml_pt;
    pt         = vt_pt;
    q          = vt_q;
    y_pred     = vt_y_pred;
    y_discr    = vt_y_discr;
    y_displ    = vt_y_displ;
    d0_displ   = vt_d0_displ;
    pt_displ   = vt_pt_displ;
    emtf_phi   = vt_emtf_phi;
    emtf_theta = vt_emtf_theta;
  }

  // Properties
  int16_t endcap;         // +1: positive; -1: negative
  int16_t sector;         // 1 to 6
  int16_t ipt;            // 0 to 8: prompt; 9 to 17: displaced
  int16_t ieta;           // 0 to 6: zone 0-6
  int16_t iphi;           // 0 to 159: quadstrip number
  road_hits_t hits;       // hits that belong to this road
  int16_t mode;           // 4-bit word: see create_road()
  int16_t quality;        // 0 to 9: see find_emtf_road_quality()
  int16_t sort_code;      // 10-bit word: see find_emtf_road_sort_code()
  float   xml_pt;         // track pt, before scaling to 90% eff WP.
  float   pt;             // track pt, after scaling to 90% eff WP.
  int16_t q;              // track charge.
  float   y_pred;         // track curvature q/pt (from NN).
  float   y_discr;        // track PU discr (from NN).
  float   y_displ;        // track curvature q/pt without vertex constraint (from NN).
  float   d0_displ;       // track transverse impact parameter (from NN).
  float   pt_displ;       // track pt without vertex constraint after scaling to ?? eff WP
  int32_t emtf_phi;       // 13-bit integer (0 to 8191): median phi of the hits. Same as Road::phi_median.
  int32_t emtf_theta;     // 7-bit integer (0 to 127): median theta of the hits. Same as Road::theta_median.
};

// A 'Feature' holds 36 values
using Feature = std::array<float, NFEATURES>;

// A 'Prediction' holds 2 values
using Prediction = std::array<float, NPREDICTIONS>;


// _____________________________________________________________________________
// Specific functions
// (adapted from rootpy_trackbuilding9.py)

template<typename Container, typename Predicate>
Container my_filter(Predicate pred, const Container& input) {
  Container output;
  std::copy_if(input.cbegin(), input.cend(), std::back_inserter(output), pred);
  return output;
}

template<typename ForwardIt>
ForwardIt my_remove(ForwardIt first, ForwardIt last, const std::vector<bool>& mask) {
  if (first == last)
    return last;

  assert(std::distance(mask.cbegin(), mask.cend()) == std::distance(first, last));
  std::size_t i = 0;

  for (; first != last; ++first, ++i)
    if (mask[i] == false)  // to be removed
      break;

  if (first == last)
    return last;

  ForwardIt result = first;
  ++first;
  ++i;

  for (; first != last; ++first, ++i)
    if (mask[i] == true)  // to be kept
      *result++ = std::move(*first);
  return result;
}

template<typename Container, typename Predicate>
void my_inplace_filter(Predicate pred, Container& input) {
  std::vector<bool> mask(input.size(), false);
  std::size_t i = 0;

  for (auto it = input.cbegin(); it != input.cend(); ++it, ++i) {
    if (pred(*it)) {
      mask[i] = true;
    }
  }
  input.erase(my_remove(input.cbegin(), input.cend(), mask), input.end());
}

template<typename T>
std::vector<size_t> my_argsort(const std::vector<T>& v, bool reverse=false) {
  std::vector<size_t> indices(v.size());
  std::iota(indices.begin(), indices.end(), 0);
  auto sort_f = [&v](size_t i, size_t j) { return v[i] < v[j]; };
  std::stable_sort(indices.begin(), indices.end(), sort_f);
  if (reverse) {
    std::reverse(indices.begin(), indices.end());
  }
  return indices;
}

template<typename T>
T my_median_sorted(const std::vector<T>& vec) {
  std::size_t middle = (vec.size() == 0) ? 0 : (vec.size() - 1)/2;
  return vec[middle];
}

template<typename T>
T my_median_unsorted(std::vector<T>& vec) {
  std::size_t middle = (vec.size() == 0) ? 0 : (vec.size() - 1)/2;
  std::nth_element(vec.begin(), vec.begin() + middle, vec.end());  // input vec will be partially sorted while finding median
  return vec[middle];
}


// _____________________________________________________________________________
// Specific modules
// (adapted from rootpy_trackbuilding9.py)

#include "utility.icc"

class Utility {
public:
  // Constructor
  constexpr Utility() {
    // Initialize 3-D array
    for (size_t i=0; i<find_emtf_layer_lut.size(); i++) {
      for (size_t j=0; j<find_emtf_layer_lut[i].size(); j++) {
        for (size_t k=0; k<find_emtf_layer_lut[i][j].size(); k++) {
          find_emtf_layer_lut[i][j][k] = std::move(_find_emtf_layer_lut[i][j][k]);
        }
      }
    }

    // Initialize 5-D array
    for (size_t i=0; i<find_emtf_zones_lut.size(); i++) {
      for (size_t j=0; j<find_emtf_zones_lut[i].size(); j++) {
        for (size_t k=0; k<find_emtf_zones_lut[i][j].size(); k++) {
          for (size_t l=0; l<find_emtf_zones_lut[i][j][k].size(); l++) {
            for (size_t m=0; m<find_emtf_zones_lut[i][j][k][l].size(); m++) {
              find_emtf_zones_lut[i][j][k][l][m] = std::move(_find_emtf_zones_lut[i][j][k][l][m]);
            }
          }
        }
      }
    }
  }  // end constructor

  bool isFront_detail(int subsystem, int station, int ring, int chamber, int subsector) const {
    bool result = false;

    if (subsystem == TriggerPrimitive::kCSC) {
      bool isOverlapping = !(station == 1 && ring == 3);
      // not overlapping means back
      if (isOverlapping)
      {
        bool isEven = (chamber % 2 == 0);
        // odd chambers are bolted to the iron, which faces
        // forward in 1&2, backward in 3&4, so...
        result = (station < 3) ? isEven : !isEven;
      }
    } else if (subsystem == TriggerPrimitive::kRPC) {
      //// 10 degree rings have even subsectors in front
      //// 20 degree rings have odd subsectors in front
      //bool is_10degree = !((station == 3 || station == 4) && (ring == 1));
      //bool isEven = (subsector % 2 == 0);
      //result = (is_10degree) ? isEven : !isEven;

      // Use the equivalent CSC chamber F/R
      bool isEven = (chamber % 2 == 0);
      result = (station < 3) ? isEven : !isEven;
    } else if (subsystem == TriggerPrimitive::kGEM) {
      //
      result = (chamber % 2 == 0);
    } else if (subsystem == TriggerPrimitive::kME0) {
      //
      result = (chamber % 2 == 0);
    } else if (subsystem == TriggerPrimitive::kDT) {
      //
      result = (chamber % 2 == 0);
    }
    return result;
  }

  bool find_fr(const EMTFHit& conv_hit) const {
    return isFront_detail(conv_hit.Subsystem(), conv_hit.Station(), conv_hit.Ring(), conv_hit.Chamber(),
                          (conv_hit.Subsystem() == TriggerPrimitive::kRPC ? conv_hit.Subsector_RPC() : conv_hit.Subsector()));
  }

  int32_t find_endsec(int32_t endcap, int32_t sector) const {
    return (endcap == 1) ? (sector - 1) : (sector - 1 + 6);
  }

  int32_t find_endsec(const EMTFHit& conv_hit) const {
    int32_t endcap     = conv_hit.Endcap();
    int32_t sector     = conv_hit.PC_sector();
    return find_endsec(endcap, sector);
  }

  // A coarse-graining operation
  // divide by 'quadstrip' unit (4 * 8), and adjust for rounding
  int32_t find_pattern_x(int32_t emtf_phi) const {
    return (emtf_phi+16)/32;
  }

  // Undo the coarse-graining operation
  // multiply by 'quadstrip' unit (4 * 8)
  int32_t find_pattern_x_inverse(int32_t x) const {
    return (x*32);
  }

  // Calculate transverse impact parameter, d0
  double calculate_d0(double invPt, double phi, double xv, double yv, double B=3.811) const {
    double _invPt = (std::abs(invPt) < 1./10000) ? (invPt < 0 ? -1./10000 : +1./10000) : invPt;
    double _R = -1.0 / (0.003 * B * _invPt);                          // R = -pT/(0.003 q B)  [cm]
    double _xc = xv - (_R * std::sin(phi));                           // xc = xv - R sin(phi)
    double _yc = yv + (_R * std::cos(phi));                           // yc = yv + R cos(phi)
    double _d0 = _R - ((_R < 0 ? -1. : +1.) * std::hypot(_xc, _yc));  // d0 = R - sign(R) * sqrt(xc^2 + yc^2)
    return _d0;
  }

  // Decide EMTF hit layer number
  int32_t find_emtf_layer(const EMTFHit& conv_hit) const {
    int32_t type       = conv_hit.Subsystem();
    int32_t station    = conv_hit.Station();
    int32_t ring       = conv_hit.Ring();

    int32_t emtf_layer = find_emtf_layer_lut[type][station][ring];
    return emtf_layer;
  }

  // Decide EMTF hit zones
  std::vector<int32_t> find_emtf_zones(const EMTFHit& conv_hit) const {
    std::vector<int32_t> zones;

    int32_t emtf_theta = conv_hit.Theta_fp();
    int32_t type       = conv_hit.Subsystem();
    int32_t station    = conv_hit.Station();
    int32_t ring       = conv_hit.Ring();

    for (size_t zone=0; zone<find_emtf_zones_lut[type][station][ring].size(); zone++) {
      int32_t low  = find_emtf_zones_lut[type][station][ring][zone][0];
      int32_t high = find_emtf_zones_lut[type][station][ring][zone][1];
      if ((low <= emtf_theta) && (emtf_theta <= high)) {
        zones.push_back(zone);
      }
    }
    return zones;
  }

  std::vector<int32_t> find_emtf_zones(const Hit& hit) const {
    std::vector<int32_t> zones;

    int32_t emtf_theta = hit.emtf_theta;
    int32_t type       = hit.type;
    int32_t station    = hit.station;
    int32_t ring       = hit.ring;

    for (size_t zone=0; zone<find_emtf_zones_lut[type][station][ring].size(); zone++) {
      int32_t low  = find_emtf_zones_lut[type][station][ring][zone][0];
      int32_t high = find_emtf_zones_lut[type][station][ring][zone][1];
      if ((low <= emtf_theta) && (emtf_theta <= high)) {
        zones.push_back(zone);
      }
    }
    return zones;
  }

  // Decide EMTF hit bend
  int32_t find_emtf_bend(const EMTFHit& conv_hit) const {
    int32_t emtf_bend  = conv_hit.Bend();
    int32_t type       = conv_hit.Subsystem();
    int32_t station    = conv_hit.Station();
    int32_t ring       = conv_hit.Ring();
    int32_t endcap     = conv_hit.Endcap();
    int32_t quality    = conv_hit.Quality();

    if (type == TriggerPrimitive::kCSC) {
      // Special case for ME1/1a
      // rescale the bend to the same scale as ME1/1b
      if ((station == 1) && (ring == 4)) {
        emtf_bend = static_cast<int32_t>(std::round(static_cast<float>(emtf_bend) * 0.026331/0.014264));
        emtf_bend = std::clamp(emtf_bend, -32, 31);
      }
      emtf_bend *= endcap;
      emtf_bend = static_cast<int32_t>(std::round(static_cast<float>(emtf_bend) * 0.5));  // from 1/32-strip unit to 1/16-strip unit
      emtf_bend = std::clamp(emtf_bend, -16, 15);

    } else if (type == TriggerPrimitive::kME0) {
      emtf_bend = static_cast<int32_t>(std::round(static_cast<float>(emtf_bend) * 0.5));  // from 1/4-strip unit to 1/2-strip unit
      emtf_bend = std::clamp(emtf_bend, -64, 63);

    } else if (type == TriggerPrimitive::kDT) {
      if (quality >= 4) {
        emtf_bend = std::clamp(emtf_bend, -512, 511);
      } else {
        //emtf_bend = 0;
        emtf_bend = std::clamp(emtf_bend, -512, 511);
      }

    } else {  // (type == TriggerPrimitive::kRPC) || (type == TriggerPrimitive::kGEM)
      emtf_bend = 0;
    }
    return emtf_bend;
  }

  // Decide EMTF hit bend (old version)
  int32_t find_emtf_old_bend(const EMTFHit& conv_hit) const {
    static const int32_t lut[11] = {5, -5, 4, -4, 3, -3, 2, -2, 1, -1, 0};

    int32_t emtf_bend  = conv_hit.Bend();
    int32_t type       = conv_hit.Subsystem();
    int32_t endcap     = conv_hit.Endcap();
    int32_t pattern    = conv_hit.Pattern();  // CLCT pattern

    if (type == TriggerPrimitive::kCSC) {
      assert(0 <= pattern && pattern <= 10);
      emtf_bend = lut[pattern];
      emtf_bend *= endcap;

    } else if (type == TriggerPrimitive::kME0) {
      // do nothing

    } else if (type == TriggerPrimitive::kDT) {
      // do nothing

    } else {  // (type == TriggerPrimitive::kRPC) || (type == TriggerPrimitive::kGEM)
      emtf_bend = 0;
    }
    return emtf_bend;
  }

  // Decide EMTF hit phi (integer unit)
  int32_t find_emtf_phi(const EMTFHit& conv_hit) const {
    int32_t emtf_phi   = conv_hit.Phi_fp();
    int32_t type       = conv_hit.Subsystem();
    int32_t station    = conv_hit.Station();
    int32_t ring       = conv_hit.Ring();
    int32_t endcap     = conv_hit.Endcap();
    int32_t bend       = conv_hit.Bend();
    int32_t fr         = find_fr(conv_hit);

    if (type == TriggerPrimitive::kCSC) {
      if (station == 1) {
        float bend_corr = 0.;
        if (ring == 1) {
          bend_corr = ((static_cast<float>(1-fr) * -2.0832) + (static_cast<float>(fr) * 2.0497));  // ME1/1b (r,f)
        } else if (ring == 4) {
          bend_corr = ((static_cast<float>(1-fr) * -2.4640) + (static_cast<float>(fr) * 2.3886));  // ME1/1a (r,f)
        } else if (ring == 2) {
          bend_corr = ((static_cast<float>(1-fr) * -1.3774) + (static_cast<float>(fr) * 1.2447));  // ME1/2 (r,f)
        } else {
          bend_corr = 0.;  // ME1/3 (r,f): no correction
        }
        bend_corr *= bend;
        bend_corr *= endcap;
        emtf_phi += static_cast<int32_t>(std::round(bend_corr));
      } else {
        // do nothing
      }
    } else {
      // do nothing
    }
    return emtf_phi;
  }

  // Decide EMTF hit phi (integer unit) (old version)
  int32_t find_emtf_old_phi(const EMTFHit& conv_hit) const {
    static const int32_t ph_pattern_corr_lut[11] = {0, 0, 5, 5, 5, 5, 2, 2, 2, 2, 0};
    static const int32_t ph_pattern_corr_sign_lut[11] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0};
    static const int32_t ph_init_lut[2*6*61] = {
        1324,1926,2527,1332,1932,2532,1415,2015,2615,1332,1932,2532,726,732,815,732,3128,3725,4325,3132,3732,4332,3214,3815,4415,3131,3732,4331,1316,2516,
        3716,1332,1932,2532,3132,3732,4332,116,732,2580,3781,4980,1964,2564,3164,3764,4364,4964,1380,1364,2580,3780,4980,1964,2564,3164,3764,4364,4964,1380,
        1364,1326,1923,2527,1332,1932,2532,1415,2015,2614,1331,1931,2531,725,732,815,731,3123,3724,4326,3132,3732,4332,3214,3814,4415,3131,3731,4331,1316,
        2516,3715,1332,1932,2532,3132,3731,4331,116,732,2580,3780,4980,1964,2564,3164,3763,4363,4964,1380,1364,2580,3780,4979,1964,2564,3164,3763,4363,4963,
        1380,1363,1323,1926,2525,1331,1932,2532,1415,2015,2615,1331,1931,2531,726,732,816,732,3124,3727,4325,3132,3732,4332,3215,3815,4415,3131,3731,4332,
        1316,2516,3716,1332,1932,2532,3131,3731,4332,116,732,2580,3780,4980,1964,2564,3164,3764,4364,4964,1381,1364,2580,3780,4980,1964,2564,3164,3764,4364,
        4964,1380,1364,1324,1929,2531,1332,1932,2532,1416,2015,2615,1331,1932,2531,725,732,815,732,3123,3728,4327,3132,3732,4332,3215,3815,4416,3132,3731,
        4332,1316,2516,3716,1332,1932,2532,3132,3733,4332,116,732,2580,3781,4980,1964,2564,3165,3765,4365,4964,1380,1364,2580,3781,4981,1964,2564,3164,3765,
        4365,4965,1380,1364,1325,1925,2524,1332,1932,2532,1415,2015,2615,1331,1931,2531,727,732,815,731,3124,3726,4325,3132,3732,4332,3215,3815,4415,3132,
        3731,4331,1316,2516,3716,1332,1932,2532,3132,3732,4332,116,732,2580,3780,4980,1964,2564,3164,3764,4364,4964,1380,1364,2580,3780,4980,1964,2564,3164,
        3764,4364,4964,1380,1364,1321,1927,2524,1332,1932,2532,1415,2015,2615,1331,1931,2532,725,732,815,731,3128,3727,4326,3133,3732,4332,3215,3815,4415,
        3131,3731,4332,1316,2516,3716,1332,1932,2532,3132,3732,4332,116,732,2580,3780,4980,1964,2564,3164,3764,4364,4964,1380,1364,2580,3780,4980,1964,2564,
        3164,3764,4364,4964,1380,1364,1979,2578,3178,1972,2572,3172,1890,2489,3090,1973,2573,3173,1380,1372,1289,1375,3779,4380,4978,3772,4372,4972,3689,4289,
        4889,3772,4373,4972,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,2524,3724,1340,
        1940,2540,3140,3740,4340,124,740,1979,2578,3179,1972,2572,3172,1889,2489,3089,1973,2573,3173,1378,1372,1289,1372,3778,4380,4982,3772,4372,4972,3689,
        4289,4890,3773,4373,4973,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,2524,3724,
        1340,1940,2540,3140,3740,4340,124,740,1977,2580,3179,1972,2572,3172,1889,2489,3089,1975,2572,3173,1382,1372,1289,1372,3779,4379,4979,3772,4372,4972,
        3688,4289,4889,3773,4373,4973,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,2524,
        3724,1340,1940,2540,3140,3740,4340,124,740,1979,2577,3180,1972,2572,3172,1889,2489,3089,1973,2573,3173,1379,1372,1289,1373,3780,4378,4979,3772,4372,
        4972,3689,4289,4889,3773,4373,4973,2588,3788,4988,1972,2572,3172,3772,4372,4972,1388,1372,1324,2524,3724,1340,1940,2540,3140,3740,4340,124,740,1324,
        2524,3724,1340,1940,2541,3141,3741,4341,124,740,1978,2580,3179,1972,2572,3172,1889,2489,3089,1973,2573,3173,1380,1372,1290,1373,3780,4378,4981,3772,
        4372,4972,3689,4290,4889,3775,4372,4977,2588,3787,4987,1972,2572,3172,3772,4372,4972,1388,1372,1324,2523,3723,1340,1940,2540,3140,3740,4340,124,740,
        1324,2523,3724,1341,1941,2541,3141,3741,4341,124,741,1979,2581,3178,1973,2573,3173,1890,2490,3090,1973,2575,3173,1382,1373,1290,1377,3779,4380,4981,
        3773,4373,4973,3690,4290,4890,3774,4374,4976,2589,3789,4989,1973,2573,3173,3773,4373,4973,1388,1373,1325,2525,3725,1341,1941,2541,3141,3741,4341,124,
        741,1325,2525,3725,1341,1941,2541,3141,3741,4341,124,742
    };

    int32_t emtf_phi   = conv_hit.Phi_fp();
    int32_t type       = conv_hit.Subsystem();

    if (type == TriggerPrimitive::kCSC) {
      const bool is_neighbor = conv_hit.Neighbor();
      const int fw_endcap = (conv_hit.Endcap() == 1) ? 0 : 1;
      const int fw_sector = conv_hit.PC_sector()-1;
      const int fw_station = (conv_hit.Station() == 1) ? (is_neighbor ? 0 : (conv_hit.Subsector()-1)) : conv_hit.Station();
      const int fw_cscid   = (conv_hit.CSC_nID()-1);
      const int fw_strip = conv_hit.Strip();
      const int pc_station = conv_hit.PC_station();
      const int pc_chamber = conv_hit.PC_chamber();
      const bool is_me11a = (conv_hit.Station() == 1 && conv_hit.Ring() == 4);
      const bool is_me11b = (conv_hit.Station() == 1 && conv_hit.Ring() == 1);
      const bool is_me13  = (conv_hit.Station() == 1 && conv_hit.Ring() == 3);

      // Is this chamber mounted in reverse direction?
      // (i.e., phi vs. strip number is reversed)
      bool ph_reverse = false;
      if ((fw_endcap == 0 && fw_station >= 3) || (fw_endcap == 1 && fw_station < 3))  // ME+3, ME+4, ME-1, ME-2
        ph_reverse = true;

      // Is this 10-deg or 20-deg chamber?
      bool is_10degree = false;
      if (
          (fw_station <= 1) || // ME1
          (fw_station >= 2 && ((fw_cscid >= 3 && fw_cscid <= 8) || fw_cscid == 10))  // ME2,3,4/2
      ) {
        is_10degree = true;
      }

      // LUT index
      // There are 54 CSC chambers including the neighbors in a sector, but 61 LUT indices
      // This comes from dividing the 6 chambers + 1 neighbor in ME1/1 into ME1/1a and ME1/1b
      int pc_lut_id = pc_chamber;
      if (pc_station == 0) {         // ME1 sub 1: 0 - 11
        pc_lut_id = is_me11a ? pc_lut_id + 9 : pc_lut_id;
      } else if (pc_station == 1) {  // ME1 sub 2: 16 - 27
        pc_lut_id += 16;
        pc_lut_id = is_me11a ? pc_lut_id + 9 : pc_lut_id;
      } else if (pc_station == 2) {  // ME2: 28 - 36
        pc_lut_id += 28;
      } else if (pc_station == 3) {  // ME3: 39 - 47
        pc_lut_id += 39;
      } else if (pc_station == 4) {  // ME4 : 50 - 58
        pc_lut_id += 50;
      } else if (pc_station == 5 && pc_chamber < 3) {  // neighbor ME1: 12 - 15
        pc_lut_id = is_me11a ? pc_lut_id + 15 : pc_lut_id + 12;
      } else if (pc_station == 5 && pc_chamber < 5) {  // neighbor ME2: 37 - 38
        pc_lut_id += 28 + 9 - 3;
      } else if (pc_station == 5 && pc_chamber < 7) {  // neighbor ME3: 48 - 49
        pc_lut_id += 39 + 9 - 5;
      } else if (pc_station == 5 && pc_chamber < 9) {  // neighbor ME4: 59 - 60
        pc_lut_id += 50 + 9 - 7;
      }
      assert(pc_lut_id < 61);

      // Convert half-strip into 1/8-strip
      int eighth_strip = 0;

      // Apply phi correction from CLCT pattern number (from src/SectorProcessorLUT.cc)
      int clct_pat_corr = ph_pattern_corr_lut[conv_hit.Pattern()];
      int clct_pat_corr_sign = (ph_pattern_corr_sign_lut[conv_hit.Pattern()] == 0) ? 1 : -1;
      if (fw_strip == 0 && clct_pat_corr_sign == -1)
        clct_pat_corr = 0;

      if (is_10degree) {
        eighth_strip = fw_strip << 2;  // full precision, uses only 2 bits of pattern correction
        eighth_strip += clct_pat_corr_sign * (clct_pat_corr >> 1);
      } else {
        eighth_strip = fw_strip << 3;  // multiply by 2, uses all 3 bits of pattern correction
        eighth_strip += clct_pat_corr_sign * (clct_pat_corr >> 0);
      }
      assert(eighth_strip >= 0);

      // Multiplicative factor for eighth_strip
      int factor = 1024;
      if (is_me11a)
        factor = 1707;  // ME1/1a
      else if (is_me11b)
        factor = 1301;  // ME1/1b
      else if (is_me13)
        factor = 947;   // ME1/3

      int ph_tmp = (eighth_strip * factor) >> 10;
      int ph_tmp_sign = (ph_reverse == 0) ? 1 : -1;

      int endsec_pc_lut_id = (fw_endcap * 6 + fw_sector) * 61 + pc_lut_id;
      int fph = ph_init_lut[endsec_pc_lut_id];
      fph = fph + ph_tmp_sign * ph_tmp;
      assert(0 <= fph && fph < 5000);
      //
      emtf_phi = fph;

    } else {
      // do nothing
    }

    return emtf_phi;
  }

  // Decide EMTF hit theta (integer unit)
  int32_t find_emtf_theta(const EMTFHit& conv_hit) const {
    int32_t emtf_theta = conv_hit.Theta_fp();
    int32_t type       = conv_hit.Subsystem();
    int32_t station    = conv_hit.Station();
    int32_t wire       = conv_hit.Wire();
    int32_t quality    = conv_hit.Quality();

    if (type == TriggerPrimitive::kDT) {
      // wire -1 means no theta SL
      // quality 0&1 are RPC digis
      if ((wire == -1) || (quality < 2)) {
        if (station == 1) {
          emtf_theta = 112;
        } else if (station == 2) {
          emtf_theta = 122;
        } else if (station == 3) {
          emtf_theta = 131;
        }
      } else {
        // do nothing
      }
    } else {
      // do nothing
    }
    return emtf_theta;
  }

  // Decide EMTF hit z-position (floating-point)
  // Not implemented
  float   find_emtf_zee(const EMTFHit& conv_hit) const { return 0.; }

  // Decide EMTF hit quality
  int32_t find_emtf_qual(const EMTFHit& conv_hit) const {
    int32_t emtf_qual  = conv_hit.Quality();
    int32_t type       = conv_hit.Subsystem();

    int32_t fr         = find_fr(conv_hit);

    if ((type == TriggerPrimitive::kCSC) || (type == TriggerPrimitive::kME0)) {
      // front chamber  -> +1
      // rear chamber   -> -1
      if (fr == 1) {
        emtf_qual *= +1;
      } else {
        emtf_qual *= -1;
      }
    } else if ((type == TriggerPrimitive::kRPC) || (type == TriggerPrimitive::kGEM)) {
      emtf_qual = 0;
    } else {  // type == TriggerPrimitive::kDT
      // do nothing
    }
    return emtf_qual;
  }

  // Decide EMTF hit time (integer unit)
  int32_t find_emtf_time(const EMTFHit& conv_hit) const {
    //int32_t emtf_time  = static_cast<int32_t>(std::round(conv_hit.Time() * 16./25));  // integer unit is 25ns/16 (4-bit)
    int32_t emtf_time  = conv_hit.BX();
    return emtf_time;
  }

  // Decide EMTF road quality (by pattern straightness)
  int32_t find_emtf_road_quality(int32_t ipt) const {
    // First 9 patterns for prompt muons  : -1/2 <= q/pT <= +1/2
    // Next 9 patterns for displaced muons: -1/14 <= q/pT <= +1/14, -120 <= d0 <= 120
    // Total is 18 patterns.
    // ipt   0  1  2  3  4  5  6  7  8
    // strg  5  6  7  8  9  8  7  6  5
    // ipt   9 10 11 12 13 14 15 16 17
    // strg  0  1  2  3  4  3  2  1  0
    static const int32_t lut[PATTERN_BANK_NPT] = {5,6,7,8,9,8,7,6,5,0,1,2,3,4,3,2,1,0};
    assert((0 <= ipt) && (ipt < PATTERN_BANK_NPT));
    return lut[ipt];
  }

  // Decide EMTF road sort code (by hit composition)
  int32_t find_emtf_road_sort_code(int32_t road_quality, const std::vector<int32_t>& road_hits_layers) const {
    // 12   11     10     9      8      7      6      5      4      3..0
    //      ME1/1  ME1/2  ME2    ME3    ME4                         qual
    //                                                RE1&2  RE3&4
    // ME0                                     GE1/1  GE2/1
    // MB1  MB2                                MB3&4
    static const int32_t lut[NLAYERS] = {11,10,9,8,7,5,5,4,4,6,5,12,12,11,6,6};

    int32_t sort_code = 0;
    for (const auto& hit_lay : road_hits_layers) {
      int32_t mlayer = lut[hit_lay];
      sort_code |= (1 << mlayer);
    }
    assert((0 <= road_quality) && (road_quality < 16));
    sort_code |= road_quality;
    return sort_code;
  }

  // Decide EMTF road mode
  int32_t find_emtf_road_mode(const std::vector<Hit>& road_hits) const {
    // 'road_mode' is a 4-bit word where each bit indicates whether a hit was
    // found in one of the 4 stations
    // |bit| 3 | 2 | 1 | 0 |
    // |---|---|---|---|---|
    // |st | 1 | 2 | 3 | 4 |
    int32_t road_mode = 0;
    for (const auto& hit : road_hits) {
      int32_t station = hit.station;
      road_mode |= (1 << (4 - station));
    }
    return road_mode;
  }

  // Decide EMTF road accept
  bool find_emtf_road_accept(int32_t road_zone, const std::vector<Hit>& road_hits) const {
    // The other road modes are used to add specific rules for different zones.
    int road_mode          = 0;
    int road_mode_csc      = 0;
    int road_mode_me0      = 0;  // zones 0,1
    int road_mode_me12     = 0;  // zone 4
    int road_mode_csc_me12 = 0;  // zone 4
    int road_mode_mb1      = 0;  // zone 6
    int road_mode_mb2      = 0;  // zone 6
    int road_mode_me13     = 0;  // zone 6
    //int road_mode_me22     = 0;  // zone 6

    for (const auto& hit : road_hits) {
      int32_t type    = hit.type;
      int32_t station = hit.station;
      int32_t ring    = hit.ring;
      int32_t bx      = hit.bx;
      road_mode |= (1 << (4 - station));

      if ((type == TriggerPrimitive::kCSC) || (type == TriggerPrimitive::kME0)) {
        road_mode_csc |= (1 << (4 - station));
      }

      if ((type == TriggerPrimitive::kME0) && (bx == 0)) {
        road_mode_me0 |= (1 << 1);
      } else if ((type == TriggerPrimitive::kCSC) && (station == 1) && ((ring == 1) || (ring == 4)) && (bx == 0)) {
        road_mode_me0 |= (1 << 0);
      }

      if ((type == TriggerPrimitive::kCSC) && (station == 1) && ((ring == 2) || (ring == 3))) {  // pretend as station 2
        road_mode_me12 |= (1 << (4 - 2));
      } else if ((type == TriggerPrimitive::kRPC) && (station == 1) && ((ring == 2) || (ring == 3))) {  // pretend as station 2
        road_mode_me12 |= (1 << (4 - 2));
      } else {
        road_mode_me12 |= (1 << (4 - station));
      }

      if ((type == TriggerPrimitive::kCSC) && (station == 1) && ((ring == 2) || (ring == 3))) {  // pretend as station 2
        road_mode_csc_me12 |= (1 << (4 - 2));
      } else if (type == TriggerPrimitive::kCSC) {
        road_mode_csc_me12 |= (1 << (4 - station));
      }

      if ((type == TriggerPrimitive::kDT) && (station == 1)) {
        road_mode_mb1 |= (1 << 1);
      } else if ((type == TriggerPrimitive::kDT) && (station >= 2)) {
        road_mode_mb1 |= (1 << 0);
      } else if ((type == TriggerPrimitive::kCSC) && (station >= 1) && ((ring == 2) || (ring == 3))) {
        road_mode_mb1 |= (1 << 0);
      } else if ((type == TriggerPrimitive::kRPC) && (station >= 1) && ((ring == 2) || (ring == 3))) {
        road_mode_mb1 |= (1 << 0);
      }

      if ((type == TriggerPrimitive::kDT) && (station == 2)) {
        road_mode_mb2 |= (1 << 1);
      } else if ((type == TriggerPrimitive::kDT) && (station >= 3)) {
        road_mode_mb2 |= (1 << 0);
      } else if ((type == TriggerPrimitive::kCSC) && (station >= 1) && ((ring == 2) || (ring == 3))) {
        road_mode_mb2 |= (1 << 0);
      } else if ((type == TriggerPrimitive::kRPC) && (station >= 1) && ((ring == 2) || (ring == 3))) {
        road_mode_mb2 |= (1 << 0);
      }

      if ((type == TriggerPrimitive::kCSC) && (station == 1) && ((ring == 2) || (ring == 3))) {
        road_mode_me13 |= (1 << 1);
      } else if ((type == TriggerPrimitive::kCSC) && (station >= 2) && ((ring == 2) || (ring == 3))) {
        road_mode_me13 |= (1 << 0);
      } else if ((type == TriggerPrimitive::kRPC) && (station == 1) && ((ring == 2) || (ring == 3))) {
        road_mode_me13 |= (1 << 1);
      } else if ((type == TriggerPrimitive::kRPC) && (station >= 2) && ((ring == 2) || (ring == 3))) {
        road_mode_me13 |= (1 << 0);
      }

      //if ((type == TriggerPrimitive::kCSC) && (station == 2) && ((ring == 2) || (ring == 3))) {
      //  road_mode_me22 |= (1 << 1);
      //} else if ((type == TriggerPrimitive::kCSC) && (station >= 3) && ((ring == 2) || (ring == 3))) {
      //  road_mode_me22 |= (1 << 0);
      //} else if ((type == TriggerPrimitive::kRPC) && (station == 2) && ((ring == 2) || (ring == 3))) {
      //  road_mode_me22 |= (1 << 1);
      //} else if ((type == TriggerPrimitive::kRPC) && (station >= 3) && ((ring == 2) || (ring == 3))) {
      //  road_mode_me22 |= (1 << 0);
      //}
    }  // end loop over road_hits

    // Apply SingleMu requirement
    // + (zones 0,1) any road with ME0 and ME1
    // + (zone 4) any road with ME1/1, ME1/2 + one more station
    // + (zone 5) any road with 2 stations
    // + (zone 6) any road with MB1+MB2, MB1+MB3, MB1+ME1/3, MB1+ME2/2, MB2+MB3, MB2+ME1/3, MB2+ME2/2, ME1/3+ME2/2
    bool accept = ((is_emtf_singlemu(road_mode) && is_emtf_muopen(road_mode_csc)) ||
        (((road_zone == 0) || (road_zone == 1)) && (road_mode_me0 == 3)) ||
        ((road_zone == 4) && is_emtf_singlemu(road_mode_me12) && is_emtf_muopen(road_mode_csc_me12)) ||
        ((road_zone == 5) && is_emtf_doublemu(road_mode) && is_emtf_muopen(road_mode_csc)) ||
        ((road_zone == 6) && ((road_mode_mb1 == 3) || (road_mode_mb2 == 3) || (road_mode_me13 == 3))) );
    return accept;
  }

  bool is_emtf_singlemu(int mode) const {
    static const std::set<int> s {11,13,14,15};
    return (s.find(mode) != s.end());  // s.contains(mode);
  }

  bool is_emtf_doublemu(int mode) const {
    //static const std::set<int> s {7,10,12,11,13,14,15};
    static const std::set<int> s {9,10,12,11,13,14,15};  // replace 2-3-4 with 1-4
    return (s.find(mode) != s.end());  // s.contains(mode);
  }

  bool is_emtf_muopen(int mode) const {
    //static const std::set<int> s {3,5,6,9,7,10,12,11,13,14,15};
    static const std::set<int> s {5,6,9,7,10,12,11,13,14,15};  // remove 3-4
    return (s.find(mode) != s.end());  // s.contains(mode);
  }

  bool is_emtf_singlehit(int mode) const {
    return bool(mode & (1 << 3));
  }

  bool is_emtf_singlehit_me2(int mode) const {
    return bool(mode & (1 << 2));
  }

  // For now, only consider BX=0
  bool is_emtf_legit_hit_check_bx(const EMTFHit& conv_hit) const {
    int32_t type       = conv_hit.Subsystem();
    int32_t bx         = conv_hit.BX();

    if (type == TriggerPrimitive::kCSC) {
      return (bx == -1) || (bx == 0);
    } else if (type == TriggerPrimitive::kDT) {
      return (bx == -1) || (bx == 0);
    }
    return (bx == 0);
  }

  bool is_emtf_legit_hit_check_phi(const EMTFHit& conv_hit) const {
    int32_t type       = conv_hit.Subsystem();
    int32_t emtf_phi   = conv_hit.Phi_fp();

    if (type == TriggerPrimitive::kME0) {
      return (emtf_phi > 0);
    } else if (type == TriggerPrimitive::kDT) {
      return (emtf_phi > 0);
    }
    return true;
  }

  bool is_emtf_legit_hit(const EMTFHit& conv_hit) const {
    return is_emtf_legit_hit_check_bx(conv_hit) && is_emtf_legit_hit_check_phi(conv_hit);
  }

  int find_pt_bin(float x) const {
    static const std::vector<float> v = {-0.49376795, -0.38895044, -0.288812, -0.19121648, -0.0810074, 0.0810074, 0.19121648, 0.288812, 0.38895044, 0.49376795};  // bin edges

    x = std::clamp(x, v.front(), static_cast<float>(v.back() - 1e-5));
    unsigned ind = std::upper_bound(v.begin(), v.end(), x) - v.begin() - 1;
    assert(ind < v.size());
    return ind;
  }

  int find_eta_bin(float x) const {
    static const std::vector<float> v = {0.8, 1.24, 1.56, 1.7, 1.8, 1.98, 2.16, 2.4};  // bin edges

    x = std::abs(x);  // abs(eta)
    x = std::clamp(x, v.front(), static_cast<float>(v.back() - 1e-5));
    unsigned ind = std::upper_bound(v.begin(), v.end(), x) - v.begin() - 1;
    assert(ind < v.size());
    ind = (v.size()-1) - ind;  // zone 0 starts at highest eta
    return ind;
  }

private:
  // 3-D array of size [# types][# stations][# rings]
  using lut_5_5_5_t = std::array<std::array<std::array<int32_t, 5>, 5>, 5>;
  lut_5_5_5_t find_emtf_layer_lut {};

  // 5-D array of size [# types][# stations][# rings][# zones][low, high]
  using lut_5_5_5_7_2_t = std::array<std::array<std::array<std::array<std::array<int32_t, 2>, 7>, 5>, 5>, 5>;
  lut_5_5_5_7_2_t find_emtf_zones_lut {};
};

constexpr Utility util;

#include "patternbank.icc"

class PatternBank {
public:
  // Constructor
  constexpr PatternBank() {
    // Initialize 4-D array
    for (size_t i=0; i<x_array.size(); i++) {
      for (size_t j=0; j<x_array[i].size(); j++) {
        for (size_t k=0; k<x_array[i][j].size(); k++) {
          for (size_t l=0; l<x_array[i][j][k].size(); l++) {
            x_array[i][j][k][l] = std::move(_patternbank[i][j][k][l]);
          }
        }
      }
    }
  }  // end constructor

  // 4-D array of size [NLAYERS][NETA][NVARS][NPT]
  // Note: rearranged for cache-friendliness. In the original python script,
  // it's arranged as [NPT][NETA][NLAYERS][NVARS]
  using patternbank_t = std::array<std::array<std::array<std::array<int32_t, PATTERN_BANK_NPT>,
      PATTERN_BANK_NVARS>, PATTERN_BANK_NETA>, PATTERN_BANK_NLAYERS>;

  patternbank_t x_array {};
};

constexpr PatternBank bank;


// _____________________________________________________________________________
// PatternRecognition class matches hits to pre-defined patterns.
// Before the pattern matching, it also converts the EMTFHitCollection into a
// vector<Hit>, where Hit is a simple data struct. A set of 18 patterns is
// used for each zone and for each 'quadstrip'. The pattern matching is done by
// comparing the phi value of each hit to the window encoded for the station of
// the hit in a pattern. When a pattern fires, a road is produced. A road
// contains information about the pattern that fires and the hits that
// belong to the road. The output of this class is a vector<Road>, which
// contains all the roads.
// In this C++ version, the pattern matching is done with some trick to speed
// up software processing. It is not the logic meant to be implemented in
// firmware, but it should give the same results.

class PatternRecognition {
public:
  void run(int32_t endcap, int32_t sector, const EMTFHitCollection& conv_hits,
           std::vector<Hit>& sector_hits, std::vector<Road>& sector_roads) const {

    // Optimize for CPU processing?
    bool optimize_for_cpu = true;

    // Use endcap = +1 or -1
    if (endcap == 2)
      endcap = -1;

    // Convert all the hits again and apply the filter to get the legit hits
    int32_t sector_mode = 0;

    for (size_t ihit = 0; ihit < conv_hits.size(); ++ihit) {
      const EMTFHit& conv_hit = conv_hits.at(ihit);

      int32_t dummy_sim_tp = -1;

      if (util.is_emtf_legit_hit(conv_hit)) {
        //Hit(int16_t vh_type, int16_t vh_station, int16_t vh_ring,
        //    int16_t vh_endsec, int16_t vh_fr, int16_t vh_bx,
        //    int32_t vh_emtf_layer, int32_t vh_emtf_phi, int32_t vh_emtf_theta,
        //    int32_t vh_emtf_bend, int32_t vh_emtf_qual, int32_t vh_emtf_time,
        //    int32_t vh_old_emtf_phi, int32_t vh_old_emtf_bend,
        //    int32_t vh_sim_tp, int32_t vh_ref)
        sector_hits.emplace_back(conv_hit.Subsystem(), conv_hit.Station(), conv_hit.Ring(),
            util.find_endsec(conv_hit), util.find_fr(conv_hit), conv_hit.BX(),
            util.find_emtf_layer(conv_hit), util.find_emtf_phi(conv_hit), util.find_emtf_theta(conv_hit),
            util.find_emtf_bend(conv_hit), util.find_emtf_qual(conv_hit), util.find_emtf_time(conv_hit),
            util.find_emtf_old_phi(conv_hit), util.find_emtf_old_bend(conv_hit),
            dummy_sim_tp, ihit);

        // Set sector_mode
        const Hit& hit = sector_hits.back();
        assert(0 <= hit.endsec && hit.endsec <= 11);
        assert(hit.emtf_layer != -99);

        if (optimize_for_cpu) {
          if (hit.type == TriggerPrimitive::kCSC) {
            sector_mode |= (1 << (4 - hit.station));
          } else if (hit.type == TriggerPrimitive::kME0) {
            sector_mode |= (1 << (4 - 1));
          } else if (hit.type == TriggerPrimitive::kDT) {
            sector_mode |= (1 << (4 - 1));
          }
        }
      }
    }  // end loop over conv_hits

    // Provide early exit if no hit in stations 1&2 (check CSC, ME0, DT)
    if (optimize_for_cpu) {
      if (!util.is_emtf_singlehit(sector_mode) && !util.is_emtf_singlehit_me2(sector_mode)) {
        return;
      }
    }

    // Apply patterns to the sector hits
    if (optimize_for_cpu) {
      apply_patterns(endcap, sector, sector_hits, sector_roads);
    } else {
      apply_patterns_unoptimized(endcap, sector, sector_hits, sector_roads);
    }

    // Sort the roads according to the road_id
    constexpr auto sort_roads_f = [](const Road& lhs, const Road& rhs) {
      return lhs.id() < rhs.id();
    };
    std::sort(sector_roads.begin(), sector_roads.end(), sort_roads_f);
    return;
  }

private:
  void create_road(const Road::road_id_t road_id, const Road::road_hits_t road_hits, std::vector<Road>& sector_roads) const {
    int32_t ipt  = road_id[2];  // road_id = (endcap, sector, ipt, ieta, iphi)
    int32_t ieta = road_id[3];

    // Check whether the road is OK
    bool accept = util.find_emtf_road_accept(ieta, road_hits);
    if (accept) {
      std::vector<int32_t> road_hits_layers;
      std::transform(road_hits.begin(), road_hits.end(), std::back_inserter(road_hits_layers),
          [](const auto& hit) -> int32_t { return hit.emtf_layer; });

      int32_t road_mode = util.find_emtf_road_mode(road_hits);
      int32_t road_quality = util.find_emtf_road_quality(ipt);
      int32_t road_sort_code = util.find_emtf_road_sort_code(road_quality, road_hits_layers);
      int32_t road_phi_median = 0;   // to be determined later
      int32_t road_theta_median = 0; // to be determined later

      //Road(int16_t vr_endcap, int16_t vr_sector, int16_t vr_ipt, int16_t vr_ieta, int16_t vr_iphi,
      //     const road_hits_t& vr_hits, int16_t vr_mode, int16_t vr_quality,
      //     int16_t vr_sort_code, int32_t vr_theta_median)
      sector_roads.emplace_back(road_id[0], road_id[1], road_id[2], road_id[3], road_id[4],
                                road_hits, road_mode, road_quality,
                                road_sort_code, road_phi_median, road_theta_median);
    }
    return;
  }

  void apply_patterns_in_zone(int32_t hit_zone, int32_t hit_lay,
                              std::vector<std::pair<int32_t, int32_t> >& result) const {
    result.clear();

    // Given zone & lay & quadstrip, only have to check against straightness
    //const auto& patterns_x0 = bank.x_array[hit_lay][hit_zone][0];
    //const auto& patterns_x1 = bank.x_array[hit_lay][hit_zone][2];
    //assert(patterns_x0.size() == patterns_x1.size());

    auto x0_iter = bank.x_array[hit_lay][hit_zone][0].begin();  // zero-copy op when using the iterators
    auto x0_end  = bank.x_array[hit_lay][hit_zone][0].end();
    auto x1_iter = bank.x_array[hit_lay][hit_zone][2].begin();
    auto x1_end  = bank.x_array[hit_lay][hit_zone][2].end();

    int32_t ipt  = 0;
    int32_t iphi = 0;
    for (; (x0_iter != x0_end) && (x1_iter != x1_end); ++x0_iter, ++x1_iter) {
      auto x0 = (*x0_iter);
      auto x1 = (*x1_iter);
      for (iphi = x0; iphi != (x1+1); ++iphi) {
        result.emplace_back(ipt, iphi);
      }
      ++ipt;
    }
    return;
  }

  void apply_patterns(int32_t endcap, int32_t sector,
                      const std::vector<Hit>& sector_hits, std::vector<Road>& sector_roads) const {

    // Create a map of road_id -> road_hits
    std::unordered_map<Road::road_id_t, Road::road_hits_t, Road::Hasher> amap;

    // Stores the results from pattern recognition (pairs of (ipt, iphi)-indices).
    std::vector<std::pair<int32_t, int32_t> > result;

    // Loop over hits
    for (const auto& hit : sector_hits) {
      int32_t hit_lay = hit.emtf_layer;
      int32_t hit_x   = util.find_pattern_x(hit.emtf_phi);
      const auto& hit_zones = util.find_emtf_zones(hit);

      // Loop over the zones that the hit is belong to
      for (const auto& hit_zone : hit_zones) {
        if (hit_zone == 6) {  // For now, ignore zone 6
          continue;
        }

        // Pattern recognition
        apply_patterns_in_zone(hit_zone, hit_lay, result);

        // Loop over the results from pattern recognition
        for (const auto& index : result) {
          int32_t ipt  = index.first;
          int32_t iphi = index.second;
          iphi         = (hit_x - iphi);
          int32_t ieta = hit_zone;

          // 'x' is the unit used in the patterns
          // Full range is 0 <= iphi <= 160. but a reduced range is sufficient (27% saving on patterns)
          if ((PATTERN_X_SEARCH_MIN <= iphi) && (iphi <= PATTERN_X_SEARCH_MAX)) {
            Road::road_id_t road_id {{endcap, sector, ipt, ieta, iphi}};
            amap[road_id].push_back(hit);
          }
        }
      }  // end loop over hit_zones
    }  // end loop over sector_hits

    // Create roads
    for (const auto& kv : amap) {
      const Road::road_id_t&   road_id   = kv.first;
      const Road::road_hits_t& road_hits = kv.second;
      create_road(road_id, road_hits, sector_roads);  // only valid roads are being appended to sector_roads
    }
    return;
  }

  void apply_patterns_unoptimized(int32_t endcap, int32_t sector,
                                  const std::vector<Hit>& sector_hits, std::vector<Road>& sector_roads) const {

    // Loop over all zones
    for (int32_t ieta = 0; ieta != PATTERN_BANK_NETA; ++ieta) {
      if (ieta == 6) {  // For now, ignore zone 6
        continue;
      }

      // Loop over all hits, find the ones that belong to this station and this zone
      std::vector<Hit> zone_hits;
      for (const auto& hit : sector_hits) {
        //int32_t hit_lay = hit.emtf_layer;
        //int32_t hit_x   = util.find_pattern_x(hit.emtf_phi);
        const auto& hit_zones = util.find_emtf_zones(hit);

        for (const auto& hit_zone : hit_zones) {
          if (hit_zone == ieta) {
            zone_hits.push_back(hit);
          }
        }  // end loop over hit_zones
      }  // end loop over sector_hits

      // Now loop over all the different shapes (straightness)
      for (int32_t ipt = 0; ipt != PATTERN_BANK_NPT; ++ipt) {

        // 'x' is the unit used in the patterns
        // Full range is 0 <= iphi <= 160. but a reduced range is sufficient (27% saving on patterns)
        for (int32_t iphi = PATTERN_X_SEARCH_MIN; iphi != (PATTERN_X_SEARCH_MAX+1); ++iphi) {
          Road::road_id_t road_id {{endcap, sector, ipt, ieta, iphi}};
          Road::road_hits_t road_hits;

          for (const auto& hit : zone_hits) {
            int32_t hit_lay = hit.emtf_layer;
            int32_t hit_x   = util.find_pattern_x(hit.emtf_phi);

            int32_t iphi_low  = bank.x_array[hit_lay][ieta][0][ipt];
            int32_t iphi_high = bank.x_array[hit_lay][ieta][2][ipt];

            if ((iphi + iphi_low <= hit_x) && (hit_x <= iphi + iphi_high)) {
              road_hits.push_back(hit);
            }
          }

          if (!road_hits.empty()) {
            create_road(road_id, road_hits, sector_roads);  // only valid roads are being appended to sector_roads
          }
        }  // end loop over x
      }  // end loop over all the different shapes (straightness)
    }  // end loop over all zones

    return;
  }
};

// RoadCleaning class removes ghost roads.
// During pattern matching, a set of hits can fire multiple patterns and create
// ghost roads. These ghost roads typically appear adjacent to each other in
// phi, so they appear to be clustered. We want to pick only one road out of
// the cluster. The roads are ranked by a sort code, which consists of the hit
// composition and the pattern straightness. The ghost cleaning should be
// aggressive enough, but not too aggressive such that di-muon efficiency is
// affected. At the end, a check of consistency with BX=0 is also applied.
// The output of this classis the subset of roads that are not identified as
// ghosts.
// In this C++ version, a more complicated algorithm is implemented, which uses
// clustering to do aggressive cleaning. But it might not be implementable
// in firmware. A simple local maximum finding algorithm should also work
// although the results are slightly different (not fully tested).

class RoadCleaning {
public:
  void run(const std::vector<Road>& roads, std::vector<Road>& clean_roads) const {

    // Optimize for CPU processing?
    bool optimize_for_cpu = true;

    if (optimize_for_cpu) {
      apply_cleaning(roads, clean_roads);
    } else {
      apply_cleaning_unoptimized(roads, clean_roads);
    }

    apply_additional_cleaning(clean_roads);
    return;
  }

private:
  void apply_cleaning(const std::vector<Road>& roads, std::vector<Road>& clean_roads) const {
    // Skip if no roads
    if (roads.empty()) {
      return;
    }

    // Create a map of road_id -> road
    using RoadPtr = const Road*;
    std::unordered_map<Road::road_id_t, RoadPtr, Road::Hasher> amap;

    // and a (sorted) vector of road_id's
    std::vector<Road::road_id_t> road_ids;

    for (const auto& road : roads) {
      Road::road_id_t road_id = road.id();
      amap[road_id] = &road;
      road_ids.push_back(road_id);
    }

    std::sort(road_ids.begin(), road_ids.end());

    constexpr auto make_row_splits = [](auto first, auto last) {
      // assume the input vector is sorted

      constexpr auto is_adjacent = [](const Road::road_id_t& prev, const Road::road_id_t& curr) {
        // adjacent if (x,y,z') == (x,y,z+1)
        return ((prev[0] == curr[0]) &&
                (prev[1] == curr[1]) &&
                (prev[2] == curr[2]) &&
                (prev[3] == curr[3]) &&
                ((prev[4]+1) == curr[4]));
      };

      std::vector<std::size_t> row_splits;
      row_splits.push_back(0);
      if (first == last) {
        return row_splits;
      }

      auto prev = first;
      auto curr = first;
      std::size_t i = 0;

      ++curr;
      ++i;

      for (; curr != last; ++prev, ++curr, ++i) {
        if (!is_adjacent(*prev, *curr)) {
          row_splits.push_back(i);
        }
      }
      row_splits.push_back(i);
      return row_splits;
    };

    // Make road clusters (groups)
    const std::vector<std::size_t>& splits = make_row_splits(road_ids.begin(), road_ids.end());
    assert(splits.size() >= 2);

    // Loop over the groups, pick the road with best sort code in each group
    using int32_t_pair = std::pair<int32_t, int32_t>;
    std::vector<Road> tmp_clean_roads;                    // the "best" roads in each group
    std::vector<int32_t> tmp_clean_roads_sortcode;        // keep track of the sort code of each group
    std::vector<int32_t_pair> tmp_clean_roads_groupinfo;  // keep track of the iphi range of each group

    std::vector<Road::road_id_t> group; // a group of road_id's
    int32_t best_sort_code = -1;        // keeps max sort code
    std::vector<RoadPtr> best_roads;    // keeps all the roads sharing the max sort code

    for (size_t igroup=0; igroup<(splits.size()-1); ++igroup) {
      group.clear();
      for (size_t i=splits[igroup]; i<splits[igroup+1]; ++i) {
        assert(i < road_ids.size());
        group.push_back(road_ids[i]);
      }

      best_sort_code = -1;
      for (const auto& road_id : group) {
        auto road_ptr = amap[road_id];
        if (best_sort_code < road_ptr->sort_code) {
          best_sort_code = road_ptr->sort_code;
        }
      }

      best_roads.clear();
      for (const auto& road_id : group) {
        auto road_ptr = amap[road_id];
        if (best_sort_code == road_ptr->sort_code) {
          best_roads.push_back(road_ptr);
        }
      }

      RoadPtr best_road = my_median_sorted(best_roads);
      tmp_clean_roads.push_back(*best_road);
      tmp_clean_roads_sortcode.push_back(best_sort_code);

      RoadPtr first_road = amap[group.front()];  // iphi range
      RoadPtr last_road = amap[group.back()];    // iphi range
      tmp_clean_roads_groupinfo.emplace_back(first_road->iphi, last_road->iphi);
    }  // end loop over groups

    if (tmp_clean_roads.empty())
      return;

    // Sort by 'sort code'
    const std::vector<size_t>& ind = my_argsort(tmp_clean_roads_sortcode, true);  // sort reverse

    // Loop over the sorted roads, kill the siblings
    for (size_t i=0; i<tmp_clean_roads.size(); ++i) {
      bool keep = true;

      // Check for intersection in the iphi range
      for (size_t j=0; j<i; ++j) {
        if (tmp_clean_roads[ind[i]].ieta == tmp_clean_roads[ind[j]].ieta) {  // same zone
          const auto& group_i = tmp_clean_roads_groupinfo[ind[i]];
          const auto& group_j = tmp_clean_roads_groupinfo[ind[j]];
          int32_t x1 = group_i.first;
          int32_t x2 = group_i.second;
          int32_t y1 = group_j.first;
          int32_t y2 = group_j.second;

          // No intersect between two ranges (x1, x2), (y1, y2): (x2 < y1) || (x1 > y2)
          // Intersect: !((x2 < y1) || (x1 > y2)) = (x2 >= y1) and (x1 <= y2)
          // Allow +/-2 due to extrapolation-to-EMTF error
          if (((x2+2) >= y1) && ((x1-2) <= y2)) {
            keep = false;
            break;
          }
        }
      }  // end inner loop over tmp_clean_roads[:i]

      // Do not share ME1/1, ME1/2, RE1/2, GE1/1, ME0, MB1, MB2
      if (keep) {
        using int32_t_pair = std::pair<int32_t, int32_t>;  // emtf_layer, emtf_phi

        constexpr auto make_hit_set = [](const auto& hits) {
          std::set<int32_t_pair> s;
          for (const auto& hit : hits) {
            if ((hit.emtf_layer == 0) ||
                (hit.emtf_layer == 1) ||
                (hit.emtf_layer == 5) ||
                (hit.emtf_layer == 9) ||
                (hit.emtf_layer == 11) ||
                (hit.emtf_layer == 12) ||
                (hit.emtf_layer == 13) ) {
              s.insert(std::make_pair(hit.endsec*100 + hit.emtf_layer, hit.emtf_phi));
            }
          }
          return s;
        };

        const auto& road_i = tmp_clean_roads[ind[i]];
        const std::set<int32_t_pair>& s1 = make_hit_set(road_i.hits);
        for (size_t j=0; j<i; ++j) {
          const auto& road_j = tmp_clean_roads[ind[j]];
          const std::set<int32_t_pair>& s2 = make_hit_set(road_j.hits);

          std::vector<int32_t_pair> v_intersection;
          std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v_intersection));
          if (!v_intersection.empty()) {  // has sharing
            keep = false;
            break;
          }
        }  // end inner loop over tmp_clean_roads[:i]
      }

      if (keep) {
        const auto& road_i = tmp_clean_roads[ind[i]];
        clean_roads.push_back(road_i);
      }
    }  // end loop over tmp_clean_roads.size()
    return;
  }

  void apply_cleaning_unoptimized(const std::vector<Road>& roads, std::vector<Road>& clean_roads) const {
    std::vector<bool> roads_mask(roads.size(), true);  // true: keep the road

    for (int32_t ieta = 0; ieta != PATTERN_BANK_NETA; ++ieta) {
      if (ieta == 6) {  // For now, ignore zone 6
        continue;
      }

      // Fill the quality codes
      std::array<int, PATTERN_X_SEARCH_MAX+1> quality_codes;
      quality_codes.fill(0);

      for (const auto& road : roads) {
        if (road.ieta == ieta) {
          int this_code = road.sort_code;
          if (quality_codes.at(road.iphi) < this_code) {
            quality_codes.at(road.iphi) = this_code;
          }
        }
      }

      // Check if this quality code is the (local) maximum
      size_t iroad = 0;
      for (const auto& road : roads) {
        if (road.ieta == ieta) {
          int this_code = road.sort_code;

          // Center quality is the current one
          int qc = quality_codes.at(road.iphi);
          // Left and right qualities are the neighbors
          // Protect against the right end and left end special cases
          int qr = (road.iphi == PATTERN_X_SEARCH_MAX) ? 0 : quality_codes.at(road.iphi+1);
          int ql = (road.iphi == 0) ? 0 : quality_codes.at(road.iphi-1);

          // Cancellation conditions
          if ((this_code <= ql) || (this_code < qr) || (this_code < qc)) {  // this pattern is lower quality than neighbors
            roads_mask.at(iroad) = false;  // cancel
          }
        }
        ++iroad;
      }
    }

    // Do the cancellation
    std::vector<Road> tmp_clean_roads;
    std::vector<int32_t> tmp_clean_roads_sortcode;
    {
      size_t iroad = 0;
      for (const auto& road : roads) {
        if (roads_mask.at(iroad) == true) {
          tmp_clean_roads.push_back(road);
          tmp_clean_roads_sortcode.push_back(road.sort_code);
        }
        ++iroad;
      }
    }

    if (tmp_clean_roads.empty())
      return;

    // Sort by 'sort code'
    const std::vector<size_t>& ind = my_argsort(tmp_clean_roads_sortcode, true);  // sort reverse

    // Loop over the sorted roads, kill the siblings
    for (size_t i=0; i<tmp_clean_roads.size(); ++i) {
      bool keep = true;

      // Do not share ME1/1, ME1/2, RE1/2, GE1/1, ME0, MB1, MB2
      if (keep) {
        using int32_t_pair = std::pair<int32_t, int32_t>;  // emtf_layer, emtf_phi

        constexpr auto make_hit_set = [](const auto& hits) {
          std::set<int32_t_pair> s;
          for (const auto& hit : hits) {
            if ((hit.emtf_layer == 0) ||
                (hit.emtf_layer == 1) ||
                (hit.emtf_layer == 5) ||
                (hit.emtf_layer == 9) ||
                (hit.emtf_layer == 11) ||
                (hit.emtf_layer == 12) ||
                (hit.emtf_layer == 13) ) {
              s.insert(std::make_pair(hit.endsec*100 + hit.emtf_layer, hit.emtf_phi));
            }
          }
          return s;
        };

        const auto& road_i = tmp_clean_roads[ind[i]];
        const std::set<int32_t_pair>& s1 = make_hit_set(road_i.hits);
        for (size_t j=0; j<i; ++j) {
          const auto& road_j = tmp_clean_roads[ind[j]];
          const std::set<int32_t_pair>& s2 = make_hit_set(road_j.hits);

          std::vector<int32_t_pair> v_intersection;
          std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v_intersection));
          if (!v_intersection.empty()) {  // has sharing
            keep = false;
            break;
          }
        }  // end inner loop over tmp_clean_roads[:i]
      }

      if (keep) {
        const auto& road_i = tmp_clean_roads[ind[i]];
        clean_roads.push_back(road_i);
      }
    }  // end loop over tmp_clean_roads.size()
    return;
  }

  bool select_bx_zero(const Road& road) const {
    int bx_counter1 = 0;  // count hits with BX <= -1
    int bx_counter2 = 0;  // count hits with BX == 0
    int bx_counter3 = 0;  // count hits with BX >= +1

    std::set<int32_t> s;  // check if layer has been used

    for (const auto& hit : road.hits) {
      if (s.find(hit.emtf_layer) == s.end()) {  // !s.contains(hit.emtf_layer)
        s.insert(hit.emtf_layer);
        if (hit.bx <= -1) {
          ++bx_counter1;
        } else if (hit.bx == 0) {
          ++bx_counter2;
        } else if (hit.bx >= +1) {
          ++bx_counter3;
        }
      }
    }

    //bool ret = (bx_counter1 < 2) && (bx_counter2 >= 2);
    bool ret = (bx_counter1 <= 2) && (bx_counter2 >= 2) && (bx_counter3 <= 2);
    return ret;
  }

  bool select_theta_aligned(Road& road) const {  // pass by reference to modify road.hits
    std::vector<int32_t> road_hits_thetas;
    std::transform(road.hits.begin(), road.hits.end(), std::back_inserter(road_hits_thetas),
        [](const auto& hit) -> int32_t { return hit.emtf_theta; });
    int32_t road_theta_median = my_median_unsorted(road_hits_thetas);

    std::vector<Hit> tmp_road_hits;
    for (const auto& hit : road.hits) {
      int32_t type    = hit.type;
      int32_t station = hit.station;
      int32_t ring    = hit.ring;

      int32_t dtheta = std::abs(hit.emtf_theta - road_theta_median);
      int32_t cut = 12;
      if (type == TriggerPrimitive::kCSC) {
        if (station == 1) {
          cut = 4;
        } else {
          cut = 2;
        }
      } else if (type == TriggerPrimitive::kRPC) {
        if (ring == 1) {
          cut = 2;
        } else {
          cut = 8;
        }
      } else if (type == TriggerPrimitive::kGEM) {
        cut = 6;
      } else if (type == TriggerPrimitive::kME0) {
        cut = 4;
      } else {
        cut = 12;
      }

      if (road.quality >= 5) {
        if (dtheta <= cut) {
          tmp_road_hits.push_back(hit);
        }
      } else {
        if (dtheta <= (cut*2)) {
          tmp_road_hits.push_back(hit);
        }
      }
    }

    // Overwrite road hits
    std::swap(road.hits, tmp_road_hits);

    // Overwrite road mode
    road.mode = util.find_emtf_road_mode(road.hits);

    // Check whether the road is OK, after the road hits are overwritten
    bool accept = util.find_emtf_road_accept(road.ieta, road.hits);
    return accept;
  }

  void apply_additional_cleaning(std::vector<Road>& clean_roads) const {
    // Finally, check consistency with BX=0
    std::vector<Road> tmp_clean_roads;
    std::swap(tmp_clean_roads, clean_roads);

    for (auto&& road_i : tmp_clean_roads) {
      if (select_bx_zero(road_i) && select_theta_aligned(road_i)) {
        clean_roads.push_back(road_i);
      }
    }
    return;
  }
};

// RoadSlimming class selects the unique hit for each station in a given road.
// Multiple hits in the same station can fire the same pattern and belong to
// the same road. We want to pick only one hit in each station. The hits are
// ranked by (max qual, min dtheta, min dphi), where dtheta is the difference
// between the hit theta and the road theta, and dphi is the difference
// between the hit phi and the road phi + an offset term. The offset terms
// are pre-defined according to the pattern straightness. The output of this
// class is the same road collection but the roads have at most one hit in
// each station.
// In this C++ version, it might look straight forward, thanks to being able
// to keep a collection of hits in each road. But this would not work in
// firmware. The 'primitive matching' step is also known to be the most
// resource-expensive operation in the firmware. So this might be tricky.

class RoadSlimming {
public:
  void run(const std::vector<Road>& clean_roads, std::vector<Road>& slim_roads) const {

    // Loop over roads
    for (const auto& road : clean_roads) {
      const int32_t ipt  = road.ipt;
      const int32_t ieta = road.ieta;

      // Retrieve the offset terms for each emtf_layer
      std::array<int32_t, NLAYERS> patterns_xc;
      for (size_t i=0; i<patterns_xc.size(); ++i) {
        int32_t xc = bank.x_array[i][ieta][1][ipt];
        patterns_xc[i] = util.find_pattern_x_inverse(xc);
      }

      // Find the median phi and theta
      // Note: they do not have to be exact. An approximation is good enough, provided that it is stable against outliers.
      std::vector<int32_t> road_hits_phis;
      std::transform(road.hits.begin(), road.hits.end(), std::back_inserter(road_hits_phis),
          [&patterns_xc](const auto& hit) -> int32_t { return (hit.emtf_phi - patterns_xc[hit.emtf_layer]); });
      int32_t road_phi_median = my_median_unsorted(road_hits_phis);

      std::vector<int32_t> road_hits_thetas;
      std::transform(road.hits.begin(), road.hits.end(), std::back_inserter(road_hits_thetas),
          [](const auto& hit) -> int32_t { return hit.emtf_theta; });
      int32_t road_theta_median = my_median_unsorted(road_hits_thetas);

      // Loop over all the emtf_layer's, select unique hit for each emtf_layer
      std::vector<Hit> slim_road_hits;

      using int32_t_tuple = std::tuple<int32_t, int32_t, int32_t, int32_t>;  // ihit, dphi, dtheta, neg_qual
      std::vector<int32_t_tuple> sort_criteria;  // for sorting hits

      constexpr auto sort_criteria_f = [](const int32_t_tuple& lhs, const int32_t_tuple& rhs) {
        // (max qual, min dtheta, min dphi, min ihit) is better
        auto [lhs0, lhs1, lhs2, lhs3] = lhs;
        auto [rhs0, rhs1, rhs2, rhs3] = rhs;
        return std::tie(lhs3, lhs2, lhs1, lhs0) < std::tie(rhs3, rhs2, rhs1, rhs0);
      };

      for (size_t i=0; i<patterns_xc.size(); ++i) {
        sort_criteria.clear();

        int32_t hit_lay = i;
        int32_t phi_offset = patterns_xc[hit_lay];

        int32_t ihit = 0;

        for (const auto& hit : road.hits) {
          if (hit_lay == hit.emtf_layer) {
            int32_t dphi     = std::abs(hit.emtf_phi - (road_phi_median + phi_offset));
            int32_t dtheta   = std::abs(hit.emtf_theta - road_theta_median);
            int32_t neg_qual = -std::abs(hit.emtf_qual);
            sort_criteria.emplace_back(ihit, dphi, dtheta, neg_qual);
          }
          ++ihit;
        }

        // Find the best hit, which is (max qual, min dtheta, min dphi)
        if (!sort_criteria.empty()) {
          std::nth_element(sort_criteria.begin(), sort_criteria.begin() + 1, sort_criteria.end(), sort_criteria_f);  // only care about the best
          int32_t best_ihit = std::get<0>(sort_criteria.front());
          slim_road_hits.emplace_back(road.hits[best_ihit]);
        }
      }

      //Road(int16_t vr_endcap, int16_t vr_sector, int16_t vr_ipt, int16_t vr_ieta, int16_t vr_iphi,
      //     const road_hits_t& vr_hits, int16_t vr_mode, int16_t vr_quality,
      //     int16_t vr_sort_code, int32_t vr_theta_median)
      slim_roads.emplace_back(road.endcap, road.sector, road.ipt, road.ieta, road.iphi,
                              slim_road_hits, road.mode, road.quality,
                              road.sort_code, road_phi_median, road_theta_median);
    }  // end loop over clean_roads
    return;
  }
};

// PtAssignment class assigns 2 parameters: pT and PU discr
// Currently we take 36 variables for each road, send them to the NN (by
// calling Tensorflow lib) and get the 2 parameters. The NN is stored in a
// Tensorflow 'protobuf' file. The outputs of this class are the input and
// output of the NN.

class PtAssignment {
public:
  explicit PtAssignment() {
    std::string cmssw_base = std::getenv("CMSSW_BASE");
    pbFileName = "/src/L1Trigger/L1TMuonEndCap/data/emtfpp_tf_graphs/model_graph.29.pb";
    pbFileName = cmssw_base + pbFileName;
    inputName = "batch_normalization_1_input";
    outputNames = {"dense_4/BiasAdd"};

    graphDef = tensorflow::loadGraphDef(pbFileName);
    assert(graphDef != nullptr);
    session = tensorflow::createSession(graphDef);
    assert(session != nullptr);

    // Add 2nd NN dedicated for displaced muons.
    pbFileName2 = "/src/L1Trigger/L1TMuonEndCap/data/emtfpp_tf_graphs/model_graph.displ.3.pb";
    pbFileName2 = cmssw_base + pbFileName2;
    inputName2 = "batch_normalization_1_input";
    outputNames2 = {"dense_5/BiasAdd"};

    graphDef2 = tensorflow::loadGraphDef(pbFileName2);
    assert(graphDef2 != nullptr);
    session2 = tensorflow::createSession(graphDef2);
    assert(session2 != nullptr);
  }

  // Destructor
  ~PtAssignment() {
    tensorflow::closeSession(session);
    delete graphDef;

    tensorflow::closeSession(session2);
    delete graphDef2;
  }

  // Copy constructor
  PtAssignment(const PtAssignment& other) {
    graphDef = other.graphDef;
    session = other.session;
    pbFileName = other.pbFileName;
    inputName = other.inputName;
    outputNames = other.outputNames;

    graphDef2 = other.graphDef2;
    session2 = other.session2;
    pbFileName2 = other.pbFileName2;
    inputName2 = other.inputName2;
    outputNames2 = other.outputNames2;
  }

  // Copy assignment
  PtAssignment& operator=(const PtAssignment& other) {
    if (this != &other) {
      graphDef = other.graphDef;
      session = other.session;
      pbFileName = other.pbFileName;
      inputName = other.inputName;
      outputNames = other.outputNames;

      graphDef2 = other.graphDef2;
      session2 = other.session2;
      pbFileName2 = other.pbFileName2;
      inputName2 = other.inputName2;
      outputNames2 = other.outputNames2;
    }
    return *this;
  }

  void run(const std::vector<Road>& slim_roads,
           std::vector<Feature>& features, std::vector<Prediction>& predictions) const {

    // Loop over roads
    for (const auto& road : slim_roads) {
      Feature feature;
      Prediction prediction;
      feature.fill(0);
      prediction.fill(0);

      predict(road, feature, prediction);
      features.push_back(feature);
      predictions.push_back(prediction);
    }  // end loop over slim_roads

    assert(slim_roads.size() == features.size());
    assert(slim_roads.size() == predictions.size());
    return;
  }

  void predict(const Road& road, Feature& feature, Prediction& prediction) const {
    preprocessing(road, feature);
    call_tensorflow(feature, prediction);
    postprocessing(road, feature, prediction);

    // Add 2nd NN dedicated for displaced muons.
    // It overwrites 'prediction', but not 'feature'.
    bool add_displ = true;
    if (add_displ) {
      Feature feature_displ;
      Prediction prediction_displ;
      feature_displ.fill(0);
      prediction_displ.fill(0);

      preprocessing_displ(road, feature_displ);
      call_tensorflow_displ(feature_displ, prediction_displ);
      postprocessing_displ(road, feature_displ, prediction_displ);

      prediction.at(3) = prediction_displ.at(0);  // y_displ
      prediction.at(6) = prediction_displ.at(1);  // d0_displ
    }
    return;
  }

private:
  void preprocessing(const Road& road, Feature& feature) const {
    static std::array<float, NLAYERS> x_phi;   // delta-phis = (raw phis - road_phi_median)
    static std::array<float, NLAYERS> x_theta; // raw thetas
    static std::array<float, NLAYERS> x_bend;
    static std::array<float, NLAYERS> x_qual;
    static std::array<float, NLAYERS> x_time;

    // Initialize to zeros
    x_phi.fill(0);
    x_theta.fill(0);
    x_bend.fill(0);
    x_qual.fill(0);
    x_time.fill(0);

    // Set the values
    for (const auto& hit : road.hits) {
      int32_t hit_lay = hit.emtf_layer;
      assert(std::abs(x_phi.at(hit_lay)) < 1e-7);   // sanity check
      x_phi[hit_lay] = (hit.emtf_phi - road.phi_median);
      assert(std::abs(x_theta.at(hit_lay)) < 1e-7); // sanity check
      x_theta[hit_lay] = hit.emtf_theta;
      assert(std::abs(x_bend.at(hit_lay)) < 1e-7);  // sanity check
      x_bend[hit_lay] = hit.emtf_bend;
      assert(std::abs(x_qual.at(hit_lay)) < 1e-7);  // sanity check
      x_qual[hit_lay] = hit.emtf_qual;
      assert(std::abs(x_time.at(hit_lay)) < 1e-7);  // sanity check
      x_time[hit_lay] = hit.emtf_time;
    }

    // Pack the 36 features
    // 20 (CSC) + 8 (RPC) + 4 (GEM) + 4 (ME0)
    feature = {{
        x_phi  [0], x_phi  [1], x_phi  [2], x_phi  [3], x_phi  [4] , x_phi  [5] ,
        x_phi  [6], x_phi  [7], x_phi  [8], x_phi  [9], x_phi  [10], x_phi  [11],
        x_theta[0], x_theta[1], x_theta[2], x_theta[3], x_theta[4] , x_theta[5] ,
        x_theta[6], x_theta[7], x_theta[8], x_theta[9], x_theta[10], x_theta[11],
        x_bend [0], x_bend [1], x_bend [2], x_bend [3], x_bend [4] , x_bend [11],
        x_qual [0], x_qual [1], x_qual [2], x_qual [3], x_qual [4] , x_qual [11]
    }};
    return;
  }

  void call_tensorflow(const Feature& feature, Prediction& prediction) const {
    static tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, NFEATURES });
    static std::vector<tensorflow::Tensor> outputs;
    assert(feature.size() == NFEATURES);

    float* d = input.flat<float>().data();
    std::copy(feature.begin(), feature.end(), d);
    tensorflow::run(session, { { inputName, input } }, outputNames, &outputs);
    assert(outputs.size() == 1);
    assert(prediction.size() == NPREDICTIONS);

    const float reg_pt_scale = 100.;  // a scale factor applied to regression during training
    const float reg_dxy_scale = 0.4;  // a scale factor applied to regression during training

    auto relu = [&](double x) {
      // ReLU(x) = max(0, x)
      return (x >= 0.) ? x : 0;
    };
    auto softplus = [&](double x) {
      // Softplus f(x) = log(1+exp(x))
      return relu(x) + std::log1p(std::exp(-std::abs(x)));
    };
    auto get_loc = [&](double loc) {
      double c = std::log(std::expm1(0.5 * reg_pt_scale));  // shifted to 2 GeV
      loc = loc + c;
      loc = (loc < 1./500 * reg_pt_scale) ? (1./500 * reg_pt_scale) : loc;
      return loc;
    };
    auto get_sign = [&](double sign) {
      return (sign >= 0.) ? 1 : -1;
    };
    auto get_loc_dxy = [&](double loc) {
      return loc;
    };
    auto get_sign_dxy = [&](double sign) {
      return 1e0;
    };
    auto get_scale = [&](double scale) {
      return 1e-5 + softplus(0.01 * scale);
    };

    prediction.at(0) = get_loc(outputs[0].matrix<float>()(0, 0));
    prediction.at(1) = get_sign(outputs[0].matrix<float>()(0, 1));
    prediction.at(2) = get_scale(outputs[0].matrix<float>()(0, 2));
    prediction.at(3) = get_loc(outputs[0].matrix<float>()(0, 3));
    prediction.at(4) = get_sign(outputs[0].matrix<float>()(0, 4));
    prediction.at(5) = get_scale(outputs[0].matrix<float>()(0, 5));
    prediction.at(6) = get_loc_dxy(outputs[0].matrix<float>()(0, 6));
    prediction.at(7) = get_sign_dxy(outputs[0].matrix<float>()(0, 7));
    prediction.at(8) = get_scale(outputs[0].matrix<float>()(0, 8));

    // Remove scale factor used during training
    prediction.at(0) /= reg_pt_scale;
    prediction.at(3) /= reg_pt_scale;
    prediction.at(6) /= reg_dxy_scale;

    // Include sign
    prediction.at(0) *= prediction.at(1);
    prediction.at(3) *= prediction.at(4);
    prediction.at(6) *= prediction.at(7);
    return;
  }

  void postprocessing(const Road& road, const Feature& feature, Prediction& prediction) const {
    // Demote tracks with large d0
    bool demote = false;
    float y_pred = prediction.at(0);
    float d0_pred = prediction.at(6);

    float discr_pt_cut_low = 4.;
    float discr_pt_cut_med = 8.;
    float discr_pt_cut_high = 14.;
    if (std::abs(1.0/y_pred) > discr_pt_cut_high) {       // >14 GeV
      demote = (std::abs(d0_pred) > 20.);
    } else if (std::abs(1.0/y_pred) > discr_pt_cut_med) { // 8-14 GeV
      demote = (std::abs(d0_pred) > 25.);
    } else if (std::abs(1.0/y_pred) > discr_pt_cut_low) { // 4-8 GeV
      demote = (std::abs(d0_pred) > 30.);
    }

    if (demote) {
      prediction.at(0) = 0.5;  // demote to 2 GeV
    }
    return;
  }

  void preprocessing_displ(const Road& road, Feature& feature) const {
    static std::array<float, NLAYERS> x_phi;   // delta-phis = (raw phis - road_phi_median)
    static std::array<float, NLAYERS> x_theta; // raw thetas
    static std::array<float, NLAYERS> x_bend;
    static std::array<float, NLAYERS> x_qual;
    static std::array<float, NLAYERS> x_time;

    // Initialize to zeros
    x_phi.fill(0);
    x_theta.fill(0);
    x_bend.fill(0);
    x_qual.fill(0);
    x_time.fill(0);

    // Set the values
    for (const auto& hit : road.hits) {
      int32_t hit_lay = hit.emtf_layer;

      // Drop iRPC
      if ((hit_lay == 7 || hit_lay == 8) && hit.ring == 1)
        continue;
      // Drop GE1/1, GE2/1, ME0, DT
      if (hit_lay >= 9)
        continue;

      assert(std::abs(x_phi.at(hit_lay)) < 1e-7);   // sanity check
      x_phi[hit_lay] = hit.old_emtf_phi;  // uses old_emtf_phi
      assert(std::abs(x_theta.at(hit_lay)) < 1e-7); // sanity check
      x_theta[hit_lay] = hit.emtf_theta;
      assert(std::abs(x_bend.at(hit_lay)) < 1e-7);  // sanity check
      x_bend[hit_lay] = hit.old_emtf_bend;  // uses old_emtf_bend
      assert(std::abs(x_qual.at(hit_lay)) < 1e-7);  // sanity check
      x_qual[hit_lay] = hit.emtf_qual;
      assert(std::abs(x_time.at(hit_lay)) < 1e-7);  // sanity check
      x_time[hit_lay] = hit.emtf_time;
    }

    // Mimic Phase-1 EMTF input calculations
    // 6 delta Phis: S1-S2, S1-S3, S1-S4, S2-S3, S2-S4, S3-S4
    // 6 delta Thetas: S1-S2, S1-S3, S1-S4, S2-S3, S2-S4, S3-S4
    // 4 bends : set to zero if no CSC hit and thus RPC hit is used
    // 1 FR bit: for ME1 only
    // 1 Ring bit: for ME1 only
    // 1 track Theta taken from stub coordinate in ME2, ME3, ME4 (in this priority)
    // 4 RPC bits indicating if ME or RE hit was used in each station (S1, S2, S3, S4)
    // Total: 23 variables
    static std::array<float, 6> x_dphi;
    static std::array<float, 6> x_dtheta;

    static std::array<float, 4> x_phi_emtf;   // temporary
    static std::array<float, 4> x_theta_emtf; // temporary
    static std::array<float, 4> x_bend_emtf;
    static std::array<float, 1> x_fr_emtf;
    static std::array<float, 1> x_trk_theta;
    static std::array<float, 1> x_me11ring;
    static std::array<float, 4> x_rpcbit;

    // Initialize to zeros
    x_dphi.fill(0);
    x_dtheta.fill(0);
    //
    x_phi_emtf.fill(0);
    x_theta_emtf.fill(0);
    x_bend_emtf.fill(0);
    x_fr_emtf.fill(0);
    x_trk_theta.fill(0);
    x_me11ring.fill(0);
    x_rpcbit.fill(0);

    // Station 1
    if (x_theta[0] > 1e-7) {  // ME1/1
      x_phi_emtf[0]   = x_phi[0];
      x_theta_emtf[0] = x_theta[0];
      x_bend_emtf[0]  = x_bend[0];
      x_fr_emtf[0]    = (x_qual[0] > 0) ? 1 : -1;
    } else if (x_theta[1] > 1e-7) {  // ME1/2
      x_phi_emtf[0]   = x_phi[1];
      x_theta_emtf[0] = x_theta[1];
      x_bend_emtf[0]  = x_bend[1];
      x_fr_emtf[0]    = (x_qual[1] > 0) ? 1 : -1;
      x_me11ring[0]   = 1;
    } else if (x_theta[5] > 1e-7) {  // RE1
      x_phi_emtf[0]   = x_phi[5];
      x_theta_emtf[0] = x_theta[5];
      x_bend_emtf[0]  = 0;
      x_rpcbit[0]     = 1;
    }

    // Station 2
    if (x_theta[2] > 1e-7) {  // ME2
      x_phi_emtf[1]   = x_phi[2];
      x_theta_emtf[1] = x_theta[2];
      x_bend_emtf[1]  = x_bend[2];
    } else if (x_theta[6] > 1e-7) {  // RE2
      x_phi_emtf[1]   = x_phi[6];
      x_theta_emtf[1] = x_theta[6];
      x_bend_emtf[1]  = 0;
      x_rpcbit[1]     = 1;
    }

    // Station 3
    if (x_theta[3] > 1e-7) {  // ME3
      x_phi_emtf[2]   = x_phi[3];
      x_theta_emtf[2] = x_theta[3];
      x_bend_emtf[2]  = x_bend[3];
    } else if (x_theta[7] > 1e-7) {  // RE3
      x_phi_emtf[2]   = x_phi[7];
      x_theta_emtf[2] = x_theta[7];
      x_bend_emtf[2]  = 0;
      x_rpcbit[2]     = 1;
    }

    // Station 4
    if (x_theta[4] > 1e-7) {  // ME4
      x_phi_emtf[3]   = x_phi[4];
      x_theta_emtf[3] = x_theta[4];
      x_bend_emtf[3]  = x_bend[4];
    } else if (x_theta[8] > 1e-7) {  // RE4
      x_phi_emtf[3]   = x_phi[8];
      x_theta_emtf[3] = x_theta[8];
      x_bend_emtf[3]  = 0;
      x_rpcbit[3]     = 1;
    }

    // Set x_trk_theta
    if (x_theta[2] > 1e-7) {
      x_trk_theta[0] = x_theta[2];
    } else if (x_theta[3] > 1e-7) {
      x_trk_theta[0] = x_theta[3];
    } else if (x_theta[4] > 1e-7) {
      x_trk_theta[0] = x_theta[4];
    }

    // Set x_dphi, x_dtheta
    auto calc_delta = [](float a, float b) {
      if ((a > 1e-7) && (b > 1e-7)) {
        return static_cast<float>(a-b);
      }
      return 0.f;
    };

    x_dphi[0] = calc_delta(x_phi_emtf[0], x_phi_emtf[1]);
    x_dphi[1] = calc_delta(x_phi_emtf[0], x_phi_emtf[2]);
    x_dphi[2] = calc_delta(x_phi_emtf[0], x_phi_emtf[3]);
    x_dphi[3] = calc_delta(x_phi_emtf[1], x_phi_emtf[2]);
    x_dphi[4] = calc_delta(x_phi_emtf[1], x_phi_emtf[3]);
    x_dphi[5] = calc_delta(x_phi_emtf[2], x_phi_emtf[3]);

    x_dtheta[0] = calc_delta(x_theta_emtf[0], x_theta_emtf[1]);
    x_dtheta[1] = calc_delta(x_theta_emtf[0], x_theta_emtf[2]);
    x_dtheta[2] = calc_delta(x_theta_emtf[0], x_theta_emtf[3]);
    x_dtheta[3] = calc_delta(x_theta_emtf[1], x_theta_emtf[2]);
    x_dtheta[4] = calc_delta(x_theta_emtf[1], x_theta_emtf[3]);
    x_dtheta[5] = calc_delta(x_theta_emtf[2], x_theta_emtf[3]);

    // Pack the 23 features + 13 zeros used for padding
    feature = {{
      x_dphi     [0], x_dphi     [1], x_dphi     [2], x_dphi     [3], x_dphi     [4], x_dphi     [5],
      x_dtheta   [0], x_dtheta   [1], x_dtheta   [2], x_dtheta   [3], x_dtheta   [4], x_dtheta   [5],
      x_bend_emtf[0], x_bend_emtf[1], x_bend_emtf[2], x_bend_emtf[3], x_fr_emtf  [0], x_trk_theta[0],
      x_me11ring [0], x_rpcbit   [0], x_rpcbit   [1], x_rpcbit   [2], x_rpcbit   [3], 0             ,
      0             , 0             , 0             , 0             , 0             , 0             ,
      0             , 0             , 0             , 0             , 0             , 0
    }};
    return;
  }

  void call_tensorflow_displ(const Feature& feature, Prediction& prediction) const {
    const int nfeatures_displ = 23;    // 23 features
    static tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, nfeatures_displ });
    static std::vector<tensorflow::Tensor> outputs;
    //assert(feature.size() == NFEATURES);

    float* d = input.flat<float>().data();
    //std::copy(feature.begin(), feature.end(), d);
    std::copy(feature.begin(), feature.begin() + nfeatures_displ, d);
    tensorflow::run(session2, { { inputName2, input } }, outputNames2, &outputs);
    assert(outputs.size() == 1);
    //assert(prediction.size() == NPREDICTIONS);

    const float reg_pt_scale = 100.;  // a scale factor applied to regression during training
    const float reg_dxy_scale = 1.0;  // a scale factor applied to regression during training

    prediction.at(0) = outputs[0].matrix<float>()(0, 0);
    prediction.at(1) = outputs[0].matrix<float>()(0, 1);

    // Remove scale factor used during training
    prediction.at(0) /= reg_pt_scale;
    prediction.at(1) /= reg_dxy_scale;
    return;
  }

  void postprocessing_displ(const Road& road, const Feature& feature, Prediction& prediction) const {
    // Demote tracks with less than 3 stations (CSC/RPC)
    bool demote = false;
    int road_mode = 0;

    for (const auto& hit : road.hits) {
      int32_t hit_lay = hit.emtf_layer;
      if (hit_lay == 0 || hit_lay == 1 || hit_lay == 5) {
        road_mode |= (1 << 3);
      } else if (hit_lay == 2 || hit_lay == 6) {
        road_mode |= (1 << 2);
      } else if (hit_lay == 3 || hit_lay == 7) {
        road_mode |= (1 << 1);
      } else if (hit_lay == 4 || hit_lay == 8) {
        road_mode |= (1 << 0);
      }
    }

    if (!(road_mode == 11 || road_mode == 13 || road_mode == 14 || road_mode == 15)) {
      demote = true;
    }
    if (demote) {
      prediction.at(0) = 0.5;  // demote to 2 GeV
    }
    return;
  }

  // TensorFlow components
  tensorflow::GraphDef* graphDef;
  tensorflow::Session* session;
  std::string pbFileName;
  std::string inputName;
  std::vector<std::string> outputNames;

  tensorflow::GraphDef* graphDef2;
  tensorflow::Session* session2;
  std::string pbFileName2;
  std::string inputName2;
  std::vector<std::string> outputNames2;
};

// TrackProducer class does 2 things: apply scaling (or calibration) to the NN
// pT, and apply cut on the PU discr.
// The pT from the NN needs to be scaled or calibrated so that when a cut is
// applied at L1, the efficiency at a given pT is about 90%. The PU discr cut
// is applied for tracks >8 GeV. A track with the scaled pT is created if it
// passes the PU discr cut. The output of this class is vector<Track> which
// contains all the tracks.

class TrackProducer {
public:
  constexpr TrackProducer() {
    discr_pt_cut_low = 4.;
    discr_pt_cut_med = 8.;
    discr_pt_cut_high = 14.;

    s_min   = 0.;
    s_max   = 60.;
    s_nbins = 120;
    s_step  = (s_max - s_min)/float(s_nbins);
    s_lut   = {{  1.8219,  1.5725,  1.6237,  1.8672,  2.2298,  2.6692,  3.1724,  3.7226,
                  4.3005,  4.8945,  5.4972,  6.1065,  6.7217,  7.3443,  7.9797,  8.6255,
                  9.2779,  9.9339, 10.5868, 11.2340, 11.8745, 12.5056, 13.1352, 13.7707,
                 14.3985, 15.0216, 15.6519, 16.2946, 16.9297, 17.5598, 18.1981, 18.8567,
                 19.5345, 20.2317, 20.9849, 21.7932, 22.5650, 23.2764, 23.9233, 24.5326,
                 25.1879, 25.9589, 26.8144, 27.6406, 28.3182, 28.9110, 29.4817, 30.0894,
                 30.8001, 31.5674, 32.3055, 33.0457, 33.8479, 34.6975, 35.4941, 36.2179,
                 36.9157, 37.6592, 38.5602, 39.6237, 40.7733, 41.9798, 43.2775, 44.6862,
                 45.9872, 46.8917, 47.5905, 48.2057, 48.8099, 49.4649, 50.1705, 50.8610,
                 51.5614, 52.2918, 53.0282, 53.7657, 54.5035, 55.2414, 55.9793, 56.7173,
                 57.4553, 58.1933, 58.9314, 59.6694, 60.4074, 61.1454, 61.8834, 62.6214,
                 63.3594, 64.0973, 64.8353, 65.5733, 66.3113, 67.0493, 67.7873, 68.5253,
                 69.2633, 70.0012, 70.7392, 71.4772, 72.2152, 72.9532, 73.6912, 74.4292,
                 75.1671, 75.9051, 76.6431, 77.3811, 78.1191, 78.8571, 79.5950, 80.3330,
                 81.0710, 81.8090, 82.5470, 83.2849, 84.0229, 84.7609, 85.4989, 86.2369}};
    s_lut_wp50 = {{1.8124,  1.6471,  1.6755,  1.8367,  2.0737,  2.3461,  2.6370,  2.9560,
                  3.3248,  3.7365,  4.1760,  4.6336,  5.1040,  5.5834,  6.0695,  6.5604,
                  7.0554,  7.5540,  8.0544,  8.5545,  9.0529,  9.5514, 10.0488, 10.5407,
                 11.0263, 11.5075, 11.9870, 12.4668, 12.9474, 13.4297, 13.9161, 14.4090,
                 14.9068, 15.4037, 15.8966, 16.3903, 16.8852, 17.3796, 17.8709, 18.3599,
                 18.8473, 19.3375, 19.8375, 20.3540, 20.8927, 21.4490, 21.9967, 22.5160,
                 23.0021, 23.4527, 23.8652, 24.2528, 24.6402, 25.0503, 25.4903, 25.9606,
                 26.4660, 27.0031, 27.5589, 28.1126, 28.6454, 29.1493, 29.6322, 30.1029,
                 30.5670, 31.0276, 31.4843, 31.9233, 32.3456, 32.7724, 33.2167, 33.6778,
                 34.1510, 34.6287, 35.1127, 35.6217, 36.1572, 36.7039, 37.2606, 37.8230,
                 38.3763, 38.9074, 39.4167, 39.9213, 40.4378, 40.9845, 41.5990, 42.2614,
                 42.9157, 43.5348, 44.1085, 44.6446, 45.1498, 45.6289, 46.0819, 46.5207,
                 46.9573, 47.3828, 47.7878, 48.1767, 48.5567, 48.9351, 49.3208, 49.7180,
                 50.1278, 50.5593, 51.0135, 51.4887, 51.9777, 52.4705, 52.9646, 53.4593,
                 53.9542, 54.4493, 54.9446, 55.4399, 55.9353, 56.4307, 56.9261, 57.4215}};
    assert(s_lut.size() == (size_t) s_nbins);
    assert(s_lut_wp50.size() == (size_t) s_nbins);
  }

  int digitize(float x) const {
    x = std::clamp(x, s_min, static_cast<float>(s_max - 1e-5));
    x = (x - s_min) / (s_max - s_min) * float(s_nbins);  // convert to bin number
    int binx = static_cast<int>(x);
    binx = (binx == s_nbins-1) ? (binx-1) : binx;  // avoid boundary
    return binx;
  }

  float interpolate(float x, float x0, float x1, float y0, float y1) const {
    float y = (x - x0) / (x1 - x0) * (y1 - y0) + y0;
    return y;
  }

  float get_trigger_pt(float y_pred) const {
    float xml_pt = std::abs(1.0/y_pred);
    if (xml_pt <= 2.) {  // do not use the LUT if below 2 GeV
      return xml_pt;
    }

    int binx = digitize(xml_pt);
    float x0 = float(binx) * s_step;
    float x1 = float(binx+1) * s_step;
    float y0 = s_lut.at(binx);
    float y1 = s_lut.at(binx+1);
    float trg_pt = interpolate(xml_pt, x0, x1, y0, y1);
    return trg_pt;
  }

  float get_trigger_pt_wp50(float y_pred) const {
    float xml_pt = std::abs(1.0/y_pred);
    if (xml_pt <= 2.) {  // do not use the LUT if below 2 GeV
      return xml_pt;
    }

    int binx = digitize(xml_pt);
    float x0 = float(binx) * s_step;
    float x1 = float(binx+1) * s_step;
    float y0 = s_lut_wp50.at(binx);
    float y1 = s_lut_wp50.at(binx+1);
    float trg_pt = interpolate(xml_pt, x0, x1, y0, y1);
    return trg_pt;
  }

  bool pass_trigger(int ndof, int mode, int strg, int zone, int theta_median, float y_pred, float y_discr, float d0_pred) const {
    //int ipt1 = strg;
    //int ipt2 = util.find_pt_bin(y_pred);
    //int quality1 = util.find_emtf_road_quality(ipt1);
    //int quality2 = util.find_emtf_road_quality(ipt2);
    //bool strg_ok = (quality2 <= (quality1+1));
    //float xml_pt = std::abs(1.0/y_pred);
    bool trigger = true;

    //// OBSOLETE since v4
    //// Apply cuts on d0
    //if (xml_pt > discr_pt_cut_high) {       // >14 GeV
    //  trigger = (std::abs(d0_pred) < 20.);
    //} else if (xml_pt > discr_pt_cut_med) { // 8-14 GeV
    //  trigger = (std::abs(d0_pred) < 25.);
    //} else if (xml_pt > discr_pt_cut_low) { // 4-8 GeV
    //  trigger = (std::abs(d0_pred) < 30.);
    //} else {
    //  trigger = (y_discr >= 0.);
    //}

    //// OBSOLETE since v3
    //// Apply cuts
    //bool trigger = false;
    //if (xml_pt > discr_pt_cut_high) {       // >14 GeV (98.5% coverage)
    //  trigger = (y_discr > 0.9600);
    //} else if (xml_pt > discr_pt_cut_med) { // 8-14 GeV (98.5% coverage)
    //  trigger = (y_discr > 0.8932);
    //} else if (xml_pt > discr_pt_cut_low) { // 4-8 GeV (99.0% coverage)
    //  trigger = (y_discr > 0.2000);
    //} else {
    //  trigger = (y_discr >= 0.) && strg_ok;
    //}
    return trigger;
  }

  float theta_to_eta_f(int theta_int) const {
    static const std::vector<float> theta_to_eta_lut = {
      2.599, 2.566, 2.534, 2.503, 2.473, 2.444, 2.415, 2.388, 2.361, 2.334,
      2.309, 2.284, 2.259, 2.236, 2.212, 2.190, 2.167, 2.145, 2.124, 2.103,
      2.083, 2.063, 2.043, 2.024, 2.005, 1.986, 1.968, 1.950, 1.932, 1.915,
      1.898, 1.881, 1.864, 1.848, 1.832, 1.816, 1.800, 1.785, 1.770, 1.755,
      1.740, 1.726, 1.711, 1.697, 1.683, 1.670, 1.656, 1.642, 1.629, 1.616,
      1.603, 1.590, 1.578, 1.565, 1.553, 1.541, 1.529, 1.517, 1.505, 1.493,
      1.482, 1.470, 1.459, 1.448, 1.436, 1.425, 1.415, 1.404, 1.393, 1.382,
      1.372, 1.362, 1.351, 1.341, 1.331, 1.321, 1.311, 1.301, 1.291, 1.282,
      1.272, 1.262, 1.253, 1.244, 1.234, 1.225, 1.216, 1.207, 1.198, 1.189,
      1.180, 1.171, 1.162, 1.154, 1.145, 1.136, 1.128, 1.119, 1.111, 1.103,
      1.094, 1.086, 1.078, 1.070, 1.062, 1.054, 1.046, 1.038, 1.030, 1.022,
      1.014, 1.007, 0.999, 0.991, 0.984, 0.976, 0.969, 0.961, 0.954, 0.946,
      0.939, 0.932, 0.924, 0.917, 0.910, 0.903, 0.896, 0.888, 0.881, 0.874,
      0.867, 0.860, 0.853, 0.847, 0.840, 0.833, 0.826, 0.819, 0.813, 0.806,
      0.799, 0.793, 0.786, 0.779, 0.773, 0.766, 0.760, 0.753, 0.747, 0.741
    };
    return theta_to_eta_lut.at(theta_int);
  }

  void run(const std::vector<Road>& slim_roads, const std::vector<Prediction>& predictions,
           std::vector<Track>& tracks) const {

    // Loop over roads & predictions
    auto predictions_it = predictions.begin();

    for (const auto& road : slim_roads) {
      const auto& prediction = *predictions_it++;

      float y_pred     = prediction.at(0);  // vtx-constrained pT
      float y_discr    = 1.0;               // PU discr (obsolete)
      float d1_pred    = prediction.at(3);  // vtx-unconstrained pT
      float d0_pred    = prediction.at(6);  // d0
      int ndof         = road.hits.size();
      int mode         = road.mode;
      int strg         = road.ipt;
      int zone         = road.ieta;
      int phi_median   = road.phi_median;
      int theta_median = road.theta_median;
      float eta        = theta_to_eta_f(theta_median);  // absolute eta

      bool passed = pass_trigger(ndof, mode, strg, zone, theta_median, y_pred, y_discr, d0_pred);

      if (passed) {
        float xml_pt = std::abs(1.0/y_pred);
        float pt = get_trigger_pt(y_pred);
        if ((2.15 <= eta) && (eta <= 2.25)) {
          pt = get_trigger_pt_wp50(y_pred);
        }

        float pt_displ = get_trigger_pt(d1_pred);

        int trk_q = (y_pred < 0) ? -1 : +1;
        //Track(int16_t vt_endcap, int16_t vt_sector, int16_t vt_ipt, int16_t vt_ieta, int16_t vt_iphi,
        //      const road_hits_t& vt_hits, int16_t vt_mode, int16_t vt_quality, int16_t vt_sort_code,
        //      float vt_xml_pt, float vt_pt, int16_t vt_q, float vt_y_pred, float vt_y_discr,
        //      float vt_y_displ, float vt_d0_displ, float vt_pt_displ, int32_t vt_emtf_phi, int32_t vt_emtf_theta)
        tracks.emplace_back(road.endcap, road.sector, road.ipt, road.ieta, road.iphi,
                            road.hits, mode, road.quality, road.sort_code,
                            xml_pt, pt, trk_q, y_pred, y_discr,
                            d1_pred, d0_pred, pt_displ, phi_median, theta_median);
      }
    }  // end loop over slim_roads, predictions
    return;
  }

private:
  // Used for pass_trigger()
  float discr_pt_cut_low {0.};
  float discr_pt_cut_med {0.};
  float discr_pt_cut_high {0.};

  // Used for get_trigger_pt()
  float s_min {0.};
  float s_max {0.};
  int   s_nbins {0};
  float s_step {0.};
  std::array<float, 120> s_lut {};
  std::array<float, 120> s_lut_wp50 {};
};

// GhostBusting class remove ghost tracks.
// This is very similar to the RoadCleaning, but now it is done on the tracks
// from all the sectors.

class GhostBusting {
public:
  void run(EMTFTrackCollection& tracks) const {

    EMTFTrackCollection tracks_after_gb;

    // Sort by 'sort code'
    constexpr auto sort_tracks_f = [](const EMTFTrack& lhs, const EMTFTrack& rhs) {
      return lhs.Winner() > rhs.Winner();  // EMTFTrack::Winner() was hacked to store Track::sort_code
    };
    std::sort(tracks.begin(), tracks.end(), sort_tracks_f);

    // Loop over the sorted tracks and remove duplicates (ghosts)
    for (size_t i=0; i<tracks.size(); ++i) {
      bool keep = true;

      // Do not share ME1/1, ME1/2, RE1/2, GE1/1, ME0, MB1, MB2
      // Need to check for neighbor sector hits
      if (keep) {
        using int32_t_pair = std::pair<int32_t, int32_t>;  // emtf_layer, emtf_phi

        constexpr auto make_hit_set = [](const auto& conv_hits) {
          std::set<int32_t_pair> s;
          for (const auto& conv_hit : conv_hits) {
            int hit_emtf_layer = util.find_emtf_layer(conv_hit);
            if ((hit_emtf_layer == 0) ||
                (hit_emtf_layer == 1) ||
                (hit_emtf_layer == 5) ||
                (hit_emtf_layer == 9) ||
                (hit_emtf_layer == 11) ||
                (hit_emtf_layer == 12) ||
                (hit_emtf_layer == 13) ) {

              int32_t hit_endsec = util.find_endsec(conv_hit);
              int32_t hit_emtf_phi = util.find_emtf_phi(conv_hit);
              int32_t tmp_endsec = hit_endsec;
              int32_t tmp_emtf_phi = hit_emtf_phi;
              if (hit_emtf_phi < (22*60)) {  // is a neighbor hit
                if ((hit_endsec == 0) || (hit_endsec == 6)) {
                  tmp_endsec += 5;
                } else if ((1 <= hit_endsec && hit_endsec <= 5) || (7 <= hit_endsec && hit_endsec <= 11)) {
                  tmp_endsec -= 1;
                }
                tmp_emtf_phi += (60*60);
              }
              s.insert(std::make_pair(tmp_endsec*100 + hit_emtf_layer, tmp_emtf_phi));
            }
          }
          return s;
        };

        const auto& track_i = tracks[i];
        const std::set<int32_t_pair>& s1 = make_hit_set(track_i.Hits());
        for (size_t j=0; j<i; ++j) {
          const auto& track_j = tracks[j];
          const std::set<int32_t_pair>& s2 = make_hit_set(track_j.Hits());

          std::vector<int32_t_pair> v_intersection;
          std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v_intersection));
          if (!v_intersection.empty()) {  // has sharing
            keep = false;
            break;
          }
        }  // end inner loop over tracks[:i]
      }

      if (keep) {
        const auto& track_i = tracks[i];
        tracks_after_gb.push_back(track_i);
      }
    }  // end loop over tracks

    std::swap(tracks, tracks_after_gb);
    return;
  }
};

// TrackConverter class converts the internal Track object into EMTFTrack.
// The EMTFTrackCollection can be used by the rest of CMSSW.

class TrackConverter {
public:
  void run(const PtAssignmentEngineAux& the_aux, const EMTFHitCollection& conv_hits,
           const std::vector<Track>& best_tracks, EMTFTrackCollection& best_emtf_tracks) const {

    // Loop over tracks
    for (const auto& track : best_tracks) {

      // Create PtAssignment::aux()
      auto aux = [the_aux](){ return the_aux; };

      // Create the hit collection
      EMTFHitCollection emtf_track_hits;
      emtf_track_hits.reserve(track.hits.size());
      for (const auto& hit : track.hits) {
        const EMTFHit& emtf_hit = conv_hits.at(hit.ref);
        emtf_track_hits.push_back(emtf_hit);
      }

      // Create the PTLUT data
      //CUIDADO: not filled
      EMTFPtLUT ptlut_data = {};

      // Create a track
      EMTFTrack emtf_track;

      // Setters
      // Part 1: from src/PrimitiveMatching.cc
      emtf_track.set_endcap     ( track.endcap );
      emtf_track.set_sector     ( track.sector );
      emtf_track.set_sector_idx ( (track.endcap == 1) ? (track.sector - 1) : (track.sector + 5) );
      emtf_track.set_bx         ( 0 );
      emtf_track.set_zone       ( track.ieta );
      //emtf_track.set_ph_num     ( road.Key_zhit() );
      //emtf_track.set_ph_q       ( road.Quality_code() );
      //emtf_track.set_rank       ( road.Quality_code() );
      //emtf_track.set_winner     ( road.Winner() );
      emtf_track.set_winner     ( track.sort_code ); // hacked to store Track::sort_code
      emtf_track.clear_Hits();
      emtf_track.set_Hits( emtf_track_hits );

      // Part 2: from src/AngleCalculation.cc
      //emtf_track.set_rank     ( rank );
      emtf_track.set_mode     ( track.mode );
      //emtf_track.set_mode_inv ( mode_inv );
      emtf_track.set_phi_fp   ( track.emtf_phi );
      emtf_track.set_theta_fp ( track.emtf_theta );
      emtf_track.set_PtLUT    ( ptlut_data );
      emtf_track.set_phi_loc  ( emtf::calc_phi_loc_deg(emtf_track.Phi_fp()) );
      emtf_track.set_phi_glob ( emtf::calc_phi_glob_deg(emtf_track.Phi_loc(), emtf_track.Sector()) );
      emtf_track.set_theta    ( emtf::calc_theta_deg_from_int(emtf_track.Theta_fp()) );
      emtf_track.set_eta      ( emtf::calc_eta_from_theta_deg(emtf_track.Theta(), emtf_track.Endcap()) );
      //emtf_track.clear_Hits();
      //emtf_track.set_Hits( tmp_hits );
      //emtf_track.set_first_bx  ( first_bx );
      //emtf_track.set_second_bx ( second_bx );

      // Part 3: from src/PtAssignment.cc
      //emtf_track.set_PtLUT    ( ptlut_data );
      emtf_track.set_pt_XML ( track.xml_pt );
      emtf_track.set_pt     ( track.pt );
      emtf_track.set_pt_dxy ( track.pt_displ );
      emtf_track.set_dxy    ( track.d0_displ );
      emtf_track.set_charge ( track.q );
      emtf_track.set_invpt_prompt( track.y_pred );
      emtf_track.set_invpt_displ ( track.y_displ );
      //
      int gmt_pt = aux().getGMTPt(emtf_track.Pt());
      int gmt_pt_dxy = aux().getGMTPtDxy(emtf_track.Pt_dxy());
      int gmt_dxy = aux().getGMTDxy(emtf_track.Dxy());
      int gmt_phi = aux().getGMTPhiV2(emtf_track.Phi_fp());
      int gmt_eta = aux().getGMTEta(emtf_track.Theta_fp(), emtf_track.Endcap());
      bool promoteMode7 = false;
      int modeQualVer = 2;
      int gmt_quality = aux().getGMTQuality(emtf_track.Mode(), emtf_track.Theta_fp(), promoteMode7, modeQualVer);
      int charge = 0;
      if (emtf_track.Charge() == 1)
        charge = 1;
      int charge_valid = 1;
      if (emtf_track.Charge() == 0)
        charge_valid = 0;
      std::pair<int, int> gmt_charge = std::make_pair(charge, charge_valid);
      emtf_track.set_gmt_pt           ( gmt_pt );
      emtf_track.set_gmt_pt_dxy       ( gmt_pt_dxy );
      emtf_track.set_gmt_dxy          ( gmt_dxy );
      emtf_track.set_gmt_phi          ( gmt_phi );
      emtf_track.set_gmt_eta          ( gmt_eta );
      emtf_track.set_gmt_quality      ( gmt_quality );
      emtf_track.set_gmt_charge       ( gmt_charge.first );
      emtf_track.set_gmt_charge_valid ( gmt_charge.second );

      // Part 4: from src/BestTrackSelection.cc
      emtf_track.set_track_num ( best_emtf_tracks.size() );
      //emtf_track.set_winner ( o );
      //emtf_track.set_bx ( second_bx );

      // Finally
      best_emtf_tracks.push_back(emtf_track);
    }  // end loop over tracks
    return;
  }
};

constexpr PatternRecognition recog;
constexpr RoadCleaning clean;
constexpr RoadSlimming slim;
static const PtAssignment assig;  // cannot use 'constexpr'
constexpr TrackProducer trkprod;
constexpr GhostBusting ghost;
constexpr TrackConverter trkconv;


// _____________________________________________________________________________
void Phase2SectorProcessor::build_tracks(
    // Input
    const EMTFHitCollection& conv_hits,
    // Output
    std::vector<Track>& best_tracks
) const {
  // Containers for each sector
  std::vector<Hit> hits;
  std::vector<Road> roads, clean_roads, slim_roads;
  std::vector<Feature> features;
  std::vector<Prediction> predictions;
  std::vector<Track> tracks;

  // Run the algorithms
  recog.run(endcap_, sector_, conv_hits, hits, roads);
  clean.run(roads, clean_roads);
  slim.run(clean_roads, slim_roads);
  assig.run(slim_roads, features, predictions);
  trkprod.run(slim_roads, predictions, tracks);

  best_tracks.insert(best_tracks.end(), tracks.begin(), tracks.end());

  // Debug
  bool debug = false;
  if (debug) {
    debug_tracks(hits, roads, clean_roads, slim_roads, tracks);
  }
  return;
}

// _____________________________________________________________________________
void Phase2SectorProcessor::convert_tracks(
    // Input
    const EMTFHitCollection& conv_hits,
    const std::vector<Track>& best_tracks,
    // Output
    EMTFTrackCollection& best_emtf_tracks
) const {
  // Run the algorithms
  trkconv.run(pt_assign_engine_->aux(), conv_hits, best_tracks, best_emtf_tracks);
  return;
}

// _____________________________________________________________________________
void Phase2SectorProcessor::debug_tracks(
    // Input
    const std::vector<Hit>& hits,
    const std::vector<Road>& roads,
    const std::vector<Road>& clean_roads,
    const std::vector<Road>& slim_roads,
    const std::vector<Track>& tracks
) const {
  size_t i = 0;
  size_t j = 0;

  std::cout << "SP e:" << endcap_ << " s:" << sector_ << " has "
      << hits.size() << " hits, " << roads.size() << " roads, "
      << clean_roads.size() << " clean roads, " << tracks.size() << " tracks"
      << std::endl;

  i = 0;
  for (const auto& hit : hits) {
    const auto& id = hit.id();
    std::cout << ".. hit " << i++ << " id: (" << id[0] << ", "
        << id[1] << ", " << id[2] << ", " << id[3] << ", "
        << id[4] << ", " << id[5] << ") lay: " << hit.emtf_layer
        << " ph: " << hit.emtf_phi << " (" << util.find_pattern_x(hit.emtf_phi)
        << ") th: " << hit.emtf_theta << " bd: " << hit.emtf_bend
        << " old_ph: " << hit.old_emtf_phi << " old_bd: " << hit.old_emtf_bend
        << " ql: " << hit.emtf_qual << " tp: " << hit.sim_tp << std::endl;
  }

  i = 0;
  for (const auto& road : roads) {
    const auto& id = road.id();
    std::cout << ".. road " << i++ << " id: (" << id[0] << ", "
        << id[1] << ", " << id[2] << ", " << id[3] << ", "
        << id[4] << ") nhits: " << road.hits.size()
        << " mode: " << road.mode << " qual: " << road.quality
        << " sort: " << road.sort_code << std::endl;
  }

  i = 0;
  for (const auto& road : clean_roads) {
    const auto& id = road.id();
    std::cout << ".. croad " << i++ << " id: (" << id[0] << ", "
        << id[1] << ", " << id[2] << ", " << id[3] << ", "
        << id[4] << ") nhits: " << road.hits.size()
        << " mode: " << road.mode << " qual: " << road.quality
        << " sort: " << road.sort_code << std::endl;

    j = 0;
    for (const auto& hit : road.hits) {
      const auto& hit_id = hit.id();
      std::cout << ".. .. hit " << j++ << " id: (" << hit_id[0] << ", "
          << hit_id[1] << ", " << hit_id[2] << ", " << hit_id[3] << ", "
          << hit_id[4] << ", " << hit_id[5] << ") lay: " << hit.emtf_layer
          << " ph: " << hit.emtf_phi << " th: " << hit.emtf_theta << std::endl;
    }
  }

  i = 0;
  for (const auto& road : slim_roads) {
    const auto& id = road.id();
    std::cout << ".. sroad " << i++ << " id: (" << id[0] << ", "
        << id[1] << ", " << id[2] << ", " << id[3] << ", "
        << id[4] << ") nhits: " << road.hits.size()
        << " mode: " << road.mode << " qual: " << road.quality
        << " sort: " << road.sort_code << std::endl;

    j = 0;
    for (const auto& hit : road.hits) {
      const auto& hit_id = hit.id();
      std::cout << ".. .. hit " << j++ << " id: (" << hit_id[0] << ", "
          << hit_id[1] << ", " << hit_id[2] << ", " << hit_id[3] << ", "
          << hit_id[4] << ", " << hit_id[5] << ") lay: " << hit.emtf_layer
          << " ph: " << hit.emtf_phi << " th: " << hit.emtf_theta << std::endl;
    }
  }

  i = 0;
  for (const auto& trk : tracks) {
    const auto& id = trk.id();
    std::cout << ".. trk " << i++ << " id: (" << id[0] << ", "
        << id[1] << ", " << id[2] << ", " << id[3] << ", "
        << id[4] << ") nhits: " << trk.hits.size()
        << " mode: " << trk.mode << " pt: " << trk.pt
        << " y_pred: " << trk.y_pred << " y_discr: " << trk.y_discr
        << " y_displ: " << trk.y_displ << " d0_displ: " << trk.d0_displ
        << " pt_displ: " << trk.pt_displ << std::endl;

    j = 0;
    for (const auto& hit : trk.hits) {
      const auto& hit_id = hit.id();
      std::cout << ".. .. hit " << j++ << " id: (" << hit_id[0] << ", "
          << hit_id[1] << ", " << hit_id[2] << ", " << hit_id[3] << ", "
          << hit_id[4] << ", " << hit_id[5] << ") lay: " << hit.emtf_layer
          << " ph: " << hit.emtf_phi << " th: " << hit.emtf_theta << std::endl;
    }
  }
}

// _____________________________________________________________________________
void Phase2SectorProcessor::ghost_busting(
    // Input & output
    EMTFTrackCollection& best_emtf_tracks
) const {
  // Run the algorithms
  ghost.run(best_emtf_tracks);
  return;
}

}  // namespace experimental
