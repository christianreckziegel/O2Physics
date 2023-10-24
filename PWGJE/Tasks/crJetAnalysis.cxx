// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// Personal analysis for 
// -plotting jet standard quantities: pt, eta, phi
// -comparing leading track within jet with jet
// -comparing leading track within collision with leading jet within collision
// -investigate tracks with higher pT than jet itself
// jetFinder task, inclusive jets, R=0.4
//
// Author: Christian Reckziegel
//

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/DataModel/EMCALClusters.h"

#include "Framework/runDataProcessing.h"

#include "Common/DataModel/TrackSelectionTables.h" // for DCA data

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;


/*class MyTrack{
	public:
		MyTrack();
		~MyTrack();
		float pt(){ return pt; }
		float eta(){ return eta; }
		float phi(){ return phi; }
		void setPt(float newPt){ pt = newPt}
		void setEta(float newEta){ eta = newEta}
		void setPhi(float newPhi){ phi = newPhi}
	private:
		float pt;
		float eta;
		float phi;
}
MyTrack::MyTrack(){
	pt = -1.;
	eta = 0.;
	phi = 0.;
} // constructor
MyTrack::~MyTrack(){} // destructor
*/

struct crJetAnalysis{
	Configurable<float> leadPtHistMax{"leadPtHistMax", 200., "maximum leadTrack pT histogram extreme"};
	Configurable<int> leadPtBins{"leadPtBins", 200, "lead track pT histogram bins"};
	Configurable<float> jetPtHistMax{"jetPtHistMax", 200., "maximum jet pT histogram extreme"};
	Configurable<int> hist_2D_Bins{"hist_2D_Bins", 100, "2D pT histogram bins"};
	
	HistogramRegistry registry{
		"registry",
		{
			{"h_leadTrack_pt", "Leading track within jet pT;#it{p}_{T,lead track}(GeV/#it{c});entries", {HistType::kTH1F, {{leadPtBins, 0., 50.}}}}, //0-200 GeV
			{"h_leadTrack_eta", "Leading track within jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{200, -10., 10.}}}},
			{"h_leadTrack_phi", "Leading track within jet #phi;#phi_{jet} (rads);entries", {HistType::kTH1F, {{140, -7.0, 7.0}}}},
			//{"h_jet_vs_track_pt", "Jet p_{T} in function of leading track from jet p_{T};p_{T,lead track} (GeV/c);p_{T,jet} (GeV/c)", {HistType::kTH2F, {{hist_2D_Bins, 0., 50.},{hist_2D_Bins,0.,50.}}}},
			// 2D histograms
			{"2DHistograms/h_Collision_jet_vs_track_pt", "Collision leading p_{T,jet} in function of collision leading p_{T,track};p_{T,lead track} (GeV/c);p_{T,lead jet} (GeV/c)", {HistType::kTH2F, {{hist_2D_Bins, 0., 50.},{hist_2D_Bins, 0.,50.}}}},
			{"2DHistograms/h_Collision_jet_vs_unbehavedTrack_pt", "Collision leading p_{T,jet} in function of collision unbehaved track p_{T,track};p_{T,track}^{unbehaved} (GeV/c);p_{T,lead jet} (GeV/c)", {HistType::kTH2F, {{hist_2D_Bins, 0., 50.},{hist_2D_Bins, 0.,50.}}}},
			{"2DHistograms/h_Collision_jet_vs_track_eta", "Collision leading #eta_{jet} in function of collision leading #eta_{track};#eta_{lead track};#eta_{lead jet}", {HistType::kTH2F, {{200, -10., 10.},{200, -10.,10.}}}},
			{"2DHistograms/h_Collision_jet_vs_unbehavedTrack_eta", "Collision leading #eta_{jet} in function of collision unbehaved track #eta_{track};#eta_{track}^{unbehaved};#eta_{lead jet}", {HistType::kTH2F, {{200, -10., 10.},{200, -10.,10.}}}},
			{"2DHistograms/h_Collision_jet_vs_track_phi", "Collision leading #phi_{jet} in function of collision leading #phi_{track};#phi_{lead track};#phi_{lead jet}", {HistType::kTH2F, {{80, -1., 7.},{80, -1.,7.}}}},
			{"2DHistograms/h_Collision_jet_vs_unbehavedTrack_phi", "Collision leading #phi_{jet} in function of collision unbehaved track #phi_{track};#phi_{track}^{unbehaved};#phi_{lead jet}", {HistType::kTH2F, {{80, -1., 7.},{80, -1.,7.}}}},
			{"2DHistograms/h_Collision_jet_vs_track_deltas", "|#Delta#phi| in function of |#Delta#eta| between collision leading jet and track;|#Delta#eta|;|#Delta#phi|", {HistType::kTH2F, {{100, 0., 1.4},{100, 0.,7.}}}},
			{"2DHistograms/h_Collision_track_eta_vs_phi", "Collision leading track #phi in function of #eta;#eta_{lead track};#phi_{lead track}", {HistType::kTH2F, {{200, -1., 1.},{200, 0.,7.}}}},
			{"2DHistograms/h_Collision_jet_eta_vs_phi", "Collision leading jet #phi in function of #eta;#eta_{lead jet};#phi_{lead jet}", {HistType::kTH2F, {{200, -1., 1.},{200, 0.,7.}}}},
			// unbehaved tracks data
			{"Unbehaved/h_pt", "Unbehaved tracks within collision pT;#it{p}_{T,lead track}(GeV/#it{c});entries", {HistType::kTH1F, {{leadPtBins, 0., 50.}}}}, //0-200 GeV
			{"Unbehaved/h_eta", "Unbehaved tracks within collision #eta;#eta_{track};entries", {HistType::kTH1F, {{200, -10., 10.}}}},
			{"Unbehaved/h_phi", "Unbehaved tracks within collision #phi;#phi_{track} (rads);entries", {HistType::kTH1F, {{140, -7.0, 7.0}}}},
			{"Unbehaved/h_dcaxy", "Unbehaved tracks within collision DCA_{xy};DCA_{xy} (cm?);entries", {HistType::kTH1F, {{1000, 0., 40.}}}},
			{"Unbehaved/h_dcaz", "Unbehaved tracks within collision DCA_{z};DCA_{z} (cm?);entries", {HistType::kTH1F, {{200, 0., 120.}}}},
			{"Unbehaved/h_time", "Unbehaved tracks within collision ET;t (ns);entries", {HistType::kTH1F, {{200, -1000., 6000.}}}},
			{"Unbehaved/h_TPCNClsFound", "Unbehaved tracks' number of found TPC clusters;nº of clusters;entries", {HistType::kTH1I, {{200, 0., 200.}}}},
			{"Unbehaved/h_TPCChi2NCl", "Unbehaved tracks' #Chi^{2}/cluster for the TPC track segment;#frac{#Chi^{2}}{cluster};entries", {HistType::kTH1F, {{100, 0., 10.}}}},
			{"Unbehaved/h_ITSNCls", "Unbehaved tracks' number of ITS clusters;nº of clusters;entries", {HistType::kTH1I, {{200, 0., 200.}}}},
			{"Unbehaved/h_ITSChi2NCl", "Unbehaved tracks' #Chi^{2}/cluster for the TPC track segment;#frac{#Chi^{2}}{cluster};entries", {HistType::kTH1F, {{100, 0., 10.}}}},
			{"Unbehaved/h_constNum", "Number of constituents of leading jet when leading track wasn't included;nº of tracks;entries", {HistType::kTH1I, {{20, 0., 20.}}}},
			// well behaved tracks data
			{"Behaved/h_pt", "Behaved tracks within collision pT;#it{p}_{T,lead track}(GeV/#it{c});entries", {HistType::kTH1F, {{leadPtBins, 0., 50.}}}}, //0-200 GeV
			{"Behaved/h_eta", "Behaved tracks within collision #eta;#eta_{track};entries", {HistType::kTH1F, {{200, -10., 10.}}}},
			{"Behaved/h_phi", "Behaved tracks within collision #phi;#phi_{track} (rads);entries", {HistType::kTH1F, {{140, -7.0, 7.0}}}},
			{"Behaved/h_dcaxy", "Behaved tracks within collision DCA_{xy};DCA_{xy} (cm?);entries", {HistType::kTH1F, {{1000, 0., 40.}}}},
			{"Behaved/h_dcaz", "Behaved tracks within collision DCA_{z};DCA_{z} (cm?);entries", {HistType::kTH1F, {{200, 0., 120.}}}},
			{"Behaved/h_time", "Behaved tracks within collision ET;t (ns);entries", {HistType::kTH1F, {{200, -1000., 6000.}}}},
			{"Behaved/h_TPCNClsFound", "Behaved tracks' number of found TPC clusters;nº of clusters;entries", {HistType::kTH1I, {{200, 0., 200.}}}},
			{"Behaved/h_TPCChi2NCl", "Behaved tracks' #Chi^{2}/cluster for the TPC track segment;#frac{#Chi^{2}}{cluster};entries", {HistType::kTH1F, {{100, 0., 10.}}}},
			{"Behaved/h_ITSNCls", "Behaved tracks' number of ITS clusters;nº of clusters;entries", {HistType::kTH1I, {{200, 0., 200.}}}},
			{"Behaved/h_ITSChi2NCl", "Behaved tracks' #Chi^{2}/cluster for the TPC track segment;#frac{#Chi^{2}}{cluster};entries", {HistType::kTH1F, {{100, 0., 10.}}}},
			{"Behaved/h_constNum", "Number of constituents of leading jet when leading track was included;nº of tracks;entries", {HistType::kTH1I, {{20, 0., 20.}}}},
			{"h_constituentsNumber", "Number of constituents on each jet;nº of tracks;entries", {HistType::kTH1I, {{20, 0., 20.}}}},
			{"h_DeltaR", "Distance between collision leading track and jet;#DeltaR;entries", {HistType::kTH1I, {{210, 0., 7.}}}}
		}
		
	};
	
	Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
	Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
	
	void init(InitContext const&){}
	
	float vertexZCut = 10.0f;
	float trackPtMin = 0.15;
	float trackPtMax = 1000.;
	float trackEtaMin = -0.9;
	float trackEtaMax = 0.9;
	
	Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
	Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax); // do we need eta cut both here and in globalselection?
	//Filter globalTracks = aod::track::IsGlobalTrack;
	Filter jetCuts = aod::jet::pt > jetPtMin && aod::jet::r == nround(jetR.node() * 100.0f);
	
	int unbe_hasITS = 0;
	int unbe_hasTPC = 0;
	int beha_hasITS = 0;
	int beha_hasTPC = 0;
	int collisionNumber = 0;
	void processDataCharged(aod::Collision const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& jets, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>> const& tracks){ 
		
		
		
		//
		// looping through tracks within the collision
		//
		// initializing values for comparison for leading track within collision
		float collisionLeadPtTrack = -1.; // lead track within collision
		float collisionLeadEtaTrack = 10.;
		float collisionLeadPhiTrack = 10.;
		float collisionLeadDCAxyTrack = -1.;
		float collisionLeadDCAzTrack = -1.;
		float collisionLeadTimeTrack = -1.;
		int collisionLeadTPCNClsTrack = -1.;
		int collisionLeadTPCChi2NClTrack = -1.;
		int collisionLeadITSNClsTrack = -1.;
		int collisionLeadITSChi2NClTrack = -1.;
		bool collisionTrack_hasITS = false;
		bool collisionTrack_hasTPC = false;
		for(auto &track : tracks){ // collecting all tracks from *collision*
			// selecting highest pT track within the collision
			if(collisionLeadPtTrack < track.pt() && track.isGlobalTrackWoPtEta()){
				collisionLeadPtTrack = track.pt();
				collisionLeadEtaTrack = track.eta();
				collisionLeadPhiTrack = track.phi();
				collisionLeadDCAxyTrack = track.dcaXY();
				collisionLeadDCAzTrack = track.dcaZ();
				collisionLeadTimeTrack = track.trackTime();
				collisionLeadTPCNClsTrack = track.tpcNClsFound();
				collisionLeadTPCChi2NClTrack = track.tpcChi2NCl();
				collisionLeadITSNClsTrack = track.itsNCls();
				collisionLeadITSChi2NClTrack = track.itsChi2NCl();
				collisionTrack_hasITS = track.hasITS();
				collisionTrack_hasTPC = track.hasTPC();
			}
		}
		
		// initializing values for comparison
		float jetLeadPtTrack; // lead track within jet
		float jetLeadEtaTrack;
		float jetLeadPhiTrack;
		float leadPtJet = -1.; // lead jet within collision
		float leadEtaJet = -50.;
		float leadPhiJet = -50.;
		int leadConstJet = 0;
		//
		// looping through tracks within each jet
		//
		bool jetsInCollision = false; // were there jets in collision? (avoiding redundant countings)
		for(auto& jet : jets){
			jetsInCollision = true; // if there were jets, jet's loop is read
			// number of current jet constituents
			int constNum = 0;
			cout << "Jets were found.\n";
			
			// initializing values for comparison
			jetLeadPtTrack = -1.; // lead track within jet
			jetLeadEtaTrack = 10.;
			jetLeadPhiTrack = 10.;
			
			// searching for leading track within jet
			for (auto& constituent : jet.tracks_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>()){ // collecting all tracks from *jet*
				constNum++;
				// search for leading track
				if(jetLeadPtTrack < constituent.pt()){
					jetLeadPtTrack = constituent.pt();
					jetLeadEtaTrack = constituent.eta();
					jetLeadPhiTrack = constituent.phi();
				}
			}
			
			// selecting highest pT jet within the collision
			if(leadPtJet < jet.pt()){
				leadPtJet = jet.pt();
				leadEtaJet = jet.eta();
				leadPhiJet = jet.phi();
				leadConstJet = constNum;
			}
			
			cout << "Finishing jet contituents loop...\n";
			registry.fill(HIST("h_constituentsNumber"), constNum);
			registry.fill(HIST("h_leadTrack_pt"), jetLeadPtTrack);
			registry.fill(HIST("h_leadTrack_eta"), jetLeadEtaTrack);
			registry.fill(HIST("h_leadTrack_phi"), jetLeadPhiTrack);
			//registry.fill(HIST("h_jet_vs_track_pt"), jetLeadPtTrack, jet.pt()); // lead track within current jet
		} // end of jets loop
		
		
		
		
		// if jets were found in collision, fill histograms
		if(jetsInCollision){
			registry.fill(HIST("2DHistograms/h_Collision_jet_vs_track_pt"), collisionLeadPtTrack, leadPtJet); // lead collision jet vs lead collision track - pT
			registry.fill(HIST("2DHistograms/h_Collision_jet_vs_track_eta"), collisionLeadEtaTrack, leadEtaJet); // lead collision jet vs lead collision track - eta
			registry.fill(HIST("2DHistograms/h_Collision_jet_vs_track_phi"), collisionLeadPhiTrack, leadPhiJet); // lead collision jet vs lead collision track - phi
			registry.fill(HIST("2DHistograms/h_Collision_track_eta_vs_phi"), collisionLeadEtaTrack, collisionLeadPhiTrack); // lead collision track - eta vs phi
			registry.fill(HIST("2DHistograms/h_Collision_jet_eta_vs_phi"), leadEtaJet, leadPhiJet); // lead collision jet - eta vs phi
			registry.fill(HIST("2DHistograms/h_Collision_jet_vs_track_deltas"), fabs(collisionLeadEtaTrack-leadEtaJet), fabs(collisionLeadPhiTrack-leadPhiJet)); // lead collision jet vs lead collision track - distance
			// calculating distance
			float deltaR = sqrt(pow(collisionLeadEtaTrack-leadEtaJet,2) + pow(collisionLeadPhiTrack-leadPhiJet,2));
			registry.fill(HIST("h_DeltaR"), deltaR);
			
			
			// if it's an unbehaved track
			if(collisionLeadPtTrack > leadPtJet){ //good opportunity for learning to use partitions?
			//if(collisionLeadPtTrack > leadPtJet*1.1){ //at least 10% bigger than the leading jet pT
				registry.fill(HIST("Unbehaved/h_pt"), collisionLeadPtTrack);
				registry.fill(HIST("Unbehaved/h_eta"), collisionLeadEtaTrack);
				registry.fill(HIST("Unbehaved/h_phi"), collisionLeadPhiTrack);
				registry.fill(HIST("Unbehaved/h_dcaxy"), fabs(collisionLeadDCAxyTrack));
				registry.fill(HIST("Unbehaved/h_dcaz"), fabs(collisionLeadDCAzTrack));
				registry.fill(HIST("Unbehaved/h_time"), collisionLeadTimeTrack);
				registry.fill(HIST("Unbehaved/h_TPCNClsFound"), collisionLeadTPCNClsTrack);
				registry.fill(HIST("Unbehaved/h_TPCChi2NCl"), collisionLeadTPCChi2NClTrack);
				registry.fill(HIST("Unbehaved/h_ITSNCls"), collisionLeadITSNClsTrack);
				registry.fill(HIST("Unbehaved/h_ITSChi2NCl"), collisionLeadITSChi2NClTrack);
				registry.fill(HIST("2DHistograms/h_Collision_jet_vs_unbehavedTrack_pt"), collisionLeadPtTrack, leadPtJet); // lead collision jet vs unbehaved collision track - pt
				registry.fill(HIST("2DHistograms/h_Collision_jet_vs_unbehavedTrack_eta"), collisionLeadEtaTrack, leadEtaJet); // lead collision jet vs unbehaved collision track - eta
				registry.fill(HIST("2DHistograms/h_Collision_jet_vs_unbehavedTrack_phi"), collisionLeadPhiTrack, leadPhiJet); // lead collision jet vs unbehaved collision track - phi
				registry.fill(HIST("Unbehaved/h_constNum"), leadConstJet);
				if(collisionTrack_hasITS){ unbe_hasITS++; }
				if(collisionTrack_hasTPC){ unbe_hasTPC++; }
			} else{
				registry.fill(HIST("Behaved/h_pt"), collisionLeadPtTrack);
				registry.fill(HIST("Behaved/h_eta"), collisionLeadEtaTrack);
				registry.fill(HIST("Behaved/h_phi"), collisionLeadPhiTrack);
				registry.fill(HIST("Behaved/h_dcaxy"), fabs(collisionLeadDCAxyTrack));
				registry.fill(HIST("Behaved/h_dcaz"), fabs(collisionLeadDCAzTrack));
				registry.fill(HIST("Behaved/h_time"), collisionLeadTimeTrack);
				registry.fill(HIST("Behaved/h_TPCNClsFound"), collisionLeadTPCNClsTrack);
				registry.fill(HIST("Behaved/h_TPCChi2NCl"), collisionLeadTPCChi2NClTrack);
				registry.fill(HIST("Behaved/h_ITSNCls"), collisionLeadITSNClsTrack);
				registry.fill(HIST("Behaved/h_ITSChi2NCl"), collisionLeadITSChi2NClTrack);
				registry.fill(HIST("Behaved/h_constNum"), leadConstJet);
				if(collisionTrack_hasITS){ beha_hasITS++; }
				if(collisionTrack_hasTPC){ beha_hasTPC++; }
			}
		} // end of jets in collision conditional
		cout << "Collision number: " << collisionNumber << endl;
		collisionNumber++;
		std::cout << "behaved_hasITS = " << beha_hasITS << endl;
		std::cout << "behaved_hasTPC = " << beha_hasTPC << endl;
		std::cout << "unbehaved_hasITS = " << unbe_hasITS << endl;
		std::cout << "unbehaved_hasTPC = " << unbe_hasTPC << endl << endl;
		
	}
	PROCESS_SWITCH(crJetAnalysis, processDataCharged, "jets data", true);
	
	
};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { 
	return WorkflowSpec{adaptAnalysisTask<crJetAnalysis>(cfgc, TaskName{"creckzie-jet-analysis"})}; 
}



