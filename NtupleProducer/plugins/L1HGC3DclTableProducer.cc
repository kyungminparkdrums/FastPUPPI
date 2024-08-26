// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterEgID.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterIdentificationBase.h"

#include <algorithm>

class L1HGC3DclTableProducer : public edm::stream::EDProducer<>  {
    public:
        explicit L1HGC3DclTableProducer(const edm::ParameterSet&);
        ~L1HGC3DclTableProducer();

    private:
        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;


        std::string name_;
        edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> clusters_;
        StringCutObjectSelector<l1t::HGCalMulticluster> sel_;
        l1tpf::HGC3DClusterEgID emVsPionID_;
        l1tpf::HGC3DClusterEgID emVsPUID_;
        l1tpf::HGC3DClusterID multiClassPID_;
   	std::unique_ptr<HGCalTriggerClusterIdentificationBase> id_;
};

L1HGC3DclTableProducer::L1HGC3DclTableProducer(const edm::ParameterSet& iConfig) :
    name_(iConfig.getParameter<std::string>("name")),
    clusters_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    sel_(iConfig.getParameter<std::string>("cut"), true),
    emVsPionID_(iConfig.getParameter<edm::ParameterSet>("emVsPionID")),
    emVsPUID_(iConfig.getParameter<edm::ParameterSet>("emVsPUID")),
    multiClassPID_(iConfig.getParameter<edm::ParameterSet>("multiClassPID"))	
{
    produces<nanoaod::FlatTable>();

    if (!emVsPionID_.method().empty()) {
        emVsPionID_.prepareTMVA();
    }
    if (!emVsPUID_.method().empty()) {
        emVsPUID_.prepareTMVA();
    }

    id_ = std::unique_ptr<HGCalTriggerClusterIdentificationBase>{
        HGCalTriggerClusterIdentificationFactory::get()->create("HGCalTriggerClusterIdentificationBDT")};
    id_->initialize(iConfig.getParameter<edm::ParameterSet>("EGIdentification"));


}

L1HGC3DclTableProducer::~L1HGC3DclTableProducer() { }

// ------------ method called for each event  ------------
    void
L1HGC3DclTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    edm::Handle<l1t::HGCalMulticlusterBxCollection> clusters;
    iEvent.getByToken(clusters_, clusters);

    std::vector<const l1t::HGCalMulticluster *> selected;


    for (const l1t::HGCalMulticluster &cl : *clusters)
        if(sel_(cl)) selected.push_back(&cl);


    // create the table
    unsigned int ncands = selected.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

    // Binary classifier legacy BDTs
    std::vector<float> vals_puid, vals_pfemid, vals_egemid;
    std::vector<bool> pass_puid, pass_pfemid; 
    vals_puid.resize(ncands);
    vals_pfemid.resize(ncands);
    vals_egemid.resize(ncands);
    pass_puid.resize(ncands); 
    pass_pfemid.resize(ncands);

    // Multi-classifier BDT
    std::vector<float> vals_multiClassMaxScore, vals_multiClassPUid, vals_multiClassPionid, vals_multiClassEmid;

    vals_multiClassMaxScore.resize(ncands);
    vals_multiClassPUid.resize(ncands);
    vals_multiClassPionid.resize(ncands);
    vals_multiClassEmid.resize(ncands);

    for (unsigned int i = 0; i < ncands; ++i) {
        auto cl3d = *selected[i];

        l1t::PFCluster cluster;
	// Binary classifier Legacy BDTs
        bool passEmVsPU = false;
        bool passPFEmVsPion = false;
        if (!emVsPUID_.method().empty()) {
            passEmVsPU = emVsPUID_.passID(cl3d, cluster);
        }
        if (!emVsPionID_.method().empty()) {
            passPFEmVsPion = emVsPionID_.passID(cl3d, cluster);
        }

        vals_puid[i] = cluster.egVsPUMVAOut();
        vals_pfemid[i] = cluster.egVsPionMVAOut();
        vals_egemid[i] = id_->value(cl3d);
        pass_puid[i] = passEmVsPU;
        pass_pfemid[i] = passPFEmVsPion;

	// Multi-class BDT
	float maxScore = multiClassPID_.evaluate(cl3d, cluster);

        vals_multiClassMaxScore[i] = maxScore;
        vals_multiClassPUid[i] = cluster.puIDScore();
 	vals_multiClassPionid[i] = cluster.piIDScore();
	vals_multiClassEmid[i] = cluster.emIDScore();

    }
    out->addColumn<float>("pfPuIdScore", vals_puid, "");
    out->addColumn<float>("pfEmIdScore", vals_pfemid, "");
    out->addColumn<float>("egEmIdScore", vals_egemid, "");
    out->addColumn<bool>("pfPuIdPass", pass_puid, "");
    out->addColumn<bool>("pfEmIdPass", pass_pfemid, "");

    out->addColumn<float>("multiClassMaxScore", vals_multiClassMaxScore, ""); 
    out->addColumn<float>("multiClassPuIdScore", vals_multiClassPUid, ""); 
    out->addColumn<float>("multiClassPionIdScore", vals_multiClassPionid, ""); 
    out->addColumn<float>("multiClassEmIdScore", vals_multiClassEmid, ""); 

    // save to the event branches
    iEvent.put(std::move(out));
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1HGC3DclTableProducer);
