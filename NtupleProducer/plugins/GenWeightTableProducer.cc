// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

#include <algorithm>

class GenWeightTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit GenWeightTableProducer(const edm::ParameterSet&);
        ~GenWeightTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;


        std::string name_;
	edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;
};

GenWeightTableProducer::GenWeightTableProducer(const edm::ParameterSet& iConfig) :
    name_(iConfig.getParameter<std::string>("name")),
    genEvtInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<nanoaod::FlatTable>();
    std::cout << "Produce gen weight table!" << std::endl;
}

GenWeightTableProducer::~GenWeightTableProducer() { }

// ------------ method called for each event  ------------
    void GenWeightTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    std::cout << "Gen Weight produce!" << std::endl;
    
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(genEvtInfo_, genEvtInfo);

    // create the table
    std::vector<float> genWeight; // make a vector of size 1 ... :( not memory-friendly but voila
    genWeight.resize(1);
    genWeight[0] = genEvtInfo->weight();

    std::cout << "gen weight = " << genEvtInfo->weight() << std::endl;
    auto out = std::make_unique<nanoaod::FlatTable>(1, name_, true);
    //auto out = std::make_unique<nanoaod::FlatTable>(1, name_, false, true);

    out->addColumn<float>("genWeight", genWeight, "");

    // save to the event branches
    iEvent.put(std::move(out));

   
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenWeightTableProducer);
