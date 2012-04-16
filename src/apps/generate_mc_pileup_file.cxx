#include <iostream>

#include <a4/application.h>

#include <a4/atlas/ntup/photon/Event.pb.h>
namespace ntup = a4::atlas::ntup::photon;

using ntup::Event;

using namespace a4::process;
//using namespace a4::hist;


//#include "a4analysis.h"
#include "external.h"

class PileupProcessor : public ProcessorOf<Event> {
public: 

    virtual void process(const Event& event) {
        pileup_tool->Fill(event.run_number(), 
                         event.mc_channel_number(), 
                         event.mc_event_weight(), 
                         event.averageintperxing());
    }

    virtual ~PileupProcessor() {
        pileup_tool->WriteToFile(output_name);
    }
    
    shared<Root::TPileupReweighting> pileup_tool;
    std::string output_name;
};

class PileupConfiguration : public ConfigurationOf<PileupProcessor> {
public:
    std::string output_name;
      
    
    virtual void add_options(po::options_description_easy_init opt) {
        opt("output-name,O", po::value(&output_name), "output filename");
    }

    virtual void read_arguments(po::variables_map& arguments) {
    }
    
    virtual void setup_processor(PileupProcessor& g) {
        // Rough check against multiple runs
        static int count = 0;
        count++;
        assert(count == 1);
        
        auto pileup_tool = new Root::TPileupReweighting("pileup_reweighting");
        pileup_tool->UsePeriodConfig("MC11c");
        pileup_tool->initialize();
        g.pileup_tool.reset(pileup_tool);
        g.output_name = output_name;
    }
};

int main(int argc, const char** argv) {
    return a4_main_configuration<PileupConfiguration>(argc, argv);
}

