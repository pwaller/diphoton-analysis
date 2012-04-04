#include "event_list.h"

#include <fstream>
#include <iostream>

uint64_t make_runevent(uint32_t run, uint32_t event) {
    return (static_cast<uint64_t>(run) << 32) | event;
}

EventList::EventList(const std::string& path) {
    std::ifstream in(path);
    if (!in.good()) {
        std::cerr << "Bad event list " << path << std::endl;
        abort();
    }
    uint32_t run, event, lb;
    while (in >> run >> event >> lb) {
        const auto runevent = make_runevent(run, event);
        if (_contents.count(runevent)) {
            std::cerr << "Duplicate run/event/lb: " << run << " - " << event << " - " << lb << std::endl;
            abort();
        }
        _contents.insert(runevent);
    }
    std::cout << "Loaded " << _contents.size() << " events" << std::endl;
}

bool EventList::present(uint32_t run, uint32_t event) const {
    return _contents.count(make_runevent(run, event));
}

