#ifndef _EVENT_LIST_H_
#define _EVENT_LIST_H_

#include <unordered_set>

#include <a4/types.h>

class EventList
{
private:
    std::unordered_set<uint64_t> _contents;
    
public:
    EventList(const std::string& path);
    bool present(uint32_t run, uint32_t event) const;
};

#endif

