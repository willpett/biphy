#include "MoveSchedule.h"

MoveSchedule::MoveSchedule(const std::vector<Move*> &m ) : moves( m ) {
    
}


MoveSchedule::~MoveSchedule() {
    // we own nothing
}
