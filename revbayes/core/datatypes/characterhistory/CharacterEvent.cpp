//
//  CharacterEvent.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 8/6/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//


#include "CharacterEvent.h"
#include <iostream>

using namespace RevBayesCore;


CharacterEvent::CharacterEvent(void)
{
    
}

CharacterEvent::CharacterEvent(size_t i, unsigned int s, double t) :  index(i), state(s), time(t)
{
    
}

CharacterEvent::CharacterEvent(const CharacterEvent& c) : index(c.index), state(c.state), time(c.time)
{
    
}

CharacterEvent::~CharacterEvent(void)
{

}

CharacterEvent* CharacterEvent::clone( void ) const
{
    return new CharacterEvent( *this );
}

bool CharacterEvent::operator<(const CharacterEvent& rhs) const
{
    return time < rhs.time;
}

double CharacterEvent::getTime(void)
{
    return time;
}

size_t CharacterEvent::getIndex(void)
{
    return index;
}

unsigned int CharacterEvent::getState(void)
{
    return state;
}

void CharacterEvent::setState(unsigned int s)
{
    state = s;
}

void CharacterEvent::print(void)
{
    std::cout << index << " " << state << " " << time << "\n";
}
