//
//  TreeChangeEventHandler.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/27/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "TreeChangeEventHandler.h"
#include "TreeChangeEventListener.h"

#include <iostream>

TreeChangeEventHandler::TreeChangeEventHandler(void) :
    listeners()
{

}

TreeChangeEventHandler::TreeChangeEventHandler(const TreeChangeEventHandler &h)  :
    listeners(h.listeners)
{
    
}


TreeChangeEventHandler::~TreeChangeEventHandler(void) 
{
    if ( listeners.empty() == false )
    {
        //std::cerr << "Deleting handler while " << listeners.size() << " listener are still listening." << std::endl;
    }
}


TreeChangeEventHandler& TreeChangeEventHandler::operator=(const TreeChangeEventHandler &h)
{
    
    if ( this != &h )
    {
        listeners = h.listeners;
    }
    
    return *this;
}


void TreeChangeEventHandler::addListener(TreeChangeEventListener *l) 
{

    listeners.insert( l );

}

void TreeChangeEventHandler::clear(void) 
{

    listeners.clear();

}


void TreeChangeEventHandler::fire(const TopologyNode &n) 
{

    for (std::set<TreeChangeEventListener*>::iterator it = listeners.begin(); it != listeners.end(); ++it) 
    {
        TreeChangeEventListener *l = *it;
        l->fireTreeChangeEvent( n );
    }
}


const std::set<TreeChangeEventListener*>& TreeChangeEventHandler::getListeners( void ) const 
{

    return listeners;
}


void TreeChangeEventHandler::removeListener(TreeChangeEventListener *l) 
{

    std::set<TreeChangeEventListener*>::iterator pos = listeners.find( l );
    if ( pos != listeners.end() ) 
    {
        listeners.erase( l );
    }
    
}
