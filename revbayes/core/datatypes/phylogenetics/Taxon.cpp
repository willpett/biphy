#include "Taxon.h"

using namespace RevBayesCore;



/**
 * Constructor simply initiating the object and its members.
 *
 * \param[in]    n     The name of the taxon.
 */
Taxon::Taxon(const std::string &n, const std::string &sn) :
    date(  ),
    name( n ),
    speciesName( sn )
{
    
}


/**
 * Get the date info for this taxon.
 *
 * \return    The date.
 */
const TimeAndDate& Taxon::getDate( void ) const
{
    return date;
}


/**
 * Get the name info for this taxon.
 *
 * \return    The name.
 */
const std::string& Taxon::getName( void ) const
{
    return name;
}


/**
 * Get the species name for this taxon.
 *
 * \return    The species name.
 */
const std::string& Taxon::getSpeciesName( void ) const
{
    return speciesName;
}


/**
 * Set the date info for this taxon.
 *
 * \param[in]    d     The date.
 */
void Taxon::setDate( const TimeAndDate &d )
{
    date = d;
}


/**
 * Set the name info for this taxon.
 *
 * \param[in]    n     The name.
 */
void Taxon::setName( const std::string &n )
{
    name = n;
}


/**
 * Set the species name for this taxon.
 *
 * \param[in]    n     The species name.
 */
void Taxon::setSpeciesName( const std::string &sn )
{
    speciesName = sn;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const Taxon& x) {
    o << x.getName() << ":" << x.getSpeciesName() << ":" << x.getDate();    
    return o;
}
