#include "include/SRIMTable.h"

#include "pugixml/src/pugixml.hpp"

#include <iostream>

// Ziegler Parameters main constructor
SRIMTable::SRIMTable( int index )
{
    // Read the XML file and fill all the vectors
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file("/opt/SRIM2013.xml");

    for( pugi::xml_node tool : doc.child("SRIM2013").child("Ziegler_Data").child("Ziegler_SRIM_Tables").children( ) ){

        if( std::stoi(tool.attribute("number").value( )) != index ) continue;

        for( pugi::xml_node tool2 : tool.children( ) ){
            Energies.push_back(std::stod(tool2.attribute("Energy").value( )));
            StoppingPower.push_back(std::stod(tool2.attribute("Stopping_Power").value( )));
        }

    }

}

SRIMTable::~SRIMTable( )
{
    Energies.clear( );
    StoppingPower.clear( );
}