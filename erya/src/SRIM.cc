#include "include/SRIM.h"

#include "pugixml/src/pugixml.hpp"

#include <iostream>

// Ziegler Parameters main constructor
SRIM::SRIM( )
{
    // Read the XML file and fill all the vectors
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file("/opt/SRIM2013.xml");

    int elementNumber;
    std::string elementName;
    double atomicMass, atomicDensity, blochParam;
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12;
    for( pugi::xml_node tool : doc.child("SRIM2013").child("Ziegler_Data").child("Ziegler_Parameters").children( ) ){

        elementNumber = std::stoi(tool.attribute("number").value( ));
        elementName = tool.child("Element_Name").child_value("value");
        atomicMass = std::stod(tool.child("Atomic_Mass").child_value("value"));
        atomicDensity = std::stod(tool.child("Atomic_Density").child_value("value"));
        blochParam = std::stod(tool.child("Bloch_Parameter").child_value("value"));

        a1 = std::stod(tool.child("A-1").child_value("value"));
        a2 = std::stod(tool.child("A-2").child_value("value"));
        a3 = std::stod(tool.child("A-3").child_value("value"));
        a4 = std::stod(tool.child("A-4").child_value("value"));
        a5 = std::stod(tool.child("A-5").child_value("value"));
        a6 = std::stod(tool.child("A-6").child_value("value"));
        a7 = std::stod(tool.child("A-7").child_value("value"));
        a8 = std::stod(tool.child("A-8").child_value("value"));
        a9 = std::stod(tool.child("A-9").child_value("value"));
        a10 = std::stod(tool.child("A-10").child_value("value"));
        a11 = std::stod(tool.child("A-11").child_value("value"));
        a12 = std::stod(tool.child("A-12").child_value("value"));

        Elements.push_back(elementName);
        Values1.push_back(a1);
        Values2.push_back(a2);
        Values3.push_back(a3);
        Values4.push_back(a4);
        Values5.push_back(a5);
        Values6.push_back(a6);
        Values7.push_back(a7);
        Values8.push_back(a8);
        Values9.push_back(a9);
        Values10.push_back(a10);
        Values11.push_back(a11);
        Values12.push_back(a12);
        Mass.push_back(atomicMass);
        Density.push_back(atomicDensity);
        Bloch.push_back(blochParam);
        Tables.push_back(SRIMTable(elementNumber));

    }

}

SRIM::~SRIM( ){
    Elements.clear( );
    Values1.clear( );
    Values2.clear( );
    Values3.clear( );
    Values4.clear( );
    Values5.clear( );
    Values6.clear( );
    Values7.clear( );
    Values8.clear( );
    Values9.clear( );
    Values10.clear( );
    Values11.clear( );
    Values12.clear( );
    Mass.clear( );
    Density.clear( );
    Bloch.clear( );
    Tables.clear( );
}