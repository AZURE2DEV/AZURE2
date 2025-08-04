#include "../include/SRIMUtilities.h"
#include "../pugixml/src/pugixml.hpp"
#include <iostream>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <algorithm>

SRIMUtilities::SRIMUtilities() : dataLoaded_(false) {
    // Initialize search paths for SRIM2013.xml with cross-platform support
    srimPaths_ = {
        // Development locations
        "../erya/data/SRIM2013.xml",
        "erya/data/SRIM2013.xml",
        
        // Current directory (highest priority for portability)
        "SRIM2013.xml",
        "./SRIM2013.xml",
        
        // Unix-style installation locations
        "/usr/local/share/azure2/SRIM2013.xml",
        "/opt/azure2/SRIM2013.xml", 
        "/usr/share/azure2/SRIM2013.xml",
        
        // macOS application bundle locations
        "../Resources/SRIM2013.xml",                    // Relative to MacOS folder in bundle
        "../../Resources/SRIM2013.xml",                 // Alternative bundle structure
        "../../../Resources/SRIM2013.xml",              // Deep bundle structure
        "/Applications/AZURE2.app/Contents/Resources/SRIM2013.xml",  // Installed bundle
        
        // Windows portable locations
        "./data/SRIM2013.xml",                          // Data subfolder
        "../data/SRIM2013.xml",                         // Data subfolder one level up
        "./erya/SRIM2013.xml",                          // ERYA subfolder
        
        // Platform-specific program locations
#ifdef _WIN32
        "C:/Program Files/AZURE2/SRIM2013.xml",
        "C:/Program Files (x86)/AZURE2/SRIM2013.xml",
#endif
        
        // Homebrew location (macOS)
#ifdef __APPLE__
        "/usr/local/share/azure2/SRIM2013.xml",         // Intel Mac homebrew
        "/opt/homebrew/share/azure2/SRIM2013.xml",      // Apple Silicon homebrew
#endif
    };
    
    // Try to load the SRIM data file
    loadSRIMFile();
}

SRIMUtilities::~SRIMUtilities() {
    elementData_.clear();
}

bool SRIMUtilities::loadSRIMFile() {
    pugi::xml_document doc;
    pugi::xml_parse_result result;
    
    for (const auto& path : srimPaths_) {
        result = doc.load_file(path.c_str());
        if (result) {
            dataLoaded_ = true;
            break;
        }
    }
    
    if (!dataLoaded_) {
        std::cerr << "Warning: Failed to load SRIM2013.xml from any of the searched locations" << std::endl;
        return false;
    }
    
    // Pre-load all element data for faster access
    for (pugi::xml_node tool : doc.child("SRIM2013").child("Ziegler_Data").child("Ziegler_Parameters").children()) {
        int elementNumber = std::stoi(tool.attribute("number").value());
        
        SRIMElementData data;
        data.elementNumber = elementNumber;
        data.elementName = tool.child("Element_Name").child("value").child_value();
        
        try {
            data.atomicMass = manual_stod(tool.child("Atomic_Mass").child("value").child_value());
            data.atomicDensity = manual_stod(tool.child("Atomic_Density").child("value").child_value());
            data.blochParameter = manual_stod(tool.child("Bloch_Parameter").child("value").child_value());
            
            // Read all 12 Ziegler parameters (A-1 through A-12)
            const char* paramNames[] = {"A-1", "A-2", "A-3", "A-4", "A-5", "A-6", 
                                      "A-7", "A-8", "A-9", "A-10", "A-11", "A-12"};
            
            for (int i = 0; i < 12; i++) {
                pugi::xml_node paramNode = tool.child(paramNames[i]);
                if (paramNode) {
                    pugi::xml_node valueNode = paramNode.child("value");
                    if (valueNode) {
                        std::string valueStr = valueNode.child_value();
                        if (!valueStr.empty()) {
                            data.zieglerParams[i] = manual_stod(valueStr);
                        }
                    }
                }
            }
            
            data.isValid = true;
            elementData_[elementNumber] = data;
            
        } catch (const std::exception& e) {
            std::cerr << "Error parsing data for element " << elementNumber << ": " << e.what() << std::endl;
        }
    }
    
    return true;
}

double SRIMUtilities::manual_stod(const std::string& str) {
    size_t i = 0;
    double result = 0.0;
    bool negative = false;

    // Skip leading spaces
    while (i < str.size() && std::isspace(str[i])) {
        i++;
    }

    // Handle optional sign
    if (i < str.size() && (str[i] == '-' || str[i] == '+')) {
        negative = (str[i] == '-');
        i++;
    }

    // Parse integer part
    while (i < str.size() && std::isdigit(str[i])) {
        result = result * 10 + (str[i] - '0');
        i++;
    }

    // Parse fractional part
    if (i < str.size() && str[i] == '.') {
        i++;
        double divisor = 10.0;
        while (i < str.size() && std::isdigit(str[i])) {
            result += (str[i] - '0') / divisor;
            divisor *= 10.0;
            i++;
        }
    }

    // Parse exponential part
    if (i < str.size() && (str[i] == 'e' || str[i] == 'E')) {
        i++;
        bool negExp = false;
        if (i < str.size() && (str[i] == '-' || str[i] == '+')) {
            negExp = (str[i] == '-');
            i++;
        }
        int exponent = 0;
        while (i < str.size() && std::isdigit(str[i])) {
            exponent = exponent * 10 + (str[i] - '0');
            i++;
        }
        result *= std::pow(10.0, negExp ? -exponent : exponent);
    }

    // Skip trailing spaces
    while (i < str.size() && std::isspace(str[i])) {
        i++;
    }

    // If there's any non-space junk left, throw
    if (i != str.size()) {
        throw std::invalid_argument("Invalid character in input string");
    }

    return negative ? -result : result;
}

SRIMElementData SRIMUtilities::readSRIMDataForElement(int elementNumber) {
    SRIMElementData data;
    
    if (!dataLoaded_) {
        return data; // Return invalid data
    }
    
    auto it = elementData_.find(elementNumber);
    if (it != elementData_.end()) {
        return it->second;
    }
    
    return data; // Return invalid data if not found
}

double SRIMUtilities::calculateZieglerStoppingPower(double energy_keV, const SRIMElementData& data) {
    if (!data.isValid || energy_keV <= 0) return 0.0;
    
    // Implementation based on erya/src/Layer.cc EvaluateZiegler function
    // Use high-energy formula for all energies to match equation generation
    double stoppingLow = data.zieglerParams[0] * std::pow(energy_keV, data.zieglerParams[1]) +
                        data.zieglerParams[2] * std::pow(energy_keV, data.zieglerParams[3]);
    
    double stoppingHigh = (data.zieglerParams[4] / std::pow(energy_keV, data.zieglerParams[5])) *
                          std::log((data.zieglerParams[6] / energy_keV) + (data.zieglerParams[7] * energy_keV));
    
    // Compound formula: (Low*High)/(High+Low)
    // Apply unit conversion from SRIM units to AZURE2 units (MeV)
    return ((stoppingLow * stoppingHigh) / (stoppingHigh + stoppingLow)) / 1e+21;
}

double SRIMUtilities::calculateCompoundStoppingPower(double energy_keV, const std::vector<CompoundElement>& elements) {
    if (elements.empty()) return 0.0;
    
    double totalStoppingPower = 0.0;
    double totalStoichiometry = 0.0;
    
    // Calculate total stoichiometry for normalization
    for (const auto& element : elements) {
        totalStoichiometry += element.stoichiometry;
    }
    
    if (totalStoichiometry <= 0.0) return 0.0;
    
    // Calculate weighted stopping power
    for (const auto& element : elements) {
        SRIMElementData data = readSRIMDataForElement(element.elementNumber);
        if (data.isValid) {
            double elementStoppingPower = calculateZieglerStoppingPower(energy_keV, data);
            double weight = element.stoichiometry / totalStoichiometry;
            totalStoppingPower += weight * elementStoppingPower;
        }
    }
    
    return totalStoppingPower;
}

std::string SRIMUtilities::generateAZUREEquation(int elementNumber, std::vector<double>& parameters) {
    SRIMElementData data = readSRIMDataForElement(elementNumber);
    
    if (!data.isValid) {
        return ""; // Return empty string if data is invalid
    }
    
    parameters.clear();
    
    // Generate equation compatible with Equation.cpp
    // Note: Equation.cpp uses ^ for power and ln for natural log
    // Use only the high-energy formula (compound formulation) as it's more general
    // SRIM stopping power needs to be divided by 1e+21 to convert from SRIM units to MeV
    // Expected x is in keV, Equation expects x in MeV
    std::string stoppingLow = "(a0*(x*1e3)^a1+a2*(x*1e3)^a3)";
    std::string stoppingHigh = "(a4/((x*1e3)^a5)*ln(a6/(x*1e3)+a7*(x*1e3)))";

    std::string equation = "(" + stoppingLow + "*" + stoppingHigh + ") / (" + stoppingHigh + "+" + stoppingLow + ") / 1e+21";
    
    // Use actual Ziegler parameters from SRIM (first 8 parameters are typically used)
    for (int i = 0; i < 8; i++) {
        parameters.push_back(data.zieglerParams[i]);
    }
    
    return equation;
}

std::string SRIMUtilities::generateCompoundAZUREEquation(const std::vector<CompoundElement>& elements, std::vector<double>& parameters) {
    if (elements.empty()) {
        return "";
    }
    
    parameters.clear();
    
    // For compounds, we'll create a weighted equation
    // This is more complex and requires dynamic parameter generation
    std::ostringstream equation;
    double totalStoichiometry = 0.0;
    
    // Calculate total stoichiometry
    for (const auto& element : elements) {
        totalStoichiometry += element.stoichiometry;
    }
    
    if (totalStoichiometry <= 0.0) return "";
    
    equation << "(";
    
    int paramIndex = 0;
    bool first = true;
    
    for (const auto& element : elements) {
        SRIMElementData data = readSRIMDataForElement(element.elementNumber);
        if (!data.isValid) continue;
        
        double weight = element.stoichiometry / totalStoichiometry;
        
        if (!first) {
            equation << " + ";
        }
        first = false;
        
        // Add weighted contribution for this element
        equation << weight << "*((a" << paramIndex << "*x^a" << (paramIndex+1) 
                 << "+a" << (paramIndex+2) << "*x^a" << (paramIndex+3) << ")*(a" << (paramIndex+4) 
                 << "/x^a" << (paramIndex+5) << "*ln(a" << (paramIndex+6) << "/x+a" << (paramIndex+7) 
                 << "*x)))/((a" << (paramIndex+4) << "/x^a" << (paramIndex+5) << "*ln(a" << (paramIndex+6) 
                 << "/x+a" << (paramIndex+7) << "*x))+(a" << paramIndex << "*x^a" << (paramIndex+1) 
                 << "+a" << (paramIndex+2) << "*x^a" << (paramIndex+3) << "))";
        
        // Add parameters for this element
        for (int i = 0; i < 8; i++) {
            parameters.push_back(data.zieglerParams[i]);
        }
        
        paramIndex += 8;
    }
    
    equation << ") / 1e+21";  // Convert SRIM units to MeV
    
    return equation.str();
}

std::vector<CompoundElement> SRIMUtilities::parseCompoundFormula(const std::string& formula) {
    std::vector<CompoundElement> elements;
    
    if (formula.empty()) return elements;
    
    size_t i = 0;
    while (i < formula.length()) {
        // Parse element symbol (capital letter followed by optional lowercase letter)
        std::string symbol;
        if (i < formula.length() && std::isupper(formula[i])) {
            symbol += formula[i++];
            if (i < formula.length() && std::islower(formula[i])) {
                symbol += formula[i++];
            }
        } else {
            break; // Invalid format
        }
        
        // Parse stoichiometry number (optional)
        double stoichiometry = 1.0;
        if (i < formula.length() && std::isdigit(formula[i])) {
            std::string numStr;
            while (i < formula.length() && (std::isdigit(formula[i]) || formula[i] == '.')) {
                numStr += formula[i++];
            }
            try {
                stoichiometry = std::stod(numStr);
            } catch (...) {
                stoichiometry = 1.0;
            }
        }
        
        // Convert symbol to atomic number
        int atomicNumber = symbolToAtomicNumber(symbol);
        if (atomicNumber > 0) {
            elements.emplace_back(atomicNumber, stoichiometry);
        }
    }
    
    return elements;
}

int SRIMUtilities::symbolToAtomicNumber(const std::string& symbol) {
    // Map of element symbols to atomic numbers
    static std::map<std::string, int> symbolMap = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8},
        {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32},
        {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
        {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
        {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56},
        {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64},
        {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
        {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
        {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88},
        {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}
    };
    
    auto it = symbolMap.find(symbol);
    return (it != symbolMap.end()) ? it->second : 0;
}

std::vector<int> SRIMUtilities::getAvailableElements() const {
    std::vector<int> elements;
    for (const auto& pair : elementData_) {
        elements.push_back(pair.first);
    }
    std::sort(elements.begin(), elements.end());
    return elements;
}