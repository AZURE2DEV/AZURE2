#ifndef SRIM_UTILITIES_H
#define SRIM_UTILITIES_H

#include <string>
#include <vector>
#include <map>

/*!
 * Structure to hold SRIM element data for a single element
 */
struct SRIMElementData {
    int elementNumber;
    std::string elementName;
    double zieglerParams[12];  // A-1 through A-12
    double atomicMass;
    double atomicDensity;
    double blochParameter;
    bool isValid;
    
    SRIMElementData() : elementNumber(0), isValid(false) {
        for(int i = 0; i < 12; i++) zieglerParams[i] = 0.0;
        atomicMass = 0.0;
        atomicDensity = 0.0;
        blochParameter = 0.0;
    }
};

/*!
 * Structure to represent a compound element with its stoichiometry
 */
struct CompoundElement {
    int elementNumber;
    double stoichiometry;  // Relative proportion in the compound
    
    CompoundElement(int z, double stoich) : elementNumber(z), stoichiometry(stoich) {}
};

/*!
 * Class for SRIM stopping power utilities
 * This class provides access to SRIM stopping power data and calculations
 */
class SRIMUtilities {
public:
    /*!
     * Constructor - loads SRIM data from XML file
     */
    SRIMUtilities();
    
    /*!
     * Destructor
     */
    ~SRIMUtilities();
    
    /*!
     * Read SRIM data for a single element
     * @param elementNumber Atomic number (Z) of the element
     * @return SRIMElementData structure with element data
     */
    SRIMElementData readSRIMDataForElement(int elementNumber);
    
    /*!
     * Calculate Ziegler stopping power for a single element
     * @param energy_keV Energy in keV
     * @param data SRIM element data
     * @return Stopping power in MeV⋅cm²/mg
     */
    double calculateZieglerStoppingPower(double energy_keV, const SRIMElementData& data);
    
    /*!
     * Calculate stopping power for a compound
     * @param energy_keV Energy in keV  
     * @param elements Vector of compound elements with stoichiometry
     * @return Weighted stopping power in MeV⋅cm²/mg
     */
    double calculateCompoundStoppingPower(double energy_keV, const std::vector<CompoundElement>& elements);
    
    /*!
     * Generate AZURE2-compatible equation string for single element
     * @param elementNumber Atomic number of the element
     * @param parameters Output vector to store equation parameters
     * @return Equation string compatible with Equation.cpp
     */
    std::string generateAZUREEquation(int elementNumber, std::vector<double>& parameters);
    
    /*!
     * Generate AZURE2-compatible equation string for compound
     * @param elements Vector of compound elements with stoichiometry
     * @param parameters Output vector to store equation parameters  
     * @return Equation string compatible with Equation.cpp
     */
    std::string generateCompoundAZUREEquation(const std::vector<CompoundElement>& elements, std::vector<double>& parameters);
    
    /*!
     * Parse compound formula (e.g., "CH4", "SiO2") into compound elements
     * @param formula Chemical formula string
     * @return Vector of compound elements with parsed stoichiometry
     */
    std::vector<CompoundElement> parseCompoundFormula(const std::string& formula);
    
    /*!
     * Check if SRIM data is available and loaded
     * @return true if data is loaded, false otherwise
     */
    bool isDataLoaded() const { return dataLoaded_; }
    
    /*!
     * Get list of available elements
     * @return Vector of atomic numbers for which data is available
     */
    std::vector<int> getAvailableElements() const;

private:
    /*!
     * Find and load SRIM XML file from various locations
     * @return true if file was loaded successfully
     */
    bool loadSRIMFile();
    
    /*!
     * Manual string to double conversion to handle locale issues
     * @param str String to convert
     * @return Parsed double value
     */
    double manual_stod(const std::string& str);
    
    /*!
     * Parse element symbol to atomic number
     * @param symbol Element symbol (e.g., "C", "Si", "O")
     * @return Atomic number, or 0 if not found
     */
    int symbolToAtomicNumber(const std::string& symbol);
    
    // Member variables
    std::map<int, SRIMElementData> elementData_;  // Cache of loaded element data
    bool dataLoaded_;                             // Flag indicating if XML data was loaded
    std::vector<std::string> srimPaths_;         // Search paths for SRIM2013.xml
};

#endif // SRIM_UTILITIES_H