/*
 * GSTA (Generalized Smoothed Trajectory Analysis) Implementation for PLUMED
 *
 * The trajectory coordinates can be convolved with an arbitrary kernel function
 * This can be described as a weighted moving average of the coordinates. The kernel 
 * function can be customized via an external file.
 *
 * More details can be found in the following paper:
 * Berta, D., et al. "Nuclear Quantum Effects from the Analysis of Smoothed Trajectories: 
 * Pilot Study for Water" Journal of Chemical Theory and Computation 2020 16 (5), 3316-3334 
 * DOI: https://doi.org/10.1021/acs.jctc.9b00703
 *
 * For questions or feedback, please contact:
 * Dr. Ádám Madarász
 * Email: madarasz.adam@ttk.hu
 */
#include "../core/ActionRegister.h"
#include "../core/Colvar.h"
#include "../tools/Vector.h"
#include <fstream> // Library for file handling
#include <vector>
#include <deque>
#include <string>
#include <sstream>

using namespace PLMD;

class GSTA : public Colvar {
private:
    int windowSize; // Size of the moving average window
    std::string kernelFile = "kernel.dat"; // Default kernel file name
    std::string symFile = "symbol.xyz"; // Default file name for atom symbols
    std::string outputFile = "filtered_trajectory.xyz"; // Default output file name
    std::vector<AtomNumber> atoms; // List of selected atoms
    std::vector<std::deque<Vector>> buffers; // Buffers for each atom for the moving average
    std::vector<std::string> atomSymbols; // Atomic symbols
    std::vector<double> kernel; // Values from the kernel.dat file
    std::vector<double> weights; // Stored weights
    int outputStride = 1; // Default: write every step

    void loadKernel(const std::string& filename); // Load kernel data
    void initializeWeights(); // Initialize weights

public:
    static void registerKeywords(Keywords& keys);
    explicit GSTA(const ActionOptions&);
    void calculate() override;
};

PLUMED_REGISTER_ACTION(GSTA, "GSTA")

void GSTA::registerKeywords(Keywords& keys) {
    Colvar::registerKeywords(keys);

    keys.add("optional", "KERNEL_FILE", "Name of the kernel file to load weights from");
    keys.add("optional", "SYMBOL_FILE", "Name of the XYZ file to read the atom symbols");
    keys.add("optional", "OUTPUT_FILE", "Name of the output XYZ file");

    // Register parameters
    keys.add("atoms", "ATOMS", "The list of atoms to calculate the moving average for");
    // Register components for averaged x, y, z coordinates
    keys.addOutputComponent("dummy", "default", "dummy variable to run moving average");
    keys.add("optional", "OUTPUT_STRIDE", "Frequency of writing to the output file");
}

void GSTA::loadKernel(const std::string& filename) {
    std::ifstream kernelFile(filename);
    if (!kernelFile.is_open()) {
        error("Could not open kernel.dat file to read weights.");
    }

    kernel.clear();
    double value;
    while (kernelFile >> value) {
        kernel.push_back(value);
    }

    kernelFile.close();

    if (kernel.empty()) {
        error("Kernel file is empty.");
    }

    log.printf("Loaded kernel with %lu values.\n", kernel.size());
}

void GSTA::initializeWeights() {
    weights.resize(windowSize);
    size_t kernelSize = kernel.size();
    for (size_t j = 0; j < kernelSize; ++j) {
        weights[kernelSize - j - 1] = kernel[j];
        weights[kernelSize + j - 1] = kernel[j];
    }
    log.printf("Initialized weights with %lu values.\n", weights.size());
}

GSTA::GSTA(const ActionOptions& ao) : PLUMED_COLVAR_INIT(ao) {
    parse("KERNEL_FILE", kernelFile);
    log.printf("Kernel file set to: %s\n", kernelFile.c_str());

    parse("SYMBOL_FILE", symFile);    // Can be overridden in the plumed.dat file
    log.printf("Atom symbols are taken from: %s\n", symFile.c_str());

    parse("OUTPUT_FILE", outputFile);
    log.printf("Output file set to: %s\n", outputFile.c_str());
    parse("OUTPUT_STRIDE", outputStride);
    if (outputStride <= 0) {
        error("OUTPUT_STRIDE must be a positive integer.");
    }
    log.printf("Output stride set to %d.\n", outputStride);
    parseAtomList("ATOMS", atoms);

    // Read atomic symbols
    std::ifstream inputFile(symFile);
    if (!inputFile.is_open()) {
        error("Could not open " + symFile + " to read atom symbols.");
    }

    std::string line;
    int numAtoms = 0;
    inputFile >> numAtoms; // Read the number of atoms from the first line
    std::getline(inputFile, line); // Skip the rest of the first line
    std::getline(inputFile, line); // Skip the second line

    for (int i = 0; i < numAtoms; ++i) {
        std::getline(inputFile, line);
        std::istringstream iss(line);
        std::string symbol;
        iss >> symbol; // The first column is the atomic type
        if (!symbol.empty()) {
            atomSymbols.push_back(symbol);
        }
    }

    inputFile.close();

    log.printf("Loaded %lu atom symbols.\n", atomSymbols.size());

    // Load kernel
    loadKernel(kernelFile);

    // Set the window size
    windowSize = 2 * kernel.size() - 1;
    log.printf("Adjusted window size to %d based on kernel size.\n", windowSize);

    // Initialize weights
    initializeWeights();

    // Initialize buffers for all atoms
    buffers.resize(atoms.size());

    addComponent("dummy");
    log.printf("Registered dummy component.\n");

    // Non-periodic components
    componentIsNotPeriodic("dummy");

    // Request access to atomic data
    requestAtoms(atoms);
}

void GSTA::calculate() {
    static int step = 1 - kernel.size(); // Step counter

    std::ofstream xyzFile(outputFile, std::ios::app); // Use OUTPUT_FILE for writing
    if (!xyzFile.is_open()) {
        log.printf("Error: Could not open %s for writing.\n", outputFile.c_str());
        return;
    }

    for (size_t i = 0; i < atoms.size(); ++i) {
        Vector currentPosition = getPosition(i);

        // Update the buffer with the current position
        buffers[i].push_back(currentPosition);
        if (buffers[i].size() > static_cast<size_t>(windowSize)) {
            buffers[i].pop_front();
        }

        if (buffers[i].size() == static_cast<size_t>(windowSize)) {
          if (step % outputStride == 0) {
            // Calculate weighted average
            Vector weightedAverage(0.0, 0.0, 0.0);
            for (size_t j = 0; j < weights.size(); ++j) {
                weightedAverage += buffers[i][j] * weights[j];
            }

            if (i == 0) {
                // First line of the XYZ format: number of atoms
                xyzFile << atoms.size() << "\n";
                xyzFile << "Time step: " << step << "\n";
            }

            // Update components and write to the XYZ file
            xyzFile << atomSymbols[atoms[i].serial() - 1] << " "
                    << weightedAverage[0] << " "
                    << weightedAverage[1] << " "
                    << weightedAverage[2] << "\n";

            Value* dummyValue = getPntrToComponent("dummy");
            dummyValue->set(0.0); // Always set to 0.0
          }
        }
    }

    xyzFile.close(); // Close the file
    step++; // Increment the step counter
}
	
