
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <numeric>
using namespace std;

// C++ Mersenne Twister 19937 generator
std::random_device rd;
std::mt19937 MTgen(rd());  // Initialize MTgen with a random seed
std::uniform_real_distribution<> R01(0.0, 1.0);
std::normal_distribution<> RMB;

class Grid {
/*
The Grid class creates a discrete grid with width and height parameters,
each cell of the grid can be either

    * Vacuum: Do not interact with the particles
    * Mask: Absorbs all the particles
    * SI: Pure Silicium
    * SIF, SIF2, SIF3: Products of the reaction between Silicium and the particles
    * SiF4 final product which automatically leaves the substrat and turns this cell into vacuum
*/
public:
    enum CellType { VACUUM, MASK, SI, SIF, SIF2, SIF3, SIF4 }; // Defines the kind of material a cell can be

    Grid(int width, int height);
    void initialize(int n);
    void display() const;
    void updateCell(int x, int y, CellType newType);
    CellType getCellType(int x, int y) const;
    int getWidth() const;
    int getHeight() const;
    void exportToCSV(const string& filename) const;

private:
    int width, height;
    vector<vector<CellType>> grid;
    mutable std::mutex mtx; // Added for thread-safe access
};

Grid::Grid(int width, int height) : width(width), height(height) {
    grid.resize(height, vector<CellType>(width, VACUUM));
}

void Grid::initialize(int n) { // Initializes a new grid with either a straight, inwards, or outwards structure

    // Setting the boundaries for the mask structure
    int si_boundary = round(0.3 * height);
    int upper_boundary = round(0.6 * height);
    int lower_boundary = round(0.4 * height);
    int left_boundary = width / 3 - (upper_boundary - lower_boundary);
    int right_boundary = 2 * width / 3 + (upper_boundary - lower_boundary);

    switch (n) {
    case 0:
        // Straight case
        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                if (y < si_boundary) { // The first rows will be Si
                    grid[y][x] = SI;
                } else if (y > lower_boundary && y < upper_boundary) { // These rows will be a mask with a hole in the middle
                    grid[y][x] = (x >= width / 3 && x < 2 * width / 3) ? VACUUM : MASK; // The hole has a size of a third of the total width
                } else {
                    grid[y][x] = VACUUM; // The rest of the cells will be vacuum
                }
            }
        }
        break;

    case 1:
        // Inwards case
        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                if (y < si_boundary) {
                    grid[y][x] = SI; // The first rows will be Si
                } else if (y > lower_boundary && y < upper_boundary) { // After some vacuum the next rows will be part of the mask
                    if (x > left_boundary && x < width / 3) {
                        if (x <= y - lower_boundary + left_boundary) { // Build the triangle on the left side
                            grid[y][x] = MASK;
                        } else {
                            grid[y][x] = VACUUM;
                        }
                    } else if (x > 2 * width / 3 && x < right_boundary) { // Build the triangle on the right side
                        if (x >= -y + 2 * width / 3 + upper_boundary) {
                            grid[y][x] = MASK;
                        } else {
                            grid[y][x] = VACUUM;
                        }
                    } else if (x >= width/3 && x <= 2*width/3){ // The rest between these two boundaries is set to VACUUM
                        grid[y][x] = VACUUM;
                    } else {
                        grid[y][x] = MASK;
                    }
                } else { // The rest of the structure is set to VACUUM
                    grid[y][x] = VACUUM;
                }
            }
        }
        break;

    case 2:
        // Outwards case
        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                if (y < si_boundary) {
                    grid[y][x] = SI; // The first rows will be Si
                } else if (y > lower_boundary && y < upper_boundary) { // After some vacuum the next rows will be part of the mask
                    if (x > left_boundary && x < width / 3) {
                        if (x <= -y + lower_boundary + width/3) { // Build the triangle on the left side
                            grid[y][x] = MASK;
                        } else {
                            grid[y][x] = VACUUM;
                        }
                    } else if (x > 2 * width / 3 && x < right_boundary) { // Build the triangle on the right side
                        if (x >= y - lower_boundary + 2 * width / 3) {
                            grid[y][x] = MASK;
                        } else {
                            grid[y][x] = VACUUM;
                        }
                    } else if (x >= width/3 && x <= 2*width/3){ // The rest between these two boundaries is set to VACUUM
                        grid[y][x] = VACUUM;
                    } else { // The rest is set to MASK
                        grid[y][x] = MASK;
                    }
                } else { // The rest of the structure is Vacuum
                    grid[y][x] = VACUUM;
                }
            }
        }
        break;
    }
}


void Grid::display() const { // Show the grid for future exporting
    std::lock_guard<std::mutex> lock(mtx); // Added for thread-safe access
    for (const auto& row : grid) {
        for (const auto& cell : row) {
            switch (cell) {
                case VACUUM: cout << " "; break;
                case MASK: cout << "M"; break;
                case SI: cout << "S"; break;
                case SIF: cout << "1"; break;
                case SIF2: cout << "2"; break;
                case SIF3: cout << "3"; break;
                case SIF4: cout << "4"; break;
            }
        }
        cout << "\n";
    }
}

void Grid::updateCell(int x, int y, CellType newType) {
    std::lock_guard<std::mutex> lock(mtx); // Added for thread-safe access
    grid[y][x] = newType; // Changes the type of a cell
}

Grid::CellType Grid::getCellType(int x, int y) const {
    std::lock_guard<std::mutex> lock(mtx); // Added for thread-safe access
    return grid[y][x]; // Getter for the cell type
}

int Grid::getWidth() const { // Getter for the width
    return width;
}

int Grid::getHeight() const { // Getter for the height
    return height;
}

void Grid::exportToCSV(const std::string& filepath) const {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file for writing\n";
        return;
    }

    std::lock_guard<std::mutex> lock(mtx); // Added for thread-safe access
    for (const auto& row : grid) {
        for (size_t i = 0; i < row.size(); ++i) {
            switch (row[i]) {
                case VACUUM: file << "VACUUM"; break;
                case MASK: file << "MASK"; break;
                case SI: file << "SI"; break;
                case SIF: file << "SIF"; break;
                case SIF2: file << "SIF2"; break;
                case SIF3: file << "SIF3"; break;
                case SIF4: file << "SIF4"; break;
            }
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    std::cout << "Grid exported to " << filepath << "\n";
}


class Particle {
/*
The Particle class creates a particle which is initially created on the top of the grid and randomly moves down, left or right. 
Its parameters are the position in x and y, and its energy.

It interacts with the class Grid resulting in chemical reactions with the purpose of etching the Silicium substrat to create an etching profile. Now there is no
probability involved in the particle and the Si cells reacting, they always react. In the future this process should be controlled by a MonteCarlo method which
decides if the particle reacts or not depending on its energy.
*/
public:
    Particle(int x, int y, double energy, Grid& grid);
    void move();
    bool isDeleted() const;

private:
    int x, y;
    double energy;
    Grid& grid;
    bool deleted;

    bool moveStep();
    bool interact();
    void reflect();
    void react();
};

Particle::Particle(int x, int y, double energy, Grid& grid)
    : x(x), y(y), energy(energy), grid(grid), deleted(false) {}

bool Particle::isDeleted() const {
    return deleted;
}

void Particle::move() { // The whole functionality of the particle lies on this function, in which direction the particle moves and if it reacts with the substrat
    while (moveStep()) {
        if (interact()) {
            deleted = true; // Deletes the particle if it interacts with its sorroundings
            break;
        }
    }
}

bool Particle::moveStep() {
    // Generate a random value to decide the direction of movement
    double randomValue = R01(MTgen);

    if (randomValue < 0.3) {
        if (y > 0) {
            --y; // Move down
        }
    } else if (randomValue >= 0.3 && randomValue < 0.6) {
        x = (x + 1) % grid.getWidth(); // Move right with periodic boundary
    } else if (randomValue >= 0.6) {
        x = (x - 1 + grid.getWidth()) % grid.getWidth(); // Move left with periodic boundary
    }
    return y > 0;
}

bool Particle::interact() {
    Grid::CellType cell = grid.getCellType(x, y);
    if (cell == Grid::MASK) {
        deleted = true; // Delete the particle if it touches the mask
        return true;
    } else if (cell == Grid::SI || cell == Grid::SIF || cell == Grid::SIF2 || cell == Grid::SIF3) {
        react();
        return true;
    }
    return false;
}

void Particle::reflect() {
    // Reflection of the particle, still not in use in the main function
    y = max(0, y - 2);
}

void Particle::react() {
    // Easy example of the reaction between the particles and the grid
    Grid::CellType cell = grid.getCellType(x, y);
    if (cell == Grid::SI) {
        grid.updateCell(x, y, Grid::SIF);
    } else if (cell == Grid::SIF) {
        grid.updateCell(x, y, Grid::SIF2);
    } else if (cell == Grid::SIF2) {
        grid.updateCell(x, y, Grid::SIF3);
    } else if (cell == Grid::SIF3) {
        grid.updateCell(x, y, Grid::VACUUM);
    }
}

void create_n_particles(int width, int height, double initialEnergy, Grid& grid, vector<unique_ptr<Particle>> particles, int n) {
    particles.reserve(n); // Pre-allocate memory for particles
    for (int i = 0; i < n; ++i) {
        particles.push_back(make_unique<Particle>(i, height - 1, initialEnergy, grid));
    }
}

void particle_simulation(int start, int end, double initialEnergy, Grid& grid, std::vector<int>& distribution) {
    /*
    This function helps to preprocess the particles and improve the performance of the program by deleting the particles after interacting with the grid
    */
    for (int i = start; i < end; ++i) {
        int x = i % grid.getWidth();
        if (distribution[x] > 0) {
            Particle particle(x, grid.getHeight() - 1, initialEnergy, grid);
            particle.move();
            distribution[x]--; // Decrease the count of particles to be created in this column
        }
    }
}
//------------------------------------------------------//
//------------Distribution functions--------------------//
//------------------------------------------------------//

std::vector<int> asymmetric_distribution(int n, int amount, double portion, bool left) {
    /*
    Creates an asymmetric distribution vector of size n with amount numbers in a portion of the array.
    If the boolean left is true, then only the left side of the vector is filled, otherwise the right side is filled.
    Exmaple of assymetric_distribution(10, 5, 3, true) = {5, 5, 5, 0, 0, 0, 0, 0, 0, 0}
    */
    std::vector<int> distribution(n);
    if (left) {
        for (int i = 0; i < round(n*portion); ++i) {
            distribution[i] = amount;
        }
    } else {
        for (int i = 0; i > round(n*portion); ++i) {
            distribution[n-1-i] = amount;
        }

    }
    return distribution;
}


std::vector<int> triangle_distribution(int n) {
    /*
    Create a vector with a triangular distribution of size n. For example for a vector with 5 elements would look like {1, 2, 3, 2, 1}
    */
    std::vector<int> distribution(n);
    
    if (n % 2 == 0) {
        // For even n
        for (int i = 0; i < n / 2; ++i) {
            distribution[i] = i;
            distribution[n - 1 - i] = i;
        }
    } else {
        // For odd n
        for (int i = 0; i <= n / 2; ++i) {
            distribution[i] = i;
            distribution[n - 1 - i] = i;
        }
    }
    
    return distribution;
}



int main() {

    cout << "Enter the width of the grid : ";
    int width {};  // Initializing the width of the grid
    cin >> width;

    cout << "Enter the height of the grid : ";
    int height {}; // Initializing the height of the grid
    cin >> height;

    int p = 100; // Number of particles per column
    double initialEnergy = 100.0; // Example initial energy of the particles

    Grid grid(width, height);
    
    cout << "Select the grid configuration : \n 0 : Straight \n 1 : Inwards \n 2 : Outwards" << endl;
    int config {};
    cin >> config;

    cout << "Initializing the grid..." << endl;
    grid.initialize(config); // Create the grid

    cout << "Set the amount of particles created in each cell : ";
    int amount {};
    cin >> amount;

    cout << "initializing the particles" << endl;
    // Choose and generate a distribution
    std::vector<int> distribution = asymmetric_distribution(width, amount, 0.3, true);
    // std::vector<int> distribution = triangle_distribution(width); // Alternatively, use the triangle distribution


    int num_threads = std::thread::hardware_concurrency(); // Added for multi-threading
    vector<thread> threads; // Added for multi-threading


    int total_particles = std::accumulate(distribution.begin(), distribution.end(), 0);
    int particles_per_thread = total_particles / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        int start = i * particles_per_thread;
        int end = (i == num_threads - 1) ? total_particles : start + particles_per_thread;
        threads.emplace_back(particle_simulation, start, end, initialEnergy, std::ref(grid), std::ref(distribution));
    }

    for (auto& t : threads) { // Added for multi-threading
        t.join();
    }

    cout << "Displaying the grid" << endl;
    grid.display();

    cout << "Exporting the grid to a CSV file" << endl;
    string filepath = "C:\\Users\\Asus\\Desktop\\etching_simulations\\data\\testing_straight.csv";
    grid.exportToCSV(filepath);


    return 0;
}