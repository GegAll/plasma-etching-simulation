# Plasma-Etching-Simulation

This repository contains a parallelized C++ simulation of plasma etching using CF4 gas on a Si wafer, where a mask is placed between the wafer and the plasma. The simulation offers three different mask configurations (straight, inwards, and outwards) and two particle distribution patterns (triangular and asymmetric).

## Features

### Configurations and Distributions
* **Mask Configurations:**
  
  * **Straight:** A simple mask with a straight opening.
  * **Inwards:** A mask with inwards-sloping sides.
  * **Outwards:** A mask with outwards-sloping sides.
  
* **Particle Distributions:**
  * **Triangular Distribution:** Particles are distributed with a triangular pattern.
  * **Asymmetric Distribution:** Particles are distributed asymmetrically to simulate specific etching conditions.
  
## Classes
The code includes two main classes:

* **Grid Class:** Represents the simulation grid containing cells which can be one of the following types:
  * VACUUM
  * MASK
  * SI (Silicon)
  * SIF, SIF2, SIF3: Intermediate etching products
  * SIF4: Final product that transitions to a gas and leaves the substrate
* **Particle Class:** Represents the particles in the simulation that move randomly through the grid. Particles interact with the mask or the substrate and trigger reactions.
  
## Particle Behavior
* **Movement:** Particles move randomly in the grid until they interact with either the mask or the substrate.
* **Interactions:**
  * **Mask:** Particles are redirected downwards upon touching the mask from the sides.
  * **Substrate:** Particles react with a 100% probability, converting Si to SiF, SiF to SiF2, and SiF2 to SiF3. Interaction with SiF3 produces SiF4, which turns the cell into a VACUUM cell.
  
## Parallelization
The simulation leverages multi-threading to parallelize particle movement and interaction, improving performance on multi-core systems.

## Future Improvements
* **Reaction Probabilities:** Incorporate probabilistic reactions based on particle energy and species type.
* **Reflection Logic:** Implement reflection probabilities when particles interact with the mask.
* **Performance Optimization:** Enhance the simulation efficiency to reduce computation times.

## Customization
* **Grid Size:** The grid dimensions can be specified via terminal arguments.
* **Particle Count:** The number of particles per simulation can also be adjusted through terminal inputs.

## Output and Visualization
The simulation outputs a CSV file representing the final grid state, which can be visualized using a Jupyter Notebook. The provided notebook uses Pandas and Matplotlib to generate heatmaps, offering a clear visualization of the etching profile.

## Usage
1. **Compile the Code:** Use a C++ compiler to compile the source code.
2. **Run the Simulation:** Execute the compiled binary with appropriate arguments for grid size and particle count.
3. **Visualization:** Open the generated CSV file in the provided Jupyter Notebook to visualize the results.

## Example Visualization
### Straight configuration
![test_straight_asy03_p50000_v_03_06](https://github.com/GegAll/plasma-etching-simulation/assets/170819708/ee222e0b-35f1-4db6-a873-a7c584e4098c)
### Inwards configuration
![test_inwards_asy03_p50000_v_03_06](https://github.com/GegAll/plasma-etching-simulation/assets/170819708/240e22d8-af21-4836-8cdb-0d773be6f299)
### Outwards configuration
![test_outwards_asy03_p50000_v_03_06](https://github.com/GegAll/plasma-etching-simulation/assets/170819708/6dd75f20-cde4-4fff-a51b-2f586c19c7ee)

## Performance
Simulation times vary based on the number of particles and grid size. While parallelization improves performance, further optimizations are planned to handle larger simulations efficiently.




