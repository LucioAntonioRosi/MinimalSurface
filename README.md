# Minimal Surface Project

The Minimal Surface Project is a comprehensive exploration of minimal surfaces—geometric surfaces that locally minimize area. This project provides both a command-line interface for generating minimal surface graphs and Python scripts for interactive plotting and visualization.

## Overview

Minimal surfaces are fascinating structures with applications in mathematics, physics, and engineering. In this project you can:

- **Compute and visualize minimal surfaces** using a dedicated command-line executable.
- **Generate plots** of the minimal surface equation along with its corresponding graph.
- **Visualize the parametric boundary** of the minimal surface.
- **Automate build and cleanup processes** with a Makefile.

## Features

- **Command-Line Interface:**  
  Generate minimal surface graphs by running the executable with customizable parameters.
  
- **Python Plotting Scripts:**  
  Visualize minimal surfaces with:
  - The implicit equation and its graph.
  - A detailed view showing the parametric boundary.

- **Makefile:**  
  Easily build and clean the project with simple `make` commands.

## Requirements

- **Compilation/Execution:**
  - A Unix-like operating system with `make` installed.
  - A C/C++ compiler (e.g., `gcc`) if the project contains compiled source code.
- **Python Plotting:**
  - Python 3.x.
  - [Matplotlib](https://matplotlib.org/) for plotting.
  - [NumPy](https://numpy.org/) for numerical computations.

To install the required Python packages, run:

```bash
pip install matplotlib numpy
```

## Installation and Compilation

### Using the Makefile

The provided Makefile automates the compilation process. You can:

- **Build both projects** by running:

  ```bash
  make
  ```

  This will compile all necessary components. On the other hand, type either:

   ```bash
  make test_parametric_minimal
  ```

   ```bash
  make test_minimal_graph
  ```

  to compile the respective code.

- **Clean the project** by removing compiled files and build artifacts:

  ```bash
  make clean
  ```

## Running the Program

### Command-Line Execution

You can generate minimal surface for boundary given as a graph or as a parametric function using the command-line executable located at `./src/bin`.

**Usage (graph boundary):**

```bash
./src/bin/test_minimal_graph <domain> <N> <size1> <size2>
```

- `<domain>`: Specify the type of domain. Valid options are:
  - `square` — for a square or rectangular domain.
  - `circle` — for a circular domain.
- `<N>`:
  - For a square domain: The number of nodes along one side of the domain.
  - For a circle domain: The number of subcircles in the domain.
- `<size1>`:
  - For a square domain: The length along the x-axis.
  - For a circle domain: The radius of the circle (default is 1).
- `<size2>`: Only applicable for the square domain; it represents the length along the y-axis.

**Examples:**

- For a square domain:

  ```bash
  ./src/bin/test_minimal_graph square 100 1 1.5
  ```

- For a circular domain:

  ```bash
  ./src/bin/test_minimal_graph circle 50 1
  ```

**Usage (parametric boundary):**

```bash
./src/bin/test_parametric_minimal <N> 
```

- `<N>`:
  - The number of subcircles in the domain.
 

**Example:**

  ```bash
  ./src/bin/test_minimal_graph 30
  ```

### Python Plotting

To visualize the minimal surfaces, you can use the provided Python scripts:

- **Plot the Minimal Surface Equation and Graph:**

  ```bash
  python3 plot_mesh.py
  ```

  This command displays a plot of the minimal surface based on its implicit equation along with the corresponding graph.

- **Plot with Parametric Boundary:**

  ```bash
  python3 plot_mesh.py parametric
  ```

  This command displays the minimal surface with a detailed view of its parametric boundary.
  
