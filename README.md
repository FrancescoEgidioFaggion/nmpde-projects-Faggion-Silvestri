# Darcy–Stokes Coupled Solver

## Description
This project implements a solver for the coupled **Darcy–Stokes problem** on a 2D domain.  
The domain is split into two subdomains:
- **Left**: Stokes equations  
- **Right**: Darcy equations  

The coupling is enforced iteratively on the interface using a residual-based fixed-point algorithm.  
The code is written in **C++** with the [deal.II](https://www.dealii.org/) finite element library.

---

## Repository Structure
```
├── src/                # Source code
│   ├── main.cpp        # Entry point
│   ├── Darcy.hpp/.cpp  # Darcy solver
│   ├── Stokes.hpp/.cpp # Stokes solver
├── CMakeLists.txt      # Build configuration
└── README.md           # This file
```

---

## Compilation
Requirements: **CMake** and **deal.II**.

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/nmpde-projects.git
cd nmpde-projects/coupled-darcy-stokes

# Create build directory
mkdir build && cd build

# Run CMake
cmake ..

# Compile
make
```

The executable will be available in:
```
build/coupled_executable
```

---

## Run
Execute the program with:
```bash
./coupled_executable
```

Polynomial degrees and solver settings can be modified in `src/main.cpp`.

---

## Tests
Two sample tests are implemented:
1. **Manufactured solution (sin/cos forcing terms)** → checks convergence.  
2. **Inflow/Outflow setup** → tests coupling between Stokes and Darcy regions.

---

## Author
Francesco Faggion (Politecnico di Milano)

---

## Delivery Notes
- Code delivery via pull request to the original repository.  
- Report sent by email to *michele.bucelli@polimi.it*.  
