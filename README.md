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
Download the meshes folder from: https://drive.google.com/drive/folders/1b53ZrOJ64bX6KrA3baiffa9-QB8Bxg3C?usp=share_link

```bash
# Clone the repository
git clone https://github.com/FrancescoEgidioFaggion/nmpde-projects-Faggion-Silvestri
cd nmpde-projects-Faggion-Silvestri

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
**Manufactured solution with trigonomotric functions and with polynomial functions.



---

## Author
Francesco Egidio Faggion
Martina Silvestri

