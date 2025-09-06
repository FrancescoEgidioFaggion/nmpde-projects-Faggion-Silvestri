#include "Stokes.hpp"
#include "Darcy.hpp"

#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>

namespace fs = std::filesystem;

void write_pvd(const std::string &folder, const std::string &pvd_name = "solution.pvd")
{
    std::vector<std::string> vtu_files;

    // Collect all .vtu files in the folder
    for (const auto &entry : fs::directory_iterator(folder))
    {
        if (entry.path().extension() == ".vtu")
            vtu_files.push_back(entry.path().filename().string());
    }

    // Sort alphabetically to respect output_index order
    std::sort(vtu_files.begin(), vtu_files.end());

    // Open the PVD file
    std::ofstream pvd(folder + "/" + pvd_name);
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd << "  <Collection>\n";

    unsigned int timestep = 0;
    for (const auto &file : vtu_files)
    {
        pvd << "    <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\"" << file << "\"/>\n";
        ++timestep;
    }

    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";

    std::cout << "PVD file written: " << folder + "/" + pvd_name << std::endl;
}

void write_pvd_for_domain(const std::string &folder, const std::string &domain_name)
{
    std::vector<std::string> vtu_files;
    const std::string vtu_prefix = "solution_" + domain_name;
    const std::string pvd_name = vtu_prefix + ".pvd";

    // Raccogli tutti i file .vtu che iniziano con il prefisso specifico
    for (const auto &entry : fs::directory_iterator(folder))
    {
        if (entry.path().extension() == ".vtu" && entry.path().filename().string().rfind(vtu_prefix, 0) == 0)
        {
            vtu_files.push_back(entry.path().filename().string());
        }
    }

    // Ordina i file alfabeticamente per rispettare l'ordine delle iterazioni
    std::sort(vtu_files.begin(), vtu_files.end());

    // Scrivi il file PVD
    std::ofstream pvd(folder + "/" + pvd_name);
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd << "  <Collection>\n";

    unsigned int timestep = 0;
    for (const auto &file : vtu_files)
    {
        pvd << "    <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\"" << file << "\"/>\n";
        ++timestep;
    }

    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";

    std::cout << "PVD file written: " << folder + "/" + pvd_name << std::endl;
}

int main(){

    const std::string  mesh_file_name1  = "../meshes/mesh_stokes_structured2.msh";
    const std::string  mesh_file_name2  = "../meshes/mesh_darcy2.msh";
    const unsigned int degree_velocity = 2;
    const unsigned int degree_pressure = 1;
    const unsigned int dim            = 2;

Triangulation<dim> mesh_serial_1; 
Triangulation<dim> mesh_serial_2;

std::vector<std::tuple<double, double, double>> input_boundary_data; 


for (unsigned int i = 0; i < 17; ++i) {
    input_boundary_data.emplace_back(0.0, 0.0, 2.0/16.0 * i);
}

cout << "Initializing the meshes" << std::endl;

GridIn<dim> grid_in;
grid_in.attach_triangulation(mesh_serial_1);

std::ifstream grid_in_file(mesh_file_name1);
grid_in.read_msh(grid_in_file);

GridIn<dim> grid_in2;
grid_in2.attach_triangulation(mesh_serial_2);

std::ifstream grid_in_file2(mesh_file_name2);
grid_in2.read_msh(grid_in_file2);


    std::vector<std::tuple<double, double, double>> boundary_data_1 = input_boundary_data;
    std::vector<std::tuple<double, double, double>> boundary_data_2 = input_boundary_data;
    Vector<double> residual; 

    Stokes problem1(mesh_serial_1, degree_velocity, degree_pressure);
    Darcy problem2(mesh_serial_2, degree_velocity, degree_pressure);
    problem2.setup();
    problem1.setup();


    for (unsigned int cycle = 0; cycle < 20; ++cycle)
    {
        std::cout << "================ CYCLE " << cycle << " ================" << std::endl;

        
        problem2.apply_gamma_boundary_conditions(boundary_data_2);
        problem2.assemble();
        problem2.solve();
        problem2.output();
        problem2.compute_residual(residual);
        problem2.get_boundary_data(boundary_data_1);
        

        
        problem1.assemble();
        problem1.apply_gamma_boundary_conditions(residual);
        problem1.apply_boundary_conditions();
        problem1.solve();
        problem1.output();
        problem1.get_boundary_data(boundary_data_2);
        


        double difference_sum_pressure = 0.0;
        double difference_sum_velocity_x = 0.0;
        for (size_t i = 0; i < boundary_data_1.size(); ++i) {
            double diff_velocity_x = std::get<0>(boundary_data_2[i]) - std::get<0>(boundary_data_1[i]);
            double diff_pressure = std::get<1>(boundary_data_2[i]) - std::get<1>(boundary_data_1[i]);
            
            difference_sum_pressure += diff_pressure * diff_pressure;
            difference_sum_velocity_x += diff_velocity_x * diff_velocity_x;                 
        }
        difference_sum_pressure = std::sqrt(difference_sum_pressure);
        difference_sum_velocity_x = std::sqrt(difference_sum_velocity_x);
        std::cout << "Iteration " << cycle << ": Difference sum velocity_x: " 
                  << difference_sum_velocity_x << ", Difference sum pressure: " 
                  << difference_sum_pressure << std::endl;
        if (difference_sum_velocity_x < 1e-10) {
            std::cout << "The two boundary velocity data vectors are sufficiently similar. Convergence reached! :)" << std::endl;
            break;
        } else {
            std::cout << "The two boundary data vectors are different. Difference sum velocity_x: " 
                      << difference_sum_velocity_x << ", Difference sum pressure: " 
                      << difference_sum_pressure << std::endl;
        }
    }


    write_pvd("output");  // folder where your solution_*.vtu files are stored
    write_pvd_for_domain("output", "Darcy");
    write_pvd_for_domain("output", "Stokes");

}

