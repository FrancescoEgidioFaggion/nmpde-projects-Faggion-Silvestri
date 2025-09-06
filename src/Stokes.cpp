#include "Stokes.hpp"

void
Stokes::setup()
{
  // Load the mesh.
  {
   
    std::cout << "  Number of elements = " << mesh.n_active_cells() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    std::cout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity, dim,
                                         fe_scalar_pressure, 1);

    std::cout << "  Velocity degree:           = " << fe_scalar_velocity.degree << std::endl;
    std::cout << "  Pressure degree:           = " << fe_scalar_pressure.degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);
    std::cout << "  Quadrature points per cell = " << quadrature->size() << std::endl;

    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);
    std::cout << "  Quadrature points per face = " << quadrature_face->size() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  
    std::cout << "Initializing the DoF handler" << std::endl;

    
    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // Reorder DoFs so that velocity DoFs come first, then pressure.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

    n_u = dofs_per_block[0];
    n_p = dofs_per_block[1];

    total_dofs = n_u + n_p;

    std::cout << "  Number of DoFs: " << std::endl;
    std::cout << "    velocity = " << n_u << std::endl;
    std::cout << "    pressure = " << n_p << std::endl;
    std::cout << "    total    = " << total_dofs << std::endl;
  

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    std::cout << "Initializing the linear system" << std::endl;
    std::cout << "  Initializing the sparsity pattern" << std::endl;

    // Coupling table for Stokes problem
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
      for (unsigned int d = 0; d < dim + 1; ++d)
        coupling[c][d] = (c == dim && d == dim) ? DoFTools::none : DoFTools::always;

    BlockDynamicSparsityPattern dsp(2, 2);

    dsp.block(0, 0).reinit(n_u, n_u);
    dsp.block(0, 1).reinit(n_u, n_p);
    dsp.block(1, 0).reinit(n_p, n_u);
    dsp.block(1, 1).reinit(n_p, n_p);
    dsp.collect_sizes();
    
    AffineConstraints<double> constraints;
    constraints.close();

        DoFTools::make_sparsity_pattern(dof_handler,
                                coupling,
                                dsp,
                                constraints,
                                false);

    sparsity.copy_from(dsp);
    
    std::cout << "  A1" << std::endl;

    // Pressure mass matrix sparsity
    for (unsigned int c = 0; c < dim + 1; ++c)
      for (unsigned int d = 0; d < dim + 1; ++d)
        coupling[c][d] = (c == dim && d == dim) ? DoFTools::always : DoFTools::none;

    BlockDynamicSparsityPattern dsp_pressure_mass(2, 2);
    dsp_pressure_mass.block(0, 0).reinit(n_u, n_u);
    dsp_pressure_mass.block(0, 1).reinit(n_u, n_p);
    dsp_pressure_mass.block(1, 0).reinit(n_p, n_u);
    dsp_pressure_mass.block(1, 1).reinit(n_p, n_p);
    dsp_pressure_mass.collect_sizes();

    std::cout << "  A2" << std::endl;

    DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp_pressure_mass);
    sparsity_pressure_mass.copy_from(dsp_pressure_mass);

    std::cout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);


    system_rhs.reinit(2);
    system_rhs.block(0).reinit(n_u);
    system_rhs.block(1).reinit(n_p);
    system_rhs.collect_sizes();

    solution.reinit(2);
    solution.block(0).reinit(n_u);
    solution.block(1).reinit(n_p);
    solution.collect_sizes();
  }
}

void
Stokes::assemble()
{
  std::cout << "===============================================" << std::endl;
  std::cout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  

  FEValues<dim>     fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_normal_vectors | update_quadrature_points |
                                     update_JxW_values);

  FullMatrix<double> cell_A_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_B_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;
  pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);

      cell_A_matrix             = 0.0;
      cell_rhs                  = 0.0;
      cell_B_matrix = 0.0;

      
      for (unsigned int q = 0; q < n_q; ++q)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Viscosity term.
                  cell_A_matrix(i, j) +=
                    mu *
                    scalar_product(fe_values[velocity].gradient(i, q),
                                   fe_values[velocity].gradient(j, q)) *
                    fe_values.JxW(q);

                  // Pressure term in the momentum equation.
                    cell_A_matrix(i, j) -= fe_values[pressure].value(j, q) *
                               trace(fe_values[velocity].gradient(i, q)) *
                               fe_values.JxW(q);
                    cell_A_matrix(i, j) -= fe_values[pressure].value(i, q) *
                               trace(fe_values[velocity].gradient(j, q)) *
                               fe_values.JxW(q);         

                  // Pressure mass matrix.
                  cell_B_matrix(i, j) +=
                    fe_values[pressure].value(i, q) *
                    fe_values[pressure].value(j, q) * fe_values.JxW(q);
                }

              // Forcing term scalar velocity test function.
                for (unsigned int d = 0; d < dim; ++d)
                {
                cell_rhs(i) -= 
                  forcing_term.value(fe_values.quadrature_point(q), d) *
                  fe_values[velocity].value(i, q)[d] *
                  fe_values.JxW(q);
                }
            }
        }
    

      for (unsigned int i=0; i< dofs_per_cell; ++i){
        for (unsigned int j=0; j< dofs_per_cell; ++j){
          system_matrix.add(dof_indices[i],
                            dof_indices[j],
                            cell_A_matrix(i, j));
        }
      }

      system_rhs.add(dof_indices, cell_rhs);


      // Pressure mass matrix.
      for (unsigned int i=0; i< dofs_per_cell; ++i){
        for (unsigned int j=0; j< dofs_per_cell; ++j){
          pressure_mass.add(dof_indices[i],
                            dof_indices[j],
                            cell_B_matrix(i, j));
        }
      }

    }

}

void Stokes::apply_boundary_conditions(){

 std::map<types::global_dof_index, double> boundary_values;
 std::map<types::boundary_id, const Function<dim> *> boundary_functions;

boundary_functions.clear();
Functions::ZeroFunction<dim> zero_function(dim);  // Zero velocity
//boundary_functions[2] = &zero_function;  // Boundary ID 2, velocity = 0
boundary_functions[2] = &ExactSolution;


// Apply zero velocity condition on boundary ID 2
VectorTools::interpolate_boundary_values(dof_handler,
                     boundary_functions,
                     boundary_values,
                     ComponentMask({true, true, false}));

 MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs, false);

}

void
Stokes::solve()
{
  std::cout << "===============================================" << std::endl;

  std::cout << " norm of the rhs = " << system_rhs.l2_norm() << std::endl;

  SolverControl solver_control(20000, 1e-6 /* * system_rhs.l2_norm()*/);


  SolverGMRES<BlockVector<double>> solver(solver_control);

   //PreconditionBlockDiagonal preconditioner;
   //preconditioner.initialize(system_matrix.block(0, 0),
   //                          pressure_mass.block(1, 1));
  
  PreconditionBlockTriangular preconditioner;
  preconditioner.initialize(system_matrix.block(0, 0),
                         pressure_mass.block(1, 1),
                        system_matrix.block(1, 0));


  std::cout << "Solving the linear system" << std::endl;
  solver.solve(system_matrix, solution, system_rhs, preconditioner);
  std::cout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;
}

void Stokes::output()
{
  static unsigned int output_index = 0; // contatore che rimane tra chiamate

  std::cout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  // Interpret the first 'dim' components as part of a vector (velocity),
  // the last one as a scalar (pressure)
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  // Field names
  std::vector<std::string> names;
  for (unsigned int d = 0; d < dim; ++d)
    names.push_back("velocity");
  names.push_back("pressure");

  // Add the solution
  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           data_component_interpretation);

  data_out.build_patches();

  // File name with progressive index
  const std::string filename = "output/solution_Stokes" + std::to_string(output_index) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  std::cout << "Output written to file: " << filename << std::endl;
  std::cout << "===============================================" << std::endl;

  ++output_index; 
}



void Stokes::get_boundary_data(std::vector<std::tuple<double, double, double>> &boundary_data) const {
  boundary_data.clear();
  std::set<unsigned int> processed_vertex_indices;

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (cell->at_boundary()) {
      for (unsigned int f = 0; f < cell->n_faces(); ++f) {
        if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == 1) {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
            const unsigned int vertex_index = cell->face(f)->vertex_index(v);
            if (processed_vertex_indices.find(vertex_index) != processed_vertex_indices.end())
              continue;
            processed_vertex_indices.insert(vertex_index);

            const Point<dim> &vertex = cell->face(f)->vertex(v);
            const unsigned int dof_index_x = cell->vertex_dof_index(v, 0);
            const unsigned int dof_index_p = cell->vertex_dof_index(v, dim);

            const double velocity_x = solution[dof_index_x];
            const double pressure = solution[dof_index_p];
            const double y_coord = vertex[1];

            boundary_data.emplace_back(velocity_x, pressure, y_coord);
          }
        }
      }
    }
  }

  // Sort the boundary_data vector based on y-coordinates
    // Sort the vector of tuples by the third element (y_coord) in ascending order
    std::sort(boundary_data.begin(), boundary_data.end(),
    [](const std::tuple<double, double, double> &a,
       const std::tuple<double, double, double> &b) {
      return std::get<2>(a) < std::get<2>(b);
    });


  // Compute the average of the elements in position 0 of the tuples
  double sum = 0.0;
  for (const auto &[velocity_x, pressure, y_coord] : boundary_data)
    sum += velocity_x;

  sum -= 0.5 * std::get<0>(boundary_data[0]); // Subtract the first element
  sum -= 0.5 * std::get<0>(boundary_data[boundary_data.size() - 1]); // Subtract the last element
  double average = sum / (boundary_data.size()-1);

  // Subtract the average from each element in position 0
  for (auto &[velocity_x, pressure, y_coord] : boundary_data)
    velocity_x -= average;
}


void Stokes::apply_gamma_boundary_conditions(Vector<double> &boundary_data)
{
  AssertDimension(system_rhs.size(), boundary_data.size());
  for (unsigned int i = 0; i < boundary_data.size(); ++i)
    system_rhs[i] += boundary_data[i];
}

