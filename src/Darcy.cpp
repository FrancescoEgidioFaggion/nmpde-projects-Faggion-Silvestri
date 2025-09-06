#include "Darcy.hpp"

void
Darcy::setup()
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

    // Coupling table for Darcy problem
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


    DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp_pressure_mass);
    sparsity_pressure_mass.copy_from(dsp_pressure_mass);

    std::cout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);

    std::cout << "  Initializing the system right-hand side" << std::endl;

    const unsigned int n_dofs = dof_handler.n_dofs();
    std::cout << "    Number of DoFs = " << n_dofs << std::endl;

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
Darcy::assemble()
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
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  FullMatrix<double> cell_A_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_B_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs = 0.0;
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
                    scalar_product(fe_values[velocity].value(i, q),
                                   fe_values[velocity].value(j, q)) *
                    fe_values.JxW(q);

                  // Pressure term in the momentum equation.
                  cell_A_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                       fe_values[pressure].value(j, q) *
                                       fe_values.JxW(q);
                  cell_A_matrix(i,j) -= fe_values[velocity].divergence(j,q) *
                                       fe_values[pressure].value(i, q) *
                                       fe_values.JxW(q); 

                  // Pressure mass matrix.
                  cell_B_matrix(i, j) +=
                    fe_values[pressure].value(i, q) *
                    fe_values[pressure].value(j, q) * fe_values.JxW(q);
                }

              // Forcing term.
              cell_rhs(i) -= 
                K*forcing_term.value(fe_values.quadrature_point(q)) *
                fe_values[pressure].value(i, q) *
                             fe_values.JxW(q);
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

  std::cout << "matrices assembled!" << std::endl;

 std::map<types::global_dof_index, double> boundary_values;
 std::map<types::boundary_id, const Function<dim> *> boundary_functions;

 
 // 1. Impose velocity flux = 0 everywhere on the boundary (no flow across boundaries)
 boundary_functions.clear();
 Functions::ZeroFunction<dim> zero_function(dim);  // Zero velocity
 //boundary_functions[3] = &zero_function;
 boundary_functions[3] = &ExactSolution;

 

 VectorTools::interpolate_boundary_values(dof_handler,
                                          boundary_functions,
                                          boundary_values,
                                          ComponentMask({true, true, false}));


  boundary_functions.clear();
  boundary_functions[1] = &FunctionH;
  
  VectorTools::interpolate_boundary_values(dof_handler,
                                           boundary_functions,
                                           boundary_values,
                                           ComponentMask({true, false, false}));  // Apply to pressure (last component)
 
 MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs, false);
 
  std::cout << "Boundary conditions applied!" << std::endl;
}

void
Darcy::solve()
{
  std::cout << "===============================================" << std::endl;

  SolverControl solver_control(20000, 1e-6 /* * system_rhs.l2_norm()*/);

  SolverGMRES<BlockVector<double>> solver(solver_control);

   //PreconditionBlockDiagonal preconditioner;
   //preconditioner.initialize(system_matrix.block(0, 0),
                             //pressure_mass.block(1, 1));

  PreconditionBlockTriangular preconditioner;
  preconditioner.initialize(system_matrix.block(0, 0),
                           pressure_mass.block(1, 1),
                            system_matrix.block(1, 0));


  std::cout << "Solving the linear system" << std::endl;
  solver.solve(system_matrix, solution, system_rhs, preconditioner);
  std::cout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

}


void Darcy::output()
{
  static unsigned int output_index = 0; // counter

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

  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           data_component_interpretation);

  data_out.build_patches();

  // File name with progressive index
  const std::string filename = "output/solution_Darcy" + std::to_string(output_index) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  std::cout << "Output written to file: " << filename << std::endl;
  std::cout << "===============================================" << std::endl;

  ++output_index; 
}





void Darcy::get_boundary_data(std::vector<std::tuple<double, double, double>> &boundary_data) const {
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

  // Sort the vector of tuples by the third element (y_coord) in ascending order
  std::sort(boundary_data.begin(), boundary_data.end(),
            [](const std::tuple<double, double, double> &a,
               const std::tuple<double, double, double> &b) {
              return std::get<2>(a) < std::get<2>(b);
            });

  // Compute the average of the elements in position 1 of the tuples
  double sum = 0.0;
  for (const auto &[velocity_x, pressure, y_coord] : boundary_data)
    sum += pressure;

  sum -= 0.5 * std::get<1>(boundary_data[0]); // Subtract the first element
  sum -= 0.5 * std::get<1>(boundary_data[boundary_data.size() - 1]); // Subtract the last element
  double average = sum / (boundary_data.size()-1);

  // Subtract the average from each element in position 0
  for (auto &tuple : boundary_data)
    std::get<1>(tuple) -= average;

}

void Darcy::compute_residual(Vector<double> &residual) const
{
  std::cout << "Computing residual..." << std::endl;

  residual.reinit(solution.size());

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  Vector<double>     cell_residual(dofs_per_cell);

  FEValues<dim> fe_values(*fe, *quadrature, update_values | update_gradients | update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {

      fe_values.reinit(cell);
      cell_residual = 0.0; // Reset cell_residual for each cell

      std::vector<Tensor<1, dim>> velocity_values(n_q);
      std::vector<double> pressure_values(n_q);

      fe_values[velocity].get_function_values(solution, velocity_values);
      fe_values[pressure].get_function_values(solution, pressure_values);

      for (unsigned int q = 0; q < n_q; ++q)
        {

            const Tensor<1, dim> velocity_q = velocity_values[q];
            const double pressure_q = pressure_values[q];
        
            // Extractors for velocity and pressure
            FEValuesExtractors::Vector velocity(0);
            FEValuesExtractors::Scalar pressure(dim);

            // Get global DoF indices for this cell
            std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
            cell->get_dof_indices(dof_indices);

            // Loop over test functions (basis functions) i
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                double local_residual = 0.0;

                // Compute mu * ∫ (velocity_solution ⋅ basis_i) dx
                local_residual += mu *
                  scalar_product(
                    fe_values[velocity].value(i, q),
                    // velocity_solution at this quadrature point
                    // Interpolate solution vector to quadrature point
                    // For simplicity, use FEValues::get_function_values
                    velocity_q) *
                   fe_values.JxW(q);

                // Compute -K * ∫ (pressure_solution * div(basis_i)) dx
                local_residual -= K *
                  pressure_q *
                  fe_values[velocity].divergence(i, q) *
                  fe_values.JxW(q);

                cell_residual(i) += local_residual;
              }

            // Add local cell residual to global residual vector
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              residual[dof_indices[i]] += cell_residual(i);
          
        }
    }
}

void Darcy::apply_gamma_boundary_conditions(std::vector<std::tuple<double, double, double>> &input_boundary_data)
{
  FunctionH.set_values(input_boundary_data);
}


