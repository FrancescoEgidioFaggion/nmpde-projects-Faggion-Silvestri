#ifndef STOKES_HPP
#define STOKES_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>


#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>



#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class Stokes
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  // Trigonometric functions example:

  

class ExactSolution : public Function<dim>
{
public:
  ExactSolution() : Function<dim>(2) {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    const double x = p[0];
    const double y = p[1];
    const double pi = M_PI;

    const double sin_px = std::sin(pi * x);
    const double cos_px = std::cos(pi * x);
    const double sin_py = std::sin(pi * y);
    const double cos_py = std::cos(pi * y);

    // Componente x
    values[0] = pi * sin_px * cos_py;

    // Componente y
    values[1] = -pi * cos_px * sin_py;
  }

  virtual double
  value(const Point<dim> &p,
        const unsigned int component = 0) const override
  {
    const double x = p[0];
    const double y = p[1];
    const double pi = M_PI;

    const double sin_px = std::sin(pi * x);
    const double cos_px = std::cos(pi * x);
    const double sin_py = std::sin(pi * y);
    const double cos_py = std::cos(pi * y);

    if (component == 0)
      return pi * sin_px * cos_py;
    else if (component == 1)
      return -pi * cos_px * sin_py;
    else
      return 0.0;
  }
};

  class ForcingTerm : public Function<dim>
{
public:
  ForcingTerm() : Function<dim>(2) {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    const double x = p[0];
    const double y = p[1];
    const double pi = M_PI;

    // Componente x
    values[0] = 2.0 * std::pow(pi, 3) * std::sin(pi * x) * std::cos(pi * y)
                - pi * std::sin(pi * x);

    // Componente y
    values[1] = -2.0 * std::pow(pi, 3) * std::cos(pi * x) * std::sin(pi * y)
                - pi * std::cos(pi * y);
  }
  

  virtual double
  value(const Point<dim> &p,
        const unsigned int component = 0) const override
  {
    const double x = p[0];
    const double y = p[1];
    const double pi = M_PI;

    if (component == 0)
      return 2.0 * std::pow(pi, 3) * std::sin(pi * x) * std::cos(pi * y)
             - pi * std::sin(pi * x);
    else if (component == 1)
      return -2.0 * std::pow(pi, 3) * std::cos(pi * x) * std::sin(pi * y)
             - pi * std::cos(pi * y);
    else
      return 0.0;
  }
};

// Polynomial functions example:

/*

class ExactSolution : public Function<dim>
{
public:
  ExactSolution() : Function<dim>(2) {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    const double x = p[0];
    const double y = p[1];

    // Componente x
    values[0] = 2.0 * x * x * y;

    // Componente y
    values[1] = -2.0 * x * y * y;
  }

  virtual double
  value(const Point<dim> &p,
        const unsigned int component = 0) const override
  {
    const double x = p[0];
    const double y = p[1];

    if (component == 0)
      return 2.0 * x * x * y;
    else if (component == 1)
      return -2.0 * x * y * y;
    else
      return 0.0;
  }
};

class ForcingTerm : public Function<dim>
{
public:
  ForcingTerm() : Function<dim>(2) {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    const double x = p[0];
    const double y = p[1];

    values[0] = 1.0-4.0*y;
    values[1] = -2.0+4.0*x;
  }

  virtual double
  value(const Point<dim> &p,
        const unsigned int component = 0) const override
  {
    const double x = p[0];
    const double y = p[1];

    if (component == 0)
      return 1.0 - 4.0 * y;
    else if (component == 1)
      return -2.0 + 4.0 * x;
    else
      return 0.0;
  }
};

*/


  // Identity preconditioner.
  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(BlockVector<double>       &dst,
          const BlockVector<double> &src) const
    {
      dst = src;
    }

  protected:
  };


  class PreconditionBlockTriangular
{
public:
  void
  initialize(const SparseMatrix<double> &velocity_stiffness_,
             const SparseMatrix<double> &pressure_mass_,
             const SparseMatrix<double> &B_)
  {
    velocity_stiffness = &velocity_stiffness_;
    pressure_mass      = &pressure_mass_;
    B                  = &B_;

    preconditioner_velocity.initialize(velocity_stiffness_);
    preconditioner_pressure.initialize(pressure_mass_);
  }

  void
  vmult(BlockVector<double>       &dst,
        const BlockVector<double> &src) const
  {
    SolverControl solver_control_velocity(1000, 1e-12 * src.block(0).l2_norm());
    SolverCG<Vector<double>> solver_cg_velocity(solver_control_velocity);

    solver_cg_velocity.solve(*velocity_stiffness,
                             dst.block(0),
                             src.block(0),
                             preconditioner_velocity);

    Vector<double> tmp(src.block(1).size());
    B->vmult(tmp, dst.block(0));
    tmp.sadd(-1.0, src.block(1));

    SolverControl solver_control_pressure(1000, 1e-12 * src.block(1).l2_norm());
    SolverCG<Vector<double>> solver_cg_pressure(solver_control_pressure);

    solver_cg_pressure.solve(*pressure_mass,
                             dst.block(1),
                             tmp,
                             preconditioner_pressure);
  }

private:
  const SparseMatrix<double> *velocity_stiffness;
  const SparseMatrix<double> *pressure_mass;
  const SparseMatrix<double> *B;

  PreconditionSSOR<SparseMatrix<double>> preconditioner_velocity;
  PreconditionSSOR<SparseMatrix<double>> preconditioner_pressure;
};

  // Constructor.
  Stokes(Triangulation<dim> &mesh_stokes_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_)
    : mesh(mesh_stokes_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)

  {}

  // Setup system.
  void
  setup();

  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void
  assemble();

  // Solve system.
  void
  solve();

  // Output results.
  void
  output();

  void
  get_boundary_data(std::vector<std::tuple<double, double, double>> &boundary_data) const;

  void
  apply_gamma_boundary_conditions(Vector<double> &boundary_data);

  void 
  apply_boundary_conditions();


protected:

  const double mu = 1.0;

  // Forcing term.
  ForcingTerm forcing_term;

  ExactSolution ExactSolution;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh.
  
  Triangulation<dim> &mesh;


  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  unsigned int total_dofs;
  
  unsigned int n_u;
  unsigned int n_p;

  BlockSparsityPattern sparsity;
  BlockSparsityPattern sparsity_pressure_mass;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;
  

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for face integrals.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // System matrix.
 BlockSparseMatrix<double> system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
 BlockSparseMatrix<double> pressure_mass;

  // Right-hand side vector in the linear system.
  BlockVector<double> system_rhs;

  // System solution 
  BlockVector<double> solution;
};




#endif