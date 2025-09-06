#ifndef DARCY_HPP
#define DARCY_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>



#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

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

// Class implementing a solver for the Darcy problem.
class Darcy
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  // Function for the forcing term.
   class ForcingTerm : public Function<dim>
   {
   public:
     // Constructor.
     ForcingTerm()
     {}
 
     // Evaluation.
     virtual double
     value(const Point<dim> & /*p*/,
           const unsigned int /*component*/ = 0) const override
     {
       return 0.0;
     }
    };

    //for the example with trigonometric functions:

        class ExactSolution : public Function<dim>
    {
    public:
      ExactSolution()
      {}

      virtual void
  vector_value(const Point<dim> &/*p*/, Vector<double> &values) const override
  {
    values[0] = 0.0;
    values[1] = -1.0;
  }

  virtual double
  value(const Point<dim> & /*p*/,
        const unsigned int component = 0) const override
  {
    if (component == 0)
      return 0.0;
    else if (component == 1)
      return -1.0;
    else
      return 0.0;
  }
  };

  //polynomial example:

    /*

    class ExactSolution : public Function<dim>
    {
    public:
      ExactSolution()
      {}

      virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    values[0] = -2.0* (p[0]-1.0);
    values[1] = 2.0*p[1];
  }

  virtual double
  value(const Point<dim> &p,
        const unsigned int component = 0) const override
  {
    if (component == 0)
      return -2.0* (p[0]-1.0);
    else if (component == 1)
      return 2.0*p[1];
    else
      return 0.0;
  }
  };

*/



   class FunctionH : public Function<dim>
   {
   public:
     // Constructor.
     FunctionH()
     {}
 
     // Evaluation.
     virtual double
     value(const Point<dim> & p,
           const unsigned int /*component*/ = 0) const override
     {
        int temp = (int)(p[1] / (2.0 / (nodes_on_boundary - 1)));
        if (temp != nodes_on_boundary - 1)
        {
          return (values[temp] + (values[temp + 1] - values[temp]) / (2.0 / (nodes_on_boundary - 1)) * (p[1] - temp * (2.0 / (nodes_on_boundary - 1))));
        }
        else
        {
          return (values[temp]);
        }
     }

     void set_values(std::vector<std::tuple<double, double, double>> &values_)
     {
       values.clear();
       for (const auto &value : values_)
         {
           values.push_back(std::get<0>(value));
         }
 
       nodes_on_boundary = values_.size();

    }
    protected:
      std::vector<double> values;
      int nodes_on_boundary;
   };                                                                                  


 


  // Since we're working with block matrices, we need to make our own
  // preconditioner class.

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
  Darcy(Triangulation<dim> &mesh_darcy_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_)
    : mesh(mesh_darcy_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
  {}

  // Setup system.
  void
  setup();

  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner.
  void
  assemble();

  // Solve system.
  void
  solve();

  // Output results.
  void
  output();

  // Get boundary data.
  void
  get_boundary_data(std::vector<std::tuple<double, double, double>> &boundary_data) const;
  // set the boundary data on gamma:
  void
  apply_gamma_boundary_conditions(std::vector<std::tuple<double, double, double>> &boundary_data);

    void
  compute_residual(Vector<double> &residual) const;


protected:

  // Problem definition. ///////////////////////////////////////////////////////
  
  const double K = 1.0;
  const double mu = 1.0;

  // Outlet pressure [Pa].
  const double p_out = 10;

    // Mesh.
  Triangulation<dim> &mesh;

  // Forcing term.
  ForcingTerm forcing_term;

  // Dirichlet boundary values.
  FunctionH FunctionH;

  // Discretization. ///////////////////////////////////////////////////////////

  ExactSolution ExactSolution;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for face integrals.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  BlockSparsityPattern sparsity;
  BlockSparsityPattern sparsity_pressure_mass;

  unsigned int total_dofs;

  unsigned int n_u;
  unsigned int n_p;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // System matrix.
 BlockSparseMatrix<double> system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
 BlockSparseMatrix<double> pressure_mass;

  // Right-hand side vector in the linear system.
  BlockVector<double> system_rhs;

  BlockVector<double> solution;
};

#endif