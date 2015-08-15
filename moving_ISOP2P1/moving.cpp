#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"

#define DIM 2

void ISOP2P1::syncMesh()
{
    int n_point = this->n_point();

    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();

    for (int i = 0; i < mesh_p.n_geometry(0); ++i)
    {
    	(*mesh_p.h_geometry<0>(i))[0] = this->point(i)[0];
    	(*mesh_p.h_geometry<0>(i))[1] = this->point(i)[1];
    }

    /// 更新p网格的中点. 学习这个过程, 似乎可以去掉单元对应? 
    for (int j = 0; j < mesh_p.n_geometry(1); ++j)
    {
    	GeometryBM &bnd = mesh_p.geometry(1, j);
    	(*mesh_p.h_geometry<1>(bnd.index())->child[1]->vertex[0])[0]
    	    = 0.5 * ((*mesh_p.h_geometry<0>(bnd.vertex(0)))[0] + 
    		     (*mesh_p.h_geometry<0>(bnd.vertex(1)))[0]);
    	(*mesh_p.h_geometry<1>(bnd.index())->child[1]->vertex[0])[1]
    	    = 0.5 * ((*mesh_p.h_geometry<0>(bnd.vertex(0)))[1] + 
    		     (*mesh_p.h_geometry<0>(bnd.vertex(1)))[1]);
    }

    irregular_mesh_p->semiregularize();
    irregular_mesh_p->regularize(false);
    irregular_mesh_v->semiregularize();
    irregular_mesh_v->regularize(false);
};

void ISOP2P1::getMonitor()
{
    syncMesh();
    int i, l;
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space_p.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space_p.endElement();
    for (i = 0; the_element != end_element; ++the_element) 
    {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(1);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
	std::vector<std::vector<double> > u_h_gradient = p_h.gradient(q_point, *the_element);
	float d = 0, area = 0;
	for (l = 0; l < n_quadrature_point; l++) 
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    area += Jxw;
	    d += Jxw * innerProduct(u_h_gradient[l], u_h_gradient[l]);
	}
	monitor(i++) = d / area;
    }
    std::cout << "max monitor=" << *std::max_element(monitor().begin(), monitor().end())
	      << "\tmin monitor=" << *std::min_element(monitor().begin(), monitor().end())
	      << std::endl;
    smoothMonitor(2);
    for (i = 0; i < n_geometry(2); i++)
	monitor(i) = 1.0 / sqrt(1. +  monitor(i));

};

void ISOP2P1::updateSolution()
{
    fem_space_p.updateDofInterpPoint();

    int i, j, l;
    FEMFunction<double, DIM> _u_h(p_h);
    const double& msl = moveStepLength();
    MassMatrix<DIM, double> matrix(fem_space_p);
    matrix.algebricAccuracy() = 2;
    matrix.build();
    for (i = 1; i > 0; i--) 
    {
	Vector<double> rhs(fem_space_p.n_dof());
	FEMSpace<double, DIM>::ElementIterator the_element = fem_space_p.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element = fem_space_p.endElement();
	for (; the_element != end_element; ++the_element) 
	{
	    double volume = the_element->templateElement().volume();
	    const QuadratureInfo<DIM> &quad_info = the_element->findQuadratureInfo(2);
	    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	    int n_quadrature_point = quad_info.n_quadraturePoint();
	    std::vector<AFEPack::Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	    std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
	    std::vector<double> _u_h_value = _u_h.value(q_point, *the_element);
	    std::vector<std::vector<double> > u_h_gradient = p_h.gradient(q_point, *the_element);
	    std::vector<std::vector<double> > move_vector = moveDirection(q_point, the_element->index());
	    int n_element_dof = the_element->n_dof();
	    const std::vector<int>& element_dof = the_element->dof();
	    for (l = 0; l < n_quadrature_point; l++) 
	    {
		double Jxw = quad_info.weight(l) * jacobian[l] * volume;
		for (j = 0; j < n_element_dof; j++) 
		{
		    rhs(element_dof[j]) += Jxw * basis_value[j][l] * (_u_h_value[l]
								      + (1.0 / i) * msl * innerProduct(move_vector[l], u_h_gradient[l]));
		}
	    }
	}

	const BurgV u(viscosity, t);
	BoundaryFunction<double, DIM> boundary1(BoundaryConditionInfo::DIRICHLET, 1, *(dynamic_cast<const Function<double> *>(&u)));
	BoundaryFunction<double, DIM> boundary2(BoundaryConditionInfo::DIRICHLET, 2, *(dynamic_cast<const Function<double> *>(&u)));
	BoundaryFunction<double, DIM> boundary3(BoundaryConditionInfo::DIRICHLET, 3, *(dynamic_cast<const Function<double> *>(&u)));
	BoundaryFunction<double, DIM> boundary4(BoundaryConditionInfo::DIRICHLET, 4, *(dynamic_cast<const Function<double> *>(&u)));

	BoundaryConditionAdmin<double,2> boundary_admin(fem_space_p);

	boundary_admin.add(boundary1);
	boundary_admin.add(boundary2);
	boundary_admin.add(boundary3);
	boundary_admin.add(boundary4);

	boundary_admin.apply(matrix, p_h, rhs);

	AMGSolver solver(matrix);

	solver.solve(p_h, rhs);
    };

}; 

void ISOP2P1::outputSolution()
{
	outputPhysicalMesh("E");
	p_h.writeOpenDXData("u_h.dx");
};


void ISOP2P1::Matrix::getElementMatrix(
		const Element<double,2>& element0,
		const Element<double,2>& element1,
		const ActiveElementPairIterator<2>::State state)
{
	int n_element_dof0 = elementDof0().size();
	int n_element_dof1 = elementDof1().size();
	double volume = element0.templateElement().volume();
	const QuadratureInfo<2>& quad_info = element0.findQuadratureInfo(algebricAccuracy());
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<2> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<double> > basis_value = element0.basis_function_value(q_point);
	std::vector<std::vector<std::vector<double> > > basis_gradient = element0.basis_function_gradient(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
		double Jxw = quad_info.weight(l)*jacobian[l]*volume;
		for (int j = 0;j < n_element_dof0;j ++) {
			for (int k = 0;k < n_element_dof1;k ++) {
				elementMatrix(j,k) += Jxw*((1/dt)*basis_value[j][l]*basis_value[k][l]
						       + a*innerProduct(basis_gradient[j][l], basis_gradient[k][l]));
			}
		}
	}
};


void ISOP2P1::stepForward()
{
    std::cout << "step" << std::endl;
    std::cout << "dt = " << dt << std::endl;
    outputSolution();
    int i, j, k, l;
    FEMFunction<double,2> _u_h(p_h);
    Matrix matrix(fem_space_p, dt, viscosity);
    matrix.algebricAccuracy() = 2;
    matrix.build();
    Vector<double> rhs(fem_space_p.n_dof());
    FEMSpace<double,2>::ElementIterator the_element = fem_space_p.beginElement();
    FEMSpace<double,2>::ElementIterator end_element = fem_space_p.endElement();
    for (; the_element != end_element; ++the_element) 
    {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(2);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<2> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
	std::vector<double> u_h_value = p_h.value(q_point, *the_element);
	std::vector<std::vector<double> > u_h_gradient = p_h.gradient(q_point, *the_element);
	int n_element_dof = the_element->n_dof();
	const std::vector<int>& element_dof = the_element->dof();
	for (l = 0; l < n_quadrature_point; l++) 
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    for (j = 0; j < n_element_dof; j++) 
	    {
		rhs(element_dof[j]) += Jxw * (u_h_value[l] * basis_value[j][l] / dt
					      - u_h_value[l] * (u_h_gradient[l][0] + u_h_gradient[l][1]) * basis_value[j][l]);
	    }
	}
    }
    
    const BurgV u(viscosity, t + dt);
    BoundaryFunction<double, DIM> boundary1(BoundaryConditionInfo::DIRICHLET, 1, *(dynamic_cast<const Function<double> *>(&u)));
    BoundaryFunction<double, DIM> boundary2(BoundaryConditionInfo::DIRICHLET, 2, *(dynamic_cast<const Function<double> *>(&u)));
    BoundaryFunction<double, DIM> boundary3(BoundaryConditionInfo::DIRICHLET, 3, *(dynamic_cast<const Function<double> *>(&u)));
    BoundaryFunction<double, DIM> boundary4(BoundaryConditionInfo::DIRICHLET, 4, *(dynamic_cast<const Function<double> *>(&u)));

    BoundaryConditionAdmin<double,2> boundary_admin(fem_space_p);

    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);

    boundary_admin.apply(matrix, p_h, rhs);


    AMGSolver solver(matrix);
    solver.solve(p_h, rhs);
    t += dt;
};


void ISOP2P1::movingMesh()
{
    moveMesh();
    syncMesh();
    fem_space_p.updateDofInterpPoint();
    fem_space_v.updateDofInterpPoint();
};

