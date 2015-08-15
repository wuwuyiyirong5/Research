#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::solveStokes()
{
    buildStokesSys();
    int n_dof_v = fem_space_v.n_dof();
    int n_dof_p = fem_space_p.n_dof();
    int n_total_dof = 2 * n_dof_v + n_dof_p;

    /// 构建系数矩阵和右端项.
    /// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
    Vector<double> x(n_total_dof);
    /// 将数值解合并一个向量便于边界处理.
    for (int i = 0; i < n_dof_v; ++i)
    {
    	x(i) = v_h[0](i);
    	x(n_dof_v + i) = v_h[1](i);
    }

    for (int i = 0; i < n_dof_p; ++i)
    	x(2 * n_dof_v + i) = p_h(i);

    rhs.reinit(n_total_dof);
    /// 边界条件一起处理了.
    boundaryValueStokes(x);

    /// 矩阵求解. 
	SparseMatrix<double> mat_BTx(sp_pvx);
	SparseMatrix<double> mat_BTy(sp_pvy);
	SparseMatrix<double> mat_Bx(sp_vxp);
	SparseMatrix<double> mat_By(sp_vyp);
	SparseMatrix<double> mat_Ax(sp_vxvx);
	SparseMatrix<double> mat_Ay(sp_vyvy);

	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
	    mat_Ax.global_entry(i) = matrix.global_entry(index_vxvx[i]);
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
	    mat_Ay.global_entry(i) = matrix.global_entry(index_vyvy[i]);
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
	    mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
	    mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
	    mat_Bx.global_entry(i) = matrix.global_entry(index_vxp[i]);
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
	    mat_By.global_entry(i) = matrix.global_entry(index_vyp[i]);

    Vector<double> tmp1(n_dof_v);
    Vector<double> tmp2(n_dof_v);
    Vector<double> rhs_vx(n_dof_v);
    Vector<double> rhs_vy(n_dof_v);
    Vector<double> rhs_p(n_dof_p);

    for (int i = 0; i < n_dof_v; ++i)
    {
    	rhs_vx(i) = rhs(i);
    	v_h[0](i) = x(i);
    	rhs_vy(i) = rhs(n_dof_v + i);
    	v_h[1](i) = x(n_dof_v + i);
    }
    for (int i = 0; i < n_dof_p; ++i)
    {
    	rhs_p(i) = rhs(2 * n_dof_v + i);
    	p_h(i) = x(2 * n_dof_v + i);
    }

    Vector<double> schur_rhs (n_dof_p);
//    AMGSolver solverQ(mat_Ax, 1.0e-12, 3, 50, 0.382, 0.25);
    AMGSolver solverQ(mat_Ax);
    InverseMatrix M(mat_Ax, solverQ);
    M.vmult (tmp1, rhs_vx);
    M.vmult (tmp2, rhs_vy);
    mat_Bx.vmult(schur_rhs, tmp1);
    mat_By.vmult_add(schur_rhs, tmp2);
    schur_rhs -= rhs_p;

    SchurComplement schur_complement(mat_BTx, mat_BTy, mat_Bx, mat_By, M, M);

    SolverControl solver_control_cg (n_dof_p * 2,
				  1e-12*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control_cg);

    std::ofstream output("output.m");
    output.setf(std::ios::fixed);
    output.precision(20);
    output << "A =[" << std::endl;
	for (int i = 0; i < n_dof_p; ++i)
	{
		for (int j = 0; j < n_dof_v; ++j)
			output << mat_Bx.el(i,j) << "\t";
		output << std::endl;
	}
	output << "];" << std::endl;
	output.close();
	std::cout << "stokes matrix outputed!" << std::endl;
	getchar();
    AMGSolver AQ(mat_p_mass);
    ApproxSchurComplement asc(mat_p_mass, AQ);
    cg.solve (schur_complement, p_h, schur_rhs, asc);
    std::cout << solver_control_cg.last_step()
              << " CG Schur complement iterations to obtain convergence."
              << std::endl;

    mat_BTx.vmult(tmp1, *dynamic_cast<const Vector<double>* >(&p_h));
    mat_BTy.vmult(tmp2, *dynamic_cast<const Vector<double>* >(&p_h));
    tmp1 *= -1;
    tmp2 *= -1;
    tmp1 += rhs_vx;
    tmp2 += rhs_vy;

    M.vmult(v_h[0], tmp1);
    M.vmult(v_h[1], tmp2);
};

void ISOP2P1::buildStokesSys()
{
    int n_dof_v = fem_space_v.n_dof();
    int n_dof_p = fem_space_p.n_dof();
    int n_total_dof = 2 * n_dof_v + n_dof_p;
    
    matrix.reinit(sp_stokes);

    /// (0, 0)
    for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
	matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i);
    /// (1, 1)
    for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
	matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i);
    /// (0, 2)
    for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
	matrix.global_entry(index_pvx[i]) = mat_pvx_divT.global_entry(i);

    /// (1, 2)
    for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
	matrix.global_entry(index_pvy[i]) = mat_pvy_divT.global_entry(i);

    /// (2, 0)
    for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
	matrix.global_entry(index_vxp[i]) =  mat_vxp_div.global_entry(i);
	
    /// (2, 1)
    for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
	matrix.global_entry(index_vyp[i]) =  mat_vyp_div.global_entry(i);

};

#undef DIM
