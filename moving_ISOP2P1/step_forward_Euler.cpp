#include "ISOP2P1.h"
#include "preconditioner.h"
#define DIM 2

void ISOP2P1::stepForwardEuler()
{
    int n_dof_v = fem_space_v.n_dof();
    int n_dof_p = fem_space_p.n_dof();
    int n_total_dof = 2 * n_dof_v + n_dof_p;

    /// 系数矩阵直接使用 Stokes 矩阵结构.
    matrix.reinit(sp_stokes);

    /// (0, 0) 
    for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
	matrix.global_entry(index_vxvx[i]) = dt * viscosity * mat_v_stiff.global_entry(i) 
	    + mat_v_mass.global_entry(i);
    /// (1, 1) 这两个对角块对应扩散算子和质量算子, dt 直接乘上, 以避免
    /// 矩阵系数过大.
    for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
	matrix.global_entry(index_vyvy[i]) = dt * viscosity * mat_v_stiff.global_entry(i)
	    + mat_v_mass.global_entry(i);
    /// (0, 2) 这个不是方阵. 在矩阵结构定义的时候已经直接排除了对角元优
    /// 先.
    for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
     	matrix.global_entry(index_pvx[i]) = dt * mat_pvx_divT.global_entry(i);

    /// (1, 2)
    for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
     	matrix.global_entry(index_pvy[i]) = dt * mat_pvy_divT.global_entry(i);

    /// (2, 0)
    for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
     	matrix.global_entry(index_vxp[i]) = dt * mat_vxp_div.global_entry(i);
	
    /// (2, 1) 这四块直接复制散度矩阵. 
    for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
     	matrix.global_entry(index_vyp[i]) = dt * mat_vyp_div.global_entry(i);

    /// 问题右端项.
    rhs.reinit(n_total_dof);

    FEMSpace<double,2>::ElementIterator the_element_v = fem_space_v.beginElement();
    FEMSpace<double,2>::ElementIterator end_element_v = fem_space_v.endElement();
    FEMSpace<double,2>::ElementIterator the_element_p = fem_space_p.beginElement();
    FEMSpace<double,2>::ElementIterator end_element_p = fem_space_p.endElement();

    /// 遍历速度单元, 拼装相关系数矩阵和右端项.
    for (the_element_v = fem_space_v.beginElement(); 
	 the_element_v != end_element_v; ++the_element_v) 
    {
	/// 当前单元信息.
	double volume = the_element_v->templateElement().volume();
	/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
	/// 装时积分精度不用超过 1 次. (验证一下!)
	const QuadratureInfo<2>& quad_info = the_element_v->findQuadratureInfo(4);
	std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<2> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
	/// 速度单元信息.
	std::vector<std::vector<double> > basis_value_v = the_element_v->basis_function_value(q_point);
	std::vector<std::vector<std::vector<double> > > basis_gradient_v = the_element_v->basis_function_gradient(q_point);
	std::vector<double> vx_value = v_h[0].value(q_point, *the_element_v);
	std::vector<double> vy_value = v_h[1].value(q_point, *the_element_v);
	std::vector<double> fx_value = source_v[0].value(q_point, *the_element_v);
	std::vector<double> fy_value = source_v[1].value(q_point, *the_element_v);
	std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
	std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);
	const std::vector<int>& element_dof_v = the_element_v->dof();
	int n_element_dof_v = the_element_v->n_dof();
	Element<double, 2> &p_element = fem_space_p.element(index_v2p[the_element_v->index()]);
	const std::vector<int>& element_dof_p = p_element.dof();
	std::vector<std::vector<std::vector<double> > > basis_gradient_p = p_element.basis_function_gradient(q_point);
	std::vector<std::vector<double> >  basis_value_p = p_element.basis_function_value(q_point);
	int n_element_dof_p = p_element.n_dof();
	std::vector<double> p_value = p_h.value(q_point, p_element);
	for (int l = 0; l < n_quadrature_point; ++l)
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    for (int i = 0; i < n_element_dof_v; ++i)
	    {
		double rhs_cont = fx_value[l] * basis_value_v[i][l] + vx_value[l] * basis_value_v[i][l];
		rhs_cont -= dt * (vx_value[l] * vx_gradient[l][0] + 
				  vy_value[l] * vx_gradient[l][1]) * basis_value_v[i][l];
		rhs_cont *= Jxw;
		rhs(element_dof_v[i]) += rhs_cont;

		rhs_cont = fy_value[l] * basis_value_v[i][l] + vy_value[l] * basis_value_v[i][l];
		rhs_cont -= dt * (vx_value[l] * vy_gradient[l][0] + 
				  vy_value[l] * vy_gradient[l][1]) * basis_value_v[i][l];
		rhs_cont *= Jxw;
		rhs(n_dof_v + element_dof_v[i]) += rhs_cont;
	    }
	}
    }

    /// 构建系数矩阵和右端项.
    /// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
    Vector<double> x(n_total_dof);

    for (int i = 0; i < n_dof_v; ++i)
    {
    	x(i) = v_h[0](i);
    	x(i + n_dof_v) = v_h[1](i);
    }

    for (int i = 0; i < n_dof_p; ++i)
    	x(i + 2 * n_dof_v) = p_h(i);

    /// 边界条件一起处理了. 这里需要传递 x 因为 x 是临时的. 这里似乎应
    /// 该把 v_h, p_h 和 x 统一起来, 避免冗余错误.
    boundaryValueStokes(x);

    clock_t t_cost = clock();
    dealii::SolverControl solver_control (4000000, l_Euler_tol, check);

    /// 对于时间问题, 在没有解决好基础光滑算子前, 所有的预处理效果都不
    /// 好, 甚至 MIC 也失效. 奇怪.
    // StokesPreconditioner preconditioner;
    // /// 预处理矩阵.
    // SparseMatrix<double> matrix_vxvx(sp_vxvx);
    // SparseMatrix<double> matrix_vyvy(sp_vyvy);
    // /// 这里从 Stokes 取是因为加了边界条件.
    // for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
    // 	matrix_vxvx.global_entry(i) = matrix.global_entry(index_vxvx[i]);
    // for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
    // 	matrix_vyvy.global_entry(i) = matrix.global_entry(index_vyvy[i]);
    // preconditioner.initialize(matrix_vxvx, matrix_vyvy, mat_p_mass);
    // SparseILU<double> preconditioner;
    // preconditioner.initialize(matrix);//, SparseMIC<double>::AdditionalData());
    // SparseILU<double> precondition_ilu;
    // precondition_ilu.initialize(stiff_matrix);
    // 矩阵求解. 
    //	dealii::SolverControl solver_control (40000, G_tol);
    // SolverGMRES<Vector<double> > gmres (solver_control);
    // gmres.solve (matrix,  x,  rhs, preconditioner);
    // minres.solve (matrix, x, rhs, preconditioner);
    /// 不完全LU分解.     
    // dealii::SparseILU <double> preconditioner;
    // preconditioner.initialize(matrix);

    /// 因为时间步长小, 预处理做的不好还不如不做.
    SolverMinRes<Vector<double> > minres(solver_control);
    minres.solve (matrix, x, rhs, PreconditionIdentity());
    t_cost = clock() - t_cost;

    std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;

    /// 将整体数值解分割成速度和压力.
    for (int i = 0; i < n_dof_v; ++i)
    {
	v_h[0](i) = x(i);
	v_h[1](i) = x(i + n_dof_v);
    }
    for (int i = 0; i < n_dof_p; ++i)
	p_h(i) = x(i + 2 * n_dof_v);
};

#undef DIM
