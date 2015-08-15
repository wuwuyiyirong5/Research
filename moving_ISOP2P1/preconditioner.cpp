#include "preconditioner.h"
#define DIM 2

/** 
 * 实际估值 dst = M^{-1}src. 
 * 
 * @param dst 
 * @param src 
 */
void NavierStokesPreconditioner::vmult (Vector<double> &dst,
					const Vector<double> &src) const
{
  int n_dof_v = Axx->n();
    int n_dof_p = Q->n();
    Vector<double> d0(n_dof_v);
    Vector<double> d1(n_dof_v);
    Vector<double> d2(n_dof_p);
    Vector<double> s0(n_dof_v);
    Vector<double> s1(n_dof_v);
    Vector<double> s2(n_dof_p);

    
    Vector<double> t(n_dof_p);
    
    for (int i = 0; i < n_dof_v; ++i)
	s0(i) = src(i);
    for (int i = 0; i < n_dof_v; ++i)
	s1(i) = src(n_dof_v + i);
    for (int i = 0; i < n_dof_p; ++i)
	s2(i) = src(2 * n_dof_v + i);

   
    AMGSolver solverQ(*Q);
    solverQ.solve(d2, s2, 1e-8);
//    std::cout << "AMG 1" << std::endl;

//    s2 = t;
    
//    Fp->vmult(s2, t);

	
//    AMGSolver solverAp(*Ap);
//    solverAp.solve(d2, s2, 1e-8);
//    std::cout << "AMG 2" << std::endl;


    BTx->vmult(d0, d2);
    BTy->vmult(d1, d2);

    s0 -= d0;
    s1 -= d1;

//    AMGSolver solverA(*Axx, 1e-3, 3, 50, 0.382, 0.25);
    AMGSolver solverAxx(*Axx);
    solverAxx.solve(d0, s0, 1e-8, 1, 1);	
    AMGSolver solverAyy(*Ayy);
    solverAxx.solve(d1, s1, 1e-8, 1, 1);	


    for (int i = 0; i < n_dof_v; ++i)
	dst(i) = d0(i);
    for (int i = 0; i < n_dof_v; ++i)
	dst(n_dof_v + i) = d1(i);
    for (int i = 0; i < n_dof_p; ++i)
	dst(2 * n_dof_v + i) = d2(i);

};

/** 
 * 实际估值 dst = M^{-1}src. 
 * 
 * @param dst 
 * @param src 
 */
void NSPreconditioner::vmult (Vector<double> &dst,
			      const Vector<double> &src) const
{
    int n_dof_v = Axx->n();
    int n_dof_p = Q->n();
    Vector<double> d0(n_dof_v);
    Vector<double> d1(n_dof_v);
    Vector<double> d2(n_dof_p);
    Vector<double> s0(n_dof_v);
    Vector<double> s1(n_dof_v);
    Vector<double> s2(n_dof_p);

    for (int i = 0; i < n_dof_v; ++i)
	s0(i) = src(i);
    for (int i = 0; i < n_dof_v; ++i)
	s1(i) = src(n_dof_v + i);
    for (int i = 0; i < n_dof_p; ++i)
	s2(i) = src(2 * n_dof_v + i);
   
    AMGSolver solverQ(*Q);
    solverQ.solve(d2, s2, 1e-8);

    BTx->vmult(d0, d2);
    BTy->vmult(d1, d2);

    s0 -= d0;
    s1 -= d1;

    AMGSolver solverAxx(*Axx);
    solverAxx.solve(d0, s0, 1e-8, 1, 1);	
    AMGSolver solverAyy(*Ayy);
    solverAxx.solve(d1, s1, 1e-8, 1, 1);	

    for (int i = 0; i < n_dof_v; ++i)
	dst(i) = d0(i);
    for (int i = 0; i < n_dof_v; ++i)
	dst(n_dof_v + i) = d1(i);
    for (int i = 0; i < n_dof_p; ++i)
	dst(2 * n_dof_v + i) = d2(i);

};


void StokesPreconditioner::vmult (Vector<double> &dst,
				  const Vector<double> &src) const
{
    int n_dof_v = Ax->n();
    int n_dof_p = Q->n();
    Vector<double> d0(n_dof_v);
    Vector<double> d1(n_dof_v);
    Vector<double> s0(n_dof_v);
    Vector<double> s1(n_dof_v);

    for (int i = 0; i < n_dof_v; ++i)
	s0(i) = src(i);
    for (int i = 0; i < n_dof_v; ++i)
	s1(i) = src(n_dof_v + i);
    for (int i = 0; i < n_dof_p; ++i)
	dst(2 * n_dof_v + i) = src(2 * n_dof_v + i) / (*Q).diag_element(i);

    AMGSolver solver0(*Ax);
    solver0.solve(d0, s0, 1e-8, 1, 1);	

    AMGSolver solver1(*Ay);
    solver0.solve(d1, s1, 1e-8, 1, 1);	

    for (int i = 0; i < n_dof_v; ++i)
	dst(i) = d0(i);
    for (int i = 0; i < n_dof_v; ++i)
	dst(n_dof_v + i) = d1(i);
};

#undef DIM

