#include <AFEPack/AMGSolver.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <vector>
#define DIM 2

/**
 * 对应的预处理.
 * 
 */
class StokesPreconditioner
{
private:
    const SparseMatrix<double> *Ax; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *Ay;
    const SparseMatrix<double> *Q;

public:
    StokesPreconditioner()
    {};

    ~StokesPreconditioner()
    {};

    /** 
     * 预处理子初始化.
     * 
     * @param _stiff_vx vx 空间的刚度矩阵. 
     * @param _stiff_vy vy 空间的刚度矩阵.
     * @param _mass_p_diag p 空间的质量矩阵的对角元. 
     */
    void initialize (const SparseMatrix<double> &_stiff_vx, 
		     const SparseMatrix<double> &_stiff_vy, 
		     const SparseMatrix<double> &_mass_p_diag) 
    {
	Ax = &_stiff_vx;
	Ay = &_stiff_vy;
	Q = &_mass_p_diag;
    };
    /** 
     * 实际估值 dst = M^{-1}src. 
     * 
     * @param dst 
     * @param src 
     */
    void vmult (Vector<double> &dst,
		const Vector<double> &src) const;
};

class NavierStokesPreconditioner
{
private:
    const SparseMatrix<double> *Axx; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *Ayy; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *Wxy; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *Wyx; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *BTx;
    const SparseMatrix<double> *BTy;
    const SparseMatrix<double> *Bx;
    const SparseMatrix<double> *By;
    const SparseMatrix<double> *Q;
    const SparseMatrix<double> *Ap;
    const SparseMatrix<double> *Fp;

public:
    NavierStokesPreconditioner()
    {};

    ~NavierStokesPreconditioner()
    {};

    /** 
     * 预处理子初始化.
     * 
     * @param _stiff_vx vx 空间的刚度矩阵. 
     * @param _stiff_vy vy 空间的刚度矩阵.
     * @param _mass_p_diag p 空间的质量矩阵的对角元. 
     */
    void initialize (const SparseMatrix<double> &_vis_xx, 
		     const SparseMatrix<double> &_vis_yy, 
		     const SparseMatrix<double> &_J_xy, 
		     const SparseMatrix<double> &_J_yx, 
		     const SparseMatrix<double> &_divTx, 
		     const SparseMatrix<double> &_divTy,
 		     const SparseMatrix<double> &_mass_p, 
		     const SparseMatrix<double> &_stiff_p, 
		     const SparseMatrix<double> &_div_p)
    {
	Axx = &_vis_xx;
	Ayy = &_vis_yy;
	Wxy = &_J_xy;
	Wyx = &_J_yx;
	BTx = &_divTx;
	BTy = &_divTy;
	Q = &_mass_p;
	Ap = &_stiff_p;
	Fp = &_div_p;
    };

    /** 
     * 实际估值 dst = M^{-1}src. 
     * 
     * @param dst 
     * @param src 
     */
    void vmult (Vector<double> &dst,
		const Vector<double> &src) const;
};

class NSPreconditioner
{
private:
    const SparseMatrix<double> *Axx; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *Ayy; /**< 预处理矩阵各分块. */
    const SparseMatrix<double> *BTx;
    const SparseMatrix<double> *BTy;
    const SparseMatrix<double> *Q;

public:
    NSPreconditioner()
    {};

    ~NSPreconditioner()
    {};

    /** 
     * 预处理子初始化.
     * 
     * @param _stiff_vx vx 空间的刚度矩阵. 
     * @param _stiff_vy vy 空间的刚度矩阵.
     * @param _mass_p_diag p 空间的质量矩阵的对角元. 
     */
    void initialize (const SparseMatrix<double> &_vis_xx, 
		     const SparseMatrix<double> &_vis_yy, 
		     const SparseMatrix<double> &_divTx, 
		     const SparseMatrix<double> &_divTy,
 		     const SparseMatrix<double> &_mass_p)
    {
	Axx = &_vis_xx;
	Ayy = &_vis_yy;
	BTx = &_divTx;
	BTy = &_divTy;
	Q = &_mass_p;
    };

    /** 
     * 实际估值 dst = M^{-1}src. 
     * 
     * @param dst 
     * @param src 
     */
    void vmult (Vector<double> &dst,
		const Vector<double> &src) const;
};

#undef DIM
