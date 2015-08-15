/**
 * @file   ISOP2P1.cpp
 * @author Heyu Wang <scshw@cslin107.csunix.comp.leeds.ac.uk>
 * @date   Tue Nov  4 14:33:17 2014
 * 
 * @brief ISO P2P1 有限元的网格, 空间, 矩阵信息集中封装类的实现.
 * 
 * 
 */

#include "ISOP2P1.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::initialize()
{
    /// 读取配置文件.
    config("config");
    /// 构建网格.
    buildMesh();
    /// 构建混合有限元 ISOP1P2 空间. 
    buildFEMSpace();
    buildMatrixStruct();
  
};

void ISOP2P1::run()
{
    initialize();
    t = 0.05;
    dt = 0.002;
    if (isMoving)
    {
	std::cout << "Initialize mesh ... " << std::endl;
	double scale, scale_step = 0.2;
	scale = scale_step;
	do {
	    initialValue();
	    p_h.scale(scale);
	    movingMesh();
	    std::cout << "\r\tscale = " << scale << std::endl;
	    scale += scale_step;
	} while (scale <= 1.0);
    }
    initialValue();
    outputSolution();
    outputTecplotP("P0");
    outputTecplot("initial_value0");
    getchar();
    do {
	stepForward();
	if (isMoving)
	    movingMesh();
	outputSolution();
	std::cout << "t  = " << t << std::endl;
    } while (t < 1.0);
    outputTecplotP("P1");
    outputTecplot("V1");
};

void ISOP2P1::initialValue()
{
    BurgV u(viscosity, t);
    Operator::L2Project(u, p_h, Operator::LOCAL_LEAST_SQUARE, 3);
};

#undef DIM
