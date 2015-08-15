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
    /// 构建矩阵结构.
    buildMatrixStruct();

    buildMatrix();
};

void ISOP2P1::run()
{
	initialize();
	solveStokes();
	outputTecplotP("P0");
	outputTecplot("V0");
};

#undef DIM
