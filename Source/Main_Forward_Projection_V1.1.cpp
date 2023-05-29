// Forward_Projection_V1.1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// Version 1.1
// Author Du


/* ********** Log ************ 
* 2023.5.16
* 针对单材质的正投影计算进行优化。有栅情形，对于密度图像的积分，优化为只进行一次，因其不随能量改变而变化；对于有栅添加无探元响应的情况
* 2023.5.24
* 无栅情形，对于密度图像的积分，优化为只进行一次；对于无栅添加无探元响应的情况
* 
*/



#include <iostream>

#include "ForwardProjection.h"

int main()
{
    ForwardProjection mForwardProjection;
    //mForwardProjection.forwardPolyProjGrid();
    //mForwardProjection.forwardPolyProjNoGrid();
    //mForwardProjection.forwardSinMatPolyProjGrid();
    mForwardProjection.forwardSinMatPolyProjGridNoResponse();
    //mForwardProjection.forwardSinMatPolyProjNoGrid();
    //mForwardProjection.forwardSinMatPolyProjNoGridNoResponse();
    
    //mForwardProjection.forwarProjNoGridTestI0();
}


