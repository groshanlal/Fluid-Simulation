#ifndef __FS_COMMON_H__
#define __FS_COMMON_H__
#pragma warning(disable : 4996)

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
#include <mkl_blas.h>
#include <time.h>
#ifdef _WIN32
#include <conio.h>
#endif

using namespace std;

#define EPSILON 0.001f
#define SQR(x) ((x)*(x))
#define ROUND(x) ((int)floorf(x+0.5f))
#ifdef _DEBUG
#define SHOW_DEBUG_OUTPUT 1
#define DO_DEBUG_TESTS 1
#endif
#define frand	((float)rand()/ RAND_MAX)
#define frand2	(-0.5f + frand)

#endif