#include "Grid3f.h"
#include <mkl_blas.h>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <assert.h>
#include <algorithm>

using namespace std;

Grid3f::Grid3f(void) : origin_()
{
#ifdef _USE_VECTORS
	data_.clear();
#else
	data_ = 0;
#endif
	m_ = n_ = o_ = mno_ = 0;
	dx_ = dy_ = dz_ = 1;
}

Grid3f::~Grid3f(void)
{
#ifdef _USE_VECTORS
	data_.clear();
#else
	if(data_)
		delete[] data_;
#endif
}

Grid3f::Grid3f(int m, int n, int o) {
	m_ = m;
	n_ = n;
	o_ = o;
	mno_ = m_*n_*o_;
#ifdef _USE_VECTORS
	data_.resize(mno_);
#else
	data_ = new float[mno_];
#endif
	dx_ = dy_ = dz_ = 1;

	/// Origin is already init to (0,0,0)
}

Grid3f::Grid3f(int m, int n, int o, const Vector3f &origin) : origin_(origin)
{
	m_ = m;
	n_ = n;
	o_ = o;
	mno_ = m_*n_*o_;
#ifdef _USE_VECTORS
	data_.resize(mno_);
#else
	data_ = new float[mno_];
#endif
	dx_ = dy_ = dz_ = 1;
	//fprintf(stderr,"Grid: set to (%f,%f,%f)\n",origin_.x(),origin_.y(), origin_.z());
}

Grid3f::Grid3f(const Grid3f &g) {
	operator=(g);
}

void Grid3f::operator =(const Grid3f &g) {
#ifdef _USE_VECTORS
	data_ = g.data_;
#else
	if(data_)
		delete[] data_;
#endif
	m_ = g.m_;
	n_ = g.n_;
	o_ = g.o_;
	dx_ = g.dx_;
	dy_ = g.dy_;
	dz_ = g.dz_;
	mno_ = m_*n_*o_;
	origin_ = g.origin_;
#ifndef _USE_VECTORS
	data_ = new float[mno_];
	memcpy (data_, g.data_, (mno_ * sizeof(float)));
#endif
}

void Grid3f::reallocate() {
	mno_ = m_ * n_ * o_;
#ifdef _USE_VECTORS
	data_.resize(mno_);
#else
	if(data_)
		delete[] data_;
	data_ = new float[mno_];
#endif
}

void Grid3f::setAllOne() {
	setAll(1);
}

void Grid3f::setAllZero() {
	setAll(0);
}

void Grid3f::setAll(float value) {
	for(int i=0 ; i<mno_ ; i++)
		data_[i] = value;
}

bool Grid3f::isEmpty() {
#ifdef _USE_VECTORS
	return (data_.size() == 0);
#else
	return (data_ == 0);
#endif
}

float* Grid3f::data() {
#ifdef _USE_VECTORS
	return &data_[0];
#else
	return data_;
#endif
}

int& Grid3f::m() {
	return m_;
}

int& Grid3f::n() {
	return n_;
}

int& Grid3f::o() {
	return o_;
}

int Grid3f::dimension(int i) {
	int *p=&m_;
	return *(p+i);
}

float& Grid3f::dx() {
	return dx_;
}

float& Grid3f::dy() {
	return dy_;
}

float& Grid3f::dz() {
	return dz_;
}

int Grid3f::mno() {
	return mno_;
}

Vector3f& Grid3f::origin() {
	return origin_;
}

void Grid3f::getMinMax(float &vmin, float &vmax) {
	const float EPSILON = 0.0001f;
	vmin = 1 / EPSILON;
	vmax = 0;
	for ( int i = 0; i < mno_; i++ ) {
		vmin = min<float>(vmin, data_[i]);
		vmax = max<float>(vmax, data_[i]);
	}
}

void Grid3f::writeMatlab(char *filename, char *variable) {
	FILE *file = fopen(filename, "w");
	for ( int z = 0; z < o_; z++ ) {
		fprintf(file,"%s(:,:,%d) = [",variable, z + 1);
		for ( int x = 0; x < m_; x++ ) {
			for ( int y = 0; y < n_; y++ ) {
				fprintf(file, "%f ",operator()(x, y, z));
			}
			if(x < m_ - 1) fprintf(file, ";");
		}
		fprintf(file,"];\n");
	}
	fclose(file);
}