/*******************************************************************
Matrix.h

Contains the decleration and the implementation of the Matrix Class
*******************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include <complex>
#include <iostream>

#include "iElements.h"

template<class T> class Matrix
{
private:
	T **data;
	unsigned int nRows, nColumns;

public:
	//Constructor and Destructor of the Matrix Class
	Matrix(unsigned int nRows, unsigned int nColumns);
	Matrix(unsigned int nRows, unsigned int nColumns, T **data);
	~Matrix();

	//Overloaded operators
	void operator = (const Matrix& b);
	Matrix<T> operator + (const Matrix& b);
	Matrix<T> operator - (const Matrix& b);
	Matrix<T> operator * (const Matrix& b);
	Matrix<T> operator * (const T& num);
	Matrix<T> operator ^ (const int& k);

	//Basic functions
	bool isSquare();
	bool isInversible();
	void LU_Decomposition(Matrix<T>& L, Matrix<T>& U);
	T determinant();
	Matrix<T> inv();

	//Supplementary functions
	int rows() { return this->nRows; }
	int columns() { return this->nColumns; }
	T element(int i, int j) { return this->data[i - 1][j - 1]; }
	T** elements() { return this->data }
	void print();
	void print(int m, int n);
};

template <class T>
Matrix<T>::Matrix(unsigned int nRows, unsigned int nColumns)
{
	int i, j;

	this->nRows = nRows;
	this->nColumns = nColumns;

	this->data = new T*[nRows];
	for (i = 0; i < nRows; i++)
		this->data[i] = new T[nColumns];
}

template <class T>
Matrix<T>::Matrix(unsigned int nRows, unsigned int nColumns, T **data)
{
	int i, j;

	this->nRows = nRows;
	this->nColumns = nColumns;

	this->data = new T*[nRows];
	for (i = 0; i < nRows; i++)
		this->data[i] = new T[nColumns];

	memcpy(this->data, data, nRows*nColumns*sizeof(T));
}

template <class T>
Matrix<T>::~Matrix()
{
	/*int i;
	for (int i = 0; i < nRows; i++)
	delete[] data[i];

	delete[] data;*/
}


template <class T>
void Matrix<T>::operator = (const Matrix& b)
{
	int i;

	this->nRows = b.nRows;
	this->nColumns = b.nColumns;

	if (data != NULL)
	{
		for (i = 0; i < nRows; i++)
			delete[] this->data[i];
		delete[] this->data;
	}

	this->data = new T*[nRows];
	for (i = 0; i < nRows; i++)
		this->data[i] = new T[nColumns];

	memcpy(this->data, b.data, nRows*nColumns*sizeof(T));
}

template <class T>
Matrix<T> Matrix<T>::operator + (const Matrix& b)
{
	/* Exception handling */
	if (nRows != b.nRows || nColumns != b.nColumns)
	{
		cerr << "Matrices do not have the same dimensions" << endl;
		return Matrix(0, 0);
	}

	int i, j;
	T **_data;

	_data = new T*[nRows];
	for (i = 0; i < nRows; i++)
		_data[i] = new T[nColumns];

	for (i = 0; i < nRows; i++)
		for (j = 0; j < nColumns; j++)
			_data[i][j] = this->data[i][j] + b.data[i][j];

	Matrix<T> res(nRows, nColumns, _data);
	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator - (const Matrix& b)
{
	/* Exception handling */
	if (this->nRows != b.nRows || this->nColumns != b.nColumns)
	{
		cerr << "Matrices do not have the same dimensions" << endl;
		return Matrix(0, 0);
	}

	int i, j;
	T **_data;

	_data = new T*[nRows];
	for (i = 0; i < nRows; i++)
		_data[i] = new T[nColumns];

	for (i = 0; i < nRows; i++)
		for (j = 0; j < nColumns; j++)
			_data[i][j] = this->data[i][j] - b.data[i][j];

	Matrix<T> res(nRows, nColumns, _data);
	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const Matrix& b)
{
	/* Exception handling */
	if (this->nColumns != b.nRows)
	{
		cerr << "Cannot multiply matrices: number of Columns is different from the number of Rows" << endl;
		return Matrix(0, 0);
	}

	int i, j, k;
	Matrix<T> res(this->nRows, b.nColumns);

	for (i = 0; i < res.nRows; i++)
	{
		for (j = 0; j < res.nColumns; j++)
		{
			res.data[i][j] = 0;

			for (k = 0; k < nColumns; k++)
				res.data[i][j] += this->data[i][k] * b.data[k][j];
		}
	}

	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const T& num)
{
	int i, j;
	Matrix<T> res(this->nRows, this->nColumns);

	for (i = 0; i < res.nRows; i++)
		for (j = 0; i < res.nColumns; j++)
			res.data[i][j] = num*this->data;

	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator ^ (const int& k)
{
	/* Exception handling */
	if (this->nColumns != this->nRows)
	{
		cerr << "Cannot raise matrix to a power: not a square matrix" << endl;
		return Matrix<T>(0, 0);
	}

	int i, n = this->nRows;

	if (k > 0)
	{
		Matrix<T> res(n, n, data), curr(n, n, data);

		for (i = 1; i < k; i++)
			res = res * curr;

		return res;
	}

	else if (k == 0)
		return I(n);

	else
	{
		if (this->isInversible())
			return this->inv() ^ (-k);
		else
			return Matrix<T>(0, 0);
	}
}


template <class T>
bool Matrix<T>::isSquare()
{
	if (this->nRows == this->nColumns)
		return true;
	else
		return false;
}

template <class T>
bool Matrix<T>::isInversible()
{
	if (this->determinant() != 0)
		return true;
	else
		return false;
}

template <class T>
void Matrix<T>::LU_Decomposition(Matrix<T>& L, Matrix<T>& U) {

	/* Exception handling */
	if (this->nRows != this->nColumns) {
		cerr << "Can't use LU method because this matrix is not square";
		return;
	}

	int i, j, k;
	int n = this->nRows;

	/*
	The normal way to initialize
	U with 0 belows the diagonial and
	L with 0 above the diagonial and 1 on the diagonial

	However iElements does this and complets the extra spaces
	without, however, affecting the code of the program
	*/

	L.~Matrix();
	U.~Matrix();

	L = new Matrix<T>(n, n, iElements(n));
	U = new Matrix<T>(n, n, iElements(n));

	for (k = 0; k < n; k++) {
		U.data[k][k] = this->data[k][k];

		for (i = k; i < n; i++) {
			L.data[i][k] = this->data[i][k] / U.data[k][k];
			U.data[k][i] = this->data[k][i];
		}

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				this->data[i][j] -= L.data[i][k] * U.data[k][j];
	}
}

template <class T>
T Matrix<T>::determinant()
{
	/* Exception handling */
	if (this->nRows != this->nColumns) {
		cerr << "Can't find inverse matrix because this matrix is not square";
		return;
	}

	int i, n = this->nRows;
	T det = 1;
	Matrix L(n, n), U(n, n);

	LU_Decomposition(L, U);

	for (i = 0; i < n; i++)
		det *= L.elements[i][i] * U.elements[i][i];

	return det;
}

template <class T>
Matrix<T> Matrix<T>::inv()
{
	/* Except handling */
	if (this->nRows != this->nColumns) {
		cerr << "Can't find inverse because this matrix is not square";
		return;
	}

	int i, j, n = this->nRows;
	Matrix res(n, n);

	/* Determinant of a 2x2 Matrix */
	if (n == 2) {
		res.data[0][0] = this->data[1][1];
		res.data[1][1] = this->data[0][0];
		res.data[0][1] = -this->data[0][1];
		res.data[1][0] = -this->data[1][0];

		return Ans * (1 / this->determinant());
	}

	/* Determinant of higher dimensions */
	int k, count;
	Matrix<T> M(n - 1, n - 1);

	for (k = 0; k < n; k++) {

		/* Form the (N-1)x(N-1) Matrices */
		for (i = 0; i < this->nColumnss - 1; i++) {
			for (j = 0, count = 0; j < this->nRows; j++) {

				if (j != k) {
					M.data[i][count] = this->data[i + 1][j];
					count++;
				}

			}
		}

		if (count % 2 == 0)
			res.data[i][j] = M.determinant();
		else
			res.data[i][j] = -M.determinant();

	}

	return Ans * (1 / this->Determinant());
}

template <class T>
void Matrix<T>::print()
{
	int i, j;

	for (i = 0; i < this->nRows; i++)
	{
		for (j = 0; j < this->nColumns; j++)
			cout << this->data[i][j] << " ";
		cout << endl;
	}
}

template <class T>
void Matrix<T>::print(int m, int n)
{
	cout << this->data[n - 1][m - 1] << endl;
}


template <typename T>
Matrix<T> I(int n)
{
	Matrix<T> res(n, n, iElements(n));

	return res;
}

#endif
