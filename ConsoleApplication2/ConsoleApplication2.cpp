// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <ppl.h>
#include <ctime>
using namespace std;
const unsigned int n = 50;
const unsigned int m = 2500;

vector<vector<double>> matrixComposition(vector<vector<double>> first, vector<vector<double>> second) {
	if (first[0].size() != second.size())
		throw new exception("Матрицы нельзя перемножить");
	int rows = first.size();
	int columns = second[0].size();
	vector<vector<double>> result(rows, vector<double>(columns));
	Concurrency::parallel_for(0, rows, [&](int i) {
		for (int j = 0; j < columns; j++) {
			for (int g = 0; g < first[0].size(); g++)
				result[i][j] += first[i][g] * second[g][j];
		}
	});

	return result;
}
vector<vector<double>> transposeMatrix(vector<vector<double>> matrix) {
	double localVariable;
	if (matrix.size() != matrix[0].size())
		throw new exception("Матрицу не возможно транспонировать");
	for (int i =0; i<matrix.size(); i++)
		for (int j = i; j < matrix[0].size(); j++) {
			localVariable = matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = localVariable;
		}
	return matrix;
}
void showMatrix(vector<vector<double>> matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[0].size(); j++)
			cout << matrix[i][j] << ' ';
		cout << endl;
	}
}
double getNorm(vector<vector<double>> vec) {
	double norm = 0;
	for (int i = 0; i < vec.size(); i++)
		for (int j = 0; j < vec[0].size(); j++)
			norm += vec[i][j] * vec[i][j];
	return sqrt(norm);
}
vector<vector<double>> matrixSubstraction(vector<vector<double>> first, vector<vector<double>> second) {
	if (first.size() != second.size() || second[0].size() != first[0].size())
		throw new exception("Вычитание не возмжно, размеры матриц различны.");
	vector<vector<double>> result(first.size(), vector <double>(first[0].size()));

	for (int i = 0; i < result.size(); i++)
		for (int j = 0; j < result[0].size(); j++)
			result[i][j] = first[i][j] - second[i][j];
	return result;
}
vector<vector<double>> multiplyByConst(double value, vector<vector<double>> matrix) {
	for (int i = 0; i < matrix.size(); i++)
		for (int j = 0; j < matrix[0].size(); j++)
			matrix[i][j] *= value;
	return matrix;
}
double topNorm(vector<vector<double>> matrix) {
	double norm = 0;
	for (int i = 0; i < matrix.size(); i++)
		for (int j = i; j < matrix[0].size(); j++)
			norm += matrix[i][j] * matrix[i][j];
	return sqrt(norm);
}
double botNorm(vector<vector<double>> matrix) {
	double norm = 0;
	for (int i = 0; i < matrix.size(); i++)
		for (int j = 0; j <= i; j++)
			norm += matrix[i][j] * matrix[i][j];
	return sqrt(norm);
}
int main()
{
	double epsilon = 0.0005;
	vector<double> x(n);
	vector<double> y(n);
	vector<vector<double>> h1(n, vector<double>(n));
	vector<vector<double>> h2(n, vector<double>(n));
	vector<vector<double>> B(n, vector<double>(n));
	vector<vector<double>> vectorB(m, vector<double>(1));
	vector<vector<double>> A(m, vector<double>(m));
	unsigned int startTime = clock(); // начальное время

	ifstream hh1File("hh1.dat");
	ifstream hh12File("hh2.dat");
	ifstream f1SqFile("f_1sq.dat");
	double local1, local2;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			hh1File >> x[i]; hh1File >> y[j]; hh1File >> h1[i][j];
			hh12File >> x[i]; hh12File >> y[j]; hh12File >> h2[i][j];
			f1SqFile >> x[j]; f1SqFile >> y[i]; f1SqFile >> B[j][i];
		}
	}
	hh1File.close();
	hh12File.close();
	y = x;
	f1SqFile.close();
	double deltaX = x[1] - x[0];
	double deltaY = y[1] - y[0];

	double f = 0.00667;
	int s = 0;
	for (int k = 0; k < n; k++) {
		for (int l = 0; l < n; l++) {
			int t = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					A[50 * (k) + l][50 * (i) + j] = f*(1.0 / sqrt((x[k] - x[i])*(x[k] - x[i]) + (y[l] - y[j])*(y[l] - y[j]) + h1[i][j]*h1[i][j])
						- 1.0 / sqrt((x[k] - x[i])*(x[k] - x[i]) + (y[l] - y[j])*(y[l] - y[j]) + h2[i][j]*h2[i][j]))*
						deltaX*deltaY;
					if (k == 0 && l == 0)
						vectorB[t][0] = B[i][j];
					t += 1;
				}
			}
			s += 1;
		}
	}
	cout << "Top: " << topNorm(A) << endl << "Bot: " << botNorm(A) << endl;
	vector<vector<double>> transposeA = transposeMatrix(A); //траспонированная А
	vector<vector<double>> AonTransposeA = matrixComposition(A, transposeA); // A*Atransp

	cout << "Matrix was found" << endl;
	double check1 = 0;
	vector<vector<double>> z(m, vector<double>(1));
	int iterationCount = 0;
	while (true) {
		vector<vector<double>> zNew(m, vector<double>(1));
		vector<vector<double>> valueInBrackets = matrixSubstraction(matrixComposition(A, z), vectorB);
		check1 = getNorm(valueInBrackets) / getNorm(vectorB);
		cout << iterationCount << " " << check1 << endl;
		if (check1 < epsilon) {
			break;
		}

		iterationCount++;
		vector<vector<double>> topValue = matrixComposition(transposeA, valueInBrackets);
		double norm1 = getNorm(topValue);
		norm1 = pow(norm1, 2);
		double norm2 = getNorm(matrixComposition(AonTransposeA,valueInBrackets));
		norm2 = pow(norm2, 2);
		zNew = matrixSubstraction(z, multiplyByConst(norm1 / norm2, topValue));
		z = zNew;
	}

	ofstream myFile;
	string fileName = "out1.dat";
	myFile.open(fileName.c_str());
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			myFile << x[i] << " " << y[j] << " " << z[i + j*n][0] << "\n";
	myFile.close();
	unsigned int endTime = clock();
	unsigned int searchTime = endTime - startTime;
	int seconds = searchTime / 1000;
	int minutes = seconds / 60;
	cout << "Elapsed time: " << minutes << " minutes " << seconds % 60 << " seconds" << endl;
	return 0;
}

