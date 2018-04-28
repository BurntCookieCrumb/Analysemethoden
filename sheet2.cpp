
/*-------------------------------------------------------------------------------------------------------*/
// Analysemethoden der Hochenergiephysik
// Übungsblatt 3

// bearbeitet von Kristina Schmitt
/*-------------------------------------------------------------------------------------------------------*/

// Includes

#include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <cmath>

// Namespace

using namespace std;

// Arrays x and y

#define arrSize 10

double x[arrSize] = { 0.62, 1.25, 1.91, 3.01, 4.22, 5.89, 6.12, 6.86, 8.07, 9.43 };
double y[arrSize] = { 0.572, 0.663, 1.32, 1.79, 2.71, 3.52, 3.49, 4.0, 4.84, 5.58 };

double x2[arrSize] = { 1., 3., 7., 8., 2., 5., 4., 6., 10., 9. };
double y2[arrSize] = { 1., 3., 7., 8., 2., 5., 4., 6., 10., 9. };

double x3[arrSize] = { 0.62, 6.86, 8.07, 9.43, 1.25, 1.91, 3.01, 4.22, 5.89, 6.12 };
double y3[arrSize] = { 0.572, 4.0, 4.84, 5.58, 0.663, 1.32, 1.79, 2.71, 3.52, 3.49 };

/*-------------------------------------------------------------------------------------------------------*/
//		Exercise 1
/*-------------------------------------------------------------------------------------------------------*/

// Function to print (x,y)

void printXY(double *x,  double *y, int numElements) {

	cout << " Value Pairs: " << endl << endl;

	for (int i = 0; i < numElements; i++) {

		cout << " (" << left << setw(5) << setprecision(3) << x[i] << " | " 
			         << left << setw(5) << setprecision(3) << y[i] << ")" << endl;

	}

}

// Print Mean and Standard Deviation

void printMeanSigmaYbyX(double *x, double *y, int numElements) {

	double YbyX;

	double sum = 0;
	double sumSquared = 0;

	double mean;
	double meanSquared;
	double sigma;

	for (int i = 0; i < numElements; i++) {

		YbyX = y[i] / x[i];

		sum += YbyX;
		sumSquared += pow(YbyX, 2);

	}

	mean = sum / numElements;
	meanSquared = sumSquared / numElements;
	sigma = sqrt(fabs(meanSquared - pow(mean,2)));

	cout << endl << "Mean of y/x: " << mean << " and Standard Deviation: " << sigma << endl << endl;

}

void Exercise_2_1(void) {

	cout << endl  << " -- Exercise 1: ---------------------- " << endl << endl;
	printXY(x, y, arrSize);

	printMeanSigmaYbyX(x, y, arrSize);
	

}

/*-------------------------------------------------------------------------------------------------------*/
//		Alternative Exercise
/*-------------------------------------------------------------------------------------------------------*/

// Route

typedef struct pointRoute {

	int pointOrder[arrSize];
	double length;

} route;

// help functions

bool contains(int *array, int numElements, int elem) {

	for (int i = 0; i < numElements; i++) {

		if (array[i] == elem) return true;

	}

	return false;

}

void printRoute(route Route, double* x, double* y) {

	cout << "Route: " << endl << endl;

	int elem;

	for (int i = 0; i < arrSize; i++) {

		elem = Route.pointOrder[i];

		cout << " (" << left << setw(5) << setprecision(3) << x[elem] << " | "
			<< left << setw(5) << setprecision(3) << y[elem] << ")" << endl;

	}

	cout << endl << "Length of Route: " << setprecision(5) << Route.length << endl << endl << endl;
}

// Calculate shortest route

route getShortestRoute(double *x, double *y, int numElements) {

	route shortestRoute;

	double shortestLength = 0;
	double currentLength;

	double currentDist;
	double shortestDist;

	double currentX;
	double currentY;

	double nextX;
	double nextY;

	int currentPointOrder[arrSize];
	int pointOrder[arrSize];

	// Test all points as starting point

	for (int i = 0; i < numElements; i++){

		// reset variables		

		currentLength = 0;

		// set start point

		nextX = x[i];
		nextY = y[i];

		currentPointOrder[0] = i;

		// From start point, calculate length to all other points

		for (int point = 1; point < numElements; point++) {

			currentX = nextX;
			currentY = nextY;
			
			shortestDist = 0;
			currentPointOrder[point] = -1;

			for (int j = 0; j < numElements; j++) {

				if (contains(currentPointOrder, numElements, j)) continue;

				currentDist = sqrt(pow(currentX - x[j], 2) + pow(currentY - y[j], 2));

				// update current clostest point
				if (currentDist < shortestDist || shortestDist == 0) {
					shortestDist = currentDist;
					currentPointOrder[point] = j;
					nextX = x[j];
					nextY = y[j];

				}				
			}

			currentLength += shortestDist;

		}

		// update current route
		if (shortestLength == 0 || currentLength < shortestLength) {

			shortestLength = currentLength;
			for (int elem = 0; elem < numElements; elem++) {
				pointOrder[elem] = currentPointOrder[elem];
			}

		}

	}

	shortestRoute.length = shortestLength;
	for (int elem = 0; elem < numElements; elem++) {
		shortestRoute.pointOrder[elem] = pointOrder[elem];
	}


	return shortestRoute;

}


void Exercise_2_2(void){

	cout << endl << " -- Exercise 2: ---------------------- " << endl << endl;

	cout << "--> shortest route based on given x y pairs from exercise 1" << endl << endl;
	printRoute(getShortestRoute(x, y, arrSize), x, y);

	cout << "--> shortest route based on x y pairs on a straight line" << endl << endl;
	printRoute(getShortestRoute(x2, y2, arrSize), x2, y2);

}


/*
int main()
{

	Exercise_2_1();
	Exercise_2_2();

	return 0;
}
*/