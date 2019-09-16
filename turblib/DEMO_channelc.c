//	Copyright 2011 Johns Hopkins University
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.


#include <stdio.h>
#include <float.h>
#include <math.h>
#include "turblib.h"

/*
 * Turbulence Database sample C client code
 */
#define N 10

int main(int argc, char *argv[]) {

	char * authtoken = "edu.jhu.pha.turbulence.testing-201406";
	char * dataset = "channel";
	enum SpatialInterpolation spatialInterp = Lag6;
	enum TemporalInterpolation temporalInterp = NoTInt;

	float time = 0.364F;

	float points[N][3];    /* input of x,y,z */
	float result1[N];      /* input of x,y,z */
	float result2[N][2];   /* results of invariant */
	float result3[N][3];   /* results of x,y,z */
	float result4[N][4];   /* results of x,y,z,p */
	float result6[N][6];   /* results from Pressure Hessian queries */
	float result9[N][9];   /* results from Velocity Gradient queries */
	float result18[N][18];
	int p;

	int x_start=1, y_start=1, z_start=1, x_end=4, y_end=4, z_end=4;

	char *threshold_field = "vorticity";
	ThresholdInfo *threshold_array;      /* dynamic array for the results of Threshold queries */
	int threshold_array_size;            /* size of the threshold array */
	float threshold = 70.0f;

	/* Initialize gSOAP */
	soapinit();

	/* Enable exit on error.  See README for details. */
	turblibSetExitOnError(1);

	for (p = 0; p < N; p++) {
		points[p][0] = (float)rand() / RAND_MAX * 8 * 3.141592F;
		points[p][1] = (float)rand() / RAND_MAX * 2 - 1;
		points[p][2] = (float)rand() / RAND_MAX * 3 * 3.141592F;
	}

	printf("\nCoordinates of %d points where variables are requested:\n", N);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, points[p][0], points[p][1], points[p][2]);
	}

	printf("\nRequesting velocity at %d points...\n", N);
	getVelocity(authtoken, dataset, time, spatialInterp, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting velocity and pressure at %d points...\n", N);
	getVelocityAndPressure(authtoken, dataset, time, spatialInterp, temporalInterp, N, points, result4);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e, p=%13.6e\n", p, result4[p][0], result4[p][1], result4[p][2], result4[p][3]);
	}

	printf("\nRequesting velocity gradient at %d points...\n", N);
	getVelocityGradient(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result9);
	for (p = 0; p < N; p++) {
		printf("%d: duxdx=%13.6e, duxdy=%13.6e, duxdz=%13.6e, ", p, result9[p][0], result9[p][1], result9[p][2]);
		printf("duydx=%13.6e, duydy=%13.6e, duydz=%13.6e, ", result9[p][3], result9[p][4], result9[p][5]);
		printf("duzdx=%13.6e, duzdy=%13.6e, duzdz=%13.6e\n", result9[p][6], result9[p][7], result9[p][8]);
	}

	printf("\nRequesting velocity hessian at %d points...\n", N);
	getVelocityHessian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result18);
	for (p = 0; p < N; p++) {
		printf("%d: d2uxdxdx=%13.6e, d2uxdxdy=%13.6e, d2uxdxdz=%13.6e, ", p, result18[p][0], result18[p][1], result18[p][2]);
		printf("d2uxdydy=%13.6e, d2uxdydz=%13.6e, d2uxdzdz=%13.6e, ", result18[p][3], result18[p][4], result18[p][5]);
		printf("d2uydxdx=%13.6e, d2uydxdy=%13.6e, d2uydxdz=%13.6e, ", result18[p][6], result18[p][7], result18[p][8]);
		printf("d2uydydy=%13.6e, d2uydydz=%13.6e, d2uydzdz=%13.6e, ", result18[p][9], result18[p][10], result18[p][11]);
		printf("d2uzdxdx=%13.6e, d2uzdxdy=%13.6e, d2uzdxdz=%13.6e, ", result18[p][12], result18[p][13], result18[p][14]);
		printf("d2uzdydy=%13.6e, d2uzdydz=%13.6e, d2uzdzdz=%13.6e\n", result18[p][15], result18[p][16], result18[p][17]);
	}

	printf("\nRequesting velocity laplacian at %d points...\n", N);
	getVelocityLaplacian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: grad2ux=%13.6e, grad2uy=%13.6e, grad2uz=%13.6e\n",
			p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting pressure at %d points...\n", N);
	getPressure(authtoken, dataset, time, spatialInterp, temporalInterp, N, points, result1);
	for (p = 0; p < N; p++) {
		printf("%d: p=%13.6e\n", p, result1[p]);
	}

	printf("\nRequesting pressure gradient at %d points...\n", N);
	getPressureGradient(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: dpdx=%13.6e, dpdy=%13.6e, dpdz=%13.6e\n", p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting pressure hessian at %d points...\n", N);
	getPressureHessian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result6);
	for (p = 0; p < N; p++) {
		printf("%d: d2pdxdx=%13.6e, d2pdxdy=%13.6e, d2pdxdz=%13.6e, d2pdydy=%13.6e, d2pdydz=%13.6e, d2pdzdz=%13.6e\n", p,
			result6[p][0], result6[p][1], result6[p][2], result6[p][3], result6[p][4], result6[p][5]);

	}

	printf("\nRequesting invariant at %d points...\n", N);
	getInvariant(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result2);
	for (p = 0; p < N; p++) {
		printf("%d: S2=%13.6e, O2=%13.6e\n", p, result2[p][0], result2[p][1]);
	}

	printf("\nRequesting threshold...\n");
	//NOTE: The array storing the results is dynamically allocated inside the getThreshold function,
	//because it's size is not known. It needs to be freed after it has been used to avoid leaking the memory.
	getThreshold(authtoken, dataset, threshold_field, time, threshold, FD4NoInt, x_start, y_start, z_start, x_end, y_end, z_end,
		&threshold_array, &threshold_array_size);
	for (p = 0; p < threshold_array_size; p++) {
		printf("(%d, %d, %d): %13.6e\n", threshold_array[p].x, threshold_array[p].y, threshold_array[p].z,
			threshold_array[p].value);
	}
	// Free the threshold array after using it.
	free(threshold_array);

	/* Free gSOAP resources */
	soapdestroy();

	return 0;
}
