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
	char * dataset = "mhd1024";
	enum SpatialInterpolation spatialInterp = Lag6;
	enum TemporalInterpolation temporalInterp = NoTInt;

	float time = 0.364F;

	float points[N][3];    /* input of x,y,z */
	float result1[N];       /* results from GetPressure */
	float result2[N][2];   /* results of invariant */
	float result3[N][3];   /* results of x,y,z */
	float result4[N][4];   /* results of x,y,z velocity and pressure */
	float result6[N][6];   /* results from Pressure Hessian and SGS queries */
	float result9[N][9];   /* results from Gradient queries */
	float result18[N][18]; /* results from Hessian queries */
	int p;

	float startTime = 0.364F;
	float endTime = 0.376F;
	float lag_dt = 0.0004F;

	int x_start=1, y_start=1, z_start=1, x_end=4, y_end=4, z_end=4;

	char * field = "velocity"; /* field used for the calls to getBoxFilter, getBoxFilterSGSsymtensor
								and getBoxFilterGradient */
	char * scalar_fields = "pp"; /* two scalar fields ("p" and "p") used for the call to getBoxFilterSGSscalar */
	char * vector_scalar_fields = "up"; /* a vector and a scalar field ("u" and "p") used for the call to
										 getBoxFilterSGSvector */
	char * vector_fields = "ub"; /* two vector fields ("u" and "b") used for the call to getBoxFilterSGStensor */
	float dx = 2.0f * 3.14159265f / 1024.0f;
	float filterwidth = 7.0f * dx;
	float spacing = 4.0f*dx;

	char *threshold_field = "vorticity";
	ThresholdInfo *threshold_array;      /* dynamic array for the results of Threshold queries */
	int threshold_array_size;            /* size of the threshold array */
	float threshold = 10.0f;

	/* Initialize gSOAP */
	soapinit();

	/* Enable exit on error.  See README for details. */
	turblibSetExitOnError(1);

	/* If working with cutout files, CUTOUT_SUPPORT should be defined during compilation.
	 Make sure to run make with the "CUTOUT_SUPPORT=1" option, e.g.:
	 $ make mhdc CUTOUT_SUPPORT=1
	 Change the filename to the name of the downloaded cutout file or supply it at the command line.
	 */
#ifdef CUTOUT_SUPPORT
	if (argc > 1)
	{
		turblibAddLocalSource(argv[1]);
	}
	//turblibAddLocalSource("mhd1024.h5");
#endif

	/* The client library implements all of the server-side functionality "locally" (except
	 for particle tracking and filtering). Therefore, if an hdf5 file with cutout data is
	 available and loaded as above all queries for data that are within the region defined
	 in the file will be evaluated locally (without being sent to the server). An example
	 is provided below.

	 Please note that the use of this feature of the client library requires an hdf5
	 installation. The standard approach of simply executing queries through the server
	 does not require hdf5 or downloading any cutout data.

	 For users that make frequent calls for data in a particular region and would like to
	 download this data to their local machine the only change in their existing code
	 will be to use the turblibAddLocalSource function to load the hdf5 file that they
	 have downloaded. Each function call will determine whether the data is available
	 locally and will evaluate the query locally if it is and will make a call to the
	 Web-services on the server if it is not.

	 The steps below can be followed to download data cutouts in hdf5 format and use the
	 client library locally:
	 1) Download an hdf5 file containing a cubic region of the velocity data at timestep 0
	 with the following download link:
	 http://turbulence.pha.jhu.edu/download.aspx/[authorization Token]/mhd1024/u/0,1/0,16/0,16/0,16/
	 2) Compile this sample code with the "CUTOUT_SUPPORT=1" option, e.g.:
	 $ make mhdc CUTOUT_SUPPORT=1
	 3) Uncomment the line below to Load the hdf5 cutout file:
	 */

	 //  turblibAddLocalSource("mhd1024.h5");

	 /* 4) Uncomment the code below to restrict target locations to within the data region downloaded:
	  */

	  //  for (p = 0; p < N; p++) {
	  //    points[p][0] = (float)rand()/RAND_MAX*2*3.141592F/64.0F;
	  //    points[p][1] = (float)rand()/RAND_MAX*2*3.141592F/64.0F;
	  //    points[p][2] = (float)rand()/RAND_MAX*2*3.141592F/64.0F;
	  //  }

	  /*
	   5) Uncomment the code below to call getVelocity, which will be evaluated locally as
	   the time chosen is 0.0, which corresponds to timestep 0 and no spatial or temporal
	   interpolation is requested:
	   */

	   //  printf("\nRequesting velocity at %d points...\n", N);
	   //  getVelocity (authtoken, dataset, 0.0F, NoSInt, NoTInt, N, points, result3);
	   //  for (p = 0; p < N; p++) {
	   //    printf("%d: %13.6e, %13.6e, %13.6e\n", p, result3[p][0],  result3[p][1],  result3[p][2]);
	   //  }

	   /* In this sample code, the default time chosen is 0.364, so all of the function calls in the
		remainder of the sample code will be evaluated at the server.
		*/

	for (p = 0; p < N; p++) {
		points[p][0] = (float)rand() / RAND_MAX * 2 * 3.141592F;
		points[p][1] = (float)rand() / RAND_MAX * 2 * 3.141592F;
		points[p][2] = (float)rand() / RAND_MAX * 2 * 3.141592F;
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

	printf("\nRequesting forcing at %d points...\n", N);
	getForce(authtoken, dataset, time, spatialInterp, temporalInterp, N, points, result3);
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
		printf("%d: %13.6e\n", p, result1[p]);
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

	printf("\nRequesting magnetic field at %d points...\n", N);
	getMagneticField(authtoken, dataset, time, spatialInterp, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting magnetic field gradient at %d points...\n", N);
	getMagneticFieldGradient(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result9);
	for (p = 0; p < N; p++) {
		printf("%d: duxdx=%13.6e, duxdy=%13.6e, duxdz=%13.6e, ", p, result9[p][0], result9[p][1], result9[p][2]);
		printf("duydx=%13.6e, duydy=%13.6e, duydz=%13.6e, ", result9[p][3], result9[p][4], result9[p][5]);
		printf("duzdx=%13.6e, duzdy=%13.6e, duzdz=%13.6e\n", result9[p][6], result9[p][7], result9[p][8]);
	}

	printf("\nRequesting magnetic field hessian at %d points...\n", N);
	getMagneticFieldHessian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result18);
	for (p = 0; p < N; p++) {
		printf("%d: d2uxdxdx=%13.6e, d2uxdxdy=%13.6e, d2uxdxdz=%13.6e, ", p, result18[p][0], result18[p][1], result18[p][2]);
		printf("d2uxdydy=%13.6e, d2uxdydz=%13.6e, d2uxdzdz=%13.6e, ", result18[p][3], result18[p][4], result18[p][5]);
		printf("d2uydxdx=%13.6e, d2uydxdy=%13.6e, d2uydxdz=%13.6e, ", result18[p][6], result18[p][7], result18[p][8]);
		printf("d2uydydy=%13.6e, d2uydydz=%13.6e, d2uydzdz=%13.6e, ", result18[p][9], result18[p][10], result18[p][11]);
		printf("d2uzdxdx=%13.6e, d2uzdxdy=%13.6e, d2uzdxdz=%13.6e, ", result18[p][12], result18[p][13], result18[p][14]);
		printf("d2uzdydy=%13.6e, d2uzdydz=%13.6e, d2uzdzdz=%13.6e\n", result18[p][15], result18[p][16], result18[p][17]);
	}

	printf("\nRequesting magnetic field laplacian at %d points...\n", N);
	getMagneticFieldLaplacian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: grad2ux=%13.6e, grad2uy=%13.6e, grad2uz=%13.6e\n",
			p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting vector potential at %d points...\n", N);
	getVectorPotential(authtoken, dataset, time, spatialInterp, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting vector potential gradient at %d points...\n", N);
	getVectorPotentialGradient(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result9);
	for (p = 0; p < N; p++) {
		printf("%d: duxdx=%13.6e, duxdy=%13.6e, duxdz=%13.6e, ", p, result9[p][0], result9[p][1], result9[p][2]);
		printf("duydx=%13.6e, duydy=%13.6e, duydz=%13.6e, ", result9[p][3], result9[p][4], result9[p][5]);
		printf("duzdx=%13.6e, duzdy=%13.6e, duzdz=%13.6e\n", result9[p][6], result9[p][7], result9[p][8]);
	}

	printf("\nRequesting vector potential hessian at %d points...\n", N);
	getVectorPotentialHessian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result18);
	for (p = 0; p < N; p++) {
		printf("%d: d2uxdxdx=%13.6e, d2uxdxdy=%13.6e, d2uxdxdz=%13.6e, ", p, result18[p][0], result18[p][1], result18[p][2]);
		printf("d2uxdydy=%13.6e, d2uxdydz=%13.6e, d2uxdzdz=%13.6e, ", result18[p][3], result18[p][4], result18[p][5]);
		printf("d2uydxdx=%13.6e, d2uydxdy=%13.6e, d2uydxdz=%13.6e, ", result18[p][6], result18[p][7], result18[p][8]);
		printf("d2uydydy=%13.6e, d2uydydz=%13.6e, d2uydzdz=%13.6e, ", result18[p][9], result18[p][10], result18[p][11]);
		printf("d2uzdxdx=%13.6e, d2uzdxdy=%13.6e, d2uzdxdz=%13.6e, ", result18[p][12], result18[p][13], result18[p][14]);
		printf("d2uzdydy=%13.6e, d2uzdydz=%13.6e, d2uzdzdz=%13.6e\n", result18[p][15], result18[p][16], result18[p][17]);
	}

	printf("\nRequesting vector potential laplacian at %d points...\n", N);
	getVectorPotentialLaplacian(authtoken, dataset, time, FD4Lag4, temporalInterp, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: grad2ux=%13.6e, grad2uy=%13.6e, grad2uz=%13.6e\n",
			p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting position at %d points, starting at time %f and ending at time %f...\n", N, startTime, endTime);
	getPosition(authtoken, dataset, startTime, endTime, lag_dt, spatialInterp, N, points, result3);

	printf("\nCoordinates of 10 points at startTime:\n");
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, points[p][0], points[p][1], points[p][2]);
	}
	printf("\nCoordinates of 10 points at endTime:\n");
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting box filter of velocity at %d points...\n", N);
	getBoxFilter(authtoken, dataset, field, time, filterwidth, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p, result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting sub-grid stress symmetric tensor at %d points...\n", N);
	getBoxFilterSGSsymtensor(authtoken, dataset, field, time, filterwidth, N, points, result6);
	for (p = 0; p < N; p++) {
		printf("%d: xx=%13.6e, yy=%13.6e, zz=%13.6e, xy=%13.6e, xz=%13.6e, yz=%13.6e\n", p,
			result6[p][0], result6[p][1], result6[p][2],
			result6[p][3], result6[p][4], result6[p][5]);
	}

	printf("\nRequesting sub-grid stress of two scalar fields at %d points...\n", N);
	getBoxFilterSGSscalar(authtoken, dataset, scalar_fields, time, filterwidth, N, points, result1);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e\n", p, result1[p]);
	}

	printf("\nRequesting sub-grid stress of a vector-scalar combination at %d points...\n", N);
	getBoxFilterSGSvector(authtoken, dataset, vector_scalar_fields, time, filterwidth, N, points, result3);
	for (p = 0; p < N; p++) {
		printf("%d: %13.6e, %13.6e, %13.6e\n", p,
			result3[p][0], result3[p][1], result3[p][2]);
	}

	printf("\nRequesting sub-grid stress of two vector fields at %d points...\n", N);
	getBoxFilterSGStensor(authtoken, dataset, vector_fields, time, filterwidth, N, points, result9);
	for (p = 0; p < N; p++) {
		printf("%d: xx=%13.6e, xy=%13.6e, xz=%13.6e, "
			"yx=%13.6e, yy=%13.6e, yz=%13.6e, zx=%13.6e, zy=%13.6e, zz=%13.6e\n", p,
			result9[p][0], result9[p][1], result9[p][2],
			result9[p][3], result9[p][4], result9[p][5],
			result9[p][6], result9[p][7], result9[p][8]);
	}

	printf("\nRequesting box filter of velocity gradient at %d points...\n", N);
	getBoxFilterGradient(authtoken, dataset, field, time, filterwidth, spacing, N, points, result9);
	for (p = 0; p < N; p++) {
		printf("%d: duxdx=%13.6e, duxdy=%13.6e, duxdz=%13.6e, ", p, result9[p][0], result9[p][1], result9[p][2]);
		printf("duydx=%13.6e, duydy=%13.6e, duydz=%13.6e, ", result9[p][3], result9[p][4], result9[p][5]);
		printf("duzdx=%13.6e, duzdy=%13.6e, duzdz=%13.6e\n", result9[p][6], result9[p][7], result9[p][8]);
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
