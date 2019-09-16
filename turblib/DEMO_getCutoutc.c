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
#define N 8

int main(int argc, char *argv[])
{

	char * authtoken = "edu.jhu.pha.turbulence.testing-201406";
	char * dataset;
	char * field;

	float *result;

	int time_step = 1, x_start, y_start, z_start, x_end, y_end, z_end;
	int x_step = 1, y_step = 1, z_step = 1, filter_width = 1;
	int size, i;

	/* Initialize gSOAP */
	soapinit();

	/* Enable exit on error.  See README for details. */
	turblibSetExitOnError(1);

	printf("\n.........getAnyCutoutWeb.........\n");


	printf("\n.........channel u.........\n");
	dataset = "channel";
	x_start = 1; y_start = 10; z_start = 1; x_end = 2; y_end = 11; z_end = 2;
	field = "u";
	if (field[0] == 'p')
	{
		size = (x_end-x_start+1) * (y_end-y_start+1) * (z_end-z_start+1);
	}
	else if (field[0] == 'u')
	{
		size = (x_end-x_start+1) * (y_end-y_start+1) * (z_end-z_start+1) * 3;
	}
	result = (float *)malloc(sizeof(float)*size);
	getCutout(authtoken, dataset, field, time_step, x_start, y_start, z_start, x_end, y_end, z_end, x_step, y_step, z_step, filter_width, result);
	for (i = 0; i < size; i++) {
		printf("%d: u=%f\n", i, result[i]);
	}
	free(result);


	printf("\n.........transition_bl p.........\n");
	dataset = "transition_bl";
	x_start = 1; y_start = 1; z_start = 1; x_end = 2; y_end = 2; z_end = 2;
	field = "p";
	if (field[0] == 'p')
	{
		size = (x_end-x_start+1) * (y_end-y_start+1) * (z_end-z_start+1);
	}
	else if (field[0] == 'u')
	{
		size = (x_end-x_start+1) * (y_end-y_start+1) * (z_end-z_start+1) * 3;
	}
	result = (float *)malloc(sizeof(float)*size);
	getCutout(authtoken, dataset, field, time_step, x_start, y_start, z_start, x_end, y_end, z_end, x_step, y_step, z_step, filter_width, result);
	for (i = 0; i < size; i++) {
		printf("%d: p=%f\n", i, result[i]);
	}


	printf("\n.........isotropic1024coarse u filter.........\n");
	dataset = "isotropic1024coarse";
	x_start = 1; y_start = 1; z_start = 1; x_end = 6; y_end = 6; z_end = 6;
	x_step = 4; y_step = 4; z_step = 4; filter_width = 4;
	field = "u";
	size = 2 * 2 * 2 * 3;
	result = (float *)malloc(sizeof(float)*size);
	getCutout(authtoken, dataset, field, time_step, x_start, y_start, z_start, x_end, y_end, z_end, x_step, y_step, z_step, filter_width, result);
	for (i = 0; i < size; i++) {
		printf("%d: u=%f\n", i, result[i]);
	}
	free(result);
	/* Free gSOAP resources */
	soapdestroy();

	return 0;
}
