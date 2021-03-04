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


/* $Id: turblib.c,v 1.13 2009-10-23 17:58:57 eric Exp $ */
#include <stdio.h>
#include <math.h>

#include "soapH.h"
#include "TurbulenceServiceSoap.nsmap"
#include "turblib.h"

/* global gSOAP runtime environment
 * Temporary until we can figure out how to pass a pointer
 * to and from a function correctly in Fortran.
 */
struct soap __jhuturbsoap;

char * __version_info = "20210108";

/* Error reporting - C */
char __turblib_err[TURB_ERROR_LENGTH];
int __turblib_errno = 0;
int __turblib_exit_on_error = 1;
int __turblib_prefetching = 0;

//Linked list of all added cutout files

set_info DataSets[14] = {
	{ 0, 0, 0 },//=0
	{ 0, 0, 0 },//=1
	{ 0, 0, 0 },//=2
	{ 2.0f * 3.14159265358979f / 1024.0f, .0025f, 1024 },  //mhd1024=3
	{ 2.0f * 3.14159265358979f / 1024.0f, .002f,  1024 }, //isotropic1024coarse=4
	{ 2.0f * 3.14159265358979f / 1024.0f, .0002f, 1024 }, //isotropic1024fine=5
	{ 0, 0, 0 },//channel=6
	{ 2.0f * 3.14159265358979f / 1024.0f, .04f, 1024 }, //mixing_dataset=7
	{ 0, 0, 0 },//rmhd=8
	{ 0, 0, 0 },//=9
	{ 2.0f * 3.14159265358979f / 4096.0f, .0002f, 4096 }, //isotropic4096=10
	{ 2.0f * 3.14159265358979f / 4096.0f, 1.0f, 4096 }, //strat4096=11
	{ 0, 0, 0 },//transition_bl=12
	{ 0, 0, 0 }//channel5200=13
};

turb_fn TurbFields[6] =
{
	{ 'u', 3}, //velocity
	{ 'p', 1}, //pressure
	{ 'b', 3}, //magnetic
	{ 'a', 3}, //vector potential
	{ 'd', 1},  //density
	{ 't', 1}  //temperature
};

char * turblibGetErrorString() {
	return __turblib_err;
}

int turblibGetErrorNumber() {
	return __turblib_errno;
}

void turblibPrintError() {
	fprintf(stderr, "%d: %s\n", turblibGetErrorNumber(), turblibGetErrorString());
}

void turblibSetExitOnError(int v) {
	__turblib_exit_on_error = v;
}


/* Error reporting - Fortran */
void turblibgeterrorstring_(char *dest, int len) {
	strncpy(dest, __turblib_err, len);
}

int turblibgeterrornumber_() {
	return turblibGetErrorNumber();
}

void turblibprinterror_() {
	turblibPrintError();
}

void turblibsetexitonerror_(int *v) {
	turblibSetExitOnError(*v);
}

/* Determine appropriate error behavior */
void turblibHandleError() {
    turblibPrintError();
	if (__turblib_exit_on_error) {
		exit(1);
	}
}


/* Return the enum relating to the Fortran constant */
enum turb1__SpatialInterpolation SpatialIntToEnum(enum SpatialInterpolation spatial)
{
	switch (spatial)
	{
	case 0:
		return turb1__SpatialInterpolation__None;
	case 4:
		return turb1__SpatialInterpolation__Lag4;
	case 6:
		return turb1__SpatialInterpolation__Lag6;
	case 8:
		return turb1__SpatialInterpolation__Lag8;
	case 40:
		return turb1__SpatialInterpolation__None_USCOREFd4;
	case 44:
		return turb1__SpatialInterpolation__Fd4Lag4;
	case 60:
		return turb1__SpatialInterpolation__None_USCOREFd6;
	case 80:
		return turb1__SpatialInterpolation__None_USCOREFd8;
	case 104:
		return turb1__SpatialInterpolation__M1Q4;
	case 106:
		return turb1__SpatialInterpolation__M1Q6;
	case 108:
		return turb1__SpatialInterpolation__M1Q8;
	case 110:
		return turb1__SpatialInterpolation__M1Q10;
	case 112:
		return turb1__SpatialInterpolation__M1Q12;
	case 114:
		return turb1__SpatialInterpolation__M1Q14;
	case 204:
		return turb1__SpatialInterpolation__M2Q4;
	case 206:
		return turb1__SpatialInterpolation__M2Q6;
	case 208:
		return turb1__SpatialInterpolation__M2Q8;
	case 210:
		return turb1__SpatialInterpolation__M2Q10;
	case 212:
		return turb1__SpatialInterpolation__M2Q12;
	case 214:
		return turb1__SpatialInterpolation__M2Q14;
	case 304:
		return turb1__SpatialInterpolation__M3Q4;
	case 306:
		return turb1__SpatialInterpolation__M3Q6;
	case 308:
		return turb1__SpatialInterpolation__M3Q8;
	case 310:
		return turb1__SpatialInterpolation__M3Q10;
	case 312:
		return turb1__SpatialInterpolation__M3Q12;
	case 314:
		return turb1__SpatialInterpolation__M3Q14;
	case 404:
		return turb1__SpatialInterpolation__M4Q4;
	case 406:
		return turb1__SpatialInterpolation__M4Q6;
	case 408:
		return turb1__SpatialInterpolation__M4Q8;
	case 410:
		return turb1__SpatialInterpolation__M4Q10;
	case 412:
		return turb1__SpatialInterpolation__M4Q12;
	case 414:
		return turb1__SpatialInterpolation__M4Q14;
	default:
		return -1;
	}
	return -1;
}

/* Return the enum relating to the Fortran constant */
enum turb1__TemporalInterpolation TemporalIntToEnum(enum TemporalInterpolation temporal)
{
	switch (temporal)
	{
	case 0:
		return turb1__TemporalInterpolation__None;
	case 1:
		return turb1__TemporalInterpolation__PCHIP;
	}
	return -1;
}

/* Get turblib version */
char * getVersion_(void) {
	return getVersion();
}

char * getVersion(void) {
	return __version_info;
}

/* Intialize the gSOAP runtime environment */
void soapinit_() {
	soapinit();
}

void soapinit() {
	soap_init(&__jhuturbsoap);
}

/* Destroy the gSOAP environment */
void soapdestroy_() {
	soap_destroy(&__jhuturbsoap);
}

void soapdestroy() {
	soap_destroy(&__jhuturbsoap);
}

int getVelocity(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getVelocitySoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVelocityAndPressure(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][4])
{
		return getVelocityAndPressureSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getPressure(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[])
{
		return getPressureSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getPressureHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][6])
{
		return getPressureHessianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVelocityGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9])
{
		return getVelocityGradientSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVelocityHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18])
{
		return getVelocityHessianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVelocityLaplacian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getVelocityLaplacianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getPressureGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getPressureGradientSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getMagneticFieldGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9])
{
		return getMagneticFieldGradientSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVectorPotentialGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9])
{
		return getVectorPotentialGradientSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getMagneticFieldHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18])
{
		return getMagneticFieldHessianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getMagneticFieldLaplacian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getMagneticFieldLaplacianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getMagneticField(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getMagneticFieldSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVectorPotentialHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18])
{
		return getVectorPotentialHessianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVectorPotentialLaplacian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getVectorPotentialLaplacianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVectorPotential(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getVectorPotentialSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getDensity(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[])
{
		return getDensitySoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getDensityGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
		return getDensityGradientSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getDensityHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][6])
{
		return getDensityHessianSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getInvariant(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][2])
{
		return getInvariantSoap(authToken, dataset, time, spatial, temporal, count, datain, dataout);
}

int getVelocitySoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetVelocity input;
	struct _turb1__GetVelocityResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVelocity(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVelocityResult->Vector3,
			output.GetVelocityResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getthreshold_(char *authToken,
	char *dataset, char *field, float *time, float *threshold,
	int *spatial,
	int *x_start, int *y_start, int *z_start, int *x_end, int *y_end, int *z_end,
	ThresholdInfo** dataout, int *result_size)
{
	return getThreshold(authToken, dataset, field, *time, *threshold, *spatial,
		*x_start, *y_start, *z_start, *x_end, *y_end, *z_end, dataout, result_size);
}

void deallocate_array_(ThresholdInfo **threshold_array)
{
	free(*threshold_array);
}

int getThreshold(char *authToken,
	char *dataset, char *field, float time, float threshold,
	enum SpatialInterpolation spatial,
	int x_start, int y_start, int z_start, int x_end, int y_end, int z_end,
	ThresholdInfo **dataout, int *result_size)
{
	int rc;

	struct _turb1__GetThreshold input;
	struct _turb1__GetThresholdResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.threshold = threshold;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.x_USCOREstart = x_start;
	input.y_USCOREstart = y_start;
	input.z_USCOREstart = z_start;
	input.x_USCOREend = x_end;
	input.y_USCOREend = y_end;
	input.z_USCOREend = z_end;
	input.addr = NULL;

	rc = soap_call___turb1__GetThreshold(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		*result_size = output.GetThresholdResult->__sizeThresholdInfo;
		*dataout = (ThresholdInfo *)malloc(sizeof(ThresholdInfo) * (*result_size));
		memcpy(*dataout, output.GetThresholdResult->ThresholdInfo,
			output.GetThresholdResult->__sizeThresholdInfo * sizeof(ThresholdInfo));
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); //  detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}


int getboxfilter_(char *authToken,
	char *dataset, char *field, float *time, float *filterwidth,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getBoxFilter(authToken,
		dataset, field, *time, *filterwidth,
		*count, datain, dataout);
}

int getBoxFilter(char *authToken,
	char *dataset, char *field, float time, float filterwidth,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetBoxFilter input;
	struct _turb1__GetBoxFilterResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilter(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterResult->Vector3,
			output.GetBoxFilterResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getboxfiltersgs_(char *authToken,
	char *dataset, char *field, float *time, float *filterwidth,
	int *count, float datain[][3], float dataout[][6],
	int len_a, int len_d)
{
	return getBoxFilterSGS(authToken,
		dataset, field, *time, *filterwidth,
		*count, datain, dataout);
}

int getBoxFilterSGS(char *authToken,
	char *dataset, char *field, float time, float filterwidth,
	int count, float datain[][3], float dataout[][6])
{
	int rc;

	struct _turb1__GetBoxFilterSGS input;
	struct _turb1__GetBoxFilterSGSResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilterSGS(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterSGSResult->SGSTensor,
			output.GetBoxFilterSGSResult->__sizeSGSTensor * sizeof(float) * 6);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getboxfiltersgsscalar_(char *authToken,
	char *dataset, char *field, float *time, float *filterwidth,
	int *count, float datain[][3], float dataout[],
	int len_a, int len_d)
{
	return getBoxFilterSGSscalar(authToken,
		dataset, field, *time, *filterwidth,
		*count, datain, dataout);
}

int getBoxFilterSGSscalar(char *authToken,
	char *dataset, char *field, float time, float filterwidth,
	int count, float datain[][3], float dataout[])
{
	int rc;

	struct _turb1__GetBoxFilterSGSscalar input;
	struct _turb1__GetBoxFilterSGSscalarResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilterSGSscalar(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterSGSscalarResult->float_,
			output.GetBoxFilterSGSscalarResult->__sizefloat_ * sizeof(float));
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getboxfiltersgsvector_(char *authToken,
	char *dataset, char *field, float *time, float *filterwidth,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getBoxFilterSGSvector(authToken,
		dataset, field, *time, *filterwidth,
		*count, datain, dataout);
}

int getBoxFilterSGSvector(char *authToken,
	char *dataset, char *field, float time, float filterwidth,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetBoxFilterSGSvector input;
	struct _turb1__GetBoxFilterSGSvectorResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilterSGSvector(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterSGSvectorResult->Vector3,
			output.GetBoxFilterSGSvectorResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getboxfiltersgssymtensor_(char *authToken,
	char *dataset, char *field, float *time, float *filterwidth,
	int *count, float datain[][3], float dataout[][6],
	int len_a, int len_d)
{
	return getBoxFilterSGSsymtensor(authToken,
		dataset, field, *time, *filterwidth,
		*count, datain, dataout);
}

int getBoxFilterSGSsymtensor(char *authToken,
	char *dataset, char *field, float time, float filterwidth,
	int count, float datain[][3], float dataout[][6])
{
	int rc;

	struct _turb1__GetBoxFilterSGSsymtensor input;
	struct _turb1__GetBoxFilterSGSsymtensorResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilterSGSsymtensor(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterSGSsymtensorResult->SGSTensor,
			output.GetBoxFilterSGSsymtensorResult->__sizeSGSTensor * sizeof(float) * 6);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getboxfiltersgstensor_(char *authToken,
	char *dataset, char *field, float *time, float *filterwidth,
	int *count, float datain[][3], float dataout[][9],
	int len_a, int len_d)
{
	return getBoxFilterSGStensor(authToken,
		dataset, field, *time, *filterwidth,
		*count, datain, dataout);
}

int getBoxFilterSGStensor(char *authToken,
	char *dataset, char *field, float time, float filterwidth,
	int count, float datain[][3], float dataout[][9])
{
	int rc;

	struct _turb1__GetBoxFilterSGStensor input;
	struct _turb1__GetBoxFilterSGStensorResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilterSGStensor(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterSGStensorResult->VelocityGradient,
			output.GetBoxFilterSGStensorResult->__sizeVelocityGradient * sizeof(float) * 9);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getboxfiltergradient_(char *authToken,
	char *dataset, char* field, float *time,
	float *filterwidth, float *spacing,
	int *count, float datain[][3], float dataout[][9],
	int len_a, int len_d)
{
	return getBoxFilterGradient(authToken,
		dataset, field, *time,
		*filterwidth, *spacing,
		*count, datain, dataout);
}

int getBoxFilterGradient(char *authToken,
	char *dataset, char *field, float time,
	float filterwidth, float spacing,
	int count, float datain[][3], float dataout[][9])
{
	int rc;

	struct _turb1__GetBoxFilterGradient input;
	struct _turb1__GetBoxFilterGradientResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.time = time;
	input.filterwidth = filterwidth;
	input.spacing = spacing;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetBoxFilterGradient(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetBoxFilterGradientResult->VelocityGradient,
			output.GetBoxFilterGradientResult->__sizeVelocityGradient * sizeof(float) * 9);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVelocityAndPressureSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][4])
{
	int rc;

	struct _turb1__GetVelocityAndPressure input;
	struct _turb1__GetVelocityAndPressureResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVelocityAndPressure(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVelocityAndPressureResult->Vector3P,
			output.GetVelocityAndPressureResult->__sizeVector3P * sizeof(float) * 4);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getPressureHessianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][6])
{
	int rc;

	struct _turb1__GetPressureHessian input;
	struct _turb1__GetPressureHessianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetPressureHessian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetPressureHessianResult->PressureHessian,
			output.GetPressureHessianResult->__sizePressureHessian * sizeof(float) * 6);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVelocityGradientSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9])
{
	int rc;

	struct _turb1__GetVelocityGradient input;
	struct _turb1__GetVelocityGradientResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVelocityGradient(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVelocityGradientResult->VelocityGradient,
			output.GetVelocityGradientResult->__sizeVelocityGradient * sizeof(float) * 9);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getMagneticFieldGradientSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9])
{
	int rc;

	struct _turb1__GetMagneticFieldGradient input;
	struct _turb1__GetMagneticFieldGradientResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetMagneticFieldGradient(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetMagneticFieldGradientResult->VelocityGradient,
			output.GetMagneticFieldGradientResult->__sizeVelocityGradient * sizeof(float) * 9);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVectorPotentialGradientSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9])
{
	int rc;

	struct _turb1__GetVectorPotentialGradient input;
	struct _turb1__GetVectorPotentialGradientResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVectorPotentialGradient(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVectorPotentialGradientResult->VelocityGradient,
			output.GetVectorPotentialGradientResult->__sizeVelocityGradient * sizeof(float) * 9);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getPressureGradientSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetPressureGradient input;
	struct _turb1__GetPressureGradientResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetPressureGradient(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetPressureGradientResult->Vector3,
			output.GetPressureGradientResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVelocityHessianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18])
{
	int rc;

	struct _turb1__GetVelocityHessian input;
	struct _turb1__GetVelocityHessianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVelocityHessian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVelocityHessianResult->VelocityHessian,
			output.GetVelocityHessianResult->__sizeVelocityHessian * sizeof(float) * 18);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVelocityLaplacianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetVelocityLaplacian input;
	struct _turb1__GetVelocityLaplacianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVelocityLaplacian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVelocityLaplacianResult->Vector3,
			output.GetVelocityLaplacianResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getMagneticFieldHessianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18])
{
	int rc;

	struct _turb1__GetMagneticHessian input;
	struct _turb1__GetMagneticHessianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetMagneticHessian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetMagneticHessianResult->VelocityHessian,
			output.GetMagneticHessianResult->__sizeVelocityHessian * sizeof(float) * 18);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getMagneticFieldLaplacianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetMagneticFieldLaplacian input;
	struct _turb1__GetMagneticFieldLaplacianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetMagneticFieldLaplacian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetMagneticFieldLaplacianResult->Vector3,
			output.GetMagneticFieldLaplacianResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVectorPotentialHessianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18])
{
	int rc;

	struct _turb1__GetVectorPotentialHessian input;
	struct _turb1__GetVectorPotentialHessianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVectorPotentialHessian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVectorPotentialHessianResult->VelocityHessian,
			output.GetVectorPotentialHessianResult->__sizeVelocityHessian * sizeof(float) * 18);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getVectorPotentialLaplacianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetVectorPotentialLaplacian input;
	struct _turb1__GetVectorPotentialLaplacianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVectorPotentialLaplacian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVectorPotentialLaplacianResult->Vector3,
			output.GetVectorPotentialLaplacianResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int nullop_(char *authToken, int *count,
	float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return nullOp(authToken, *count,
		datain, dataout);
}

int nullOp(char *authToken, int count,
	float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__NullOp input;
	struct _turb1__NullOpResponse output;

	input.authToken = authToken;

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;

	rc = soap_call___turb1__NullOp(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.NullOpResult->Vector3,
			output.NullOpResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getForce(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetForce input;
	struct _turb1__GetForceResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetForce(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetForceResult->Vector3,
			output.GetForceResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getforce_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getForce(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getPosition(char *authToken,
	char *dataset, float startTime, float endTime,
	float dt,
	enum SpatialInterpolation spatial,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetPosition input;
	struct _turb1__GetPositionResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.StartTime = startTime;
	input.EndTime = endTime;
	input.dt = dt;
	input.spatialInterpolation = SpatialIntToEnum(spatial);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetPosition(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetPositionResult->Point3,
			output.GetPositionResult->__sizePoint3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	return rc;
}

int getposition_(char *authToken,
	char *dataset, float *startTime, float *endTime,
	float *dt,
	int *spatial,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getPosition(authToken,
		dataset, *startTime, *endTime,
		*dt,
		*spatial,
		*count, datain, dataout);
}


int getrawvelocity_(char *authToken, char *dataset, int *time_step,
	int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
	float dataout[])
{
	return getRawVelocity(authToken, dataset, *time_step, *X, *Y, *Z,
		*Xwidth, *Ywidth, *Zwidth, (char*)dataout);
}

int getRawVelocity(char *authToken,
	char *dataset, int time_step,
	int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[])
{
	fprintf(stderr, "%s\n", "******getRawVelocity has been deprecated. Please use getCutout instead******");
	return -999;
	int rc;

	struct _turb1__GetRawVelocity input;
	struct _turb1__GetRawVelocityResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.T = time_step;
	input.X = X;
	input.Y = Y;
	input.Z = Z;
	input.Xwidth = Xwidth;
	input.Ywidth = Ywidth;
	input.Zwidth = Zwidth;
	input.addr = NULL;

	rc = soap_call___turb1__GetRawVelocity(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetRawVelocityResult.__ptr,
			output.GetRawVelocityResult.__size);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); //  detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getMagneticFieldSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetMagneticField input;
	struct _turb1__GetMagneticFieldResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetMagneticField(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetMagneticFieldResult->Vector3,
			output.GetMagneticFieldResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  //remove deserialized data and clean up
	soap_done(&__jhuturbsoap); //detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getrawmagneticfield_(char *authToken,
	char *dataset, int *time_step,
	int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
	float dataout[])
{
	return getRawMagneticField(authToken, dataset, *time_step, *X, *Y, *Z,
		*Xwidth, *Ywidth, *Zwidth, (char*)dataout);
}

int getRawMagneticField(char *authToken,
	char *dataset, int time_step,
	int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[])
{
	fprintf(stderr, "%s\n", "******getRawMagneticField has been deprecated. Please use getCutout instead******");
	return -999;
	int rc;

	struct _turb1__GetRawMagneticField input;
	struct _turb1__GetRawMagneticFieldResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.T = time_step;
	input.X = X;
	input.Y = Y;
	input.Z = Z;
	input.Xwidth = Xwidth;
	input.Ywidth = Ywidth;
	input.Zwidth = Zwidth;
	input.addr = NULL;

	rc = soap_call___turb1__GetRawMagneticField(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetRawMagneticFieldResult.__ptr,
			output.GetRawMagneticFieldResult.__size);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getVectorPotentialSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetVectorPotential input;
	struct _turb1__GetVectorPotentialResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetVectorPotential(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetVectorPotentialResult->Vector3,
			output.GetVectorPotentialResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); //  detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getrawvectorpotential_(char *authToken,
	char *dataset, int *time_step,
	int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
	float dataout[])
{
	return getRawVectorPotential(authToken, dataset, *time_step, *X, *Y, *Z,
		*Xwidth, *Ywidth, *Zwidth, (char*)dataout);
}

int getRawVectorPotential(char *authToken,
	char *dataset, int time_step,
	int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[])
{
	fprintf(stderr, "%s\n", "******getRawVectorPotential has been deprecated. Please use getCutout instead******");
	return -999;
	int rc;

	struct _turb1__GetRawVectorPotential input;
	struct _turb1__GetRawVectorPotentialResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.T = time_step;
	input.X = X;
	input.Y = Y;
	input.Z = Z;
	input.Xwidth = Xwidth;
	input.Ywidth = Ywidth;
	input.Zwidth = Zwidth;
	input.addr = NULL;

	rc = soap_call___turb1__GetRawVectorPotential(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetRawVectorPotentialResult.__ptr,
			output.GetRawVectorPotentialResult.__size);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}


int getPressureSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[])
{
	int rc;

	struct _turb1__GetPressure input;
	struct _turb1__GetPressureResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetPressure(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetPressureResult->Pressure,
			output.GetPressureResult->__sizePressure * sizeof(float));
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getrawpressure_(char *authToken, char *dataset, int *time_step,
	int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
	float dataout[])
{
	return getRawPressure(authToken, dataset, *time_step, *X, *Y, *Z,
		*Xwidth, *Ywidth, *Zwidth, (char*)dataout);
}

int getRawPressure(char *authToken,
	char *dataset, int time_step,
	int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[])
{
	fprintf(stderr, "%s\n", "******getRawPressure has been deprecated. Please use getCutout instead******");
	return -999;
	int rc;

	struct _turb1__GetRawPressure input;
	struct _turb1__GetRawPressureResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.T = time_step;
	input.X = X;
	input.Y = Y;
	input.Z = Z;
	input.Xwidth = Xwidth;
	input.Ywidth = Ywidth;
	input.Zwidth = Zwidth;
	input.addr = NULL;

	rc = soap_call___turb1__GetRawPressure(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetRawPressureResult.__ptr,
			output.GetRawPressureResult.__size);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getDensitySoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[])
{
	int rc;

	struct _turb1__GetDensity input;
	struct _turb1__GetDensityResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetDensity(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetDensityResult->Pressure,
			output.GetDensityResult->__sizePressure * sizeof(float));
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

int getDensityGradientSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3])
{
	int rc;

	struct _turb1__GetDensityGradient input;
	struct _turb1__GetDensityGradientResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetDensityGradient(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetDensityGradientResult->Vector3,
			output.GetDensityGradientResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getDensityHessianSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][6])
{
	int rc;

	struct _turb1__GetDensityHessian input;
	struct _turb1__GetDensityHessianResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetDensityHessian(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetDensityHessianResult->PressureHessian,
			output.GetDensityHessianResult->__sizePressureHessian * sizeof(float) * 6);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;
	return rc;
}

int getInvariantSoap(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][2])
{
	int rc, i;
	float* full_InvariantOutput;
	full_InvariantOutput = malloc(sizeof(float)*count*3);
	//dataout = (float*) malloc(sizeof(float)*count);
	
	struct _turb1__GetInvariant input;
	struct _turb1__GetInvariantResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.time = time;
	input.spatialInterpolation = SpatialIntToEnum(spatial);
	input.temporalInterpolation = TemporalIntToEnum(temporal);

	struct turb1__ArrayOfPoint3 pointArray;
	pointArray.__sizePoint3 = count;
	pointArray.Point3 = (void *)datain;
	input.points = &pointArray;
	input.addr = NULL;

	rc = soap_call___turb1__GetInvariant(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(full_InvariantOutput, output.GetInvariantResult->Vector3,
			output.GetInvariantResult->__sizeVector3 * sizeof(float) * 3);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  /* remove deserialized data and clean up */
	soap_done(&__jhuturbsoap); /*  detach the gSOAP environment  */

	__turblib_errno = rc;

	for (i=0; i<count; i++ )
	{
		dataout[i][0] = *(full_InvariantOutput+i*3+0);
		dataout[i][1] = *(full_InvariantOutput+i*3+1);
	}

	return rc;
}

int getrawdensity_(char *authToken, char *dataset, int *time_step,
	int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
	float dataout[])
{
	return getRawDensity(authToken, dataset, *time_step, *X, *Y, *Z,
		*Xwidth, *Ywidth, *Zwidth, (char*)dataout);
}

int getRawDensity(char *authToken,
	char *dataset, int time_step,
	int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[])
{
	fprintf(stderr, "%s\n", "******getRawDensity has been deprecated. Please use getCutout instead******");
	return -999;
	int rc;

	struct _turb1__GetRawDensity input;
	struct _turb1__GetRawDensityResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.T = time_step;
	input.X = X;
	input.Y = Y;
	input.Z = Z;
	input.Xwidth = Xwidth;
	input.Ywidth = Ywidth;
	input.Zwidth = Zwidth;
	input.addr = NULL;

	rc = soap_call___turb1__GetRawDensity(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetRawDensityResult.__ptr,
			output.GetRawDensityResult.__size);
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}
/*
 int getAnyCutoutWeb (char *authToken,
 char *dataset, char *field, int T,
 int x0, int y0, int z0,
 int nx, int ny, int nz,
 int x_step, int y_step, int z_step, int filter_width,
 float dataout[])
 {
 return getAnyCutoutWebSoap (authToken, dataset, field, T, x0, y0, z0, nx, ny, nz, x_step, y_step, z_step, filter_width, dataout);
 }*/
int getCutout(char *authToken,
	char *dataset, char *field, int time_step, int x_start, int y_start, int z_start,
	int x_end, int y_end, int z_end,
	int x_step, int y_step, int z_step, int filter_width,
	float dataout[])
{
	int rc;

	struct _turb1__GetAnyCutoutWeb input;
	struct _turb1__GetAnyCutoutWebResponse output;

	input.authToken = authToken;
	input.dataset = dataset;
	input.field = field;
	input.T = time_step;
	input.x_USCOREstart = x_start;
	input.y_USCOREstart = y_start;
	input.z_USCOREstart = z_start;
	input.x_USCOREend = x_end;
	input.y_USCOREend = y_end;
	input.z_USCOREend = z_end;
	input.x_USCOREstep = x_step;
	input.y_USCOREstep = y_step;
	input.z_USCOREstep = z_step;
	input.filter_USCOREwidth = filter_width;
	input.addr = NULL;

	rc = soap_call___turb1__GetAnyCutoutWeb(&__jhuturbsoap, NULL, NULL, &input, &output);
	if (rc == SOAP_OK) {
		memcpy(dataout, output.GetAnyCutoutWebResult.__ptr,
			output.GetAnyCutoutWebResult.__size);
		//memcpy(dataout, output.GetAnyCutoutWebResult->float_,
		//	output.GetAnyCutoutWebResult->__sizefloat_ * sizeof(float));
		bzero(__turblib_err, TURB_ERROR_LENGTH);
	}
	else {
		soap_sprint_fault(&__jhuturbsoap, __turblib_err, TURB_ERROR_LENGTH);
		turblibHandleError();
	}

	soap_end(&__jhuturbsoap);  // remove deserialized data and clean up
	soap_done(&__jhuturbsoap); // detach the gSOAP environment

	__turblib_errno = rc;

	return rc;
}

////////////////////////////////

/* Interpolation Functions */

int lagrangianInterp(int comps, float *kernel, float position[3], int nOrder, float dx, float *result)
{
	int node[3];
	node[0] = (int)(floor(position[0] / dx));
	node[1] = (int)(floor(position[1] / dx));
	node[2] = (int)(floor(position[2] / dx));

	int x, y, z;
	float lagInt[3][8];

	for (x = 0; x < 3; x++)
	{
		float z1 = position[x] / dx - (float)node[x];
		float z2 = z1 * z1;
		float z3 = z2 * z1;
		switch (nOrder) {
		case 4:
		{
			lagInt[x][0] = (-2 * z1 + 3 * z2 - z3) / 6;
			lagInt[x][1] = (2 - z1 - 2 * z2 + z3) / 2;
			lagInt[x][2] = (2 * z1 + z2 - z3) / 2;
			lagInt[x][3] = (-z1 + z3) / 6;
			break;
		}
		case 6:
		{
			float z4 = z2 * z2;
			float z5 = z3 * z2;
			lagInt[x][0] = (6 * z1 - 5 * z2 - 5 * z3 + 5 * z4 - z5) / 120;
			lagInt[x][1] = (-12 * z1 + 16 * z2 - z3 - 4 * z4 + z5) / 24;
			lagInt[x][2] = (12 - 4 * z1 - 15 * z2 + 5 * z3 + 3 * z4 - z5) / 12;
			lagInt[x][3] = (12 * z1 + 8 * z2 - 7 * z3 - 2 * z4 + z5) / 12;
			lagInt[x][4] = (-6 * z1 - z2 + 7 * z3 + z4 - z5) / 24;
			lagInt[x][5] = (4 * z1 - 5 * z3 + z5) / 120;
			break;
		}
		case 8:
		{
			float z4 = z3 * z1;
			float z5 = z4 * z1;
			float z6 = z5 * z1;
			float z7 = z6 * z1;
			lagInt[x][0] = -z1 * (z6 - 7 * z5 + 7 * z4 + 35 * z3 - 56 * z2 - 28 * z1 + 48) / 5040;
			lagInt[x][1] = z1 * (z6 - 6 * z5 - 2 * z4 + 60 * z3 - 71 * z2 - 54 * z1 + 72) / 720;
			lagInt[x][2] = -z1 * (z6 - 5 * z5 - 9 * z4 + 65 * z3 - 16 * z2 - 180 * z1 + 144) / 240;
			lagInt[x][3] = (z7 - 4 * z6 - 14 * z5 + 56 * z4 + 49 * z3 - 196 * z2 - 36 * z1 + 144) / 144;
			lagInt[x][4] = -z1 * (z6 - 3 * z5 - 17 * z4 + 39 * z3 + 88 * z2 - 108 * z1 - 144) / 144;
			lagInt[x][5] = z1 * (z6 - 2 * z5 - 18 * z4 + 20 * z3 + 89 * z2 - 18 * z1 - 72) / 240;
			lagInt[x][6] = -z1 * (z6 - z5 - 17 * z4 + 5 * z3 + 64 * z2 - 4 * z1 - 48) / 720;
			lagInt[x][7] = z1 * (z6 - 14 * z4 + 49 * z2 - 36) / 5040;
			break;
		}
		}
	}

	int comp;
	for (comp = 0; comp < comps; comp++)
		result[comp] = 0;
	int index = 0;
	for (z = 0; z < nOrder; z++)
	{
		for (y = 0; y < nOrder; y++)
		{
			for (x = 0; x < nOrder; x++)
			{
				for (comp = 0; comp < comps; comp++)
				{
					result[comp] += kernel[index++] * lagInt[0][x] * lagInt[1][y] * lagInt[2][z];
				}
			}
		}
	}
	return 0;
}

int lagrangianInterp2(int comps, dataKernel* kernel, float position[3], int nOrder, float dx, float *result)
{
	int node[3];
	node[0] = (int)(floor(position[0] / dx));
	node[1] = (int)(floor(position[1] / dx));
	node[2] = (int)(floor(position[2] / dx));
	int x, y, z;
	float lagInt[3][8];
	float* data = kernel->data;
	for (x = 0; x < 3; x++)
	{
		float z1 = position[x] / dx - (float)node[x];
		float z2 = z1 * z1;
		float z3 = z2 * z1;
		switch (nOrder) {
		case 4:
		{
			lagInt[x][0] = (-2 * z1 + 3 * z2 - z3) / 6;
			lagInt[x][1] = (2 - z1 - 2 * z2 + z3) / 2;
			lagInt[x][2] = (2 * z1 + z2 - z3) / 2;
			lagInt[x][3] = (-z1 + z3) / 6;
			break;
		}
		case 6:
		{
			float z4 = z2 * z2;
			float z5 = z3 * z2;
			lagInt[x][0] = (6 * z1 - 5 * z2 - 5 * z3 + 5 * z4 - z5) / 120;
			lagInt[x][1] = (-12 * z1 + 16 * z2 - z3 - 4 * z4 + z5) / 24;
			lagInt[x][2] = (12 - 4 * z1 - 15 * z2 + 5 * z3 + 3 * z4 - z5) / 12;
			lagInt[x][3] = (12 * z1 + 8 * z2 - 7 * z3 - 2 * z4 + z5) / 12;
			lagInt[x][4] = (-6 * z1 - z2 + 7 * z3 + z4 - z5) / 24;
			lagInt[x][5] = (4 * z1 - 5 * z3 + z5) / 120;
			break;
		}
		case 8:
		{
			float z4 = z3 * z1;
			float z5 = z4 * z1;
			float z6 = z5 * z1;
			float z7 = z6 * z1;
			lagInt[x][0] = -z1 * (z6 - 7 * z5 + 7 * z4 + 35 * z3 - 56 * z2 - 28 * z1 + 48) / 5040;
			lagInt[x][1] = z1 * (z6 - 6 * z5 - 2 * z4 + 60 * z3 - 71 * z2 - 54 * z1 + 72) / 720;
			lagInt[x][2] = -z1 * (z6 - 5 * z5 - 9 * z4 + 65 * z3 - 16 * z2 - 180 * z1 + 144) / 240;
			lagInt[x][3] = (z7 - 4 * z6 - 14 * z5 + 56 * z4 + 49 * z3 - 196 * z2 - 36 * z1 + 144) / 144;
			lagInt[x][4] = -z1 * (z6 - 3 * z5 - 17 * z4 + 39 * z3 + 88 * z2 - 108 * z1 - 144) / 144;
			lagInt[x][5] = z1 * (z6 - 2 * z5 - 18 * z4 + 20 * z3 + 89 * z2 - 18 * z1 - 72) / 240;
			lagInt[x][6] = -z1 * (z6 - z5 - 17 * z4 + 5 * z3 + 64 * z2 - 4 * z1 - 48) / 720;
			lagInt[x][7] = z1 * (z6 - 14 * z4 + 49 * z2 - 36) / 5040;
			break;
		}
		}
	}

	int comp;
	for (comp = 0; comp < comps; comp++)
		result[comp] = 0;
	int index = 0;
	for (z = 0; z < nOrder; z++)
	{
		for (y = 0; y < nOrder; y++)
		{
			index = (kernel->x)*comps + (y + kernel->y)*kernel->hx*comps +
				(z + kernel->z)*kernel->hx*kernel->hy*comps;
			for (x = 0; x < nOrder; x++)
			{
				for (comp = 0; comp < comps; comp++)
				{
					result[comp] += data[index++] * lagInt[0][x] * lagInt[1][y] * lagInt[2][z];
				}
			}
		}
	}
	return 0;
}

int pchipInterp(int comps, float *data[4], float time, int timestep, float dt, float *result)
{
	float times[4] = { (timestep - 1) * dt, (timestep)* dt, (timestep + 1) * dt, (timestep + 2) * dt };
	int j;
	for (j = 0; j < comps; j++)
	{
		float a, b, c, d;
		float delta = times[2] - times[1];
		float drv1 = ((data[2][j] - data[1][j]) / (times[2] - times[1]) + (data[1][j] - data[0][j]) / (times[1] - times[0]
			)) / 2;
		float drv2 = ((data[3][j] - data[2][j]) / (times[3] - times[2]) + (data[2][j] - data[1][j]) / (times[2] - times[1]
			)) / 2;

		a = data[1][j];
		b = drv1;
		c = ((data[2][j] - data[1][j]) / delta - drv1) / delta;
		d = 2 / delta / delta * ((drv1 + drv2) / 2 - (data[2][j] - data[1][j]) / (times[2] - times[1]));
		result[j] = a + b * (time - times[1]) + c * (time - times[1]) * (time - times[1]) +
			d * (time - times[1]) * (time - times[1]) * (time - times[2]);
	}
	return 0;
}

/* Differentiation Functions */

int computeGradient(dataKernel* kernel, int comps, float dx, int size, int nOrder, float *output)
{
	int x, y, z, w;
	//Differentiate each point
	int   hx = comps, hy = comps * kernel->hx, hz = comps * kernel->hx*kernel->hy;
	float* data = kernel->data;

	for (z = 0; z < size; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				//Each component
				for (w = 0; w < comps; w++)
				{
					int x_ = nOrder / 2 + x + kernel->x, y_ = nOrder / 2 + y + kernel->y, z_ = nOrder / 2 + z + kernel->z;
					int index = z * size*size*comps * 3 + y * size*comps * 3 + x * comps * 3 + w * 3;
					int center = z_ * hz + y_ * hy + x_ * hx + w;
					switch (nOrder)
					{
					case 4: {
						// dfw/dx
						output[index] =
							2.0f / 3.0f / dx * (data[center + hx] - data[center - hx]) -
							1.0f / 12.0f / dx * (data[center + 2 * hx] - data[center - 2 * hx]);

						// dfw/dy
						output[index + 1] =
							2.0f / 3.0f / dx * (data[center + hy] - data[center - hy]) -
							1.0f / 12.0f / dx * (data[center + 2 * hy] - data[center - 2 * hy]);

						// dfw/dz
						output[index + 2] =
							2.0f / 3.0f / dx * (data[center + hz] - data[center - hz]) -
							1.0f / 12.0f / dx * (data[center + 2 * hz] - data[center - 2 * hz]);
						break;
					}
					case 6: {
						// dfw/dx
						output[index] =
							3.0f / 4.0f / dx * (data[center + hx] - data[center - hx]) -
							3.0f / 20.0f / dx * (data[center + 2 * hx] - data[center - 2 * hx]) +
							1.0f / 60.0f / dx * (data[center + 3 * hx] - data[center - 3 * hx]);

						// dfw/dy
						output[index + 1] =
							3.0f / 4.0f / dx * (data[center + hy] - data[center - hy]) -
							3.0f / 20.0f / dx * (data[center + 2 * hy] - data[center - 2 * hy]) +
							1.0f / 60.0f / dx * (data[center + 3 * hy] - data[center - 3 * hy]);

						// dfw/dz
						output[index + 2] =
							3.0f / 4.0f / dx * (data[center + hz] - data[center - hz]) -
							3.0f / 20.0f / dx * (data[center + 2 * hz] - data[center - 2 * hz]) +
							1.0f / 60.0f / dx * (data[center + 3 * hz] - data[center - 3 * hz]);
						break;
					}
					case 8: {
						// dfw/dx
						output[index] =
							4.0f / 5.0f / dx * (data[center + hx] - data[center - hx]) -
							1.0f / 5.0f / dx * (data[center + 2 * hx] - data[center - 2 * hx]) +
							4.0f / 105.0f / dx * (data[center + 3 * hx] - data[center - 3 * hx]) -
							1.0f / 280.0f / dx * (data[center + 4 * hx] - data[center - 4 * hx]);

						// dfw/dy
						output[index + 1] =
							4.0f / 5.0f / dx * (data[center + hy] - data[center - hy]) -
							1.0f / 5.0f / dx * (data[center + 2 * hy] - data[center - 2 * hy]) +
							4.0f / 105.0f / dx * (data[center + 3 * hy] - data[center - 3 * hy]) -
							1.0f / 280.0f / dx * (data[center + 4 * hy] - data[center - 4 * hy]);

						// dfw/dz
						output[index + 2] =
							4.0f / 5.0f / dx * (data[center + hz] - data[center - hz]) -
							1.0f / 5.0f / dx * (data[center + 2 * hz] - data[center - 2 * hz]) +
							4.0f / 105.0f / dx * (data[center + 3 * hz] - data[center - 3 * hz]) -
							1.0f / 280.0f / dx * (data[center + 4 * hz] - data[center - 4 * hz]);
						break;
					}
					}

				}
			}
		}
	}
	return 0;
}

int computeLaplacian(dataKernel* kernel, int comps, float dx, int size, int nOrder, float *output)
{
	int x, y, z, w;
	float* data = kernel->data;
	int   hz = kernel->hy*kernel->hx*comps, hy = kernel->hx*comps, hx = comps;
	for (x = 0; x < size; x++)
	{
		for (y = 0; y < size; y++)
		{
			for (z = 0; z < size; z++)
			{
				//Each component
				for (w = 0; w < 3; w++)
				{
					int z_ = nOrder / 2 + z + kernel->z,
						y_ = nOrder / 2 + y + kernel->y,
						x_ = nOrder / 2 + x + kernel->x;
					int index = z * size*size*comps + y * size*comps + x * comps + w;
					int center = z_ * hz + y_ * hy + x_ * hx + w;
					switch (nOrder)
					{
					case 4: {
						// du2w/dxdx
						output[index] =
							SecFiniteDiff4(dx,
								data[center - 2 * hx],
								data[center - hx],
								data[center],
								data[center + hx],
								data[center + 2 * hx]) +

							SecFiniteDiff4(dx,
								data[center - 2 * hy],
								data[center - hy],
								data[center],
								data[center + hy],
								data[center + 2 * hy]) +

							SecFiniteDiff4(dx,
								data[center - 2 * hz],
								data[center - hz],
								data[center],
								data[center + hz],
								data[center + 2 * hz]);


						break;
					} case 6: {

						// du2w/dxdx
						output[index + 0] =
							SecFiniteDiff6(dx,
								data[center - 3 * hx],
								data[center - 2 * hx],
								data[center - hx],
								data[center],
								data[center + hx],
								data[center + 2 * hx],
								data[center + 3 * hx]) +

							SecFiniteDiff6(dx,
								data[center - 3 * hy],
								data[center - 2 * hy],
								data[center - hy],
								data[center],
								data[center + hy],
								data[center + 2 * hy],
								data[center + 3 * hy]) +

							SecFiniteDiff6(dx,
								data[center - 3 * hz],
								data[center - 2 * hz],
								data[center - hz],
								data[center],
								data[center + hz],
								data[center + 2 * hz],
								data[center + 3 * hz]);
						break;
					} case 8: {

						// du2w/dxdx
						output[index + 0] =
							SecFiniteDiff8(dx,
								data[center - 4 * hx],
								data[center - 3 * hx],
								data[center - 2 * hx],
								data[center - hx],
								data[center],
								data[center + hx],
								data[center + 2 * hx],
								data[center + 3 * hx],
								data[center + 4 * hx]) +

							SecFiniteDiff8(dx,
								data[center - 4 * hy],
								data[center - 3 * hy],
								data[center - 2 * hy],
								data[center - hy],
								data[center],
								data[center + hy],
								data[center + 2 * hy],
								data[center + 3 * hy],
								data[center + 4 * hy]) +

							SecFiniteDiff8(dx,
								data[center - 4 * hz],
								data[center - 3 * hz],
								data[center - 2 * hz],
								data[center - hz],
								data[center],
								data[center + hz],
								data[center + 2 * hz],
								data[center + 3 * hz],
								data[center + 4 * hz]);
						break;
					}
					}
				}

			}
		}
	}
	return 0;
}

int computeHessian(dataKernel* kernel, int comps, float dx, int size, int nOrder, float *output)
{
	int x, y, z, w;
	float* data = kernel->data;
	int hz = kernel->hy*kernel->hx*comps, hy = kernel->hx*comps, hx = comps;
	for (x = 0; x < size; x++)
	{
		for (y = 0; y < size; y++)
		{
			for (z = 0; z < size; z++)
			{
				//Each component
				for (w = 0; w < comps; w++)
				{
					int z_ = nOrder / 2 + z + kernel->z,
						y_ = nOrder / 2 + y + kernel->y,
						x_ = nOrder / 2 + x + kernel->x;
					int index = z * size*size*comps * 6 + y * size*comps * 6 + x * comps * 6 + w * 6;
					int center = z_ * hz + y_ * hy + x_ * hx + w;
					switch (nOrder)
					{
					case 4: {
						// du2w/dxdx
						output[index + 0] =
							SecFiniteDiff4(dx,
								data[center - 2 * hx],
								data[center - hx],
								data[center],
								data[center + hx],
								data[center + 2 * hx]);

						// du2w/dxdy
						output[index + 1] =
							CrossFiniteDiff4(dx,
								data[center + 2 * hy + 2 * hx],
								data[center + 2 * hy - 2 * hx],
								data[center - 2 * hy - 2 * hx],
								data[center - 2 * hy + 2 * hx],
								data[center + hy + hx],
								data[center + hy - hx],
								data[center - hy - hx],
								data[center - hy + hx]);

						// du2w/dxdz
						output[index + 2] =
							CrossFiniteDiff4(dx,
								data[center + 2 * hz + 2 * hx],
								data[center + 2 * hz - 2 * hx],
								data[center - 2 * hz - 2 * hx],
								data[center - 2 * hz + 2 * hx],
								data[center + hz + hx],
								data[center + hz - hx],
								data[center - hz - hx],
								data[center - hz + hx]);

						// du2w/dydy
						output[index + 3] =
							SecFiniteDiff4(dx,
								data[center - 2 * hy],
								data[center - hy],
								data[center],
								data[center + hy],
								data[center + 2 * hy]);

						// du2w/dydz
						output[index + 4] =
							CrossFiniteDiff4(dx,
								data[center + 2 * hz + 2 * hy],
								data[center + 2 * hz - 2 * hy],
								data[center - 2 * hz - 2 * hy],
								data[center - 2 * hz + 2 * hy],
								data[center + hz + hy],
								data[center + hz - hy],
								data[center - hz - hy],
								data[center - hz + hy]);

						// du2w/dzdz
						output[index + 5] =
							SecFiniteDiff4(dx,
								data[center - 2 * hz],
								data[center - hz],
								data[center],
								data[center + hz],
								data[center + 2 * hz]);

						break;
					} case 6: {

						// du2w/dxdx
						output[index + 0] =
							SecFiniteDiff6(dx,
								data[center - 3 * hx],
								data[center - 2 * hx],
								data[center - hx],
								data[center],
								data[center + hx],
								data[center + 2 * hx],
								data[center + 3 * hx]);

						// du2w/dxdy
						output[index + 1] =
							CrossFiniteDiff6(dx,
								data[center + 3 * hy + 3 * hx],
								data[center + 3 * hy - 3 * hx],
								data[center - 3 * hy - 3 * hx],
								data[center - 3 * hy + 3 * hx],
								data[center + 2 * hy + 2 * hx],
								data[center + 2 * hy - 2 * hx],
								data[center - 2 * hy - 2 * hx],
								data[center - 2 * hy + 2 * hx],
								data[center + hy + hx],
								data[center + hy - hx],
								data[center - hy - hx],
								data[center - hy + hx]);

						// du2w/dxdz
						output[index + 2] =
							CrossFiniteDiff6(dx,
								data[center + 3 * hz + 3 * hx],
								data[center + 3 * hz - 3 * hx],
								data[center - 3 * hz - 3 * hx],
								data[center - 3 * hz + 3 * hx],
								data[center + 2 * hz + 2 * hx],
								data[center + 2 * hz - 2 * hx],
								data[center - 2 * hz - 2 * hx],
								data[center - 2 * hz + 2 * hx],
								data[center + hz + hx],
								data[center + hz - hx],
								data[center - hz - hx],
								data[center - hz + hx]);

						// du2w/dydy
						output[index + 3] =
							SecFiniteDiff6(dx,
								data[center - 3 * hy],
								data[center - 2 * hy],
								data[center - hy],
								data[center],
								data[center + hy],
								data[center + 2 * hy],
								data[center + 3 * hy]);

						// du2w/dydz
						output[index + 4] =
							CrossFiniteDiff6(dx,
								data[center + 3 * hz + 3 * hy],
								data[center + 3 * hz - 3 * hy],
								data[center - 3 * hz - 3 * hy],
								data[center - 3 * hz + 3 * hy],
								data[center + 2 * hz + 2 * hy],
								data[center + 2 * hz - 2 * hy],
								data[center - 2 * hz - 2 * hy],
								data[center - 2 * hz + 2 * hy],
								data[center + hz + hy],
								data[center + hz - hy],
								data[center - hz - hy],
								data[center - hz + hy]);

						// du2w/dzdz
						output[index + 5] =
							SecFiniteDiff6(dx,
								data[center - 3 * hz],
								data[center - 2 * hz],
								data[center - hz],
								data[center],
								data[center + hz],
								data[center + 2 * hz],
								data[center + 3 * hz]);
						break;
					} case 8: {

						// du2w/dxdx
						output[index + 0] =
							SecFiniteDiff8(dx,
								data[center - 4 * hx],
								data[center - 3 * hx],
								data[center - 2 * hx],
								data[center - hx],
								data[center],
								data[center + hx],
								data[center + 2 * hx],
								data[center + 3 * hx],
								data[center + 4 * hx]);

						// du2w/dxdy
						output[index + 1] =
							CrossFiniteDiff8(dx,
								data[center + 4 * hy + 4 * hx],
								data[center + 4 * hy - 4 * hx],
								data[center - 4 * hy - 4 * hx],
								data[center - 4 * hy + 4 * hx],
								data[center + 3 * hy + 3 * hx],
								data[center + 3 * hy - 3 * hx],
								data[center - 3 * hy - 3 * hx],
								data[center - 3 * hy + 3 * hx],
								data[center + 2 * hy + 2 * hx],
								data[center + 2 * hy - 2 * hx],
								data[center - 2 * hy - 2 * hx],
								data[center - 2 * hy + 2 * hx],
								data[center + hy + hx],
								data[center + hy - hx],
								data[center - hy - hx],
								data[center - hy + hx]);

						// du2w/dxdz
						output[index + 2] =
							CrossFiniteDiff8(dx,
								data[center + 4 * hz + 4 * hx],
								data[center + 4 * hz - 4 * hx],
								data[center - 4 * hz - 4 * hx],
								data[center - 4 * hz + 4 * hx],
								data[center + 3 * hz + 3 * hx],
								data[center + 3 * hz - 3 * hx],
								data[center - 3 * hz - 3 * hx],
								data[center - 3 * hz + 3 * hx],
								data[center + 2 * hz + 2 * hx],
								data[center + 2 * hz - 2 * hx],
								data[center - 2 * hz - 2 * hx],
								data[center - 2 * hz + 2 * hx],
								data[center + hz + hx],
								data[center + hz - hx],
								data[center - hz - hx],
								data[center - hz + hx]);

						// du2w/dydy
						output[index + 3] =
							SecFiniteDiff8(dx,
								data[center - 4 * hy],
								data[center - 3 * hy],
								data[center - 2 * hy],
								data[center - hy],
								data[center],
								data[center + hy],
								data[center + 2 * hy],
								data[center + 3 * hy],
								data[center + 4 * hy]);

						// du2w/dydz
						output[index + 4] =
							CrossFiniteDiff8(dx,
								data[center + 4 * hz + 4 * hy],
								data[center + 4 * hz - 4 * hy],
								data[center - 4 * hz - 4 * hy],
								data[center - 4 * hz + 4 * hy],
								data[center + 3 * hz + 3 * hy],
								data[center + 3 * hz - 3 * hy],
								data[center - 3 * hz - 3 * hy],
								data[center - 3 * hz + 3 * hy],
								data[center + 2 * hz + 2 * hy],
								data[center + 2 * hz - 2 * hy],
								data[center - 2 * hz - 2 * hy],
								data[center - 2 * hz + 2 * hy],
								data[center + hz + hy],
								data[center + hz - hy],
								data[center - hz - hy],
								data[center - hz + hy]);

						// du2w/dzdz
						output[index + 5] =
							SecFiniteDiff8(dx,
								data[center - 4 * hz],
								data[center - 3 * hz],
								data[center - 2 * hz],
								data[center - hz],
								data[center],
								data[center + hz],
								data[center + 2 * hz],
								data[center + 3 * hz],
								data[center + 4 * hz]);
						break;
					}
					}
				}

			}
		}
	}
	return 0;
}

int getvelocity_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getVelocity(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvelocitygradient_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][9],
	int len_a, int len_d)
{
	return getVelocityGradient(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvelocityhessian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][18],
	int len_a, int len_d)
{
	return getVelocityHessian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvelocitylaplacian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getVelocityLaplacian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getmagneticfield_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getMagneticField(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvectorpotential_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getVectorPotential(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvelocityandpressure_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][4],
	int len_a, int len_d)
{
	return getVelocityAndPressure(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getpressurehessian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][6],
	int len_a, int len_d)
{
	return getPressureHessian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getmagneticfieldlaplacian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getMagneticFieldLaplacian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvectorpotentiallaplacian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getVectorPotentialLaplacian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getdensityhessian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][6],
	int len_a, int len_d)
{
	return getDensityHessian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getdensitygradient_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getDensityGradient(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getdensity_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[],
	int len_a, int len_d)
{
	return getDensity(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getpressure_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[],
	int len_a, int len_d)
{
	return getPressure(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvectorpotentialhessian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][18],
	int len_a, int len_d)
{
	return getVectorPotentialHessian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getpressuregradient_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][3],
	int len_a, int len_d)
{
	return getPressureGradient(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getmagneticfieldgradient_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][9],
	int len_a, int len_d)
{
	return getMagneticFieldGradient(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getvectorpotentialgradient_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][9],
	int len_a, int len_d)
{
	return getVectorPotentialGradient(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getmagneticfieldhessian_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][18],
	int len_a, int len_d)
{
	return getMagneticFieldHessian(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getinvariant_(char *authToken,
	char *dataset, float *time,
	int *spatial, int *temporal,
	int *count, float datain[][3], float dataout[][2],
	int len_a, int len_d)
{
	return getInvariant(authToken,
		dataset, *time,
		*spatial, *temporal,
		*count, datain, dataout);
}

int getcutout_(char *authToken,
	char *dataset, char *field, int *T,
	int *x_start, int *y_start, int *z_start,
	int *x_end, int *y_end, int *z_end,
	int *x_step, int *y_step, int *z_step, int *filter_width,
	float dataout[])
{
	return getCutout(authToken, dataset, field, *T, *x_start, *y_start, *z_start, *x_end, *y_end, *z_end, *x_step, *y_step, *z_step, *filter_width, dataout);
}
