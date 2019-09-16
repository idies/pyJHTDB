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


/* $Id: turblib.h,v 1.18 2009-12-06 16:01:01 eric Exp $ */
#ifndef TURBLIB_H_
#define TURBLIB_H_

#include "soapH.h"

#ifdef  __cplusplus
extern "C" {
#endif//__cplusplus

#define TURBLIB_VERSION "0.2"
#define TURB_ERROR_LENGTH 2048

	/* Global gSOAP handle
	 * TODO: Figure out how to pass this to/from Fortran
	 */
	extern struct soap __jhuturbsoap;

	/* Spatial Interpolation Options */
	enum SpatialInterpolation
	{
		/* Spatial Interpolation Flags for GetVelocity & GetVelocityAndPressure */
		NoSInt = 0, /* No spatial interpolation */
		Lag4 = 4,   /* 4th order Lagrangian interpolation in space */
		Lag6 = 6,   /* 4th order Lagrangian interpolation in space */
		Lag8 = 8,   /* 4th order Lagrangian interpolation in space */
		M1Q4 = 104,
		M1Q6 = 106,
		M1Q8 = 108,
		M1Q10 = 110,
		M1Q12 = 112,
		M1Q14 = 114,
		M2Q4 = 204,
		M2Q6 = 206,
		M2Q8 = 208,
		M2Q10 = 210,
		M2Q12 = 212,
		M2Q14 = 214,
		M3Q4 = 304,
		M3Q6 = 306,
		M3Q8 = 308,
		M3Q10 = 310,
		M3Q12 = 312,
		M3Q14 = 314,
		M4Q4 = 404,
		M4Q6 = 406,
		M4Q8 = 408,
		M4Q10 = 410,
		M4Q12 = 412,
		M4Q14 = 414,

		/* Spatial Differentiation and Interpolation Flags for GetVelocityGradient
		 * and GetPressureGradient. */
		FD4NoInt = 40, /* 4th order finite differential scheme for grid values, no spatial interpolation */
		FD6NoInt = 60, /* 6th order finite differential scheme for grid values, no spatial interpolation */
		FD8NoInt = 80, /* 8th order finite differential scheme for grid values, no spatial interpolation */
		FD4Lag4 = 44,  /* 4th order finite differential scheme for grid values, 4th order Lagrangian interpolation in space */

		/* Old names, for backward compatibility */
		NoSpatialInterpolation = 0,
		Lagrangian4thOrder = 4,
		Lagrangian6thOrder = 6,
		Lagrangian8thOrder = 8
	};

	/* Temporal Interpolation Options */
	enum TemporalInterpolation
	{
		NoTInt = 0,   /* No temporal interpolation */
		PCHIPInt = 1, /* Piecewise cubic Hermite interpolation in time */

		/* Old names, for backward compatibility */
		NoTemporalInterpolation = 0,
		PCHIPInterpolation = 1
	};

	typedef struct
	{
		int x;
		int y;
		int z;
		float value;
	} ThresholdInfo;

	typedef enum
	{
		//isotropic1024old = 1,
		//isotropic1024fine_old = 2,
		//mhd1024 = 3,
		//isotropic1024coarse = 4,
		//isotropic1024fine = 5,
		//channel = 6,
		//custom_dataset = 7,
		//mixing = 8
		mhd1024 = 3,
		isotropic1024coarse = 4,
		isotropic1024fine = 5,
		channel = 6,
		mixing = 7,
		rmhd = 8,
		isotropic4096 = 10,
		rotstrat4096 = 11,
		transition_bl = 12,
		channel5200 = 13
	} TurbDataset;

	typedef enum
	{
		turb_velocity = 0,
		turb_pressure = 1,
		turb_magnetic = 2,
		turb_potential = 3,
		turb_vp = 4,
		turb_density = 5,
		turb_temperature = 6,
	} TurbField;

	typedef struct
	{
		float dx, dt;
		int size;
	} set_info;

	typedef struct
	{
		char prefix;
		int comps;
	} turb_fn;

	typedef struct dataBlock dataBlock;

	struct dataBlock
	{
		float *data;
		TurbField function;
		TurbDataset dataset;
		int xl, yl, zl;  //Lower bounds
		int hx, hy, hz, comps; //Dimensions
		dataBlock* next;
	};

//    typedef struct cutoutFile cutoutFile;
//    struct cutoutFile
//    {
//        hid_t file;
//        int start[4], size[4];
//        int contents[5];
//        dataBlock *data[4][1024];
//        TurbDataset dataset;
//        cutoutFile *next;
//    };

	typedef struct
	{
		float* data;
		int x, y, z;
		int hx, hy, hz, comps;
		int persist;
	} dataKernel;

	/* C */
	char *getVersion();
	void soapinit(void);
	void soapdestroy(void);

	/* Fortran */
	char *getVersion_();
	void soapinit_(void);
	void soapdestroy_(void);

	/* C */
	char *turblibGetErrorString(void);
	int turblibGetErrorNumber(void);
	void turblibPrintError(void);
	void turblibSetExitOnError(int);

	/* Fortran */
	void turblibgeterrorstring_(char *dest, int len);
	int turblibgeterrornumber_(void);
	void turblibprinterror_(void);
	void turblibsetexitonerror_(int *);

	/* C */
	int getVelocitySoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getvelocity_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getThreshold(char *authToken,
		char *dataset, char *field, float time, float threshold,
		enum SpatialInterpolation spatial,
		int x_start, int y_start, int z_start, int x_end, int y_end, int z_end,
		ThresholdInfo **dataout, int *result_size);

	/* Fortran */
	int getthreshold_(char *authToken,
		char *dataset, char *field, float *time, float *threshold,
		int *spatial,
		int *x_start, int *y_start, int *z_start, int *x_end, int *y_end, int *z_end,
		ThresholdInfo **dataout, int *result_size);

	void deallocate_array_(ThresholdInfo **threshold_array);

	/* C */
	int getBoxFilter(char *authToken,
		char *dataset, char *field, float time, float filterwidth,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getboxfilter_(char *authToken,
		char *dataset, char *field, float *time, float *filterwidth,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getBoxFilterSGS(char *authToken,
		char *dataset, char *field, float time, float filterwidth,
		int count, float datain[][3], float dataout[][6]);

	/* Fortran */
	int getboxfiltersgs_(char *authToken,
		char *dataset, char *field, float *time, float *filterwidth,
		int *count, float datain[][3], float dataout[][6],
		int len_a, int len_d);

	/* C */
	int getBoxFilterSGSscalar(char *authToken,
		char *dataset, char *field, float time, float filterwidth,
		int count, float datain[][3], float dataout[]);

	/* Fortran */
	int getboxfiltersgsscalar_(char *authToken,
		char *dataset, char *field, float *time, float *filterwidth,
		int *count, float datain[][3], float dataout[],
		int len_a, int len_d);

	/* C */
	int getBoxFilterSGSvector(char *authToken,
		char *dataset, char *field, float time, float filterwidth,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getboxfiltersgsvector_(char *authToken,
		char *dataset, char *field, float *time, float *filterwidth,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getBoxFilterSGSsymtensor(char *authToken,
		char *dataset, char *field, float time, float filterwidth,
		int count, float datain[][3], float dataout[][6]);

	/* Fortran */
	int getboxfiltersgssymtensor_(char *authToken,
		char *dataset, char *field, float *time, float *filterwidth,
		int *count, float datain[][3], float dataout[][6],
		int len_a, int len_d);

	/* C */
	int getBoxFilterSGStensor(char *authToken,
		char *dataset, char *field, float time, float filterwidth,
		int count, float datain[][3], float dataout[][9]);

	/* Fortran */
	int getboxfiltersgstensor_(char *authToken,
		char *dataset, char *field, float *time, float *filterwidth,
		int *count, float datain[][3], float dataout[][9],
		int len_a, int len_d);

	/* C */
	int getBoxFilterGradient(char *authToken,
		char *dataset, char *field, float time,
		float filterwidth, float spacing,
		int count, float datain[][3], float dataout[][9]);

	/* Fortran */
	int getboxfiltergradient_(char *authToken,
		char *dataset, char*field, float *time,
		float *filterwidth, float* spacing,
		int *count, float datain[][3], float dataout[][9],
		int len_a, int len_d);

	/* C */
	int getVelocityAndPressureSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][4]);

	/* Fortran */
	int getvelocityandpressure_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][4],
		int len_a, int len_d);

	/* C */
	int getPressureHessianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][6]);

	/* Fortran */
	int getpressurehessian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][6],
		int len_a, int len_d);

	/* C */
	int getVelocityGradientSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][9]);

	/* Fortran */
	int getvelocitygradient_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][9],
		int len_a, int len_d);

	/* C */
	int getMagneticFieldGradientSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][9]);

	/* Fortran */
	int getmagneticfieldgradient_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][9],
		int len_a, int len_d);

	/* C */
	int getVectorPotentialGradientSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][9]);

	/* Fortran */
	int getvectorpotentialgradient_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][9],
		int len_a, int len_d);

	/* C */
	int getPressureGradientSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getpressuregradient_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getVelocityHessianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][18]);

	/* Fortran */
	int getvelocityhessian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][18],
		int len_a, int len_d);

	/* C */
	int getVelocityLaplacianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getvelocitylaplacian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getMagneticFieldHessianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][18]);

	/* Fortran */
	int getmagneticfieldhessian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][18],
		int len_a, int len_d);

	/* C */
	int getMagneticFieldLaplacianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getmagneticfieldlaplacian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getVectorPotentialHessianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][18]);

	/* Fortran */
	int getvectorpotentialhessian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][18],
		int len_a, int len_d);

	/* C */
	int getVectorPotentialLaplacianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getvectorpotentiallaplacian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int nullOp(char *authToken, int count,
		float datain[][3], float dataout[][3]);

	/* Fortran */
	int nullop_(char *authToken, int *count,
		float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getForce(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getforce_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/*C*/
	int getPosition(char *authToken,
		char *dataset, float startTime, float endTime,
		float dt,
		enum SpatialInterpolation spatial,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getposition_(char *authToken,
		char *dataset, float *startTime, float *endTime,
		float *dt,
		int *spatial,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getRawVelocity(char *authToken,
		char *dataset, int time_step,
		int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[]);

	/* Fortran */
	int getrawvelocity_(char *authToken,
		char *dataset, int *time_step,
		int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
		float dataout[]);

	/* C */
	int getMagneticFieldSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getmagneticfield_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getRawMagneticField(char *authToken,
		char *dataset, int time_step,
		int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[]);

	/* Fortran */
	int getrawmagneticfield_(char *authToken,
		char *dataset, int *time_step,
		int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
		float dataout[]);

	/* C */
	int getVectorPotentialSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getvectorpotential_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getRawVectorPotential(char *authToken,
		char *dataset, int time_step,
		int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[]);

	/* Fortran */
	int getrawvectorpotential_(char *authToken,
		char *dataset, int *time_step,
		int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth,
		float dataout[]);

	/* C */
	int getPressureSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[]);

	/* Fortran */
	int getpressure_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[],
		int len_a, int len_d);

	/* C */
	int getRawPressure(char *authToken,
		char *dataset, int time_step,
		int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[]);

	/* Fortran */
	int getrawpressure_(char *authToken,
		char *dataset, int *time_step,
		int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth, float dataout[]);

	/* C */
	int getDensitySoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[]);

	/* Fortran */
	int getdensity_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[],
		int len_a, int len_d);

	/* C */
	int getDensityGradientSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	/* Fortran */
	int getdensitygradient_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][3],
		int len_a, int len_d);

	/* C */
	int getDensityHessianSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][6]);

	/* Fortran */
	int getdensityhessian_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][6],
		int len_a, int len_d);

	/* C */
	int getInvariantSoap(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][2]);

	/* C */
	int getInvariant(char *authToken,
		char *dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][2]);

	/* Fortran */
	int getinvariant_(char *authToken,
		char *dataset, float *time,
		int *spatial, int *temporal,
		int *count, float datain[][3], float dataout[][2],
		int len_a, int len_d);

	/* C */
	int getRawDensity(char *authToken,
		char *dataset, int time_step,
		int X, int Y, int Z, int Xwidth, int Ywidth, int Zwidth, char dataout[]);

	/* Fortran */
	int getrawdensity_(char *authToken,
		char *dataset, int *time_step,
		int *X, int *Y, int *Z, int *Xwidth, int *Ywidth, int *Zwidth, float dataout[]);

	/* C */
	int getCutout(char *authToken,
		char *dataset, char *field, int time_step,
		int x_start, int y_start, int z_start,
		int x_end, int y_end, int z_end,
		int x_step, int y_step, int z_step, int filter_width,
		float dataout[]);

	/* Fortran */
	int getcutout_(char *authToken,
		char *dataset, char *field, int *time_step,
		int *x_start, int *y_start, int *z_start,
		int *x_end, int *y_end, int *z_end,
		int *x_step, int *y_step, int *z_step, int *filter_width,
		float dataout[]);


	/* Local vs Soap functions */

	int getVelocityGradientLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float input[][3], float output[][9]);

	int getVelocityLaplacianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float input[][3], float output[][3]);

	int getVelocityHessianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float input[][3], float output[][18]);

	int getVelocityAndPressureLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][4]);

	int getMagneticFieldLaplacianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	int getVectorPotentialLaplacianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	int getPressureHessianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][6]);

	int getMagneticFieldHessianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][18]);

	int getVectorPotentialHessianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][18]);

	int getPressureGradientLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	int getMagneticFieldGradientLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][9]);

	int getVectorPotentialGradientLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][9]);

	int getDensityGradientLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][3]);

	int getDensityHessianLocal(TurbDataset dataset, float time,
		enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
		int count, float datain[][3], float dataout[][6]);

	/* Cutout functions */

	int lagrangianInterp(int comps, float *kernel, float position[3], int nOrder, float dx, float *result);
	int lagrangianInterp2(int comps, dataKernel *kernel, float position[3], int nOrder, float dx, float *result);
	int pchipInterp(int comps, float *data[4], float time, int timestep, float dt, float *result);

	int computeGradient(dataKernel* kernel, int comps, float dx, int size, int nOrder, float *output);
	int computeHessian(dataKernel* kernel, int comps, float dx, int size, int nOrder, float *output);
	int computeLaplacian(dataKernel* kernel, int comps, float dx, int size, int nOrder, float *output);

#define SecFiniteDiff4(dx, x1, x2, x3, x4, x5) 4.0f / 3.0f / dx / dx * (x2 + x4 - 2.0f * x3) - 1.0f / 12.0f / dx / dx * (x1 + x5 - 2.0f * x3)

#define SecFiniteDiff6(dx, x1, x2, x3, x4, x5, x6, x7) 3.0f / 2.0f / dx / dx * (x3 + x5 - 2.0f * x4) - 3.0f / 20.0f / dx / dx * (x2 + x6 - 2.0f * x4) + 1.0f / 90.0f / dx / dx * (x1 + x7 - 2.0f * x4)

#define SecFiniteDiff8(dx, x1, x2, x3, x4, x5, x6, x7, x8, x9) 792.0f / 591.0f / dx / dx * (x4 + x6 - 2.0f * x5) - 207.0f / 2955.0f / dx / dx * (x3 + x7 - 2.0f * x5) - 104.0f / 8865.0f / dx / dx * (x2 + x8 - 2.0f * x5) + 9.0f / 3152.0f / dx / dx * (x1 + x9 - 2.0f * x5)

#define CrossFiniteDiff4(dx, x1, x2, x3, x4, x5, x6, x7, x8) 1.0f / 3.0f / dx / dx * (x5 + x7 - x6 - x8) - 1.0f / 48.0f / dx / dx * (x1 + x3 - x2 - x4)

#define CrossFiniteDiff6(dx, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)  3.0f / 8.0f / dx / dx * (x9 + x11 - x10 - x12) - 3.0f / 80.0f / dx / dx * (x5 + x7 - x6 - x8) + 1.0f / 360.0f / dx / dx * (x1 + x3 - x2 - x4)

#define CrossFiniteDiff8(dx, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16) 14.0f / 35.0f / dx / dx * (x13 + x15 - x14 - x16) - 1.0f / 20.0f / dx / dx * (x9 + x11 - x10 - x12) + 2.0f / 315.0f / dx / dx * (x5 + x7 - x6 - x8) - 1.0f / 2240.0f / dx / dx * (x1 + x3 - x2 - x4)

#ifdef  __cplusplus
}
#endif//__cplusplus

int getVelocity(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getVelocityAndPressure(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][4]);

int getVelocityGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9]);

int getVelocityHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18]);

int getVelocityLaplacian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getPressure(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[]);

int getPressureGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getPressureHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][6]);

int getMagneticField(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getMagneticFieldGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9]);

int getMagneticFieldHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18]);

int getMagneticFieldLaplacian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getVectorPotential(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getVectorPotentialGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][9]);

int getVectorPotentialHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][18]);

int getVectorPotentialLaplacian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getDensity(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[]);

int getDensityGradient(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][3]);

int getDensityHessian(char *authToken,
	char *dataset, float time,
	enum SpatialInterpolation spatial, enum TemporalInterpolation temporal,
	int count, float datain[][3], float dataout[][6]);

#endif//TURBLIB_H_

