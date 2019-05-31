/******************************************************************************
*
*    Copyright 2014 Johns Hopkins University
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*
*   Contact: turbulence@pha.jhu.edu
*   Website: http://turbulence.pha.jhu.edu/
*
******************************************************************************/
/*
 *  These subroutines are meant to be called from a pyJHTDB.libTDB object.
 *  You're on your own otherwise.
 *
 * */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <turblib.h>

// relative error of single precision numbers
float float_error = 1e-6;

int free_threshold_array(ThresholdInfo *data)
{
    free(data);
    return 0;
}

int getBline(
        char *authToken,
        char *dataset,
        float time,
        int ssteps,
        float ds,
        enum SpatialInterpolation spatial,
        enum TemporalInterpolation temporal,
        int count,
        float x[][3])
{
    int p, s;
    float bsize;
	float (*y)[3] = malloc(count*sizeof(float[3]));
	float (*bfield0)[3] = malloc(count*sizeof(float[3]));
	float (*bfield1)[3] = malloc(count*sizeof(float[3]));

    for (s = 1; s <= ssteps; s++)
    {
        // get bhat0 and compute y
        getMagneticField(authToken, dataset, time, spatial, temporal, count, x, bfield0);
        for (p = 0; p < count; p++)
        {
            bsize = sqrt(bfield0[p][0]*bfield0[p][0] + bfield0[p][1]*bfield0[p][1] + bfield0[p][2]*bfield0[p][2]);
            bfield0[p][0] /= bsize;
            bfield0[p][1] /= bsize;
            bfield0[p][2] /= bsize;
            y[p][0] = x[p][0] + ds*bfield0[p][0];
            y[p][1] = x[p][1] + ds*bfield0[p][1];
            y[p][2] = x[p][2] + ds*bfield0[p][2];
        }
        // get bhat1 and compute next position
        getMagneticField(authToken, dataset, time, spatial, temporal, count, y, bfield1);
        for (p = 0; p < count; p++)
        {
            bsize = sqrt(bfield1[p][0]*bfield1[p][0] + bfield1[p][1]*bfield1[p][1] + bfield1[p][2]*bfield1[p][2]);
            bfield1[p][0] /= bsize;
            bfield1[p][1] /= bsize;
            bfield1[p][2] /= bsize;
            x[p+count][0] = x[p][0] + .5*ds*(bfield0[p][0] + bfield1[p][0]);
            x[p+count][1] = x[p][1] + .5*ds*(bfield0[p][1] + bfield1[p][1]);
            x[p+count][2] = x[p][2] + .5*ds*(bfield0[p][2] + bfield1[p][2]);
        }
        x += count;
    }
    free(bfield0);
    free(bfield1);
    free(y);
    return 0;
}

int getRectangularBoundedBline(
        char *authToken,
        char *dataset,
        float time,
        int maxsteps,
        float ds,
        enum SpatialInterpolation spatial,
        enum TemporalInterpolation temporal,
        int count,
        float traj[][3],
        int traj_length[],
        float xmin, float xmax,
        float ymin, float ymax,
        float zmin, float zmax)
{
    // traj will contain all the trajectories, and traj_length will contain the lengths of all trajectories
    // the shape of traj is assumed to be (number of particles, allocated_traj_length, 3)
    int p, s;
    float bsize;
    float (*x)[3];
	float y[1][3];
	float bfield0[1][3];
	float bfield1[1][3];

    // keep a little way away from the actual edge...
    xmin *= 1 + (xmin > 0 ? -1 : +1)*float_error;
    ymin *= 1 + (ymin > 0 ? -1 : +1)*float_error;
    zmin *= 1 + (zmin > 0 ? -1 : +1)*float_error;
    xmax *= 1 + (xmax > 0 ? +1 : -1)*float_error;
    ymax *= 1 + (ymax > 0 ? +1 : -1)*float_error;
    zmax *= 1 + (zmax > 0 ? +1 : -1)*float_error;

    //loop after particles
    for (p = 0; p < count; p++)
    {
        // x is at the start of the current trajectory
        x = traj + p*(maxsteps+1);
        for (s = 1;
             (s <= maxsteps) && !(((*x)[0] < xmin || (*x)[0] > xmax)
                               || ((*x)[1] < ymin || (*x)[1] > ymax)
                               || ((*x)[2] < zmin || (*x)[2] > zmax)); s++)
        {
            // get bhat0 and compute y
            getMagneticField(authToken, dataset, time, spatial, temporal, 1, x, bfield0);
            bsize = sqrt(bfield0[0][0]*bfield0[0][0]
                       + bfield0[0][1]*bfield0[0][1]
                       + bfield0[0][2]*bfield0[0][2]);
            bfield0[0][0] /= bsize;
            bfield0[0][1] /= bsize;
            bfield0[0][2] /= bsize;
            y[0][0] = x[0][0] + ds*bfield0[0][0];
            y[0][1] = x[0][1] + ds*bfield0[0][1];
            y[0][2] = x[0][2] + ds*bfield0[0][2];
            // get bhat1 and compute next position
            getMagneticField(authToken, dataset, time, spatial, temporal, 1, y, bfield1);
            bsize = sqrt(bfield1[0][0]*bfield1[0][0]
                       + bfield1[0][1]*bfield1[0][1]
                       + bfield1[0][2]*bfield1[0][2]);
            bfield1[0][0] /= bsize;
            bfield1[0][1] /= bsize;
            bfield1[0][2] /= bsize;
            x[1][0] = x[0][0] + .5*ds*(bfield0[0][0] + bfield1[0][0]);
            x[1][1] = x[0][1] + .5*ds*(bfield0[0][1] + bfield1[0][1]);
            x[1][2] = x[0][2] + .5*ds*(bfield0[0][2] + bfield1[0][2]);
            x += 1;
        }
        traj_length[p] = s;
    }
    return 0;
}

int getSphericalBoundedBline(
        char *authToken,
        char *dataset,
        float time,
        int maxsteps,
        float ds,
        enum SpatialInterpolation spatial,
        enum TemporalInterpolation temporal,
        int count,
        float traj[][3],
        int traj_length[],
        float ox,
        float oy,
        float oz,
        float radius)
{
    // traj will contain all the trajectories, and traj_length will contain the lengths of all trajectories
    // the shape of traj is assumed to be (number of particles, allocated_traj_length, 3)
    int p, s;
    float bsize, rad0, rad_one, stmp, xtmp, ytmp, ztmp, valtmp;
    float (*x)[3];
	float y[1][3];
	float bfield0[1][3];
	float bfield1[1][3];

    radius *= 1 + float_error;

    //loop after particles
    for (p = 0; p < count; p++)
    {
        // x is at the start of the current trajectory
        x = traj + p*(maxsteps+1);
        for (s = 1; s <= maxsteps; s++)
        {
            // get bhat0 and compute y
            getMagneticField(authToken, dataset, time, spatial, temporal, 1, x, bfield0);
            bsize = sqrt(bfield0[0][0]*bfield0[0][0]
                       + bfield0[0][1]*bfield0[0][1]
                       + bfield0[0][2]*bfield0[0][2]);
            bfield0[0][0] /= bsize;
            bfield0[0][1] /= bsize;
            bfield0[0][2] /= bsize;
            y[0][0] = x[0][0] + ds*bfield0[0][0];
            y[0][1] = x[0][1] + ds*bfield0[0][1];
            y[0][2] = x[0][2] + ds*bfield0[0][2];
            // get bhat1 and compute next position
            getMagneticField(authToken, dataset, time, spatial, temporal, 1, y, bfield1);
            bsize = sqrt(bfield1[0][0]*bfield1[0][0]
                       + bfield1[0][1]*bfield1[0][1]
                       + bfield1[0][2]*bfield1[0][2]);
            bfield1[0][0] /= bsize;
            bfield1[0][1] /= bsize;
            bfield1[0][2] /= bsize;
            x[1][0] = x[0][0] + .5*ds*(bfield0[0][0] + bfield1[0][0]);
            x[1][1] = x[0][1] + .5*ds*(bfield0[0][1] + bfield1[0][1]);
            x[1][2] = x[0][2] + .5*ds*(bfield0[0][2] + bfield1[0][2]);
            if ((rad_one = sqrt((x[1][0] - ox)*(x[1][0] - ox)
                           + (x[1][1] - oy)*(x[1][1] - oy)
                           + (x[1][2] - oz)*(x[1][2] - oz))) > radius)
            {
                //rad0 = sqrt((x[0][0] - ox)*(x[0][0] - ox)
                //          + (x[0][1] - oy)*(x[0][1] - oy)
                //          + (x[0][2] - oz)*(x[0][2] - oz));
                //valtmp = (x[0][0] - ox)*(x[1][0] - ox)
                //       + (x[0][1] - oy)*(x[1][1] - oy)
                //       + (x[0][2] - oz)*(x[1][2] - oz);
                //stmp = (2*(1 - valtmp) + sqrt(4*(1 - 2*valtmp + valtmp*valtmp) + 4*radius*radius*(rad_one*rad_one + rad0*rad0 - 2*valtmp))) / (2 * (rad_one*rad_one + rad0*rad0 - 2*valtmp));
                //xtmp = x[0][0]*(1 - stmp) + stmp * x[1][0];
                //ytmp = x[0][1]*(1 - stmp) + stmp * x[1][1];
                //ztmp = x[0][2]*(1 - stmp) + stmp * x[1][2];
                //fprintf(stderr, "%g %g %g\n", stmp, xtmp, x[1][0]);
                //stmp = (2*(1 - valtmp) - sqrt(4*(1 - 2*valtmp + valtmp*valtmp) + 4*radius*radius*(rad_one*rad_one + rad0*rad0 - 2*valtmp))) / (2 * (rad_one*rad_one + rad0*rad0 - 2*valtmp));
                //fprintf(stderr, "%g %g %g\n", stmp, xtmp, x[1][0]);
                //x[1][0] = xtmp;
                //x[1][1] = ytmp;
                //x[1][2] = ztmp;
                break;
            }
            x += 1;
        }
        traj_length[p] = s;
    }
    return 0;
}

int getMagneticFieldDebug(
        char *authToken,
        char *dataset,
        float time,
        enum SpatialInterpolation spatial,
        enum TemporalInterpolation temporal,
        int count,
        float datain[][3],
        float dataout[][3])
{
    int p;
    for (p = 0; p < count; p++)
    {
        dataout[p][0] = 1;
        dataout[p][1] = 0;
        dataout[p][2] = 0;
    }
    return 0;
}

int getSphericalBoundedBlineDebug(
        char *authToken,
        char *dataset,
        float time,
        int maxsteps,
        float ds,
        enum SpatialInterpolation spatial,
        enum TemporalInterpolation temporal,
        int count,
        float traj[][3],
        int traj_length[],
        float ox,
        float oy,
        float oz,
        float radius)
{
    // traj will contain all the trajectories, and traj_length will contain the lengths of all trajectories
    // the shape of traj is assumed to be (number of particles, allocated_traj_length, 3)
    int p, s;
    float bsize, rad0, rad_one, stmp, xtmp, ytmp, ztmp, valtmp;
    float (*x)[3];
	float bfield0[1][3];

    ds = (ds > 0? +1 : -1) *radius / 100;

    radius *= 1 + float_error;

    //loop after particles
    for (p = 0; p < count; p++)
    {
        // x is at the start of the current trajectory
        x = traj + p*(maxsteps+1);
        for (s = 1; s <= maxsteps; s++)
        {
            rad0 = sqrt((x[0][0] - ox)*(x[0][0] - ox)
                      + (x[0][1] - oy)*(x[0][1] - oy)
                      + (x[0][2] - oz)*(x[0][2] - oz));
            getMagneticFieldDebug(authToken, dataset, time, spatial, temporal, 1, x, bfield0);
            x[1][0] = x[0][0] + ds*bfield0[0][0];
            x[1][1] = x[0][1] + ds*bfield0[0][1];
            x[1][2] = x[0][2] + ds*bfield0[0][2];
            rad_one = sqrt((x[1][0] - ox)*(x[1][0] - ox)
                      + (x[1][1] - oy)*(x[1][1] - oy)
                      + (x[1][2] - oz)*(x[1][2] - oz));
            if (rad_one > radius)
            {
                //rad0 = sqrt((x[0][0] - ox)*(x[0][0] - ox)
                //          + (x[0][1] - oy)*(x[0][1] - oy)
                //          + (x[0][2] - oz)*(x[0][2] - oz));
                //valtmp = (x[0][0] - ox)*(x[1][0] - ox)
                //       + (x[0][1] - oy)*(x[1][1] - oy)
                //       + (x[0][2] - oz)*(x[1][2] - oz);
                //stmp = (2*(1 - valtmp) + sqrt(4*(1 - 2*valtmp + valtmp*valtmp) + 4*radius*radius*(rad_one*rad_one + rad0*rad0 - 2*valtmp))) / (2 * (rad_one*rad_one + rad0*rad0 - 2*valtmp));
                //xtmp = x[0][0]*(1 - stmp) + stmp * x[1][0];
                //ytmp = x[0][1]*(1 - stmp) + stmp * x[1][1];
                //ztmp = x[0][2]*(1 - stmp) + stmp * x[1][2];
                //fprintf(stderr, "%g %g %g\n", stmp, xtmp, x[1][0]);
                //stmp = (2*(1 - valtmp) - sqrt(4*(1 - 2*valtmp + valtmp*valtmp) + 4*radius*radius*(rad_one*rad_one + rad0*rad0 - 2*valtmp))) / (2 * (rad_one*rad_one + rad0*rad0 - 2*valtmp));
                //fprintf(stderr, "%g %g %g\n", stmp, xtmp, x[1][0]);
                //x[1][0] = xtmp;
                //x[1][1] = ytmp;
                //x[1][2] = ztmp;
                break;
            }
            x += 1;
        }
        traj_length[p] = s;
    }
    return 0;
}

int set_custom_dataset_description(float dx, float dt, int size)
{
    extern set_info DataSets[8];
    DataSets[7].dx = dx;
    DataSets[7].dt = dt;
    DataSets[7].size = size;
    return 0;
}

//int read_from_file(
//        char *dataset,
//        int tindex,
//        int xindex, int xsize,
//        int yindex, int ysize,
//        int zindex, int zsize,
//        float *ufield,
//        float *bfield)
//{
//    TurbDataset dataset_ = getDataSet(dataset);
//    hsize_t dim_mem[] = {zsize, ysize, xsize, 3};
//    hid_t mspace = H5Screate_simple(4, dim_mem, NULL);
//    loadSubBlock(
//            dataset_,
//            turb_velocity,
//            tindex,
//            mspace,
//            ufield,
//            xindex, yindex, zindex,
//            xsize,
//            ysize,
//            zsize,
//            0, 0, 0);
//    loadSubBlock(
//            dataset_,
//            turb_magnetic,
//            tindex,
//            mspace,
//            bfield,
//            xindex, yindex, zindex,
//            xsize,
//            ysize,
//            zsize,
//            0, 0, 0);
//    return 0;
//}

//int getCustomPosition(
//        char *authToken,
//        char *dataset,
//        float startTime,
//        float endTime,
//        float dt,
//        enum SpatialInterpolation spatial,
//        int count,
//        float datain[][3],
//        float dataout[][3])
//{
//    TurbDataset dataset_ = getDataSet(dataset);
//    extern set_info DataSets[8];
//    float vel0[count][3], time0;
//    float vel1[count][3];
//    float vel[count][3];
//    float time, deltat;
//    int p, tcounter, nsteps;
//    deltat = DataSets[dataset_].dt;
//    nsteps = ceil((endTime - startTime) / dt);
//    dt = (endTime - startTime) / nsteps;
//    for (p = 0; p < count; p++)
//    {
//        dataout[p][0] = datain[p][0];
//        dataout[p][1] = datain[p][1];
//        dataout[p][2] = datain[p][2];
//    }
//    time = startTime;
//    for (tcounter = 0; tcounter < nsteps; tcounter++)
//    {
//        time0 = deltat * floor(time / deltat);
//        getVelocity(authToken, dataset, time         , spatial, 0, count, datain, vel0);
//        getVelocity(authToken, dataset, time + deltat, spatial, 0, count, datain, vel1);
//        for (p = 0; p < count; p++)
//        {
//            vel[p][0] = ((time - time0) / deltat) * (vel1[p][0] - vel0[p][0]) + vel0[p][0];
//            vel[p][1] = ((time - time0) / deltat) * (vel1[p][1] - vel0[p][1]) + vel0[p][1];
//            vel[p][2] = ((time - time0) / deltat) * (vel1[p][2] - vel0[p][2]) + vel0[p][2];
//            dataout[p][0] = dataout[p][0] + dt * vel[p][0];
//            dataout[p][1] = dataout[p][1] + dt * vel[p][1];
//            dataout[p][2] = dataout[p][2] + dt * vel[p][2];
//        }
//        time += dt;
//    }
//    return 0;
//}
//
//int isBLocal(
//        const char *data_set,
//        int x,
//        int y,
//        int z,
//        int time)
//{
//    TurbDataset d = getDataSet(data_set);
//    return isDataComplete(
//            d,
//            2,
//            x - 16,
//            y - 16,
//            z - 16,
//            32,
//            32,
//            32,
//            time);
//}

//int interpolateBoxFilter(
//        char *authToken,
//        char *dataset, char *field, float time, float filterwidth,
//        int count, float datain[][3], float dataout[][3])
//{
//    extern set_info DataSets[8];
//    float cell_nodes[8*count][3];
//    float cellvalues[8*count][3];
//    int p;
//    float xx, yy, zz;
//    TurbDataset dataset_ = getDataSet(dataset);
//    float dx = DataSets[dataset_].dx;
//    for (p = 0; p < count; p++)
//    {
//        // corner 000
//        cell_nodes[8*p][0] = dx*floor(datain[p][0]/dx);
//        cell_nodes[8*p][1] = dx*floor(datain[p][1]/dx);
//        cell_nodes[8*p][2] = dx*floor(datain[p][2]/dx);
//        // corner 001
//        cell_nodes[8*p + 1][0] = cell_nodes[8*p][0] + dx;
//        cell_nodes[8*p + 1][1] = cell_nodes[8*p][1]     ;
//        cell_nodes[8*p + 1][2] = cell_nodes[8*p][2]     ;
//        // corner 010
//        cell_nodes[8*p + 2][0] = cell_nodes[8*p][0]     ;
//        cell_nodes[8*p + 2][1] = cell_nodes[8*p][1] + dx;
//        cell_nodes[8*p + 2][2] = cell_nodes[8*p][2]     ;
//        // corner 011
//        cell_nodes[8*p + 3][0] = cell_nodes[8*p][0] + dx;
//        cell_nodes[8*p + 3][1] = cell_nodes[8*p][1] + dx;
//        cell_nodes[8*p + 3][2] = cell_nodes[8*p][2]     ;
//        // corner 100
//        cell_nodes[8*p + 4][0] = cell_nodes[8*p][0]     ;
//        cell_nodes[8*p + 4][1] = cell_nodes[8*p][1]     ;
//        cell_nodes[8*p + 4][2] = cell_nodes[8*p][2] + dx;
//        // corner 101
//        cell_nodes[8*p + 5][0] = cell_nodes[8*p][0] + dx;
//        cell_nodes[8*p + 5][1] = cell_nodes[8*p][1]     ;
//        cell_nodes[8*p + 5][2] = cell_nodes[8*p][2] + dx;
//        // corner 110
//        cell_nodes[8*p + 6][0] = cell_nodes[8*p][0]     ;
//        cell_nodes[8*p + 6][1] = cell_nodes[8*p][1] + dx;
//        cell_nodes[8*p + 6][2] = cell_nodes[8*p][2] + dx;
//        // corner 111
//        cell_nodes[8*p + 7][0] = cell_nodes[8*p][0] + dx;
//        cell_nodes[8*p + 7][1] = cell_nodes[8*p][1] + dx;
//        cell_nodes[8*p + 7][2] = cell_nodes[8*p][2] + dx;
//    }
//    getBoxFilter(
//            authToken,
//            dataset,
//            field,
//            time,
//            filterwidth,
//            8*count,
//            cell_nodes,
//            cellvalues);
//    for (p = 0; p < count; p++)
//    {
//        // get fractions
//        xx = (datain[p][0] - cell_nodes[8*p][0]) / dx;
//        yy = (datain[p][1] - cell_nodes[8*p][1]) / dx;
//        zz = (datain[p][2] - cell_nodes[8*p][2]) / dx;
//        // not most efficient way of writing the formula, but the most clear
//        dataout[p][0] = (((cellvalues[8*p  ][0]*(1-xx) + cellvalues[8*p+1][0]*xx)*(1-yy)
//                        + (cellvalues[8*p+2][0]*(1-xx) + cellvalues[8*p+3][0]*xx)* yy   )*(1-zz)
//                        +((cellvalues[8*p+4][0]*(1-xx) + cellvalues[8*p+5][0]*xx)*(1-yy)
//                        + (cellvalues[8*p+6][0]*(1-xx) + cellvalues[8*p+7][0]*xx)* yy   )* zz);
//        dataout[p][1] = (((cellvalues[8*p  ][1]*(1-xx) + cellvalues[8*p+1][1]*xx)*(1-yy)
//                        + (cellvalues[8*p+2][1]*(1-xx) + cellvalues[8*p+3][1]*xx)* yy   )*(1-zz)
//                        +((cellvalues[8*p+4][1]*(1-xx) + cellvalues[8*p+5][1]*xx)*(1-yy)
//                        + (cellvalues[8*p+6][1]*(1-xx) + cellvalues[8*p+7][1]*xx)* yy   )* zz);
//        dataout[p][2] = (((cellvalues[8*p  ][2]*(1-xx) + cellvalues[8*p+1][2]*xx)*(1-yy)
//                        + (cellvalues[8*p+2][2]*(1-xx) + cellvalues[8*p+3][2]*xx)* yy   )*(1-zz)
//                        +((cellvalues[8*p+4][2]*(1-xx) + cellvalues[8*p+5][2]*xx)*(1-yy)
//                        + (cellvalues[8*p+6][2]*(1-xx) + cellvalues[8*p+7][2]*xx)* yy   )* zz);
//    }
//    return 0;
//}

//int getFilteredPosition(
//        char *authToken,
//        char *dataset,
//        float startTime,
//        float endTime,
//        float dt,
//        float filterwidth,
//        int count,
//        float datain[][3],
//        float dataout[][3])
//{
//    float vel0[count][3], vel1[count][3];
//    float y[count][3];
//    float time, deltat;
//    int p, tcounter, nsteps;
//    nsteps = ceil((endTime - startTime) / dt);
//    dt = (endTime - startTime) / nsteps;
//    for (p = 0; p < count; p++)
//    {
//        dataout[p][0] = datain[p][0];
//        dataout[p][1] = datain[p][1];
//        dataout[p][2] = datain[p][2];
//    }
//    time = startTime;
//    for (tcounter = 0; tcounter < nsteps; tcounter++)
//    {
//        interpolateBoxFilter(
//                authToken,
//                dataset,
//                "velocity",
//                time,
//                filterwidth,
//                count,
//                dataout,
//                vel0);
//        for (p = 0; p < count; p++)
//        {
//            y[p][0] = dataout[p][0] + dt * vel0[p][0];
//            y[p][1] = dataout[p][1] + dt * vel0[p][1];
//            y[p][2] = dataout[p][2] + dt * vel0[p][2];
//        }
//        time += dt;
//        interpolateBoxFilter(
//                authToken,
//                dataset,
//                "velocity",
//                time,
//                filterwidth,
//                count,
//                y,
//                vel1);
//        for (p = 0; p < count; p++)
//        {
//            dataout[p][0] = dataout[p][0] + dt * .5 * (vel0[p][0] + vel1[p][0]);
//            dataout[p][1] = dataout[p][1] + dt * .5 * (vel0[p][1] + vel1[p][1]);
//            dataout[p][2] = dataout[p][2] + dt * .5 * (vel0[p][2] + vel1[p][2]);
//        }
//    }
//    return 0;
//}

