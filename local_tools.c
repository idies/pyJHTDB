/*  
 *  These subroutines are meant to be called from a pyJHTDB.libTDB object.
 *  You're on your own otherwise.
 *
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "turblib.h"

// relative error of single precision numbers
float float_error = 1e-6;

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
            x += 1;
            if (((*x)[0] < xmin || (*x)[0] > xmax)
             || ((*x)[1] < ymin || (*x)[1] > ymax)
             || ((*x)[2] < zmin || (*x)[2] > zmax))
                break;
        }
        traj_length[p] = s+1;
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
    float bsize;
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
            x += 1;
            if (sqrt((x[0][0] - ox)*(x[0][0] - ox)
                   + (x[0][1] - oy)*(x[0][1] - oy)
                   + (x[0][2] - oz)*(x[0][2] - oz)) > radius)
                break;
        }
        traj_length[p] = s+1;
    }
    return 0;
}

