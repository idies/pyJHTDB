#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "turblib.h"

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
        // get bhat1 and compute xnext
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


