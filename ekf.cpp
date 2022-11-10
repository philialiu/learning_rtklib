/*------------------------------------------------------------------------------
* ekf.c : spp ekf
*
*          Copyright (C) 2022 by Philia Liu, All rights reserved.
*
* reference :
*     [1] 
*
*-----------------------------------------------------------------------------*/
#include <iostream>
#include "ekf.h"
#include <rtklib.h>
#include <Eigen/LU>

static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (; *i < obs->n; (*i)++)
        if (obs->data[*i].rcv == rcv)
            break;
    for (n = 0; *i + n < obs->n; n++)
    {
        /* time difference (t1-t2) (s) */
        tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
        if (obs->data[*i + n].rcv != rcv || tt > DTTOL)
            break;
    }
    return n;
}

char* charcat(const char* str1, const char* str2)
{
    char* str3 = new char[1024];
    strcpy(str3, str1);
    strcat(str3, str2);
    return str3;
}

int main()
{
    gtime_t t0 = {0}, ts = {0}, te = {0};

    char filepath[] = "/home/philia/Documents/nav_data/uav_20211211/gnc01/data_20211211_150700/processed/";
    char file1[] = "eph_202112110655.nav";
    char file2[] = "ublox_20211211_14.obs";

    obs_t obs = {0};
    nav_t nav = {0};
    sta_t sta = {""};

    int stat1 = readrnxt(charcat(filepath, file1), 1, ts, te, 0.0, "", &obs, &nav, &sta);
    int stat2 = readrnxt(charcat(filepath, file2), 1, t0, t0, 0.0, "", &obs, &nav, &sta);

    if (stat1 && stat2)
    {
        std::cout << "Files are successfully opened!" << std::endl;
    }
    else
    {
        std::cout << "Fail to open files." << std::endl;
    }

    prcopt_t opt = prcopt_default;
    
    // TODO: restore wave length according to ?

    for (int i = 0; i < NSATGPS; i++)
    {
        nav.lam[i][0] = CLIGHT / FREQ1;
        nav.lam[i][1] = CLIGHT / FREQ2;
        nav.lam[i][2] = CLIGHT / FREQ5;
    }

    for (int j = NSATGPS; j < NSATGLO + NSATGPS; j++)
    {
        nav.lam[j][0] = CLIGHT / FREQ1_GLO;
        nav.lam[j][1] = CLIGHT / FREQ2_GLO;
        nav.lam[j][2] = CLIGHT / FREQ3_GLO;
    }

    // first epoch: spp
    // predict
    // residual
    // update

    return 0;
}