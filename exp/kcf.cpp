#include "ekf.h"

#define NX          (3*3+5)
#define NUM_AGENT   4

int main(int argc, char **argv)
{
    gtime_t t0 = {0}, ts = {0}, te = {0};
    int i;

    // TODO: search dir for files
    const char suffix_obs[] = ".obs", suffix_pos[] = ".pos", suffix_kcfpos[] = "_kcf.pos";
    const char *filepath [] = {"/home/philia/Documents/nav_data/uav_20211211/gnc01/data_20211211_151827/processed/",
    "/home/philia/Documents/nav_data/uav_20211211/gnc02/data_20211211_151836/processed/",
    "/home/philia/Documents/nav_data/uav_20211211/gnc04/data_20211211_151612/processed/",
    "/home/philia/Documents/nav_data/uav_20211211/gnc05/data_20211211_151900/processed/"};
    const char *navfile [] = {"eph_202112110707.nav", "eph_202112110710.nav", "eph_202112110714.nav", "eph_202112110714.nav"};
    const char *obsfile [] = {"ublox_20211211_15.obs", "ublox_20211211_15.obs", "ublox_20211211_15.obs", "ublox_20211211_15.obs"};
    
    char *navfile_[NUM_AGENT], *obsfile_[NUM_AGENT], *trufile[NUM_AGENT], *outfile[NUM_AGENT];
    for (i=0;i<NUM_AGENT;i++) navfile_[i] = charcat(filepath[i], navfile[i]);
    for (i=0;i<NUM_AGENT;i++) obsfile_[i] = charcat(filepath[i], obsfile[i]);
    for (i=0;i<NUM_AGENT;i++) trufile[i] = charcat(filepath[i], strpl(obsfile[i], suffix_obs, suffix_pos));
    for (i=0;i<NUM_AGENT;i++) outfile[i] = charcat(filepath[i], strpl(obsfile[i], suffix_obs, suffix_kcfpos));

    // time sync: compare time in obs, re-read obs by limiting ts, te

    FILE *fp[NUM_AGENT];
    for (i=0;i<NUM_AGENT;i++) fp[i]=openfile(outfile[i]);
    // infile


    // generate relative position measurement


    return 0;
}