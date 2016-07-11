/***********************************************************************************
 *                                                                                 *
 * FERNET: Fluorescence Emission Recipes and NumErical routines Toolkit            *
 * Developed by Juan F. Angiolini and Esteban Mocskos                              *
 * Facultad de Ciencias Exactas y Naturales, UBA                                   *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

 /***********************************************************************************
 * Required libraries for compiling
 ***********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <argtable2.h>
#include <libconfig.h>
#include <tiffio.h>
#include <math.h>
#include <time.h>
#include <string.h>

/***********************************************************************************
 * Program info
 ***********************************************************************************/

#define PROGNAME "FERNET"
#define VERSION "1.3"
#define DATE "July 2016"

/***********************************************************************************
 * Function protoypes
 ***********************************************************************************/

int pointRoutine(config_t, const char *);	// Point mode emission routine
int multiRoutine(config_t, const char *);	// Multi point mode emission routine
int lineRoutine(config_t, const char *);	// Linescan mode emission routine
int imageRoutine(config_t, const char *);	// Raster mode emission routine
int stackRoutine(config_t, const char *);	// 3D stack emission routine
int spimRoutine(config_t, const char *);	// SPIM emission routine
int orbitRoutine(config_t, const char *);	// Orbital scanning emission routine
int gaussPSF(double, double, double, double, double, double, double,
	     double, int, double);
int spimPSF(double, double, double, int, double);
void writeLineTIFFTags(TIFF *, int);
void writeImageTIFFtags(TIFF *, int, int);	// Write TIFFs tags
struct args parseArgs(int, char **);	// Parse arguments from console
struct commonParms parseCommon(config_t, FILE *);	// Parse common parameters from config file
struct pointParms parsePoint(config_t);	// Parse point mode parameters from config file
struct multiParms parseMulti(config_t);	// Parse multi point mode parameters from config file
struct lineParms parseLine(config_t);	// Parse linescan mode parameters from config file
struct imageParms parseImage(config_t);	// Parse raster mode parameters from config file
struct channelInfo parseChannel(config_setting_t *, float);	// Parse channel configuration block from config file
struct stackParms parseStack(config_t);	// Parse stack mode parameters from config file
struct spimParms parseSpim(config_t);	// Parse SPIM mode parameters from config file
struct orbitParms parseOrbit(config_t);	// Parse orbital scanning parameters from config file
void parseError(char *);	// Error log when parsing variables
double randn();
void printLogo();		// Print ASCII LOGO

/***********************************************************************************
 * Structures
 ***********************************************************************************/

enum fluo_modes {		// Possible emission modes 
    POINT,
    MULTI,
    LINE,
    IMAGE,
    STACK,
    SPIM,
    ORBIT
};

struct channelInfo {		// Channel information
    int status;
    int nmols;
    int *bright;
    double *q;
    char **mols;
};

struct commonParms {		// Common parameters
    double kappa;
    double w_xy, w_z;
    float simu_dt, mD;
    const char *PSFtype;
    int lambda;
    int val;
    int nevents;
    struct channelInfo sChannel[2];
};

struct pointParms {		// Point mode parameters
    double centerx;
    double centery;
    double centerz;
    const char *prefix;
};

struct multiParms {		// Multi point mode parametes
    const char *prefix;
    double centerz;
    int nPSFX, nPSFY;
    double dx, dy;
};

struct lineParms {		// Linescan mode parameters
    double centerx, centery, centerz;
    int ncolumn;
    double shift, deadtime;
    const char *tiffname;
};

struct imageParms {		// Raster mode parameters
    double pixel, centerz, deadtime;
    int width, height;
    const char *tiffname;
};

struct stackParms {		// Stack mode parameters
    double pixel, deadtime, top_z, bot_z, step;
    int width, height;
    const char *tiffname;
};

struct spimParms {		// SPIM mode parameters
    double pixel, frame_t, centerz, NA, waist, lambda;
    int width, height;
    const char *tiffname;
};

struct orbitParms {		// Orbital scanning parameters
    double centerx, centery, centerz;
    double radius, period;
    const char *tiffname;
};

struct args {			// Console arguments
    const char *filename;
    const char *mode;
    config_t cfg;
};