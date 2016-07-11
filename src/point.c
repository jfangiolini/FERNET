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

#include "fernet.h"

/***********************************************************************************
 * Point mode emission routine
 ***********************************************************************************/
int pointRoutine(config_t cfg, const char *filename)
{
    /* Parameters for simulation */
    float x, y, z, prog;
    int nphot[] = { 0, 0 };
    FILE *fileIn, *fileOut[2];
    char *outname[] =
	{ (char *) malloc(30 * sizeof(char)),
    (char *) malloc(30 * sizeof(char)) };
    char *molname = (char *) malloc(10 * sizeof(char));

    /* Open input file */
    fileIn = fopen(filename, "r");
    if (fileIn == NULL) {
	fprintf(stderr, "Error opening %s for reading.\n", filename);
	exit(1);
    }

    /* Get common parameters */
    struct commonParms cParms = parseCommon(cfg, fileIn);

    /* Get point mode parameters */
    struct pointParms pParms = parsePoint(cfg);

    /* Open output files */
    if (cParms.sChannel[0].status == 1) {
	sprintf(outname[0], "%s_c0.txt", pParms.prefix);
	fileOut[0] = fopen(outname[0], "w");
	if (fileOut[0] == NULL) {
	    fprintf(stderr, "Error opening %s for writing.\n", outname[0]);
	    fclose(fileIn);
	    exit(1);
	}
    }

    if (cParms.sChannel[1].status == 1) {
	sprintf(outname[1], "%s_c1.txt", pParms.prefix);
	fileOut[1] = fopen(outname[1], "w");
	if (fileOut[1] == NULL) {
	    fprintf(stderr, "Error opening %s for writing.\n", outname[1]);
	    fclose(fileIn);
	    exit(1);
	}
    }

    /* Info about files */
    printLogo();
    printf("\n");
    printf("%s %s starting in point mode.\n", PROGNAME, VERSION);
    printf("  Reading input file %s\n", filename);
    if (cParms.sChannel[0].status == 1) {
	printf("  Writing output file %s_c0.txt for channel 0\n",
	       pParms.prefix);
    }
    if (cParms.sChannel[1].status == 1) {
	printf("  Writing output file %s_c1.txt for channel 1\n",
	       pParms.prefix);
    }
    printf("\n");

    /* Print recovered parameter values */
    printf("Parameters recovered from config file:\n");
    printf("  Waist in XY plane (w_xy): %g um\n", cParms.w_xy);
    printf("  Waist in Z plane (w_z): %g um\n", cParms.w_z);
    printf("  Center of PSF [X Y Z]: [%0.2f %0.2f %0.2f] um\n",
	   pParms.centerx, pParms.centery, pParms.centerz);

    if (cParms.sChannel[0].status == 1) {
	printf("  Molecules emitting in channel 0: [ ");
	for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
	    printf("%s ", cParms.sChannel[0].mols[i]);
	}
	printf("]\n");
	printf
	    ("  Emission probability (q) for molecules in channel 0: [ ");
	for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
	    printf("%g ", cParms.sChannel[0].q[i]);
	}
	printf("]\n");
    }

    if (cParms.sChannel[1].status == 1) {
	printf("  Molecules emitting in channel 1: [ ");
	for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
	    printf("%s ", cParms.sChannel[1].mols[i]);
	}
	printf("]\n");
	printf
	    ("  Emission probability (q) for molecules in channel 1: [ ");
	for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
	    printf("%g ", cParms.sChannel[1].q[i]);
	}
	printf("]\n");
    }
    printf("\n");

    /* Photon emission routine */

    while (fscanf(fileIn, "%s %f %f %f\n", molname, &x, &y, &z) != EOF) {
	if (x == 100) {
	    /*Time step separator */
	    /* Restart the number of processed molecule position */
	    prog = 100 * (y / z);
	    printf("Progress: %.1f%%\r", prog);

	    if (cParms.sChannel[0].status == 1) {
		fprintf(fileOut[0], "%d\n", nphot[0]);
		nphot[0] = 0;
	    }
	    if (cParms.sChannel[1].status == 1) {
		fprintf(fileOut[1], "%d\n", nphot[1]);
		nphot[1] = 0;
	    }
	} else {
	    if (cParms.sChannel[0].status == 1) {
		for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
		    if (!strcmp(molname, cParms.sChannel[0].mols[i])) {
			nphot[0] +=
			    gaussPSF(x, y, z, cParms.w_xy, cParms.w_z,
				     pParms.centerx, pParms.centery,
				     pParms.centerz, cParms.nevents,
				     cParms.sChannel[0].q[i]);
		    }
		}
	    }

	    if (cParms.sChannel[1].status == 1) {
		for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
		    if (!strcmp(molname, cParms.sChannel[1].mols[i])) {
			nphot[1] +=
			    gaussPSF(x, y, z, cParms.w_xy, cParms.w_z,
				     pParms.centerx, pParms.centery,
				     pParms.centerz, cParms.nevents,
				     cParms.sChannel[1].q[i]);
		    }
		}
	    }
	}
    }
    printf("\n");

    /* Closing all pointers and cleaning up */
    fclose(fileIn);
    if (cParms.sChannel[0].status == 1) {
	fclose(fileOut[0]);
    }
    if (cParms.sChannel[1].status == 1) {
	fclose(fileOut[1]);
    }

    return 0;
}
