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
 * Orbital scanning mode emission routine
 ***********************************************************************************/
int orbitRoutine(config_t cfg, const char *filename)
{
    /* Parameters for simulation */
    float x, y, z, prog;
    int nphot[] = { 0, 0 };
    int pixel = 0, row = 0;
    FILE *fileIn;
    TIFF *tif[2];
    char *outname[] =
	{ (char *) malloc(30 * sizeof(char)),
    (char *) malloc(30 * sizeof(char)) };
    char *molname = (char *) malloc(10 * sizeof(char));

    /* Opening input file */
    fileIn = fopen(filename, "r");
    if (fileIn == NULL) {
	fprintf(stderr, "Error opening %s for reading.\n", filename);
	exit(1);
    }

    /* Get common parameters */
    struct commonParms cParms = parseCommon(cfg, fileIn);

    /* Get orbital scanning mode parameters */
    struct orbitParms orParms = parseOrbit(cfg);

    /* Orbit calculation */
    int n_pixels = round(orParms.period / cParms.simu_dt);
    double dtheta = 2 * M_PI / n_pixels;
    double x_o[n_pixels];
    double y_o[n_pixels];
    for (int i = 0; i < n_pixels; i++) {
	x_o[i] = orParms.radius * cos(i * dtheta);
	y_o[i] = orParms.radius * sin(i * dtheta);
    }

    /* Open output files */
    if (cParms.sChannel[0].status == 1) {
	sprintf(outname[0], "%s_c0.tif", orParms.tiffname);
	tif[0] = TIFFOpen(outname[0], "w");
	if (tif[0] == NULL) {
	    fprintf(stderr, "Error opening %s for writing.\n", outname[0]);
	    fclose(fileIn);
	    exit(1);
	}
    }

    if (cParms.sChannel[1].status == 1) {
	sprintf(outname[1], "%s_c1.tif", orParms.tiffname);
	tif[1] = TIFFOpen(outname[1], "w");
	if (tif[1] == NULL) {
	    fprintf(stderr, "Error opening %s for writing.\n", outname[1]);
	    fclose(fileIn);
	    exit(1);
	}
    }

    /* TIFF tags */
    if (cParms.sChannel[0].status == 1) {
	writeLineTIFFTags(tif[0], n_pixels);
    }
    if (cParms.sChannel[1].status == 1) {
	writeLineTIFFTags(tif[1], n_pixels);
    }

    /* Info about files */
    printLogo();
    printf("\n");
    printf("%s %s starting in orbital scanning mode.\n", PROGNAME,
	   VERSION);
    printf("  Reading input file %s\n", filename);
    if (cParms.sChannel[0].status == 1) {
	printf("  Writing output file %s for channel 0\n", outname[0]);
    }
    if (cParms.sChannel[1].status == 1) {
	printf("  Writing output file %s for channel 1\n", outname[1]);
    }
    printf("\n");

    /* Print parameters recovered from config file */
    printf("Parameters recovered from config file: \n");
    printf("  Waist in XY plane (w_xy): %g um\n", cParms.w_xy);
    printf("  Waist in Z plane (w_z): %g um\n", cParms.w_z);
    printf("  Center of orbit [X Y Z]: [%0.2f %0.2f %0.2f] um\n",
	   orParms.centerx, orParms.centery, orParms.centerz);
    printf("  Number of pixels along orbit: %d\n", n_pixels);
    printf("  Orbit radius: %0.2f um\n", orParms.radius);
    printf("  Orbit period: %0.2f ms\n", orParms.period * 1000);
    printf("\n");

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
    char buf_row[2][n_pixels];	// buffer for TIFF writing

    while (fscanf(fileIn, "%s %f %f %f\n", molname, &x, &y, &z) != EOF) {
	if (x == 100) {
	    prog = 100 * (y / z);
	    printf("Progress: %.1f%%\r", prog);
	    if (((int) y + 1) % n_pixels == 0) {
		if (cParms.sChannel[0].status == 1) {
		    TIFFWriteScanline(tif[0], buf_row[0], row, 0);
		}
		if (cParms.sChannel[1].status == 1) {
		    TIFFWriteScanline(tif[1], buf_row[1], row, 0);
		}
		pixel = 0;
		row++;
	    } else {
		if (cParms.sChannel[0].status == 1) {
		    buf_row[0][pixel] = nphot[0];
		    nphot[0] = 0;
		}
		if (cParms.sChannel[1].status == 1) {
		    buf_row[1][pixel] = nphot[1];
		    nphot[1] = 0;
		}
		pixel++;
	    }
	} else {
	    if (cParms.sChannel[0].status == 1) {
		for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
		    if (!strcmp(molname, cParms.sChannel[0].mols[i])) {
			nphot[0] +=
			    gaussPSF(x, y, z, cParms.w_xy, cParms.w_z,
				     x_o[pixel] - orParms.centerx,
				     y_o[pixel] - orParms.centery,
				     orParms.centerz, cParms.nevents,
				     cParms.sChannel[0].q[i]);
		    }
		}
	    }

	    if (cParms.sChannel[1].status == 1) {
		for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
		    if (!strcmp(molname, cParms.sChannel[1].mols[i])) {
			nphot[1] +=
			    gaussPSF(x, y, z, cParms.w_xy, cParms.w_z,
				     x_o[pixel] - orParms.centerx,
				     y_o[pixel] - orParms.centery,
				     orParms.centerz, cParms.nevents,
				     cParms.sChannel[1].q[i]);
		    }
		}
	    }
	}
    }
    printf("\n");

    /* Closing files */
    fclose(fileIn);

    if (cParms.sChannel[0].status == 1) {
	TIFFClose(tif[0]);
    }
    if (cParms.sChannel[1].status == 1) {
	TIFFClose(tif[1]);
    }

    return 0;
}
