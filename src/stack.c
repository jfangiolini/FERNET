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

/**********************************************************************************
* Stack mode emission routine
***********************************************************************************/
int stackRoutine(config_t cfg, const char *filename)
{
    /* Parameters for simulation */
    float x, y, z, prog;
    int nphot[] = { 0, 0 };
    int column = 0, row = 0, slice = 0;
    FILE *fileIn;
    TIFF *tif[2];
    char *outname[] =
	{ (char *) malloc(30 * sizeof(char)),
    (char *) malloc(30 * sizeof(char)) };
    char *molname = { (char *) malloc(10 * sizeof(char)) };

    /* Opening input file */
    fileIn = fopen(filename, "r");
    if (fileIn == NULL) {
	fprintf(stderr, "Error opening %s for reading.\n", filename);
	exit(1);
    }

    /* Get common parameters */
    struct commonParms cParms = parseCommon(cfg, fileIn);

    /* Get stack mode parameters */
    struct stackParms sParms = parseStack(cfg);

    /* Open output files */
    if (cParms.sChannel[0].status == 1) {
	sprintf(outname[0], "%s_c0.tif", sParms.tiffname);
	tif[0] = TIFFOpen(outname[0], "w");
	if (tif[0] == NULL) {
	    fprintf(stderr, "Error opening %s for writing.\n", outname[0]);
	    fclose(fileIn);
	    exit(1);
	}
    }

    if (cParms.sChannel[1].status == 1) {
	sprintf(outname[1], "%s_c1.tif", sParms.tiffname);
	tif[1] = TIFFOpen(outname[1], "w");
	if (tif[1] == NULL) {
	    fprintf(stderr, "Error opening %s for writing.\n", outname[1]);
	    fclose(fileIn);
	    exit(1);
	}
    }

    /* Write TIFF tags */
    if (cParms.sChannel[0].status == 1) {
	writeImageTIFFtags(tif[0], sParms.width, sParms.height);
    }
    if (cParms.sChannel[1].status == 1) {
	writeImageTIFFtags(tif[1], sParms.width, sParms.height);
    }

    /* Vector with center for each pixel */
    double centerx[sParms.width], centery[sParms.height];

    for (int i = 0; i < sParms.width; i++) {
	centerx[i] =
	    i * sParms.pixel - (sParms.width - 1) * sParms.pixel / 2;
    }

    for (int i = 0; i < sParms.height; i++) {
	centery[i] =
	    i * sParms.pixel - (sParms.height - 1) * sParms.pixel / 2;
    }

    /* Calculate "dummy" line length */
    int ndummy = sParms.width + round(sParms.deadtime / cParms.simu_dt);

    /* Calculate Z slices positions */

    int nslices = round((sParms.top_z - sParms.bot_z) / sParms.step);
    double zpos[nslices];

    int Niters = ndummy * sParms.height * nslices;

    for (int i = 0; i < nslices; i++) {
	zpos[i] = sParms.top_z - i * sParms.step;
    }

    /* Info about files */
    printLogo();
    printf("\n");
    printf("%s %s starting in stack mode.\n", PROGNAME, VERSION);
    printf("  Reading input file %s\n", filename);
    if (cParms.sChannel[0].status == 1) {
	printf("  Writing output file %s for channel 0\n", outname[0]);
    }
    if (cParms.sChannel[1].status == 1) {
	printf("  Writing output file %s for channel 1\n", outname[1]);
    }
    printf("\n");

    /* Print recovered parameters from config file */
    printf("Parameters recovered from config file: \n");
    printf("  Waist in XY plane (w_xy): %g um\n", cParms.w_xy);
    printf("  Waist in Z plane (w_z): %g um\n", cParms.w_z);
    printf("  Image width: %0.2f um\n", sParms.width * sParms.pixel);
    printf("  Image height: %0.2f um\n", sParms.height * sParms.pixel);
    printf("  Pixel size: %0.3f um\n", sParms.pixel);
    printf("  Line time: %0.2f ms\n", 1000 * ndummy * cParms.simu_dt);
    printf("  Frame time: %0.2f s\n",
	   ndummy * sParms.height * cParms.simu_dt);
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
    char buf_row[2][sParms.width];
    while ((fscanf(fileIn, "%s %f %f %f\n", molname, &x, &y, &z) != EOF)
	   && y < Niters) {
	if (x == 100) {
	    prog = 100 * (y / Niters);
	    printf("Progress: %0.1f%%\r", prog);
	    if ((int) y % ndummy < sParms.width) {
		if (cParms.sChannel[0].status == 1) {
		    buf_row[0][column] = nphot[0];
		    nphot[0] = 0;
		}
		if (cParms.sChannel[1].status == 1) {
		    buf_row[1][column] = nphot[1];
		    nphot[1] = 0;
		}

		if (column == sParms.width - 1) {
		    if (cParms.sChannel[0].status == 1) {
			TIFFWriteScanline(tif[0], buf_row[0], row, 0);
		    }
		    if (cParms.sChannel[1].status == 1) {
			TIFFWriteScanline(tif[1], buf_row[1], row, 0);
		    }
		    column = 0;

		    if (row == sParms.height - 1) {
			if (cParms.sChannel[0].status == 1) {
			    TIFFWriteDirectory(tif[0]);
			    writeImageTIFFtags(tif[0], sParms.width,
					       sParms.height);
			}
			if (cParms.sChannel[1].status == 1) {
			    TIFFWriteDirectory(tif[1]);
			    writeImageTIFFtags(tif[1], sParms.width,
					       sParms.height);
			}
			row = 0;
			slice++;
		    } else {
			row++;
		    }
		} else {
		    column++;
		}
	    } else {
		if (cParms.sChannel[0].status == 1) {
		    nphot[0] = 0;
		}
		if (cParms.sChannel[1].status == 1) {
		    nphot[1] = 0;
		}
	    }
	} else {
	    if (cParms.sChannel[0].status == 1) {
		for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
		    if (!strcmp(molname, cParms.sChannel[0].mols[i])) {
			nphot[0] +=
			    gaussPSF(x, y, z, cParms.w_xy, cParms.w_z,
				     centerx[column], centery[row],
				     zpos[slice], cParms.nevents,
				     cParms.sChannel[0].q[i]);
		    }
		}
	    }

	    if (cParms.sChannel[1].status == 1) {
		for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
		    if (!strcmp(molname, cParms.sChannel[1].mols[i])) {
			nphot[1] +=
			    gaussPSF(x, y, z, cParms.w_xy, cParms.w_z,
				     centerx[column], centery[row],
				     zpos[slice], cParms.nevents,
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
