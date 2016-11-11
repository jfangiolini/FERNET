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
 * Linescan mode emission routine
 ***********************************************************************************/
int lineRoutine(config_t cfg, const char *filename, gsl_rng * r)
{
	/* Parameters for simulation */
	float x, y, z, prog;
	int nphot[] = { 0, 0 };
	int column = 0, row = 0;
	FILE *fileIn;
	TIFF *tif[2];
	char *outname[] = { (char *)malloc(30 * sizeof(char)),
		(char *)malloc(30 * sizeof(char))
	};
	char *molname = (char *)malloc(10 * sizeof(char));

	/* Opening input file */
	fileIn = fopen(filename, "r");
	if (fileIn == NULL) {
		fprintf(stderr, "Error opening %s for reading.\n", filename);
		exit(1);
	}

	/* Get common parameters */
	struct commonParms cParms = parseCommon(cfg, fileIn);

	/* Get line mode parameters */
	struct lineParms lParms = parseLine(cfg);

	/* Open output files */
	if (cParms.sChannel[0].status == 1) {
		sprintf(outname[0], "%s_c0.tif", lParms.tiffname);
		tif[0] = TIFFOpen(outname[0], "w");
		if (tif[0] == NULL) {
			fprintf(stderr, "Error opening %s for writing.\n", outname[0]);
			fclose(fileIn);
			exit(1);
		}
	}

	if (cParms.sChannel[1].status == 1) {
		sprintf(outname[1], "%s_c1.tif", lParms.tiffname);
		tif[1] = TIFFOpen(outname[1], "w");
		if (tif[1] == NULL) {
			fprintf(stderr, "Error opening %s for writing.\n", outname[1]);
			fclose(fileIn);
			exit(1);
		}
	}

	/* TIFF tags */
	if (cParms.sChannel[0].status == 1) {
		writeLineTIFFTags(tif[0], lParms.ncolumn);
	}
	if (cParms.sChannel[1].status == 1) {
		writeLineTIFFTags(tif[1], lParms.ncolumn);
	}

	/* Info about files */
	printLogo();
	printf("\n");
	printf("%s %s starting in line mode.\n", PROGNAME, VERSION);
	printf("  Reading input file %s\n", filename);
	if (cParms.sChannel[0].status == 1) {
		printf("  Writing output file %s for channel 0\n", outname[0]);
	}
	if (cParms.sChannel[1].status == 1) {
		printf("  Writing output file %s for channel 1\n", outname[1]);
	}
	printf("\n");

	/* Calculate "dummy" line length */
	int ndummy = lParms.ncolumn + round(lParms.deadtime / cParms.simu_dt);

	/* Print parameters recovered from config file */
	printf("Parameters recovered from config file: \n");
	printf("  Waist in XY plane (w_xy): %g um\n", cParms.w_xy);
	printf("  Waist in Z plane (w_z): %g um\n", cParms.w_z);
	printf("  Center of line [X Y Z]: [%0.2f %0.2f %0.2f] um\n", lParms.centerx, lParms.centery, lParms.centerz);
	printf("  Number of columns: %d\n", lParms.ncolumn);
	printf("  Line time: %0.2f ms\n", 1000 * ndummy * cParms.simu_dt);
	printf("  Line length: %0.2f um\n", lParms.ncolumn * lParms.shift);
	printf("\n");

	if (cParms.sChannel[0].status == 1) {
		printf("  Molecules emitting in channel 0: [ ");
		for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
			printf("%s ", cParms.sChannel[0].mols[i]);
		}
		printf("]\n");
		printf("  Emission probability (q) for molecules in channel 0: [ ");
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
		printf("  Emission probability (q) for molecules in channel 1: [ ");
		for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
			printf("%g ", cParms.sChannel[1].q[i]);
		}
		printf("]\n");
	}
	printf("\n");

	/* Calculate pixel position of every column */
	float centros[lParms.ncolumn];
	for (int i = 0; i < lParms.ncolumn; i++) {
		centros[i] = i * lParms.shift - (lParms.ncolumn - 1) * lParms.shift / 2;
	}

	/* Photon emission routine */
	char buf_row[2][lParms.ncolumn];	// buffer for TIFF writing

	while (fscanf(fileIn, "%s %f %f %f\n", molname, &x, &y, &z) != EOF) {
		if (x == 100) {
			prog = 100 * (y / z);
			printf("Progress: %.1f%%\r", prog);
			if (((int)y % ndummy) < lParms.ncolumn) {
				if (cParms.sChannel[0].status == 1) {
					nphot[0] += noiseGenerator(nphot[0], cParms.noise, r);
					buf_row[0][column] = nphot[0];
					nphot[0] = 0;
				}
				if (cParms.sChannel[1].status == 1) {
					nphot[1] += noiseGenerator(nphot[1], cParms.noise, r);
					buf_row[1][column] = nphot[1];
					nphot[1] = 0;
				}

				if (column == lParms.ncolumn - 1) {
					if (cParms.sChannel[0].status == 1) {
						TIFFWriteScanline(tif[0], buf_row[0], row, 0);
					}
					if (cParms.sChannel[1].status == 1) {
						TIFFWriteScanline(tif[1], buf_row[1], row, 0);
					}
					column = 0;
					row++;
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
						    gaussPSF(x, y, z,
							     cParms.w_xy,
							     cParms.w_z,
							     centros[column] -
							     lParms.centerx,
							     lParms.centery,
							     lParms.centerz,
							     cParms.nevents, cParms.sChannel[0].q[i], r);
					}
				}
			}

			if (cParms.sChannel[1].status == 1) {
				for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
					if (!strcmp(molname, cParms.sChannel[1].mols[i])) {
						nphot[1] +=
						    gaussPSF(x, y, z,
							     cParms.w_xy,
							     cParms.w_z,
							     centros[column] -
							     lParms.centerx,
							     lParms.centery,
							     lParms.centerz,
							     cParms.nevents, cParms.sChannel[1].q[i], r);
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

/***********************************************************************************
 * Function to write Line TIFF tags
 ***********************************************************************************/
void writeLineTIFFTags(TIFF * tif, int width)
{
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
	//TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
	//TIFFSetField(tif, TIFFTAG_EXIFIFD,8);
	//TIFFSetField(tif, TIFFTAG_EXTRASAMPLES,     0, NULL);
	//TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 16);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
	//TIFFSetField(tif, TIFFTAG_COMPRESSION,1);
	TIFFSetField(tif, TIFFTAG_SUBFILETYPE, 0);
	TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, "carpet");
	TIFFCheckpointDirectory(tif);
}
