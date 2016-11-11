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
 * Multi point mode emission routine 
 ***********************************************************************************/
int multiRoutine(config_t cfg, const char *filename, gsl_rng * r)
{
	/* Parameters for simulation */
	float x, y, z, prog;
	int countPSF, nPSF;
	FILE *fileIn;
	char *outname = (char *)malloc(30 * sizeof(char));
	char *molname = (char *)malloc(10 * sizeof(char));

	/* Open input files */
	fileIn = fopen(filename, "r");
	if (fileIn == NULL) {
		fprintf(stderr, "Error opening %s for reading.\n", filename);
		exit(1);
	}

	/* Get common parameters */
	struct commonParms cParms = parseCommon(cfg, fileIn);

	/* Get multi mode parameters */
	struct multiParms mParms = parseMulti(cfg);

	/* Vectors with PSF centers */
	countPSF = mParms.nPSFX * mParms.nPSFY;	// total number of PSFs
	double centerx[mParms.nPSFX], centery[mParms.nPSFY], center[countPSF][2];

	for (int i = 0; i < mParms.nPSFX; i++) {
		centerx[i] = i * mParms.dx - (mParms.nPSFX - 1) * mParms.dx / 2;
	}
	for (int i = 0; i < mParms.nPSFY; i++) {
		centery[i] = i * mParms.dy - (mParms.nPSFY - 1) * mParms.dy / 2;
	}

	nPSF = 0;
	FILE *fileIdx = fopen("index.txt", "w");
	for (int i = 0; i < mParms.nPSFX; i++) {
		for (int j = 0; j < mParms.nPSFY; j++) {
			center[nPSF][0] = centerx[i];
			center[nPSF][1] = centery[j];
			fprintf(fileIdx, "%03d \t %0.3f \t %0.3f\n", nPSF, center[nPSF][0], center[nPSF][1]);
			nPSF++;
		}
	}
	fclose(fileIdx);

	/* Opening output files and initial photon number set to zero */
	FILE *fileOut[countPSF][2];
	int nphot[countPSF][2];

	if (cParms.sChannel[0].status == 1) {
		for (nPSF = 0; nPSF < countPSF; nPSF++) {
			sprintf(outname, "%s_%03d_c0.txt", mParms.prefix, nPSF);
			fileOut[nPSF][0] = fopen(outname, "w");
			if (fileOut[nPSF][0] == NULL) {
				fprintf(stderr, "Error opening output file.\n");
				fclose(fileIn);
				exit(1);
			}
			nphot[nPSF][0] = 0;
		}
	}

	if (cParms.sChannel[1].status == 1) {
		for (nPSF = 0; nPSF < countPSF; nPSF++) {
			sprintf(outname, "%s_%03d_c1.txt", mParms.prefix, nPSF);
			fileOut[nPSF][1] = fopen(outname, "w");
			if (fileOut[nPSF][1] == NULL) {
				fprintf(stderr, "Error opening output file.\n");
				fclose(fileIn);
				exit(1);
			}
			nphot[nPSF][1] = 0;
		}
	}

	/* Info about files */
	printLogo();
	printf("\n");
	printf("%s %s starting in multi mode.\n", PROGNAME, VERSION);
	printf("  Reading input file %s\n", filename);
	if (cParms.sChannel[0].status == 1) {
		printf("  Writing %d output files for channel 0\n", countPSF);
	}
	if (cParms.sChannel[1].status == 1) {
		printf("  Writing %d output files for channel 1\n", countPSF);
	}
	printf("\n");

	/* Print recovered parameter values */
	printf("Parameters recovered from config file:\n");
	printf("  Waist in XY plane (w_xy): %g um\n", cParms.w_xy);
	printf("  Waist in Z plane (w_z): %g um\n", cParms.w_z);
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

	/* Photon emission routine */
	while (fscanf(fileIn, "%s %f %f %f\n", molname, &x, &y, &z) != EOF) {
		if (x == 100) {
			prog = 100 * (y / z);
			printf("Progress: %.1f%%\r", prog);
			for (nPSF = 0; nPSF < countPSF; nPSF++) {
				if (cParms.sChannel[0].status == 1) {
					nphot[nPSF][0] += noiseGenerator(nphot[nPSF][0], cParms.noise, r);
					fprintf(fileOut[nPSF][0], "%d\n", nphot[nPSF][0]);
					nphot[nPSF][0] = 0;
				}
				if (cParms.sChannel[1].status == 1) {
					nphot[nPSF][1] += noiseGenerator(nphot[nPSF][1], cParms.noise, r);
					fprintf(fileOut[nPSF][1], "%d\n", nphot[nPSF][1]);
					nphot[nPSF][1] = 0;
				}
			}
		} else {
			for (nPSF = 0; nPSF < countPSF; nPSF++) {
				if (cParms.sChannel[0].status == 1) {
					for (int i = 0; i < cParms.sChannel[0].nmols; i++) {
						if (!strcmp(molname, cParms.sChannel[0].mols[i])) {
							nphot[nPSF][0] +=
							    gaussPSF(x, y, z,
								     cParms.w_xy,
								     cParms.w_z,
								     center
								     [nPSF][0],
								     center
								     [nPSF][1],
								     mParms.centerz,
								     cParms.nevents, cParms.sChannel[0].q[i], r);
						}
					}
				}
				if (cParms.sChannel[1].status == 1) {
					for (int i = 0; i < cParms.sChannel[1].nmols; i++) {
						if (!strcmp(molname, cParms.sChannel[1].mols[i])) {
							nphot[nPSF][1] +=
							    gaussPSF(x, y, z,
								     cParms.w_xy,
								     cParms.w_z,
								     center
								     [nPSF][0],
								     center
								     [nPSF][1],
								     mParms.centerz,
								     cParms.nevents, cParms.sChannel[1].q[i], r);
						}
					}
				}
			}
		}
	}

	printf("\n");
	/* Close and destroy file pointers */
	fclose(fileIn);
	for (nPSF = 0; nPSF < countPSF; nPSF++) {
		if (cParms.sChannel[0].status == 1) {
			fclose(fileOut[nPSF][0]);
		}
		if (cParms.sChannel[1].status == 1) {
			fclose(fileOut[nPSF][1]);
		}
	}

	return 0;
}
