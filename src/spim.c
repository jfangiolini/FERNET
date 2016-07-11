/**********************************************************************************
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
* Spim mode emission routine
***********************************************************************************/

int spimRoutine(config_t cfg, const char *filename)
{
    /* Parameters for simulation */
    float x, y, z, prog;
    FILE *fileIn;
    TIFF *tif;
    char *outname = (char *) malloc(30 * sizeof(char));
    char *molname = (char *) malloc(10 * sizeof(char));

    /* Opening input file */
    fileIn = fopen(filename, "r");
    if (fileIn == NULL) {
	fprintf(stderr, "Error opening %s for reading.\n", filename);
	exit(1);
    }
    //fprintf(stdout, "Check open input\n");

    /* Get common parameters */
    struct commonParms cParms = parseCommon(cfg, fileIn);

    /* Get spim mode parameters */
    struct spimParms spParms = parseSpim(cfg);
    //fprintf(stdout, "Check parse spim\n");

    /* Alloc memory for big files */
    //TIFF *tif = (TIFF*)_TIFFmalloc(spParms.height*spParms.width*10000*sizeof(char));

    /* Open output files */
    sprintf(outname, "%s.tif", spParms.tiffname);
    tif = TIFFOpen(outname, "w");
    if (tif == NULL) {
	fprintf(stderr, "Error opening %s for writing.\n", outname);
	fclose(fileIn);
	exit(1);
    }

    /* Write TIFF tags */
    writeImageTIFFtags(tif, spParms.width, spParms.height);
    //fprintf(stdout, "Check write tiff tags\n");

    /* CCD array allocation */
    char **CCD_buf = (char **) malloc(spParms.height * sizeof(char *));
    for (int i = 0; i < spParms.height; i++) {
	CCD_buf[i] = (char *) malloc(spParms.width * sizeof(char));
    }
    //fprintf(stdout, "Check ccd allocation\n");

    int nbin = round(spParms.frame_t / cParms.simu_dt);

    /* CCD array initialization */
    for (int i = 0; i < spParms.height; i++) {
	for (int j = 0; j < spParms.width; j++) {
	    CCD_buf[i][j] = 0;
	}
    }
    //fprintf(stdout, "Check ccd init\n");
    //fprintf(stdout, "%f\n",cParms.sChannel[0].q[0]);

    /* Position jitter */
    double R = 0.61 * (spParms.lambda / 1000) / (2 * spParms.NA);

    /* CCD bounds */
    double lx = (spParms.width * spParms.pixel) / 2;
    double ly = (spParms.height * spParms.pixel) / 2;

    /* Info about files */
    printLogo();
    printf("\n");
    printf("%s %s starting in SPIM mode.\n", PROGNAME, VERSION);
    printf("  Reading input file %s\n", filename);
    printf("  Writing output file %s for channel 0\n", outname);
    printf("\n");

    /* Print recovered parameters from config file */
    printf("Parameters recovered from config file: \n");
    printf("  Emission wavelength: %3.1f nm\n", spParms.lambda);
    printf("  Numerical aperture: %1.1f\n", spParms.NA);
    printf("  Waist in Z plane: %g um\n", spParms.waist);
    printf("  Image width: %0.2f um\n", spParms.width * spParms.pixel);
    printf("  Image height: %0.2f um\n", spParms.height * spParms.pixel);
    printf("  Pixel size: %0.3f um\n", spParms.pixel);
    printf("  Frame time: %0.3f s\n", spParms.frame_t);
    printf("  Z center of image: %0.2f um\n", spParms.centerz);
    printf("\n");


    /* SPIM routine */
    //int count = 1;
    while (fscanf(fileIn, "%s %f %f %f\n", molname, &x, &y, &z) != EOF) {
	if (x == 100) {
	    prog = 100 * (y / z);
	    printf("Progress: %.1f%%\r", prog);
	    if (((int) y + 1) % nbin == 0) {

		for (int i = 0; i < spParms.height; i++) {
		    TIFFWriteScanline(tif, CCD_buf[i], i, 0);
		}

		TIFFWriteDirectory(tif);
		//fprintf(stdout, "wrote directory %d\n",count );
		//count++;
		writeImageTIFFtags(tif, spParms.width, spParms.height);
		free(CCD_buf);
		char **CCD_buf =
		    (char **) malloc(spParms.height * sizeof(char *));
		for (int i = 0; i < spParms.height; i++) {
		    CCD_buf[i] =
			(char *) malloc(spParms.width * sizeof(char));
		}

		for (int i = 0; i < spParms.height; i++) {
		    for (int j = 0; j < spParms.width; j++) {
			CCD_buf[i][j] = 0;
		    }
		}
	    }
	} else {
	    x += R * randn();
	    y += R * randn();

	    if ((x > -lx && x < lx) && (y > -ly && y < ly)) {
		int idx_x = floor((x + lx) / spParms.pixel);
		int idx_y = floor((y + ly) / spParms.pixel);

		CCD_buf[idx_x][idx_y] +=
		    spimPSF(z, spParms.waist, spParms.centerz,
			    cParms.nevents, cParms.sChannel[0].q[0]);
	    } else {
		continue;
	    }
	}
    }
    return 0;
}


/**********************************************************************************
* Random gaussian number
***********************************************************************************/

double randn()
{
    double unif[] = { 1.0 * rand() / RAND_MAX, 1.0 * rand() / RAND_MAX };
    double out = sqrt(-2 * log(unif[0])) * cos(2 * M_PI * unif[1]);
    return out;
}
