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
 * Main function parses arguments from console and config file and calls desired
 * emission routine
 ***********************************************************************************/
int main(int argc, char *argv[])
{
	srand(time(NULL));	// randomize seed
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r, rand());

	/* Parse arguments from command line */
	struct args Args = parseArgs(argc, argv);

	/* Get desired fluorescence mode */
	enum fluo_modes desired_mode;
	if (!strcmp(Args.mode, "point")) {
		desired_mode = POINT;
	} else if (!strcmp(Args.mode, "multi")) {
		desired_mode = MULTI;
	} else if (!strcmp(Args.mode, "line")) {
		desired_mode = LINE;
	} else if (!strcmp(Args.mode, "raster")) {
		desired_mode = RASTER;
	} else if (!strcmp(Args.mode, "stack")) {
		desired_mode = STACK;
	} else if (!strcmp(Args.mode, "spim")) {
		desired_mode = SPIM;
	} else if (!strcmp(Args.mode, "orbit")) {
		desired_mode = ORBIT;
	} else {
		printf("Invalid mode. Type %s --help for usage\n", argv[0]);
		exit(1);
	}

	/* Call fluorescence routine */
	switch (desired_mode) {
	case POINT:
		pointRoutine(Args.cfg, Args.filename, r);
		break;

	case MULTI:
		multiRoutine(Args.cfg, Args.filename, r);
		break;

	case LINE:
		lineRoutine(Args.cfg, Args.filename, r);
		break;

	case RASTER:
		rasterRoutine(Args.cfg, Args.filename, r);
		break;

	case STACK:
		stackRoutine(Args.cfg, Args.filename, r);
		break;

	case SPIM:
		spimRoutine(Args.cfg, Args.filename, r);
		break;

	case ORBIT:
		orbitRoutine(Args.cfg, Args.filename, r);
	}

	/* Cleanup */
	config_destroy(&Args.cfg);
	gsl_rng_free(r);
	printf("\n");

	return 0;
}

/***********************************************************************************
 * Supplementary functions definition
 ***********************************************************************************/
void parseError(char *variable)
{
	fprintf(stderr, "Invalid type or missing '%s' parameter in configuration file.\n", variable);
	exit(1);
}

void printLogo()
{
	printf("                                                     \n");
	printf("  ███████╗███████╗██████╗ ███╗   ██╗███████╗████████╗\n");
	printf("  ██╔════╝██╔════╝██╔══██╗████╗  ██║██╔════╝╚══██╔══╝\n");
	printf("  █████╗  █████╗  ██████╔╝██╔██╗ ██║█████╗     ██║   \n");
	printf("  ██╔══╝  ██╔══╝  ██╔══██╗██║╚██╗██║██╔══╝     ██║   \n");
	printf("  ██║     ███████╗██║  ██║██║ ╚████║███████╗   ██║   \n");
	printf("  ╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝   \n");
	printf("                                                     \n");
	printf(" Fluorescence Emission Recipes and NumErical Toolkit \n");
}
