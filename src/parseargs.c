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
 * Parse console arguments
 ***********************************************************************************/
struct args parseArgs(int argc, char **argv)
{
	struct args Args;

	/* Command line options */
	struct arg_file *infile = arg_file1(NULL, NULL, "<input>", "input position file");
	struct arg_file *config = arg_file1("c", "config", "<config file>", "configuration file");
	struct arg_lit *help = arg_lit0(NULL, "help", "print this help and exit");
	struct arg_str *mode = arg_str1("m", "mode", "<point,multi,line,raster>", "sampling mode");
	struct arg_lit *version = arg_lit0(NULL, "version", "print version information and exit");
	struct arg_end *end = arg_end(20);
	int nerrors;
	void *argtable[] = { infile, mode, config, help, version, end };

	/* Verify the argtable[] entries were allocated sucessfully */
	if (arg_nullcheck(argtable) != 0) {
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n", argv[0]);
		exit(1);
	}

	/* Parse the command line as defined by argtable[] */
	nerrors = arg_parse(argc, argv, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0) {
		printf("Usage: %s", argv[0]);
		arg_print_syntax(stdout, argtable, "\n");
		printf("FERNET: Fluorescence Emission Recipes and NumErical routines Network.\n");
		printf("This program evaluates the emision of photons in different fluorescence experiments.\n");
		printf("The input file must have the positions of each molecule in each time step.\n");
		arg_print_glossary(stdout, argtable, "  %-35s %s\n");
		exit(0);
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0) {
		printf("'%s' program for generating photon emision.\n", PROGNAME);
		printf("%s, version %s\n", DATE, VERSION);
		exit(0);
	}

	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0) {
		/* Display the error details contained in the arg_end struct. */
		arg_print_errors(stdout, end, argv[0]);
		printf("Try '%s --help' for more information.\n", argv[0]);
		exit(1);
	}

	/* Reading config file */
	config_t cfg;
	config_init(&cfg);

	/* Read the file. If there is an error, report it and exit. */
	if (!config_read_file(&cfg, config->filename[0])) {
		fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
			config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		exit(1);
	}

	Args.filename = infile->filename[0];
	Args.mode = *mode->sval;
	Args.cfg = cfg;

	/* Free table */
	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

	return Args;
}
