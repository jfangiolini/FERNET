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
 * Parse common parameters from config file
 ***********************************************************************************/

struct commonParms parseCommon(config_t cfg, FILE * fileIn)
{
    struct commonParms cParms;

    /* Get common parameters */
    config_setting_t *common = config_lookup(&cfg, "common");

    /* Get event time */
    if (!config_setting_lookup_float(common, "kappa", &cParms.kappa)) {
	parseError("kappa");
    }

    /* Parse channel 0 */
    config_setting_t *pChann0 =
	config_setting_get_member(common, "channel0");
    if (pChann0 == NULL) {
	fprintf(stderr, "Missing channel0 block in configuration file.\n");
	exit(1);
    } else {
	cParms.sChannel[0] = parseChannel(pChann0, cParms.kappa);
    }

    /* Parse channel 1 */
    config_setting_t *pChann1 =
	config_setting_get_member(common, "channel1");
    if (pChann0 == NULL) {
	fprintf(stderr, "Missing channel1 block in configuration file.\n");
	exit(1);
    } else {
	cParms.sChannel[1] = parseChannel(pChann1, cParms.kappa);
    }

    /* Get w_xy waist */
    if (!config_setting_lookup_float(common, "w_xy", &cParms.w_xy)) {
	parseError("w_xy");
    }

    /* Get w_z waist */
    if (!config_setting_lookup_float(common, "w_z", &cParms.w_z)) {
	parseError("w_z");
    }

    /* Parse input file header */
    fscanf(fileIn, "%f %f\n", &cParms.simu_dt, &cParms.mD);
    cParms.nevents = round(cParms.simu_dt / cParms.kappa);
    double tauD = (cParms.w_xy * cParms.w_xy) / (4 * cParms.mD * 1e8);
    if (tauD < 10 * cParms.simu_dt) {
	fprintf(stderr,
		"MCell time step is not adequate for the maximum D simulated (%g cm^2/s). Adjust MCell time step.\n",
		cParms.mD);
	/* cerrar los files */
	exit(1);
    }

    return cParms;
}

/***********************************************************************************
 * Parse point mode parameters from config file
 ***********************************************************************************/

struct pointParms parsePoint(config_t cfg)
{
    struct pointParms pParms;

    /* Get specific point mode parameters */
    config_setting_t *point = config_lookup(&cfg, "point");

    /* Get x center */
    if (!config_setting_lookup_float(point, "centerx", &pParms.centerx)) {
	parseError("centerx");
    }

    /* Get y center */
    if (!config_setting_lookup_float(point, "centery", &pParms.centery)) {
	parseError("centery");
    }

    /* Get z center */
    if (!config_setting_lookup_float(point, "centerz", &pParms.centerz)) {
	parseError("centerz");
    }

    /* Get output name */
    if (!config_setting_lookup_string(point, "prefix", &pParms.prefix)) {
	parseError("prefix");
    }

    return pParms;
}

/***********************************************************************************
 * Parse multi point mode parameters from config file
 ***********************************************************************************/

struct multiParms parseMulti(config_t cfg)
{
    struct multiParms mParms;

    /* Get specific multi mode parameters */
    config_setting_t *multi = config_lookup(&cfg, "multi");

    /* Get z center */
    if (!config_setting_lookup_float(multi, "centerz", &mParms.centerz)) {
	parseError("centerz");
    }

    /* Get output name prefix */
    if (!config_setting_lookup_string(multi, "prefix", &mParms.prefix)) {
	parseError("prefix");
    }

    /* Get number of PSFs in each axis */
    if (!config_setting_lookup_int(multi, "nPSFX", &mParms.nPSFX)) {
	parseError("nPSFX");
    }
    if (!config_setting_lookup_int(multi, "nPSFY", &mParms.nPSFY)) {
	parseError("nPSFY");
    }

    /* Get separation between PSFs */
    if (!config_setting_lookup_float(multi, "dx", &mParms.dx)) {
	parseError("dx");
    }
    if (!config_setting_lookup_float(multi, "dy", &mParms.dy)) {
	parseError("dy");
    }

    return mParms;
}

/***********************************************************************************
 * Parse linescan mode parameters from config file
 ***********************************************************************************/

struct lineParms parseLine(config_t cfg)
{
    struct lineParms lParms;

    /* Get specific parameters for line mode */
    config_setting_t *line = config_lookup(&cfg, "line");

    /* Get x center */
    if (!config_setting_lookup_float(line, "centerx", &lParms.centerx)) {
	parseError("centerx");
    }

    /* Get y center */
    if (!config_setting_lookup_float(line, "centery", &lParms.centery)) {
	parseError("centery");
    }

    /* Get z center */
    if (!config_setting_lookup_float(line, "centerz", &lParms.centerz)) {
	parseError("centerz");
    }

    /* Get number of columns */
    if (!config_setting_lookup_int(line, "n_columns", &lParms.ncolumn)) {
	parseError("n_columns");
    }

    /* Get shift between columns */
    if (!config_setting_lookup_float(line, "shift", &lParms.shift)) {
	parseError("shift");
    }

    /* Get output name */
    if (!config_setting_lookup_string(line, "tiffname", &lParms.tiffname)) {
	parseError("tiffname");
    }

    /* Get dead time */
    if (!config_setting_lookup_float(line, "deadtime", &lParms.deadtime)) {
	parseError("deadtime");
    }

    return lParms;
}

/***********************************************************************************
 * Parse raster mode parameters from config file
 ***********************************************************************************/

struct imageParms parseImage(config_t cfg)
{
    struct imageParms iParms;

    /* Get specific parameters for image mode */
    config_setting_t *image = config_lookup(&cfg, "image");

    /* Get z center */
    if (!config_setting_lookup_float(image, "centerz", &iParms.centerz)) {
	parseError("centerz");
    }

    /* Get pixel size */
    if (!config_setting_lookup_float(image, "pixel", &iParms.pixel)) {
	parseError("pixel");
    }

    /* Get width */
    if (!config_setting_lookup_int(image, "width", &iParms.width)) {
	parseError("width");
    }

    /* Get height */
    if (!config_setting_lookup_int(image, "height", &iParms.height)) {
	parseError("height");
    }

    /* Get output name */
    if (!config_setting_lookup_string(image, "tiffname", &iParms.tiffname)) {
	parseError("tiffname");
    }

    /* Get dead time */
    if (!config_setting_lookup_float(image, "deadtime", &iParms.deadtime)) {
	parseError("linetime");
    }

    return iParms;
}

/***********************************************************************************
 * Parse stack mode parameters from config file
 ***********************************************************************************/

struct stackParms parseStack(config_t cfg)
{
    struct stackParms sParms;

    /* Get specific parameters for stack mode */
    config_setting_t *stack = config_lookup(&cfg, "stack");

    /* Get pixel size */
    if (!config_setting_lookup_float(stack, "pixel", &sParms.pixel)) {
	parseError("pixel");
    }

    /* Get width */
    if (!config_setting_lookup_int(stack, "width", &sParms.width)) {
	parseError("width");
    }

    /* Get height */
    if (!config_setting_lookup_int(stack, "height", &sParms.height)) {
	parseError("height");
    }

    /* Get output name */
    if (!config_setting_lookup_string(stack, "tiffname", &sParms.tiffname)) {
	parseError("tiffname");
    }

    /* Get dead time */
    if (!config_setting_lookup_float(stack, "deadtime", &sParms.deadtime)) {
	parseError("linetime");
    }

    /* Get top Z position */
    if (!config_setting_lookup_float(stack, "top_z", &sParms.top_z)) {
	parseError("top_z");
    }

    /* Get bottom Z position */
    if (!config_setting_lookup_float(stack, "bot_z", &sParms.bot_z)) {
	parseError("bot_z");
    }

    /* Get Z step */
    if (!config_setting_lookup_float(stack, "step", &sParms.step)) {
	parseError("step");
    }

    /* Checking if parameters are OK */

    if (sParms.bot_z > sParms.top_z) {
	fprintf(stderr,
		"Error: bot_z position is greater than top_z position.\n");
	/* cerrar files? */
	exit(1);
    }

    if (sParms.step > (sParms.top_z - sParms.bot_z)) {
	fprintf(stderr, "Error: Z step is greater than stack height.\n");
	/* cerrar files? */
	exit(1);
    }

    return sParms;
}

 /***********************************************************************************
 * Parse spim mode parameters from config file
 ***********************************************************************************/

struct spimParms parseSpim(config_t cfg)
{
    struct spimParms spParms;

    /* Get specific spim mode parameters */
    config_setting_t *spim = config_lookup(&cfg, "spim");

    /* Get numerical aperture */
    if (!config_setting_lookup_float(spim, "NA", &spParms.NA)) {
	parseError("NA");
    }

    /* Get emission wavelength */
    if (!config_setting_lookup_float(spim, "lambda", &spParms.lambda)) {
	parseError("lambda");
    }

    /* Get z center */
    if (!config_setting_lookup_float(spim, "centerz", &spParms.centerz)) {
	parseError("centerz");
    }

    /* Get waist */
    if (!config_setting_lookup_float(spim, "waist", &spParms.waist)) {
	parseError("waist");
    }

    /* Get output name */
    if (!config_setting_lookup_string(spim, "tiffname", &spParms.tiffname)) {
	parseError("tiffname");
    }

    /* Get pixel size */
    if (!config_setting_lookup_float(spim, "pixel", &spParms.pixel)) {
	parseError("pixel");
    }

    /* Get frame time */
    if (!config_setting_lookup_float(spim, "frame_t", &spParms.frame_t)) {
	parseError("frame_t");
    }

    /* Get CCD width */
    if (!config_setting_lookup_int(spim, "width", &spParms.width)) {
	parseError("width");
    }

    /* Get CCD height */
    if (!config_setting_lookup_int(spim, "height", &spParms.height)) {
	parseError("height");
    }

    return spParms;
}

/***********************************************************************************
 * Parse orbital scanning mode parameters from config file
 ***********************************************************************************/

struct orbitParms parseOrbit(config_t cfg)
{
    struct orbitParms orParms;

    /* Get specific orbital scanning mode parameters */
    config_setting_t *orbit = config_lookup(&cfg, "orbital");

    /* Get x center */
    if (!config_setting_lookup_float(orbit, "centerx", &orParms.centerx)) {
	parseError("centerx");
    }

    /* Get y center */
    if (!config_setting_lookup_float(orbit, "centery", &orParms.centery)) {
	parseError("centery");
    }

    /* Get z center */
    if (!config_setting_lookup_float(orbit, "centerz", &orParms.centerz)) {
	parseError("centerz");
    }

    /* Get orbit radius */
    if (!config_setting_lookup_float(orbit, "radius", &orParms.radius)) {
	parseError("radius");
    }

    /* Get orbit radius */
    if (!config_setting_lookup_float(orbit, "period", &orParms.period)) {
	parseError("period");
    }

    /* Get output name */
    if (!config_setting_lookup_string
	(orbit, "tiffname", &orParms.tiffname)) {
	parseError("tiffname");
    }

    return orParms;
}

/***********************************************************************************
 * Parse channel configuration block
 ***********************************************************************************/

struct channelInfo parseChannel(config_setting_t * pChan, float dt)
{
    struct channelInfo chan;
    const char *tmp;
    if (!config_setting_lookup_int(pChan, "on", &chan.status)) {
	parseError("on");
    }
    if (chan.status == 1) {
	const config_setting_t *molnames =
	    config_setting_get_member(pChan, "molec");
	if (molnames == NULL) {
	    parseError("molec");
	} else {
	    chan.nmols = config_setting_length(molnames);
	    chan.mols = (char **) malloc(chan.nmols * sizeof(char *));

	    for (int i = 0; i < chan.nmols; i++) {
		tmp = config_setting_get_string_elem(molnames, i);
		chan.mols[i] =
		    (char *) malloc((strlen(tmp) + 1) * sizeof(char));
		strcpy(chan.mols[i], tmp);
	    }

	    const config_setting_t *molbright =
		config_setting_get_member(pChan, "bright");
	    if (config_setting_length(molnames) !=
		config_setting_length(molbright)) {
		fprintf(stderr,
			"Molecules name and brightness list must be equal in size.\n");
		exit(1);
	    }

	    chan.bright = (int *) malloc(chan.nmols * sizeof(int));
	    for (int i = 0; i < chan.nmols; i++) {
		chan.bright[i] = config_setting_get_int_elem(molbright, i);
	    }

	    chan.q = (double *) malloc(chan.nmols * sizeof(double));
	    for (int i = 0; i < chan.nmols; i++) {
		chan.q[i] = (double) chan.bright[i] * dt;
	    }
	}
    }
    return chan;
}
