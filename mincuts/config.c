#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"

#define LINE_LENGTH 256

extern struct gconf gconf;

void readConfig(const char *fn)
{
    FILE *fp;
    char line[LINE_LENGTH + 1];
    char *tokenKey;
    char *tokenValue;

    if ((fp = fopen(fn, "r")) == NULL)
    {
        fprintf(stderr, "Could not open configuration file %s\n", fn);
    }
    else
    {
        while (fgets(line, sizeof(line), fp) != NULL)
        {
            tokenKey = strtok(line, "\t =\n\r");
            if (tokenKey != NULL && tokenKey[0] != '#')
            {
                tokenValue = strtok(NULL, "\t =\n\r");

                if (tokenValue == NULL) {
                    continue;
                }

                if (strcasecmp(tokenKey, "inputData") == 0)
                {
                    strcpy(gconf.inputData, tokenValue);
                }
                else if (strcasecmp(tokenKey, "nx") == 0)
                {
                    gconf.nx = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "ny") == 0)
                {
                    gconf.ny = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "nw") == 0)
                {
                    gconf.nw = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "sigma") == 0)
                {
                    gconf.sigma = strtof(tokenValue, NULL);
                }
                else if (strcasecmp(tokenKey, "evcrit") == 0)
                {
                    gconf.evcrit = strtof(tokenValue, NULL);
                }
                else if (strcasecmp(tokenKey, "t") == 0)
                {
                    gconf.t = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "maxevfact") == 0)
                {
                    gconf.maxevfact = strtof(tokenValue, NULL);
                }
                else if (strcasecmp(tokenKey, "kelbw") == 0)
                {
                    gconf.kelbw = strtof(tokenValue, NULL);
                }
                else if (strcasecmp(tokenKey, "taucardinality") == 0)
                {
                    gconf.taucardinality = strtof(tokenValue, NULL);
                }
                else if (strcasecmp(tokenKey, "krep") == 0)
                {
                    gconf.krep = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "kiter") == 0)
                {
                    gconf.kiter = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "kcent") == 0)
                {
                    gconf.kcent = strtof(tokenValue, NULL);
                }
                else if (strcasecmp(tokenKey, "verbosity") == 0)
                {
                    gconf.verbosity = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcasecmp(tokenKey, "cardmin") == 0)
                {
                    gconf.cardmin = (int) strtol(tokenValue, NULL, 10);
                }
                else
                {
                    fprintf(stderr, "Unrecognized key %s in config file, skipping\n", tokenKey);
                    continue;
                }
            }
        }
    }
}

void printConfig(void)
{
    fprintf(stderr, "Runtime configuration\n"
                    "\t inputData: %s\n"
                    "\t nx: %d\n"
                    "\t ny: %d\n"
                    "\t nw: %d\n"
                    "\t sigma: %f\n"
                    "\t evcrit: %f\n"
                    "\t maxevfact: %f\n"
                    "\t t: %d\n"
                    "\t krep: %d\n"
                    "\t kiter: %d\n"
                    "\t kcent: %d\n"
                    "\t verbosity: %d\n"
                    "\t taucardinality: %f\n"
                    "\t kelbw: %f\n"
                    "\t cardmin: %d\n\n",
            gconf.inputData, gconf.nx, gconf.ny,
            gconf.nw, gconf.sigma, gconf.evcrit,
            gconf.maxevfact, gconf.t, gconf.krep,
            gconf.kiter, gconf.kcent, gconf.verbosity,
            gconf.taucardinality, gconf.kelbw,
            gconf.cardmin);
}
