#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"

#define LINE_LENGTH 256

extern global_config_t global_config;

void readConfig(const char *fn)
{
    FILE *file_pointer;

    char line[LINE_LENGTH + 1];
    char *tokenKey;
    char *tokenValue;

    if ((file_pointer = fopen(fn, "r")) == NULL)
    {
        fprintf(stderr, "Could not open configuration file %s\n", fn);
    }
    else
    {
        while (fgets(line, sizeof(line), file_pointer) != NULL)
        {
            tokenKey = strtok(line, "\t =\n\r");
            if (tokenKey != NULL && tokenKey[0] != '#')
            {
                tokenValue = strtok(NULL, "\t =\n\r");

                if (tokenValue == NULL)
                {
                    continue;
                }

                if (strcmp(tokenKey, "inputData") == 0)
                {
                    strcpy(global_config.inputData, tokenValue);
                }
                else if (strcmp(tokenKey, "nw") == 0)
                {
                    global_config.nw = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcmp(tokenKey, "sigma") == 0)
                {
                    global_config.sigma = strtof(tokenValue, NULL);
                }
                else if (strcmp(tokenKey, "evcrit") == 0)
                {
                    global_config.evcrit = strtof(tokenValue, NULL);
                }
                else if (strcmp(tokenKey, "t") == 0)
                {
                    global_config.t = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcmp(tokenKey, "maxevfact") == 0)
                {
                    global_config.maxevfact = strtof(tokenValue, NULL);
                }
                else if (strcmp(tokenKey, "kelbw") == 0)
                {
                    global_config.kelbw = strtof(tokenValue, NULL);
                }
                else if (strcmp(tokenKey, "taucardinality") == 0)
                {
                    global_config.taucardinality = strtof(tokenValue, NULL);
                }
                else if (strcmp(tokenKey, "krep") == 0)
                {
                    global_config.krep = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcmp(tokenKey, "kiter") == 0)
                {
                    global_config.kiter = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcmp(tokenKey, "kcent") == 0)
                {
                    global_config.kcent = strtof(tokenValue, NULL);
                }
                else if (strcmp(tokenKey, "verbosity") == 0)
                {
                    global_config.verbosity = (int) strtol(tokenValue, NULL, 10);
                }
                else if (strcmp(tokenKey, "cardmin") == 0)
                {
                    global_config.cardmin = (int) strtol(tokenValue, NULL, 10);
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
            global_config.inputData, global_config.nw, global_config.sigma,
            global_config.evcrit, global_config.maxevfact, global_config.t,
            global_config.krep, global_config.kiter, global_config.kcent,
            global_config.verbosity, global_config.taucardinality,
            global_config.kelbw, global_config.cardmin);
}
