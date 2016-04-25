#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"

#define LINE_LENGTH 256

extern struct gconf gconf;

void readConfig(const char * fn)
{
  FILE * fp;
  char line[LINE_LENGTH+1];
  char *tokenKey;
  char *tokenValue;

  if((fp=fopen(fn,"r"))==NULL)
  {
    fprintf(stderr,"Could not open configuration file %s\n",fn);
  } else
  {

    while(fgets(line,sizeof(line),fp)!=NULL)
    {
      tokenKey=strtok(line,"\t =\n\r");
      if(tokenKey!=NULL && tokenKey[0]!='#')
      {
        tokenValue=strtok(NULL,"\t =\n\r");
        if(tokenValue==NULL) continue;

        if(strcasecmp(tokenKey,"inputData")==0)
        {
          strcpy(gconf.inputData,tokenValue);
        } 
        else if(strcasecmp(tokenKey,"nx")==0)
        {
          gconf.nx=strtol(tokenValue,NULL,10);
        }
        else if(strcasecmp(tokenKey,"ny")==0)
        {
          gconf.ny=strtol(tokenValue,NULL,10);
        }
        else if(strcasecmp(tokenKey,"nw")==0)
        {
          gconf.nw=strtol(tokenValue,NULL,10);
        }
        else
        {
          fprintf(stderr,"Unrecognized key %s in config file, skipping\n",tokenKey);
          continue;
        }
      }
    }
  }
}

void printConfig(void)
{
  fprintf(stderr, "Run configuration\n"
      "\t inputData: %s\n"
      "\t nx: %ld\n"
      "\t ny: %ld\n"
      "\t nw: %ld\n",
      gconf.inputData,gconf.nx,gconf.ny,gconf.nw);
}
