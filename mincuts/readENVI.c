#include <stdio.h>
#include "config.h"

void readEnvi(float *** data)
{

  FILE * fp = fopen(gconf.inputData,"r");

  if(fread(data,sizeof(float),gconf.nx*gconf.ny*gconf.nw,fp)!=
      (gconf.nx*gconf.ny*gconf.nw))
  {
    fprintf(stderr,"Could not find enough data in the file: %s",gconf.inputData);
  }

  fclose(fp);




}
