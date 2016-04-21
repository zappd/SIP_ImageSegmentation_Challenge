#include <stdio.h>

#include "config.h"
#include "readENVI.h"
#include "segmenter.h"
#include "graphCuts.h"

/* Global configuration, initialized to defaults */
struct gconf gconf = {
  .inputData="data.data",
  .nx=64,
  .ny=64,
  .nw=64
};

int main(int argc, char* argv[])
{
    readConfig("seconfig");
    /*printConfig();*/

    return 0;
}
