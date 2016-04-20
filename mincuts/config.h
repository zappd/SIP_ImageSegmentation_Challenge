#if !defined(_CONFIG_H)
#define _CONFIG_H


struct gconf
{
  char inputData[256];
  long nx;
  long ny;
  long nw;
};


void readConfig(const char *);
void printConfig(void);


#endif
