#if !defined(_CONFIG_H)
#define _CONFIG_H


struct gconf
{
  char inputData[256];
  int nx;
  int ny;
  int nw;
  float sigma;
  float evcrit;
  int t;
  float maxevfact;
  int krep;
  int kiter;
  int kcent;
};

extern struct gconf gconf;

void readConfig(const char *);
void printConfig(void);


#endif
