#if !defined(_CONFIG_H)
#define _CONFIG_H


typedef struct
{
    char  inputData[256];
    int   nw;
    float sigma;
    float evcrit;
    int   t;
    float maxevfact;
    int   krep;
    int   kiter;
    int   kcent;
    int   verbosity;
    float taucardinality;
    float kelbw;
    int   cardmin;
} global_config_t;

extern global_config_t global_config;

void readConfig(const char *);

void printConfig(void);


#endif
