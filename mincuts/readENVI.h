#include <stdint.h>

/********************************
 * Definitions
 ********************************/

#define INTERLEAVE_BSQ "bsq"
#define INTERLEAVE_BIL "bil"
#define INTERLEAVE_BIP "bip"


/********************************
 * Enums
 ********************************/

typedef enum 
{
	BSQ,
	BIL,
	BIP
} interleave_t;

typedef enum
{
	L_ENDIAN,
	B_ENDIAN
} byte_order_t;

typedef enum
{
	INT8 = 1,
	INT16 = 2,
	INT32 = 3,
	FLOAT = 4,
	DOUBLE = 5,
	COMPLEX_FLOAT = 6,
	COMPLEX_DOUBLE = 9,
	UINT16 = 12,
	UINT32 = 13,
	INT64 = 14,
	UINT64 = 15
} date_type_t;


/********************************
 * Typedef Structs
 ********************************/

typedef struct
{
	uint32_t samples;
	uint32_t lines;
	uint32_t bands;

	uint32_t header_offset;

	date_type_t data_type;

	interleave_t interleave;

	byte_order_t byte_order;
} image_info_t;


/********************************
 * Functions
 ********************************/

void freeImageCube(float ***image_cube, image_info_t *image_info);

float ***readImageCube(char const *data_file_path, char const *header_file_path, image_info_t *image_info);