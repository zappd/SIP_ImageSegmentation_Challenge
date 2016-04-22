#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "readENVI.h"

/********************************
 * Functions
 ********************************/

static int parseHeaderFile(char *header_file, image_info_t *image_info);
static int readHeaderFile(char const *header_file_path, char **string_buffer);
static int readDataFile(char const *data_file_path, float ****image_cube, image_info_t *image_info);


/********************************
 * Variables
 ********************************/

static const char *SAMPLES_KEY = "samples";
static const char *LINES_KEY = "lines";
static const char *BANDS_KEY = "bands";
static const char *HEADER_OFFSET_KEY = "header offset";
static const char *DATA_TYPE_KEY = "data type";
static const char *INTERLEAVE_KEY = "interleave";
static const char *BYTE_ORDER_KEY = "byte order";


float ***readImageCube(char const *data_file_path, char const *header_file_path, image_info_t *image_info)
{
	char *header_file;
	float ***image_cube = NULL;

	if (readHeaderFile(header_file_path, &header_file) < 0)
	{
		return NULL;
	}

	if (parseHeaderFile(header_file, image_info) < 0)
	{
		return NULL;
	}

	// check to see if we can read this file in
	if (image_info->data_type != FLOAT)
	{
		printf("Only files with FLOAT data are supported at this time");
		return NULL;
	}

	if (image_info->interleave != BSQ) 
	{
		printf("Only files with BSQ-interleaved data are supported at this time");
		return NULL;
	}

	if (image_info->header_offset > 0)
	{
		printf("Only files with 0 header offset are supported at this time");
		return NULL;
	}

	if (readDataFile(data_file_path, &image_cube, image_info) < 0)
	{
		return NULL;
	}

	return image_cube;
}

void freeImageCube(float ***image_cube, image_info_t *image_info)
{
	uint32_t sample, line;

	for (line = 0; line < image_info->lines; line++)
	{
		for (sample = 0; sample < image_info->samples; sample++)
		{
			free(image_cube[line][sample]);
		}

		free(image_cube[line]);
	}

	free(image_cube);
}

static int readDataFile(char const *data_file_path, float ****image_cube, image_info_t *image_info)
{
	FILE *file_pointer;
    long band_length;
    float *band_buffer;
	uint32_t sample, line, band;

	// initialize the cube: allocate space for each pixel
	*image_cube = (float ***)malloc(image_info->lines * sizeof(float **));

	for (line = 0; line < image_info->lines; line++)
	{
		(*image_cube)[line] = (float **)malloc(image_info->samples * sizeof(float *));

		for (sample = 0; sample < image_info->samples; sample++)
		{
			(*image_cube)[line][sample] = (float *)malloc(image_info->bands * sizeof(float));
		}
	}
	// done initializing cube

    // open the file
    if((file_pointer = fopen(data_file_path, "rb")) == NULL) 
    {
    	printf("File Not Found: %s\n", data_file_path);
        return -1;
    }

    band_length = (image_info->samples * image_info->lines);

    // allocate a buffer to hold one band of the file
    if((band_buffer = malloc(band_length * sizeof(float))) == NULL) {
    	printf("Error allocating space for the band\n");
        return -1;
    }

	fseek(file_pointer, 0, SEEK_SET);

    // read each band into the buffer
    for (band = 0; band < image_info->bands; band++)
    {
		fseek(file_pointer, (band_length*sizeof(float)*band), SEEK_SET);

	    if(fread(band_buffer, sizeof(float), band_length, file_pointer) != band_length)
	    {
	        free(band_buffer);
			printf("Error Reading Band %d\n", band);
			return -1;
	    }

	    // assign each band value to the correct pixel
		for (line = 0; line < image_info->lines; line++)
		{
			for (sample = 0; sample < image_info->samples; sample++)
			{
				(*image_cube)[line][sample][band] = band_buffer[line*image_info->samples + sample];
			}
		}

		printf("%.4f\t\n", band, (*image_cube)[400][356][band]);
    }

    // free the temporary buffer memory
    free(band_buffer);

    // close the file
    if(fclose(file_pointer) == EOF)
    {
        return -1;
    }

    return 0;
}

static int readHeaderFile(char const* header_file_path, char **string_buffer)
{
    FILE *file_pointer;
    long file_length;
    size_t file_size;

    // open the file 
    if((file_pointer = fopen(header_file_path, "rb")) == NULL) 
    {
    	printf("File Not Found: %s\n", header_file_path);
        return -1;
    }

    // seek to the end of the file
    if(fseek(file_pointer, 0, SEEK_END) != 0) 
    {
    	printf("Error seeking in file\n");
        return -1;
    }

    // find byte offset to the end of the file (size)
    if((file_length = ftell(file_pointer)) < 0) 
    {
    	printf("Error getting file length: %ld\n", file_length);
        return -1;
    }

    file_size = (size_t)file_length;

    // allocate a buffer to hold the whole file
    // add an extra space to null-terminate the string
    if((*string_buffer = malloc( file_size + 1 )) == NULL) 
    {
    	printf("Error allocating space for string\n");
        return -1;
    }

    // rewind file pointer to start of file
    rewind(file_pointer);

    // read the whole file into buffer
    if(file_size != fread(*string_buffer, sizeof(char), file_size, file_pointer)) 
    {
        free(*string_buffer);
        return -1;
    }

    // close the file
    if(fclose(file_pointer) == EOF) 
    {
        free(*string_buffer);
        return -1;
    }

    // null-terminate the string
	(*string_buffer)[file_size] = '\0';

    // return the file length
    return file_length;
}

static int parseHeaderFile(char *header_file, image_info_t *image_info)
{
	char delimiter[2] = "\r\n";
	char *token = strtok(header_file, delimiter);

	while (token != NULL)
	{
		if (memcmp(token, SAMPLES_KEY, strlen(SAMPLES_KEY)) == 0)
		{
			token += strlen(SAMPLES_KEY);
			for (; *token == ' ' || *token == '='; token++);
			image_info->samples = (uint32_t)strtoul(token, NULL, 10);
		}
		else if (memcmp(token, LINES_KEY, strlen(LINES_KEY)) == 0)
		{
			token += strlen(LINES_KEY);
			for (; *token == ' ' || *token == '='; token++);
			image_info->lines = (uint32_t)strtoul(token, NULL, 10);
		}
		else if (memcmp(token, BANDS_KEY, strlen(BANDS_KEY)) == 0)
		{
			token += strlen(BANDS_KEY);
			for (; *token == ' ' || *token == '='; token++);
			image_info->bands = (uint32_t)strtoul(token, NULL, 10);
		}
		else if (memcmp(token, HEADER_OFFSET_KEY, strlen(HEADER_OFFSET_KEY)) == 0)
		{
			token += strlen(HEADER_OFFSET_KEY);
			for (; *token == ' ' || *token == '='; token++);
			image_info->header_offset = (uint32_t)strtoul(token, NULL, 10);
		}
		else if (memcmp(token, DATA_TYPE_KEY, strlen(DATA_TYPE_KEY)) == 0)
		{
			token += strlen(DATA_TYPE_KEY);
			for (; *token == ' ' || *token == '='; token++);

			uint32_t data_type = (uint32_t)strtoul(token, NULL, 10);

			if (data_type < 1 || data_type > 15 ||
			   (data_type > 6 && data_type < 9))
			{
				printf("Error Parsing Header File: Unknown Data Type");
				return -1;
			}
			else 
			{
				image_info->data_type = data_type;
			}
 		}
		else if (memcmp(token, INTERLEAVE_KEY, strlen(INTERLEAVE_KEY)) == 0)
		{
			token += strlen(INTERLEAVE_KEY);
			for (; *token == ' ' || *token == '='; token++);
			
			if (strcmp(token, INTERLEAVE_BSQ) == 0)
			{
				image_info->interleave = BSQ;
			}
			else if (strcmp(token, INTERLEAVE_BIL) == 0)
			{
				image_info->interleave = BIL;
			}
			else if (strcmp(token, INTERLEAVE_BIP) == 0)
			{
				image_info->interleave = BIP;
			}
			else {
				printf("Error Parsing Header File: Unknown Interleave Format");
				return -1;
			}
		}
		else if (memcmp(token, BYTE_ORDER_KEY, strlen(BYTE_ORDER_KEY)) == 0)
		{
			token += strlen(BYTE_ORDER_KEY);
			for (; *token == ' ' || *token == '='; token++);

			long byte_order = strtoul(token, NULL, 10);

			if (byte_order == 0) 
			{
				image_info->byte_order = L_ENDIAN;
			}
			else if (byte_order == 1)
			{
				image_info->byte_order = B_ENDIAN;
			}
			else
			{
				printf("Error Parsing Header File: Unknown Byte Order");
				return -1;
			}
		}

		token = strtok(NULL, delimiter);
	}

	return 0;
}