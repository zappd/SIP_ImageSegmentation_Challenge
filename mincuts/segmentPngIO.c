#include <png.h>
#include <stdlib.h>
#include <stdint.h>
#include "segmentPngIO.h"

/*
 * Based Heavily on the PNG writing implementation found here:
 * http://www.lemoda.net/c/write-png/
 */


/********************************
 * Functions
 ********************************/

static png_color *pixelAt(uint32_t x, uint32_t y);

static int savePngToFile(const char *path);


/********************************
 * Global Vars
 ********************************/

static png_color *pixels;
static size_t    image_width;
static size_t    image_height;


void writeSegmentImage(uint32_t **image, char *file_path, uint32_t width, uint32_t height)
{
    image_width  = width;
    image_height = height;

    pixels = calloc(sizeof(png_color), (image_width * image_height));

    for (uint32_t y = 0; y < image_height; y++)
    {
        for (uint32_t x = 0; x < image_width; x++)
        {
            png_color *pixel = pixelAt(x, y);
            pixel->red   = (uint8_t) ((image[x][y] * 114) % 256);
            pixel->green = (uint8_t) ((image[x][y] * 19) % 256);
            pixel->blue  = (uint8_t) ((image[x][y] * 88) % 256);
        }
    }

    savePngToFile(file_path);

    free(pixels);
}

static png_color *pixelAt(uint32_t x, uint32_t y)
{
    return (pixels + (image_width * y) + x);
}

static int savePngToFile(const char *path)
{
    FILE        *file_pointer;
    png_byte    **row_pointers = NULL;
    png_structp png_ptr        = NULL;
    png_infop   info_ptr       = NULL;

    /* pixels sized for R, G, and B */
    int pixel_size = 3;

    /* with a bit depth of 8 */
    int bit_depth = 8;

    if (!(file_pointer = fopen(path, "wb")))
    {
        printf("Failed to open file '%s'\n", path);
        return -1;
    }

    if ((png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL)) == NULL)
    {
        printf("Failed to write PNG struct\n");
        fclose(file_pointer);
        return -1;
    }

    if ((info_ptr = png_create_info_struct(png_ptr)) == NULL)
    {
        printf("Failed to create PNG info\n");
        fclose(file_pointer);
        return -1;
    }

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        printf("PNG Failure");
        fclose(file_pointer);
        return -1;
    }

    /* set image attributes */
    png_set_IHDR(png_ptr, info_ptr, image_width, image_height, bit_depth,
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    /* initialize rows of PNG */
    row_pointers = png_malloc(png_ptr, image_height * sizeof(png_byte *));

    for (uint32_t y = 0; y < image_height; ++y)
    {
        row_pointers[y] = (png_byte *) png_malloc(png_ptr, sizeof(uint8_t) * image_width * pixel_size);

        for (uint32_t x = 0; x < image_width; ++x)
        {
            png_color *pixel = pixelAt(x, y);

            row_pointers[y][(pixel_size * x) + 0] = pixel->red;
            row_pointers[y][(pixel_size * x) + 1] = pixel->green;
            row_pointers[y][(pixel_size * x) + 2] = pixel->blue;
        }
    }

    /* write the image data to file_pointer */
    png_init_io(png_ptr, file_pointer);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    for (uint32_t y = 0; y < image_height; y++)
    {
        png_free(png_ptr, row_pointers[y]);
    }
    png_free(png_ptr, row_pointers);

    return 0;
}