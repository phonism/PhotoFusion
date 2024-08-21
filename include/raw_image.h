#ifndef RAW_IMAGE_H
#define RAW_IMAGE_H

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdint>

struct RawImage {
    unsigned int width;
    unsigned int height;
    unsigned short bitDepth;
    unsigned short black_level;
    unsigned short data_max;
    unsigned short maximum = 16383;
    unsigned int filters = 3031741620;
    
    unsigned short* raw_data;
    unsigned short (*image)[4];
    unsigned short (*rgb_image)[3];
    unsigned short (*yuv_image)[3];
    uint64_t (*integral)[3];
    unsigned short (*denoised_image)[3];
    float cam_mul[4] = {1.9722067118f, 0.9411969781f, 1.1376225948f, 0.9411969781f};
    unsigned short curve[0x10000];
    int (*histogram)[0x2000] = NULL;

    unsigned *oprof = NULL;
    int output_bps = 8;

    double gamm[6] = {0.45, 4.5, 0, 0, 0, 0};
    FILE *output;

    ~RawImage() {
        if (raw_data) {
            free(raw_data);
            raw_data = nullptr;
        }
        if (image) {
            free(image);
            image = nullptr;
        }
        if (rgb_image) {
            free(rgb_image);
            rgb_image = nullptr;
        }
        if (yuv_image) {
            free(yuv_image);
            yuv_image = nullptr;
        }
        if (denoised_image) {
            free(denoised_image);
            denoised_image = nullptr;
        }
        if (histogram) {
            free(histogram);
            histogram = nullptr;
        }
        if (oprof) {
            free(oprof);
            oprof = nullptr;
        }
        if (output) {
            fclose(output);
            output = nullptr;
        }
    }

    const char* to_string() {
        return "NONE";
    }
};

#endif // RAW_IMAGE_H
