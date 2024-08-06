#ifndef RAW_PROCESSOR_H
#define RAW_PROCESSOR_H

// #include <omp.h>
#include <arpa/inet.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <cmath>
#include "RawImage.h"
#include "Common.h"

struct tiff_tag {
    unsigned short tag, type;
    int count;
    union {
        char c[4];
        short s[2];
        int i;
    } val;
};

struct new_tiff_hdr {
    unsigned short t_order, magic;
    int ifd;
    unsigned short pad, ntag;
    struct tiff_tag tag[23];
    int nextifd;
    unsigned short pad2, nexif;
    struct tiff_tag exif[4];
    unsigned short pad3, ngps;
    struct tiff_tag gpst[10];
    short bps[4];
    int rat[10];
    unsigned int gps[26];
    char t_desc[512], t_make[64], t_model[64], soft[32], date[20], t_artist[64];
};

class RawProcessor {
public:
    static void raw2image(RawImage& raw_image);

    static void adjust_maximum(RawImage& raw_image);

    static void scale_colors(RawImage& raw_image);

    static void pre_interpolate(RawImage& raw_image);

    static void ahd_interpolate(RawImage& raw_image);

    static void ahd_interpolate_green_h_and_v(
            RawImage& raw_image, int top, int left, unsigned short (*out_rgb)[AHD_TILE][AHD_TILE][3]);

    static void border_interpolate(RawImage& raw_iamge, int border);

    static void ahd_interpolate_r_and_b_and_convert_to_cielab(
            RawImage& raw_image, int top, int left, unsigned short (*inout_rgb)[AHD_TILE][AHD_TILE][3],
            short (*out_lab)[AHD_TILE][AHD_TILE][3]);

    static void ahd_interpolate_r_and_b_in_rgb_and_convert_to_cielab(
            RawImage& raw_image, int top, int left, 
            unsigned short (*inout_rgb)[AHD_TILE][3], short (*out_lab)[AHD_TILE][3]);

    static void ahd_interpolate_build_homogeneity_map(
            RawImage& raw_image, int top, int left, short (*lab)[AHD_TILE][AHD_TILE][3],
            char (*out_homogeneity_map)[AHD_TILE][2]);

    static void ahd_interpolate_combine_homogeneous_pixels(
            RawImage& raw_image, int top, int left, ushort (*rgb)[AHD_TILE][AHD_TILE][3],
            char (*homogeneity_map)[AHD_TILE][2]);

    static void convert_to_rgb(RawImage& raw_image);

    static int ppm_tiff_writer(RawImage& raw_image, const char *filename);

    static void gamma_adjustment(RawImage& raw_image);

};

#endif // RAW_PROCESSOR_H

