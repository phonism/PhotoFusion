#ifndef CONST_H
#define CONST_H

#define AHD_TILE 512
#define BLOCK_SIZE 256
#define HISTOGRAM_SIZE 0x2000

#define FC(filters, row, col) (((filters) >> ((((row) << 1 & 14) | ((col) & 1)) << 1)) & 3)
#define SQR(x) ((x) * (x))
#define ABS(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define LIM(x, min, max) MAX(min, MIN(x, max))
#define ULIM(x, y, z) ((y) < (z) ? LIM(x, y, z) : LIM(x, z, y))
#define CLIP(x) LIM((int)(x), 0, 65535)

// RGB -> YUV
#define RGB2Y(R, G, B) CLIP(( (  66 * (R) + 129 * (G) +  25 * (B) + 128) >> 8) +  16)
#define RGB2U(R, G, B) CLIP(( ( -38 * (R) -  74 * (G) + 112 * (B) + 128) >> 8) + 128)
#define RGB2V(R, G, B) CLIP(( ( 112 * (R) -  94 * (G) -  18 * (B) + 128) >> 8) + 128)

// YUV -> RGB
#define C(Y) ( (Y) - 16  )
#define D(U) ( (U) - 128 )
#define E(V) ( (V) - 128 )

#define YUV2R(Y, U, V) CLIP(( 298 * C(Y)              + 409 * E(V) + 128) >> 8)
#define YUV2G(Y, U, V) CLIP(( 298 * C(Y) - 100 * D(U) - 208 * E(V) + 128) >> 8)
#define YUV2B(Y, U, V) CLIP(( 298 * C(Y) + 516 * D(U)              + 128) >> 8)


#define CLIP_16BIT(value) ((value) < 0 ? 0 : ((value) > 65535 ? 65535 : (value)))

#define RGB2Y_16BIT(R, G, B) CLIP_16BIT((((  66 * ((R) >> 8) + 129 * ((G) >> 8) +  25 * ((B) >> 8) + 128) >> 8) +  16) << 8)
#define RGB2U_16BIT(R, G, B) CLIP_16BIT(((( -38 * ((R) >> 8) -  74 * ((G) >> 8) + 112 * ((B) >> 8) + 128) >> 8) + 128) << 8)
#define RGB2V_16BIT(R, G, B) CLIP_16BIT(((( 112 * ((R) >> 8) -  94 * ((G) >> 8) -  18 * ((B) >> 8) + 128) >> 8) + 128) << 8)


#define C_16BIT(Y) ((Y) - 4096)  // 对应的偏移调整
#define D_16BIT(U) ((U) - 32768) // 对应的偏移调整
#define E_16BIT(V) ((V) - 32768) // 对应的偏移调整

#define YUV2R_16BIT(Y, U, V) CLIP_16BIT((298 * C_16BIT(Y)              + 409 * E_16BIT(V) + 32768) >> 8)
#define YUV2G_16BIT(Y, U, V) CLIP_16BIT((298 * C_16BIT(Y) - 100 * D_16BIT(U) - 208 * E_16BIT(V) + 32768) >> 8)
#define YUV2B_16BIT(Y, U, V) CLIP_16BIT((298 * C_16BIT(Y) + 516 * D_16BIT(U)              + 32768) >> 8)


// RGB -> YCbCr
#define CRGB2Y(R, G, B) CLIP((19595 * R + 38470 * G + 7471 * B ) >> 16)
#define CRGB2Cb(R, G, B) CLIP((36962 * (B - CLIP((19595 * R + 38470 * G + 7471 * B ) >> 16) ) >> 16) + 128)
#define CRGB2Cr(R, G, B) CLIP((46727 * (R - CLIP((19595 * R + 38470 * G + 7471 * B ) >> 16) ) >> 16) + 128)

// YCbCr -> RGB
#define CYCbCr2R(Y, Cb, Cr) CLIP( Y + ( 91881 * Cr >> 16 ) - 179 )
#define CYCbCr2G(Y, Cb, Cr) CLIP( Y - (( 22544 * Cb + 46793 * Cr ) >> 16) + 135)
#define CYCbCr2B(Y, Cb, Cr) CLIP( Y + (116129 * Cb >> 16 ) - 226 )

extern const double xyzd50_srgb[3][3];
extern const double rgb_rgb[3][3];
extern const double adobe_rgb[3][3];
extern const double wide_rgb[3][3];
extern const double prophoto_rgb[3][3];
extern const double aces_rgb[3][3];
extern const double dcip3d65_rgb[3][3];
extern const double rec2020_rgb[3][3];
extern const double xyz_rgb[3][3];
extern const double rgb_cam[3][4];
extern const float d65_white[3];


#ifdef USE_OPENMP
#include "omp.h"
#endif
#ifdef USE_GCD
#include <dispatch/dispatch.h>
#endif


template <typename Func>
void parallel_for_image(int height, int width, Func&& func) {
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    for (int row = 0; row < height; row += 1) {
        for (int col = 0; col < width; col += 1) {
            func(col, row);
        }
    }
#elif defined(USE_GCD)
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_apply(height, queue, ^(size_t row) {
        for (int col = 0; col < width; col += 1) {
            func(col, static_cast<int>(row));
        }
    });
#else
    // 默认单线程执行
    for (int row = 0; row < height; row += 1) {
        for (int col = 0; col < width; col += 1) {
            func(col, row);
        }
    }
#endif
}

#endif // CONST_H
