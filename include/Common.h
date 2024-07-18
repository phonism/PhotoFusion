#ifndef CONST_H
#define CONST_H

#define AHD_TILE 512
#define HISTOGRAM_SIZE 0x2000

#define FC(filters, row, col) (((filters) >> ((((row) << 1 & 14) | ((col) & 1)) << 1)) & 3)
#define SQR(x) ((x) * (x))
#define ABS(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define LIM(x, min, max) MAX(min, MIN(x, max))
#define ULIM(x, y, z) ((y) < (z) ? LIM(x, y, z) : LIM(x, z, y))
#define CLIP(x) LIM((int)(x), 0, 65535)

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

#endif // CONST_H
