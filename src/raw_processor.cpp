#include "raw_processor.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <cstring>
#include <algorithm>
#include <mutex>  // 包含 std::mutex
#ifdef USE_JPEG
#include <jpeglib.h>
#endif
#include <setjmp.h>

#ifdef USE_GCD
#include <dispatch/dispatch.h>
#endif

std::mutex log_mutex;  // 定义互斥锁

void RawProcessor::raw2image(RawImage& raw_image) {
    // raw_image.image.resize(raw_image.height * raw_image.width);

    int maxHeight = raw_image.height;
    int maxWidth = raw_image.width;

    int alloc_sz = maxHeight * maxWidth;

    // raw_image.image = (unsigned short(*)[4])realloc(raw_image.image, alloc_sz * sizeof(*raw_image.image));
    raw_image.image = (unsigned short(*)[4])malloc(alloc_sz * sizeof(*raw_image.image));
    memset(raw_image.image, 0, alloc_sz * sizeof(*raw_image.image));

#ifdef USE_GCD
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);
    dispatch_group_t group = dispatch_group_create();

    // 使用原子变量来减少临界区的开销
    __block unsigned short global_max = 0;
    dispatch_semaphore_t semaphore = dispatch_semaphore_create(1);

    int chunk_size = maxHeight / 4; // 根据实际情况调整块大小

    for (int chunk_start = 0; chunk_start < maxHeight; chunk_start += chunk_size) {
        int chunk_end = std::min(chunk_start + chunk_size, maxHeight);

        dispatch_group_async(group, queue, ^{
            unsigned short ldmax = 0;
            for (int row = chunk_start; row < chunk_end; ++row) {
                for (int col = 0; col < maxWidth; ++col) {
                    int idx = row * maxWidth + col;
                    int fc = FC(raw_image.filters, row, col); // 预先计算 FC(row, col)
                    unsigned short val = raw_image.raw_data[idx];
                    if (val > raw_image.black_level) {
                        val -= raw_image.black_level;
                        if (val > ldmax)
                            ldmax = val;
                    } else {
                        val = 0;
                    }
                    raw_image.image[idx][fc] = val;
                }
            }

            // 使用原子操作来更新 global_max
            dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
            if (ldmax > global_max) {
                global_max = ldmax;
            }
            dispatch_semaphore_signal(semaphore);
        });
    }

    dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
    raw_image.data_max = global_max;

    dispatch_release(group);
    dispatch_release(semaphore);
#else
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) default(none) shared(raw_image, maxHeight, maxWidth)
#endif 
    for (int row = 0; row < maxHeight; ++row) {
        unsigned short ldmax = 0;
        for (int col = 0; col < maxWidth; ++col) {
            int idx = row * maxWidth + col;
            int fc = FC(raw_image.filters, row, col); // 预先计算 FC(row, col)
            unsigned short val = raw_image.raw_data[idx];
            if (val > raw_image.black_level) {
                val -= raw_image.black_level;
                if (val > ldmax)
                    ldmax = val;
            } else {
                val = 0;
            }
            raw_image.image[idx][fc] = val;
        }
#ifdef USE_OPENMP
        // 临界区保护 dmax 更新
        #pragma omp critical
#endif
        {
            if (ldmax > raw_image.data_max) {
                raw_image.data_max = ldmax;
            }
        }
    }
#endif
}

// 白平衡处理
void RawProcessor::scale_colors(RawImage& raw_image) {
    float dmin = DBL_MAX, dmax = 0;
    for (const auto mul : raw_image.cam_mul) {
        if (mul < dmin) {
            dmin = mul;
        }
        if (mul > dmax) {
            dmax = mul;
        }
    }
    dmax = dmin;
    unsigned short maximum = raw_image.maximum - raw_image.black_level;

    std::vector<float> scale_mul(4);
    for (int i = 0; i < 4; ++i) {
        // TODO 这里有个问题，有点暗，当前先在这里 * 2
        scale_mul[i] = 2 * (raw_image.cam_mul[i] /= dmax) * 65535.0 / maximum;
    }
    for (int idx = 0; idx < raw_image.height * raw_image.width; ++idx) {
        auto& pix = raw_image.image[idx];
        for (int i = 0; i < 4; ++i) {
            int val = pix[i];
            val *= scale_mul[i];
            pix[i] = CLIP(val);
        }
    }
}

void RawProcessor::pre_interpolate(RawImage& raw_image) {
    for (int idx = 0; idx < raw_image.height * raw_image.width; ++idx) {
        auto& pix = raw_image.image[idx];
        if (pix[3] > 0) {
            pix[1] = pix[3];
        }
    }
    raw_image.filters &= ~((raw_image.filters & 0x55555555U) << 1);
}

void cielab(unsigned short rgb[3], short lab[3]) {
    int c, i, j, k;
    float r, xyz[3];
    static float cbrt[0x10000], xyz_cam[3][4];
    if (!rgb) {
        for (i = 0; i < 0x10000; i++) {   
            r = i / 65535.0;
            cbrt[i] = r > 0.008856 ? std::pow(r, 1.f / 3.0f) : 7.787f * r + 16.f / 116.0f;
        }   
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (xyz_cam[i][j] = k = 0; k < 3; k++) {
                    xyz_cam[i][j] += xyz_rgb[i][k] * rgb_cam[k][j] / d65_white[i];
                }
            }
        }
        return;
    }
    xyz[0] = xyz[1] = xyz[2] = 0.5;
    for (int c = 0; c < 3; ++c) {
        xyz[0] += xyz_cam[0][c] * rgb[c];
        xyz[1] += xyz_cam[1][c] * rgb[c];
        xyz[2] += xyz_cam[2][c] * rgb[c];
    }
    xyz[0] = cbrt[CLIP((int)xyz[0])];
    xyz[1] = cbrt[CLIP((int)xyz[1])];
    xyz[2] = cbrt[CLIP((int)xyz[2])];
    lab[0] = 64 * (116 * xyz[1] - 16);
    lab[1] = 64 * 500 * (xyz[0] - xyz[1]);
    lab[2] = 64 * 200 * (xyz[1] - xyz[2]);
}

void RawProcessor::ahd_interpolate_green_h_and_v(
        RawImage& raw_image, int top, int left, unsigned short (*out_rgb)[AHD_TILE][AHD_TILE][3]) {
    int width = raw_image.width;
    int height = raw_image.height;
    int row, col;
    int c, val;
    unsigned short(*pix)[4];
    const int rowlimit = MIN(top + AHD_TILE, raw_image.height - 2);
    const int collimit = MIN(left + AHD_TILE, raw_image.width - 2);

    for (row = top; row < rowlimit; row++) {
        col = left + (FC(raw_image.filters, row, left) & 1);
        for (c = FC(raw_image.filters, row, col); col < collimit; col += 2) {
            pix = raw_image.image + row * width + col;
            val = ((pix[-1][1] + pix[0][c] + pix[1][1]) * 2 - pix[-2][c] - pix[2][c]) >> 2;
            out_rgb[0][row - top][col - left][1] = ULIM(val, pix[-1][1], pix[1][1]);
            val = ((pix[-width][1] + pix[0][c] + pix[width][1]) * 2 - pix[-2 * width][c] - pix[2 * width][c]) >> 2;
            out_rgb[1][row - top][col - left][1] = ULIM(val, pix[-width][1], pix[width][1]);
        }
    }
}

/*
 * 图像边界处的像素进行插值，来平滑这些边界，以减少图像处理中的边界效应。
 */
void RawProcessor::border_interpolate(RawImage& raw_image, int border) {
    int sum[8];
    for (int row = 0; row < raw_image.height; ++row) {
        for (int col = 0; col < raw_image.width; ++col) {
            if (col == border && row >= border && row < (raw_image.height - border)) {
                col = raw_image.width - border;
            }
            memset(sum, 0, sizeof sum);
            for (int y = row - 1; y != row + 2; ++y) {
                for (int x = col - 1; x != col + 2; ++x) {
                    if (y < raw_image.height && x < raw_image.width) {
                        int f = FC(raw_image.filters, y, x);
                        sum[f] += raw_image.image[y * raw_image.width + x][f];
                        sum[f + 4]++;
                    }
                }
            }
            int f = FC(raw_image.filters, row, col);
            for (int c = 0; c < 4; ++c) {
                if (c != f && sum[c + 4]) {
                    raw_image.image[row * raw_image.width + col][c] = sum[c] / sum[c + 4];
                }
            }
        }
    }
}

void RawProcessor::ahd_interpolate_r_and_b_and_convert_to_cielab(
        RawImage& raw_image, int top, int left, unsigned short (*inout_rgb)[AHD_TILE][AHD_TILE][3],
        short (*out_lab)[AHD_TILE][AHD_TILE][3]) {
    for (int direction = 0; direction < 2; direction++) {
        ahd_interpolate_r_and_b_in_rgb_and_convert_to_cielab(raw_image, top, left, inout_rgb[direction], out_lab[direction]);
    }
}

void RawProcessor::ahd_interpolate_r_and_b_in_rgb_and_convert_to_cielab(
        RawImage& raw_image, int top, int left, 
        unsigned short (*inout_rgb)[AHD_TILE][3], short (*out_lab)[AHD_TILE][3]) {
    int height = raw_image.height;
    int width = raw_image.width;
    auto image = raw_image.image;
    unsigned row, col;
    int c, val;
    unsigned short(*pix)[4];
    unsigned short(*rix)[3];
    short(*lix)[3];
    const unsigned num_pix_per_row = 4 * width;
    const unsigned rowlimit = MIN(top + AHD_TILE - 1, height - 3);
    const unsigned collimit = MIN(left + AHD_TILE - 1, width - 3);
    unsigned short *pix_above;
    unsigned short *pix_below;
    int t1, t2;

    for (row = top + 1; row < rowlimit; row++) { 
        pix = image + row * width + left;
        rix = &inout_rgb[row - top][0];
        lix = &out_lab[row - top][0];

        for (col = left + 1; col < collimit; col++) {
            pix++; 
            pix_above = &pix[0][0] - num_pix_per_row;
            pix_below = &pix[0][0] + num_pix_per_row;
            rix++;
            lix++;

            c = 2 - FC(raw_image.filters, row, col);
            if (c == 1) {
                c = FC(raw_image.filters, row + 1, col);
                t1 = 2 - c;
                val = pix[0][1] + ((pix[-1][t1] + pix[1][t1] - rix[-1][1] - rix[1][1]) >> 1);
                rix[0][t1] = CLIP(val);
                val = pix[0][1] + ((pix_above[c] + pix_below[c] - rix[-AHD_TILE][1] - rix[AHD_TILE][1]) >> 1);
            } else {
                t1 = -4 + c; /* -4+c: pixel of color c to the left */
                t2 = 4 + c;  /* 4+c: pixel of color c to the right */
                val = rix[0][1] + ((pix_above[t1] + pix_above[t2] + pix_below[t1] + pix_below[t2] -
                            rix[-AHD_TILE - 1][1] - rix[-AHD_TILE + 1][1] -
                            rix[+AHD_TILE - 1][1] - rix[+AHD_TILE + 1][1] + 1) >> 2);
            }
            rix[0][c] = CLIP(val);
            c = FC(raw_image.filters, row, col);
            rix[0][c] = pix[0][c];
            cielab(rix[0], lix[0]);
        }
    }
}

void RawProcessor::ahd_interpolate_build_homogeneity_map(
        RawImage& raw_image, int top, int left, short (*lab)[AHD_TILE][AHD_TILE][3],
        char (*out_homogeneity_map)[AHD_TILE][2]) {
    int height = raw_image.height;
    int width = raw_image.width;
    int row, col;
    int tr;
    int direction;
    int i;
    short(*lix)[3];
    short(*lixs[2])[3];
    short *adjacent_lix;
    unsigned ldiff[2][4], abdiff[2][4], leps, abeps;
    static const int dir[4] = {-1, 1, -AHD_TILE, AHD_TILE};
    const int rowlimit = MIN(top + AHD_TILE - 2, height - 4);
    const int collimit = MIN(left + AHD_TILE - 2, width - 4);
    int homogeneity;
    char(*homogeneity_map_p)[2];

    memset(out_homogeneity_map, 0, 2 * AHD_TILE * AHD_TILE);

    for (row = top + 2; row < rowlimit; row++) {
        tr = row - top;
        homogeneity_map_p = &out_homogeneity_map[tr][1];
        for (direction = 0; direction < 2; direction++) {
            lixs[direction] = &lab[direction][tr][1];
        } 

        for (col = left + 2; col < collimit; col++) {
            homogeneity_map_p++; 
            for (direction = 0; direction < 2; direction++) {
                lix = ++lixs[direction];
                for (i = 0; i < 4; i++) {
                    adjacent_lix = lix[dir[i]];
                    ldiff[direction][i] = ABS(lix[0][0] - adjacent_lix[0]);
                    abdiff[direction][i] = SQR(lix[0][1] - adjacent_lix[1]) + SQR(lix[0][2] - adjacent_lix[2]);
                }
            }
            leps = MIN(MAX(ldiff[0][0], ldiff[0][1]), MAX(ldiff[1][2], ldiff[1][3]));
            abeps = MIN(MAX(abdiff[0][0], abdiff[0][1]), MAX(abdiff[1][2], abdiff[1][3]));
            for (direction = 0; direction < 2; direction++) {
                homogeneity = 0;
                for (i = 0; i < 4; i++) {
                    if (ldiff[direction][i] <= leps && abdiff[direction][i] <= abeps) {
                        homogeneity++;
                    }
                }
                homogeneity_map_p[0][direction] = homogeneity;
            }
        }
    }
}

void RawProcessor::ahd_interpolate_combine_homogeneous_pixels(
        RawImage& raw_image, int top, int left, ushort (*rgb)[AHD_TILE][AHD_TILE][3],
        char (*homogeneity_map)[AHD_TILE][2]) {
    int height = raw_image.height;
    int width = raw_image.width;
    auto image = raw_image.image;
    int row, col;
    int tr, tc;
    int i, j;
    int direction;
    int hm[2];
    int c;
    const int rowlimit = MIN(top + AHD_TILE - 3, height - 5);
    const int collimit = MIN(left + AHD_TILE - 3, width - 5);

    unsigned short(*pix)[4];
    unsigned short(*rix[2])[3];

    for (row = top + 3; row < rowlimit; row++) {
        tr = row - top;
        pix = &image[row * width + left + 2];
        for (direction = 0; direction < 2; direction++) {
            rix[direction] = &rgb[direction][tr][2];
        }

        for (col = left + 3; col < collimit; col++) {
            tc = col - left;
            pix++;
            for (direction = 0; direction < 2; direction++) {
                rix[direction]++;
            }
            for (direction = 0; direction < 2; direction++) {
                hm[direction] = 0;
                for (i = tr - 1; i <= tr + 1; i++) {
                    for (j = tc - 1; j <= tc + 1; j++) {
                        hm[direction] += homogeneity_map[i][j][direction];
                    }
                }
            }
            if (hm[0] != hm[1]) {
                memcpy(pix[0], rix[hm[1] > hm[0]][0], 3 * sizeof(ushort));
            } else {
                for (c = 0; c < 3; ++c) { 
                    pix[0][c] = (rix[0][0][c] + rix[1][0][c]) >> 1; 
                }
            }
        }
    }
}

void RawProcessor::ahd_interpolate(RawImage& raw_image) {
    int terminate_flag = 0;
    cielab(0, 0);
    border_interpolate(raw_image, 5);

#ifdef USE_GCD
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_apply((raw_image.height - 7) / (AHD_TILE - 6) + 1, queue, ^(size_t index) {
        int top = 2 + index * (AHD_TILE - 6);
        if (top >= raw_image.height - 5) return;
        std::vector<char> thread_local_buffer(26 * AHD_TILE * AHD_TILE);
        char* buffer = thread_local_buffer.data();
        unsigned short(*rgb)[AHD_TILE][AHD_TILE][3];
        short(*lab)[AHD_TILE][AHD_TILE][3];
        char(*homo)[AHD_TILE][2];

        rgb = (unsigned short(*)[AHD_TILE][AHD_TILE][3])buffer;
        lab = (short(*)[AHD_TILE][AHD_TILE][3])(buffer + 12 * AHD_TILE * AHD_TILE);
        homo = (char(*)[AHD_TILE][2])(buffer + 24 * AHD_TILE * AHD_TILE);

        for (int left = 2; !terminate_flag && (left < raw_image.width - 5); left += AHD_TILE - 6) {
            int block_width = std::min(AHD_TILE, (int)raw_image.width - left - 2);
            int block_height = std::min(AHD_TILE, (int)raw_image.height - top - 2);
            ahd_interpolate_green_h_and_v(raw_image, top, left, rgb);
            ahd_interpolate_r_and_b_and_convert_to_cielab(raw_image, top, left, rgb, lab);
            ahd_interpolate_build_homogeneity_map(raw_image, top, left, lab, homo);
            ahd_interpolate_combine_homogeneous_pixels(raw_image, top, left, rgb, homo);
        }
    });
#else

#ifdef USE_OPENMP
    int buffer_count = omp_get_max_threads();
#else
    int buffer_count = 1;
#endif
    size_t buffer_size = 26 * AHD_TILE * AHD_TILE;
    char** buffers = (char**)calloc(sizeof(char*), buffer_count); 
    for (int i = 0; i < buffer_count; i++) {   
        buffers[i] = (char*)calloc(buffer_size, sizeof(char));
    }   

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) default(none) shared(raw_image, terminate_flag) firstprivate(buffers)
#endif
    for (int top = 2; top < raw_image.height - 5; top += AHD_TILE - 6) {
#ifdef USE_OPENMP
        char* buffer = buffers[omp_get_thread_num()];
#else
        char* buffer = buffers[0];
#endif
        unsigned short(*rgb)[AHD_TILE][AHD_TILE][3];
        short(*lab)[AHD_TILE][AHD_TILE][3];
        char(*homo)[AHD_TILE][2];

        rgb = (unsigned short(*)[AHD_TILE][AHD_TILE][3])buffer;
        lab = (short(*)[AHD_TILE][AHD_TILE][3])(buffer + 12 * AHD_TILE * AHD_TILE);
        homo = (char(*)[AHD_TILE][2])(buffer + 24 * AHD_TILE * AHD_TILE);
        for (int left = 2; !terminate_flag && (left < raw_image.width - 5); left += AHD_TILE - 6) {
            ahd_interpolate_green_h_and_v(raw_image, top, left, rgb);
            ahd_interpolate_r_and_b_and_convert_to_cielab(raw_image, top, left, rgb, lab);
            ahd_interpolate_build_homogeneity_map(raw_image, top, left, lab, homo);
            ahd_interpolate_combine_homogeneous_pixels(raw_image, top, left, rgb, homo);
        }
    }
    for (int i = 0; i < buffer_count; i++) {
        if (buffers[i]) {
            free(buffers[i]);
        }
    }
    free(buffers);
#endif
}

void gamma_curve(RawImage& raw_image, double pwr, double ts, int mode, int imax) {
    int i;
    double g[6], bnd[2] = {0, 0}, r;

    g[0] = pwr;
    g[1] = ts;
    g[2] = g[3] = g[4] = 0;
    bnd[g[1] >= 1] = 1;
    if (g[1] && (g[1] - 1) * (g[0] - 1) <= 0) {
        for (i = 0; i < 48; i++) {
            g[2] = (bnd[0] + bnd[1]) / 2;
            if (g[0]) {
                bnd[(std::pow(g[2] / g[1], -g[0]) - 1) / g[0] - 1 / g[2] > -1] = g[2];
            } else {
                bnd[g[2] / std::exp(1 - 1 / g[2]) < g[1]] = g[2];
            }
        }
        g[3] = g[2] / g[1];
        if (g[0]) {
            g[4] = g[2] * (1 / g[0] - 1);
        }
    }
    if (g[0]) {
        g[5] = 1 / (g[1] * SQR(g[3]) / 2 - g[4] * (1 - g[3]) + (1 - std::pow(g[3], 1 + g[0])) * (1 + g[4]) / (1 + g[0])) - 1;
    } else {
        g[5] = 1 / (g[1] * SQR(g[3]) / 2 + 1 - g[2] - g[3] - g[2] * g[3] * (std::log(g[3]) - 1)) - 1;
    }
    if (!mode--) { 
        memcpy(raw_image.gamm, g, sizeof raw_image.gamm);
        return;
    }
    for (i = 0; i < 0x10000; i++) {
        raw_image.curve[i] = 0xffff;
        if ((r = (double)i / imax) < 1) {
            raw_image.curve[i] = 0x10000 * (mode ? (r < g[3] ? r * g[1] 
                        : (g[0] ? pow(r, g[0]) * (1 + g[4]) - g[4] : log(r) * g[2] + 1))
                    : (r < g[2] ? r / g[1]
                        : (g[0] ? std::pow((r + g[4]) / (1 + g[4]), 1 / g[0]) : exp((r - 1) / g[2]))));
        }
    } 
}

void pseudoinverse(double (*in)[3], double (*out)[3], int size) {
    double work[3][6], num;
    int i, j, k;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 6; j++) {
            work[i][j] = j == i + 3;
        }
        for (j = 0; j < 3; j++) {
            for (k = 0; k < size && k < 4; k++) {
                work[i][j] += in[k][i] * in[k][j]; 
            }
        }
    }

    for (i = 0; i < 3; i++) { 
        num = work[i][i];
        for (j = 0; j < 6; j++) {
            if (fabs(num) > 0.00001f) {
                work[i][j] /= num;
            }
        }
        for (k = 0; k < 3; k++) {
            if (k == i) {
                continue;
            }
            num = work[k][i];
            for (j = 0; j < 6; j++) {
                work[k][j] -= work[i][j] * num;
            }
        }   
    }
    for (i = 0; i < size && i < 4; i++) {
        for (j = 0; j < 3; j++) {
            for (out[i][j] = k = 0; k < 3; k++) {
                out[i][j] += work[j][k + 3] * in[i][k];
            }
        }
    }
}

void convert_to_rgb_loop(RawImage& raw_image, float out_cam[3][4]) {
    int row, col, c;
    float out[3];
    ushort *img;
    if (!raw_image.histogram) {   
        raw_image.histogram = (int(*)[HISTOGRAM_SIZE])malloc(sizeof(int) * HISTOGRAM_SIZE * 4); 
        if (!raw_image.histogram) {
            std::cout << "Memory allocation failed." << std::endl;
            // Handle the error, possibly by exiting the program or by not proceeding further.
        }
    }
    memset(raw_image.histogram, 0, sizeof(int) * HISTOGRAM_SIZE * 4);
    for (img = raw_image.image[0], row = 0; row < raw_image.height; row++) {
        for (col = 0; col < raw_image.width; col++, img += 4) {
            out[0] = out_cam[0][0] * img[0] + out_cam[0][1] * img[1] +
                out_cam[0][2] * img[2];
            out[1] = out_cam[1][0] * img[0] + out_cam[1][1] * img[1] +
                out_cam[1][2] * img[2];
            out[2] = out_cam[2][0] * img[0] + out_cam[2][1] * img[1] +
                out_cam[2][2] * img[2];
            img[0] = CLIP((int)out[0]);
            img[1] = CLIP((int)out[1]);
            img[2] = CLIP((int)out[2]);
            raw_image.histogram[0][img[0] >> 3]++;
            raw_image.histogram[1][img[1] >> 3]++;
            raw_image.histogram[2][img[2] >> 3]++;
        }
    }
}

void RawProcessor::convert_to_rgb(RawImage& raw_image) {
    float out_cam[3][4];
    double num, inverse[3][3];
    static const double(*out_rgb[])[3] = {
        rgb_rgb, adobe_rgb, wide_rgb, prophoto_rgb, xyz_rgb, aces_rgb, dcip3d65_rgb, rec2020_rgb
    };
    static const char *name[] = {
        "sRGB", "Adobe RGB (1998)", "WideGamut D65", "ProPhoto D65", "XYZ", "ACES", "DCI-P3 D65", "Rec. 2020"
    };
    static const unsigned phead[] = {
        1024, 0, 0x2100000,  0x6d6e7472, 0x52474220, 0x58595a20, 0,
        0,    0, 0x61637370, 0,          0,          0x6e6f6e65, 0,
        0,    0, 0,          0xf6d6,     0x10000,    0xd32d};
    unsigned pbody[] = {
        10,         0x63707274, 0,  36, 
        0x64657363, 0,          60,     
        0x77747074, 0,          20,   
        0x626b7074, 0,          20,  
        0x72545243, 0,          14, 
        0x67545243, 0,          14,
        0x62545243, 0,          14, 
        0x7258595a, 0,          20, 
        0x6758595a, 0,          20, 
        0x6258595a, 0,          20};
    static const unsigned pwhite[] = {0xf351, 0x10000, 0x116cc};
    unsigned pcurve[] = {0x63757276, 0, 1, 0x1000000};

    auto& oprof = raw_image.oprof;
    auto& gamm = raw_image.gamm;

    gamma_curve(raw_image, raw_image.gamm[0], raw_image.gamm[1], 0, 0);

    memcpy(out_cam, rgb_cam, sizeof out_cam);
	size_t prof_desc_len;
	std::vector<char> prof_desc;
    int i, j, k;
    // TODO
    int output_color = 1;

    prof_desc_len = snprintf(NULL, 0, "%s gamma %g toe slope %g", name[output_color - 1],
            floorf(1000.f / gamm[0] + .5f) / 1000.f, floorf(gamm[1] * 1000.0f + .5f) / 1000.f) + 1;
    prof_desc.resize(prof_desc_len);
    sprintf(prof_desc.data(), "%s gamma %g toe slope %g", name[output_color - 1], 
            floorf(1000.f / gamm[0] + .5f) / 1000.f,
            floorf(gamm[1] * 1000.0f + .5f) / 1000.f);

	raw_image.oprof = (unsigned *)calloc(phead[0], 1);
    memcpy(raw_image.oprof, phead, sizeof phead);
    oprof[0] = 132 + 12 * pbody[0];
    for (i = 0; i < (int)pbody[0]; i++) {
        oprof[oprof[0] / 4] = i ? (i > 1 ? 0x58595a20 : 0x64657363) : 0x74657874;
        pbody[i * 3 + 2] = oprof[0];
        oprof[0] += (pbody[i * 3 + 3] + 3) & -4;
    }
    memcpy(oprof + 32, pbody, sizeof pbody);
    oprof[pbody[5] / 4 + 2] = unsigned(prof_desc_len + 1);
    memcpy((char *)oprof + pbody[8] + 8, pwhite, sizeof pwhite);
    pcurve[3] = (short)(256 / gamm[5] + 0.5) << 16;
    for (i = 4; i < 7; i++) {
        memcpy((char *)oprof + pbody[i * 3 + 2], pcurve, sizeof pcurve);
    }
    pseudoinverse((double(*)[3])out_rgb[output_color - 1], inverse, 3);
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (num = k = 0; k < 3; k++) {
                num += xyzd50_srgb[i][k] * inverse[j][k];
            }
            oprof[pbody[j * 3 + 23] / 4 + i + 2] = num * 0x10000 + 0.5;
        }
    }
    for (i = 0; i < (int)phead[0] / 4; i++) {
        oprof[i] = htonl(oprof[i]);
    }
    strcpy((char *)oprof + pbody[2] + 8, "auto-generated by dcraw");
    if (pbody[5] + 12 + prof_desc.size() < phead[0]) {
        strcpy((char *)oprof + pbody[5] + 12, prof_desc.data());
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (out_cam[i][j] = k = 0; k < 3; k++) {
                    out_cam[i][j] += out_rgb[output_color - 1][i][k] * rgb_cam[k][j];
            }
        }
    }
    convert_to_rgb_loop(raw_image, out_cam);
}

void tiff_set(struct new_tiff_hdr *th, ushort *ntag, ushort tag, ushort type, int count, int val) {
    struct tiff_tag *tt;
    int c;

    tt = (struct tiff_tag *)(ntag + 1) + (*ntag)++;
    tt->val.i = val;
    if (type == 1 && count <= 4) {
        for (int c = 0; c < 4; ++c) {
            tt->val.c[c] = val >> (c << 3);
        }
    } else if (type == 2) {
        count = int(strnlen((char *)th + val, count - 1)) + 1;
        if (count <= 4) {
            for (int c = 0; c < 4; ++c) {
                tt->val.c[c] = ((char *)th)[val + c];
            }
        }
    } else if (type == 3 && count <= 2) {
        for (int c = 0; c < 2; ++c) {
            tt->val.s[c] = val >> (c << 4);
        }
    }
    tt->count = count;
    tt->type = type;
    tt->tag = tag;
}

#define TOFF(ptr) ((char *)(&(ptr)) - (char *)th)

void tiff_head(RawImage& raw_image, struct new_tiff_hdr *th, int full) {
    int c, psize = 0;
    struct tm *t;
    int output_bps = raw_image.output_bps;
    auto& oprof = raw_image.oprof;

    memset(th, 0, sizeof *th);
    th->t_order = htonl(0x4d4d4949) >> 16;
    th->magic = 42;
    th->ifd = 10;
    th->rat[0] = th->rat[2] = 300;
    th->rat[1] = th->rat[3] = 1;
    for (int c = 0; c < 6; ++c) {
        th->rat[4 + c] = 1000000;
    }
    th->rat[4] *= 1;
    th->rat[6] *= 1;
    th->rat[8] *= 1;
    strncpy(th->t_desc, "", 512);
    strncpy(th->t_make, "", 64);
    strncpy(th->t_model, "", 64);
    strcpy(th->soft, "dcraw v");
    sprintf(th->date, "%04d:%02d:%02d %02d:%02d:%02d", 2024, 12, 13, 21, 55, 21);
    strncpy(th->t_artist, "", 64);
    if (full) {
        tiff_set(th, &th->ntag, 254, 4, 1, 0);
        tiff_set(th, &th->ntag, 256, 4, 1, raw_image.width);
        tiff_set(th, &th->ntag, 257, 4, 1, raw_image.height);
        tiff_set(th, &th->ntag, 258, 3, 3, output_bps);
        th->tag[th->ntag - 1].val.i = TOFF(th->bps);
        for (int c = 0; c < 4; ++c) {
            th->bps[c] = output_bps;
        }
        tiff_set(th, &th->ntag, 259, 3, 1, 1);
        tiff_set(th, &th->ntag, 262, 3, 1, 2);
    }
    tiff_set(th, &th->ntag, 270, 2, 512, TOFF(th->t_desc));
    tiff_set(th, &th->ntag, 271, 2, 64, TOFF(th->t_make));
    tiff_set(th, &th->ntag, 272, 2, 64, TOFF(th->t_model));
        if (oprof) {
            // TODO
            psize = ntohl(oprof[0]);
            // psize = 496;
        }
    tiff_set(th, &th->ntag, 273, 4, 1, sizeof *th + psize);
    tiff_set(th, &th->ntag, 277, 3, 1, 3);
    tiff_set(th, &th->ntag, 278, 4, 1, raw_image.height);
    tiff_set(th, &th->ntag, 279, 4, 1, raw_image.height * raw_image.width * 3 * output_bps / 8);
    tiff_set(th, &th->ntag, 282, 5, 1, TOFF(th->rat[0]));
    tiff_set(th, &th->ntag, 283, 5, 1, TOFF(th->rat[2]));
    tiff_set(th, &th->ntag, 284, 3, 1, 1);
    tiff_set(th, &th->ntag, 296, 3, 1, 2);
    tiff_set(th, &th->ntag, 305, 2, 32, TOFF(th->soft));
    tiff_set(th, &th->ntag, 306, 2, 20, TOFF(th->date));
    tiff_set(th, &th->ntag, 315, 2, 64, TOFF(th->t_artist));
    tiff_set(th, &th->ntag, 34665, 4, 1, TOFF(th->nexif));
    if (psize) {
        tiff_set(th, &th->ntag, 34675, 7, psize, sizeof *th);
    }
    tiff_set(th, &th->nexif, 33434, 5, 1, TOFF(th->rat[4]));
    tiff_set(th, &th->nexif, 33437, 5, 1, TOFF(th->rat[6]));
    tiff_set(th, &th->nexif, 34855, 3, 1, 1);
    tiff_set(th, &th->nexif, 37386, 5, 1, TOFF(th->rat[8]));
}

void RawProcessor::gamma_adjustment(RawImage& raw_image) {
    int c;
    int perc, val, total, t_white = 0x2000;
    int height = raw_image.height;
    int width = raw_image.width;

    auto& histogram = raw_image.histogram;
    auto& gamm = raw_image.gamm;

    perc = width * height * 0.01;
    // TODO: 这里有有个auto_bright的问题,注释的是auto birght       
    /*
    for (t_white = c = 0; c < 3; c++) {
        for (val = 0x2000, total = 0; --val > 32;) {
            std::cout << val << std::endl;
            if ((total += histogram[c][val]) > perc) {
                break;
            }
        }
        if (t_white < val) {
            t_white = val;
        }
    }*/
    gamma_curve(raw_image, gamm[0], gamm[1], 2, int((t_white << 3) / 1.));

    raw_image.rgb_image = (unsigned short(*)[3])malloc(height * width * sizeof(*raw_image.rgb_image));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            int pix = row * width + col;
            for (int c = 0; c < 3; ++c) {
                raw_image.rgb_image[pix][c] = raw_image.curve[raw_image.image[pix][c]];
            }
        }
    }
    std::cout << raw_image.rgb_image[657 * 8280 + 4394][1] << std::endl;
}

void write_ppm_tiff(RawImage& raw_image) {
    try { 
        struct new_tiff_hdr th;
        ushort *ppm2;
        int c, row, col, soff, rstep, cstep;
        int height = raw_image.height;
        int width = raw_image.width;

        auto& curve = raw_image.curve;
        auto& oprof = raw_image.oprof;

        int output_bps = raw_image.output_bps;
        std::vector<unsigned char> ppm(width * 4 * output_bps / 8);
        ppm2 = (unsigned short *)ppm.data();
        tiff_head(raw_image, &th, 1);
        fwrite(&th, sizeof th, 1, raw_image.output);
        if (oprof) {
            fwrite(oprof, ntohl(oprof[0]), 1, raw_image.output);
        }
        soff = 0;
        cstep = 1;
        rstep = 0;
        for (row = 0; row < height; row++, soff += rstep) {
            for (col = 0; col < width; col++, soff += cstep) {
                for (int c = 0; c < 3; ++c) {
                    // ppm[col * 3 + c] = curve[raw_image.image[soff][c]] >> 8;
                    //ppm[col * 3 + c] = raw_image.rgb_image[soff][c];
                    if (output_bps == 8) {
                        ppm[col * 3 + c] = raw_image.rgb_image[soff][c] >> 8;
                    } else {
                        ppm2[col * 3 + c] = raw_image.rgb_image[soff][c];
                    }
                }
            }
            if (output_bps == 8) {
                fwrite(ppm.data(), 3 * output_bps / 8, width, raw_image.output);
            } else {
                fwrite(ppm2, 3 * output_bps / 8, width, raw_image.output);
            }
        }
    } catch (...) {
    }
}

int RawProcessor::ppm_tiff_writer(RawImage& raw_image, const char *filename) {
    if (!filename)
        return -1;
    FILE *f = NULL;
    if (!strcmp(filename, "-")) {
        f = stdout;
    } else {
        f = fopen(filename, "wa");
    }
    raw_image.output = f;
    write_ppm_tiff(raw_image);
    raw_image.output = NULL;
    fclose(f);
    return 0;
}

#ifdef USE_JPEG
void RawProcessor::write_jpeg(RawImage& raw_image, const char* filename) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr; 
    int width = raw_image.width;
    int height = raw_image.height;
    
    FILE* outfile;
    if ((outfile = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        return;
    } 
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE); // 质量值为75
    jpeg_start_compress(&cinfo, TRUE);
    JSAMPROW row_pointer[1];
    unsigned char* row_data = new unsigned char[width * 3];

    while (cinfo.next_scanline < cinfo.image_height) {
        // 将16位数据缩放到8位范围
        for (int i = 0; i < width; ++i) {
            row_data[i * 3 + 0] = static_cast<unsigned char>(raw_image.rgb_image[cinfo.next_scanline * raw_image.width + i][0] >> 8);
            row_data[i * 3 + 1] = static_cast<unsigned char>(raw_image.rgb_image[cinfo.next_scanline * raw_image.width + i][1] >> 8);
            row_data[i * 3 + 2] = static_cast<unsigned char>(raw_image.rgb_image[cinfo.next_scanline * raw_image.width + i][2] >> 8);
        }
        row_pointer[0] = row_data;
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    delete[] row_data; // 释放动态分配的内存
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}
#endif
