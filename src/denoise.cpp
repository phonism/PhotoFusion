#include "denoise.h"
#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>


// 高斯模糊核生成
std::vector<std::vector<double>> generate_gaussian_kernel(int kernel_size, double sigma) {
    std::vector<std::vector<double>> kernel(kernel_size, std::vector<double>(kernel_size));
    double sum = 0.0;
    int radius = kernel_size / 2;
    double s = 2.0 * sigma * sigma;

    for (int x = -radius; x <= radius; ++x) {
        for (int y = -radius; y <= radius; ++y) {
            double r = sqrt(x * x + y * y);
            kernel[x + radius][y + radius] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += kernel[x + radius][y + radius];
        }
    }
    for (int i = 0; i < kernel_size; ++i) {
        for (int j = 0; j < kernel_size; ++j) {
            kernel[i][j] /= sum;
        }
    }
    return kernel;
}

// 应用高斯模糊
void apply_gaussian_blur_channel(RawImage* raw_image, int channel, int kernel_size, double sigma) {
    int height = raw_image->height;
    int width = raw_image->width;
    auto image = raw_image->yuv_image;
    std::vector<std::vector<double>> kernel = generate_gaussian_kernel(kernel_size, sigma);
    int radius = kernel_size / 2;
    std::vector<unsigned short> temp(width * height);

#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double sum = 0.0;
            for (int ky = -radius; ky <= radius; ++ky) {
                for (int kx = -radius; kx <= radius; ++kx) {
                    int ny = std::min(std::max(y + ky, 0), height - 1);
                    int nx = std::min(std::max(x + kx, 0), width - 1);
                    sum += image[ny * width + nx][channel] * kernel[ky + radius][kx + radius];
                }
            }
            temp[y * width + x] = static_cast<unsigned short>(sum);
        }
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            raw_image->denoised_image[y * width + x][channel] = temp[y * width + x];
        }
    }
}

void apply_bilateral_filter(RawImage* raw_image, int channel, int kernel_size, double sigma_color, double sigma_space) {
    int width = raw_image->width;
    int height = raw_image->height;
    auto& yuv_image = raw_image->yuv_image;
    int radius = kernel_size / 2;
    std::vector<unsigned short> temp(width * height);

    // 预计算空间高斯核
    std::vector<std::vector<double>> spatial_kernel(kernel_size, std::vector<double>(kernel_size));
    for (int ky = -radius; ky <= radius; ++ky) {
        for (int kx = -radius; kx <= radius; ++kx) {
            spatial_kernel[ky + radius][kx + radius] = std::exp(-(ky * ky + kx * kx) / (2 * sigma_space * sigma_space));
        }
    }


    // 双边滤波
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double sum = 0.0;
            double weight_sum = 0.0;
            unsigned short center_value = yuv_image[y * width + x][channel];

            for (int ky = -radius; ky <= radius; ++ky) {
                for (int kx = -radius; kx <= radius; ++kx) {
                    int ny = std::min(std::max(y + ky, 0), height - 1);
                    int nx = std::min(std::max(x + kx, 0), width - 1);
                    unsigned short neighbor_value = yuv_image[ny * width + nx][channel];

                    double color_dist = (center_value - neighbor_value) * (center_value - neighbor_value) / (2 * sigma_color * sigma_color);
                    double spatial_weight = spatial_kernel[ky + radius][kx + radius];
                    double color_weight = std::exp(-color_dist);
                    double weight = spatial_weight * color_weight;

                    sum += neighbor_value * weight;
                    weight_sum += weight;
                }
            }
            temp[y * width + x] = static_cast<unsigned short>(sum / weight_sum);
        }
    }

    // 更新图像
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            raw_image->denoised_image[y * width + x][channel] = temp[y * width + x];
        }
    }
}

void convert_rgb_to_yuv(RawImage* raw_image) {
    int width = raw_image->width;
    int height = raw_image->height;

    raw_image->yuv_image = (unsigned short(*)[3])malloc(height * width * sizeof(*raw_image->yuv_image));

#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int row = 0; row < height; ++row) {
        for (int col = 0; col < width; ++col) { 
            int pix = row * width + col;
            unsigned short r = raw_image->rgb_image[pix][0];
            unsigned short g = raw_image->rgb_image[pix][1];
            unsigned short b = raw_image->rgb_image[pix][2];
            raw_image->yuv_image[pix][0] = RGB2Y(r, g, b);
            raw_image->yuv_image[pix][1] = RGB2U(r, g, b);
            raw_image->yuv_image[pix][2] = RGB2V(r, g, b);
        }
    }
}

void convert_yuv_to_rgb(RawImage* raw_image) {
    int width = raw_image->width;
    int height = raw_image->height;

#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int row = 0; row < height; ++row) {
        for (int col = 0; col < width; ++col) {
            int pix = row * width + col;
            unsigned short y = raw_image->denoised_image[pix][0];
            unsigned short u = raw_image->denoised_image[pix][1];
            unsigned short v = raw_image->denoised_image[pix][2];

            unsigned short r = YUV2R(y, u, v);
            unsigned short g = YUV2G(y, u, v);
            unsigned short b = YUV2B(y, u, v);

            raw_image->rgb_image[pix][0] = r;
            raw_image->rgb_image[pix][1] = g;
            raw_image->rgb_image[pix][2] = b;
        }
    }
}

void compute_integral_image(RawImage* raw_image) {
    int width = raw_image->width;
    int height = raw_image->height;
    auto integral = raw_image->integral;
    auto yuv_image = raw_image->yuv_image;
    for (int c = 0; c < 3; ++c) {
        integral[0][c] = yuv_image[0][c];
    }
    for (int i = 1; i < height; ++i) {
        for (int c = 0; c < 3; ++c) {
            integral[i * width][c] = integral[(i - 1) * width][c] + yuv_image[i * width][c];
        }
    }
    for (int i = 1; i < width; ++i) {
        for (int c = 0; c < 3; ++c) {
            integral[i][c] = integral[i - 1][c] + yuv_image[i][c];
        }
    }
    for (int i = 1; i < height; ++i) {
        for (int j = 1; j < width; ++j) {
            for (int c = 0; c < 3; ++c) {
                int cur_idx = i * width + j;
                int up_idx = (i - 1) * width + j;
                int left_idx = i * width + j - 1;
                int ul_idx = (i - 1) * width + j - 1;
                integral[cur_idx][c] = integral[up_idx][c] + integral[left_idx][c] - integral[ul_idx][c] + yuv_image[cur_idx][c];
            }
        }
    }
}

inline uint64_t get_block_sum(RawImage* raw_image, int c, int x1, int y1, int x2, int y2) {
    int width = raw_image->width;
    int height = raw_image->height;
    auto integral = raw_image->integral;
    uint64_t A = 0, B = 0, C = 0, D;

    if (x1 > 0 && y1 > 0) {
        A = integral[(y1 - 1) * width + (x1 - 1)][c];
    }
    if (y1 > 0) {
        B = integral[(y1 - 1) * width + x2][c];
    }
    if (x1 > 0) {
        C = integral[y2 * width + (x1 - 1)][c];
    }
    D = integral[y2 * width + x2][c];

    uint64_t block_sum = D - B - C + A;
    return block_sum;
}

int sum_max = 256 * 7 * 7;
int sum_min = -sum_max;
std::vector<float> gaussian_cache(sum_max * 2 + 1);

void initialize_gaussian_cache(int h) {
    for (int sum = sum_min; sum <= sum_max; sum += 1) {
        int idx = sum - sum_min;
        gaussian_cache[idx] = std::exp(-1.0 * sum * sum / (h * h));
    }
}


// 高斯权重函数
inline float gaussian(int x, int sigma) {
    if (x < sum_min || x > sum_max) {
        std::cout << "FUCK" << std::endl;
    }
    return gaussian_cache[x - sum_min];
}

// 计算两个邻域的相似性
inline float compute_similarity(
        RawImage* raw_image, int channel, 
        int x1, int y1, int x2, int y2, int patch_size, float h) {
    float sum = 0.0;
    int half_patch = patch_size / 2;
    int width = raw_image->width;
    int height = raw_image->height;
    auto yuv_image = raw_image->yuv_image;

    float diff = 0;
    for (int i = -half_patch; i <= half_patch; ++i) {
        for (int j = -half_patch; j <= half_patch; ++j) {
            int nx1 = x1 + i;
            int ny1 = y1 + j;
            int nx2 = x2 + i;
            int ny2 = y2 + j;
            if (nx1 >= 0 && nx1 < width && ny1 >= 0 && ny1 < height 
                    && nx2 >= 0 && nx2 < width && ny2 >= 0 && ny2 < height) {
                int idx1 = ny1 * width + nx1;
                int idx2 = ny2 * width + nx2;
                diff += yuv_image[idx1][channel] - yuv_image[idx2][channel];
                sum += diff * diff;
            }
        }   
    }
    return std::exp(-1.0 * sum / (h * h));
}

float fast_compute_similarity(
        RawImage* raw_image, int channel, 
        int x1, int y1, int x2, int y2, int patch_size, float h) {
    int half_patch = patch_size / 2;
    int width = static_cast<int>(raw_image->width);
    int height = static_cast<int>(raw_image->height);
    int x1_start = std::max(x1 - half_patch, 0);
    int y1_start = std::max(y1 - half_patch, 0);
    int x1_end = std::min(x1 + half_patch, width - 1);
    int y1_end = std::min(y1 + half_patch, height - 1);

    int x2_start = std::max(x2 - half_patch, 0);
    int y2_start = std::max(y2 - half_patch, 0);
    int x2_end = std::min(x2 + half_patch, width - 1);
    int y2_end = std::min(y2 + half_patch, height - 1);

    // Compute the block sums for each patch
    uint64_t sum1 = get_block_sum(raw_image, channel, x1_start, y1_start, x1_end, y1_end);
    uint64_t sum2 = get_block_sum(raw_image, channel, x2_start, y2_start, x2_end, y2_end);

    // Compute the squared difference
    //float diff = static_cast<float>(sum1 - sum2);

    // Apply the Gaussian function
    return gaussian(sum1 - sum2, h);
}

// Non-Local Means去噪算法
void nonlocal_means_denoising(
        RawImage* raw_image, int channel, 
        int template_window_size, int search_window_size, float h) {
    int half_search = search_window_size / 2;
    int width = raw_image->width;
    int height = raw_image->height;
    auto yuv_image = raw_image->yuv_image;
    auto denoised_image = raw_image->denoised_image;

    // 遍历图像中的每个像素
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int y = half_search; y < height - half_search; ++y) {
        for (int x = half_search; x < width - half_search; ++x) {
            float sum_weights = 0.0;
            float sum_pixels = 0.0;

            // 遍历搜索窗口
            for (int j = -half_search; j <= half_search; ++j) {
                for (int i = -half_search; i <= half_search; ++i) {
                    float weight = fast_compute_similarity(raw_image, channel, x, y, x + i, y + j, template_window_size, h);
                    //float weight = compute_similarity(raw_image, channel, x, y, x + i, y + j, template_window_size, h);
                    int idx = (y + j) * width + (x + i);
                    sum_weights += weight;
                    sum_pixels += weight * yuv_image[idx][channel];
                }
            }

            int idx = y * width + x;
            denoised_image[idx][channel] = static_cast<unsigned short>(sum_pixels / sum_weights);
        }
    }
}

void DenoiseProcessor::process(RawImage* raw_image, int algorithm) {
    int width = raw_image->width;
    int height = raw_image->height;
    convert_rgb_to_yuv(raw_image);
    if (!raw_image->denoised_image) {
        raw_image->denoised_image = (unsigned short(*)[3])malloc(height * width * sizeof(*raw_image->denoised_image));
    }
    for (int i = 0; i < width * height; ++i) {
        for (int c = 0; c < 3; ++c) {
            raw_image->denoised_image[i][c] = raw_image->yuv_image[i][c];
        }
    }
    if (algorithm == 1) {
        if (!raw_image->integral) {
            raw_image->integral = (uint64_t(*)[3])malloc(height * width * sizeof(*raw_image->integral));
        }
        int h = 100;
        int template_window_size = 7;
        int search_window_size = 35;
        initialize_gaussian_cache(h);
        std::cout << "initialize_gaussian_cache done!" << std::endl;
        compute_integral_image(raw_image);

        nonlocal_means_denoising(raw_image, 1, template_window_size, search_window_size, h);
        nonlocal_means_denoising(raw_image, 2, template_window_size, search_window_size, h);
    } else if (algorithm == 2) {
        apply_gaussian_blur_channel(raw_image, 1, 15, 1);
        apply_gaussian_blur_channel(raw_image, 2, 15, 1);
    } else if (algorithm == 3) {
        apply_bilateral_filter(raw_image, 1, 11, 20, 5);
        apply_bilateral_filter(raw_image, 2, 11, 20, 5);
    }

    convert_yuv_to_rgb(raw_image);
    std::cout << raw_image->yuv_image[657 * 8280 + 4394][1] << std::endl;
    std::cout << raw_image->denoised_image[657 * 8280 + 4394][1] << std::endl;
}
