#include "Stack.h"
#include <cmath>
#ifdef USE_FFTW
#include <fftw3.h>
#endif

// 计算均方误差
double mse_loss(RawImage* img1, RawImage* img2, int tx, int ty) {
    double loss = 0.0;
    int count = 0;
#ifdef USE_OPENMP
    #pragma omp parallel for reduction(+:loss, count) collapse(2)
#endif
    for (unsigned int y = 0; y < img1->height; ++y) {
        for (unsigned int x = 0; x < img1->width; ++x) {
            int x2 = x + tx;
            int y2 = y + ty;

            if (x2 >= 0 && x2 < static_cast<int>(img2->width) && y2 >= 0 && y2 < static_cast<int>(img2->height)) {
                for (int c = 0; c < 3; ++c) {
                    double diff = img1->image[y * img1->width + x][c] - img2->image[y2 * img2->width + x2][c];
                    loss += diff * diff;
                }
                ++count;
            }
        }
    }
    return loss / count;
}



void brute_force(std::vector<RawImage*>& images) {
    int height = images[0]->height;
    int width = images[0]->width;
    int num_images = images.size();

    for (int i = 1; i < num_images; ++i) {
        RawImage* img1 = images[0];
        RawImage* img2 = images[i];
        double min_loss = 1e9;
        int tx = 0, ty = 0;
#ifdef USE_OPENMP
        #pragma omp parallel for collapse(2) default(none) shared(img1, img2, min_loss, tx, ty)
#endif
        for (int i = -40; i <= 16; ++i) {
            for (int j = -20; j <= 15; ++j) {
                double loss = mse_loss(img1, img2, i, j);
#ifdef USE_OPENMP
                // 使用临界区以确保更新共享变量的安全性
                #pragma omp critical
#endif
                {
                    if (loss < min_loss) {
                        tx = i;
                        ty = j;
                        min_loss = loss;
                    }
                }
            }
        }
        std::cout << "Final Results: " << std::endl;
        std::cout << "min_loss = " << min_loss << ", tx = " << tx << ", ty = " << ty << std::endl;
        unsigned short* aligned_image = new unsigned short[img2->width * img2->height * 4]();

        for (unsigned int y = 0; y < img2->height; ++y) {
            for (unsigned int x = 0; x < img2->width; ++x) {
                int x2 = x + tx;
                int y2 = y + ty;

                if (x2 >= 0 && x2 < static_cast<int>(img2->width) && y2 >= 0 && y2 < static_cast<int>(img2->height)) {
                    for (int c = 0; c < 4; ++c) {
                        aligned_image[(y * img2->width + x) * 4 + c] = img2->image[(y2 * img2->width + x2)][c];
                    }
                } else {
                    // 边界处理，设置为0
                    for (int c = 0; c < 4; ++c) {
                        aligned_image[(y * img2->width + x) * 4 + c] = 0;
                    }
                }
            }
        }
        delete[] img2->image;
        img2->image = reinterpret_cast<unsigned short(*)[4]>(aligned_image);
    }
}

#ifdef USE_FFTW
void forward_FFT(double* input, fftw_complex* output, int height, int width) {
    fftw_plan plan = fftw_plan_dft_r2c_2d(height, width, input, output, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

// 将频域结果转换回空间域
void inverse_FFT(fftw_complex* input, double* output, int height, int width) {
    fftw_plan plan = fftw_plan_dft_c2r_2d(height, width, input, output, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

// 计算相位谱
void compute_phase_correlation(fftw_complex* img1, fftw_complex* img2, fftw_complex* result, int height, int width) {
    for (int i = 0; i < height * (width / 2 + 1); ++i) {
        double re1 = img1[i][0];
        double im1 = img1[i][1];
        double re2 = img2[i][0];
        double im2 = img2[i][1];

        // 计算共轭
        double denom = std::sqrt((re1 * re1 + im1 * im1) * (re2 * re2 + im2 * im2));
        result[i][0] = (re1 * re2 + im1 * im2) / denom;
        result[i][1] = (im1 * re2 - re1 * im2) / denom;
    }
}

// 寻找峰值位置
std::pair<int, int> find_peak(double* data, int height, int width) {
    int maxIdx = 0;
    for (int i = 1; i < height * width; ++i) {
        if (data[i] > data[maxIdx]) {
            maxIdx = i;
        }
    }
    int x = maxIdx % width;
    int y = maxIdx / width;
    return std::make_pair(x, y);
}


void phase_correlation(std::vector<RawImage*>& images) {
    // 初始化FFTW多线程支持
    fftw_init_threads();
    
    // 设置使用的线程数
    int num_threads = 128; // 你可以根据你的硬件情况设置合适的线程数
    fftw_plan_with_nthreads(num_threads);

    int num_images = images.size();
    RawImage* img2 = images[0];
    int height = img2->height;
    int width = img2->width;
    double* img1_fft = new double[height * width];
    double* img2_fft = new double[height * width];
    fftw_complex* fft1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * height * (width / 2 + 1));
    fftw_complex* fft2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * height * (width / 2 + 1));
    fftw_complex* result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * height * (width / 2 + 1));
    double* corr = new double[height * width];
    double* combined_corr = new double[height * width]();
    for (int idx = 1; idx < num_images; ++idx) {
        RawImage* img1 = images[idx];
        for (int c = 0; c < 3; ++c) {
            // 将每个通道的数据复制到 img1_fft 和 img2_fft
            for (int i = 0; i < height * width; ++i) {
                img1_fft[i] = img1->image[i][c];
                img2_fft[i] = img2->image[i][c];
            }

            forward_FFT(img1_fft, fft1, height, width);
            forward_FFT(img2_fft, fft2, height, width);

            compute_phase_correlation(fft1, fft2, result, height, width);

            inverse_FFT(result, corr, height, width);

            for (int i = 0; i < height * width; ++i) {
                combined_corr[i] += corr[i];
            }
        }

        // 归一化逆傅里叶变换结果
        for (int i = 0; i < height * width; ++i) {
            combined_corr[i] /= (height * width);
        }

        // 找到峰值位置
        auto peak = find_peak(combined_corr, height, width);
        int dx = peak.first;
        int dy = peak.second;

        // 计算平移量（考虑周期性边界条件）
        if (dx > width / 2) {
            dx -= width;
        }
        if (dy > height / 2) {
            dy -= height;
        }

        unsigned short* aligned_image = new unsigned short[img1->width * img1->height * 4]();

        for (unsigned int y = 0; y < img1->height; ++y) {
            for (unsigned int x = 0; x < img1->width; ++x) {
                int x2 = x + dx;
                int y2 = y + dy;

                if (x2 >= 0 && x2 < static_cast<int>(img1->width) && y2 >= 0 && y2 < static_cast<int>(img1->height)) {
                    for (int c = 0; c < 4; ++c) {
                        aligned_image[(y * img1->width + x) * 4 + c] = img1->image[(y2 * img1->width + x2)][c];
                    }
                } else {
                    // 边界处理，设置为0
                    for (int c = 0; c < 4; ++c) {
                        aligned_image[(y * img1->width + x) * 4 + c] = 0;
                    }
                }
            }
        }
        delete[] img1->image;
        img1->image = reinterpret_cast<unsigned short(*)[4]>(aligned_image);

        std::cout << "dx: " << dx << ", dy: " << dy << std::endl;
    }

    // 释放内存
    delete[] img1_fft;
    delete[] img2_fft;
    delete[] corr;
    delete[] combined_corr;
    fftw_free(fft1);
    fftw_free(fft2);
    fftw_free(result);
}
#endif

RawImage* StackProcessor::process(std::vector<RawImage*>& images) {
    if (images.empty()) {
        std::cerr << "[StackProcessor]: images is empty!" << std::endl;
        return new RawImage();
    }
    int height = images[0]->height;
    int width = images[0]->width;
    int num_images = images.size();

#ifdef USE_FFTW
    phase_correlation(images);
#else
    brute_force(images);
#endif

    RawImage* result = new RawImage();
    result->height = height;
    result->width = width;
    result->image = (unsigned short(*)[4])malloc(height * width * sizeof(*result->image));
    memset(result->image, 0, height * width * sizeof(*result->image));
    for (int row = 0; row < height; ++row) {
        for (int col = 0; col < width; ++col) {
            for (int c = 0; c < 3; ++c) {
                int pix = 0;
                for (auto& img : images) {
                    pix += img->image[row * width + col][c];
                    // pix = std::max(pix, (int)img->image[row * width + col][c]);
                }
                pix = pix / num_images;
                result->image[row * width + col][c] = (unsigned short)pix;
            }
        }
    }
    if (!result->histogram) {   
        result->histogram = (int(*)[HISTOGRAM_SIZE])malloc(sizeof(int) * HISTOGRAM_SIZE * 4); 
        if (!result->histogram) {
            std::cerr << "Memory allocation failed." << std::endl;
        }
    }
    memset(result->histogram, 0, sizeof(int) * HISTOGRAM_SIZE * 4);
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            auto& img = result->image[row * width + col];
            result->histogram[0][img[0] >> 3]++;
            result->histogram[1][img[1] >> 3]++;
            result->histogram[2][img[2] >> 3]++;
        }
    }
    return result;
}
