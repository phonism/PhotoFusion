#ifndef RAW_IMAGE_H
#define RAW_IMAGE_H

#include <vector>
#include <cstring>
#include <array>
#include <cstdint>
#include <cfloat>
#include <string>
#include <iostream>
#include <cmath>

struct RawImage {
    unsigned int width;
    unsigned int height;
    unsigned short bitDepth;
    unsigned short black_level;
    unsigned short data_max;
    unsigned short maximum = 16383;
    unsigned int filters = 3031741620;
    
    std::vector<unsigned short> raw_data;
    unsigned short (*image)[4];
    std::vector<float> pre_mul = {1.9722067118, 0.9411969781, 1.1376225948, 0.9411969781};
    unsigned short curve[0x10000];
    int (*histogram)[0x2000] = NULL;

    unsigned *oprof = NULL;

    double gamm[6] = {0.45, 4.5, 0, 0, 0, 0};
    FILE *output;

    ~RawImage() {
        if (image) {
            delete[] image;
            image = nullptr;
        }
        if (histogram) {
            delete[] histogram;
            histogram = nullptr;
        }
        if (oprof) {
            delete[] oprof;
            oprof = nullptr;
        }
        if (output) {
            fclose(output);
            output = nullptr;
        }
    }

    std::string to_string() {
        std::cout << "==========================================================" << std::endl;
        std::cout << "height: " << height << " width:" << width << std::endl;
        std::cout << "black_level:" << black_level << std::endl;
        std::cout << "data_max:" << data_max << std::endl;
        std::cout << "raw_image[0]:" << raw_data[0] << std::endl;
        std::cout << "raw_image[1]:" << raw_data[1] << std::endl;
        std::cout << "raw_image[2]:" << raw_data[2] << std::endl;
        std::cout << "raw_image[3]:" << raw_data[3] << std::endl;
        std::cout << "image[0]:" << image[0][0] << " " << image[0][1] << " " << image[0][2] << " " << image[0][3] << std::endl;
        std::cout << "image[1]:" << image[1][0] << " " << image[1][1] << " " << image[1][2] << " " << image[1][3] << std::endl;
        std::cout << "image[2]:" << image[width][0] << " " << image[width][1] << " " << image[width][2] << " " << image[width][3] << std::endl;
        std::cout << "image[3]:" << image[width + 1][0] << " " << image[width + 1][1] << " " << image[width + 1][2] << " " << image[width + 1][3] << std::endl;
        int idx = 100 * width + 100;
        std::cout << "image[100]:" << image[idx][0] << " " << image[idx][1] << " " << image[idx][2] << " " << image[idx][3] << std::endl;
        std::cout << "image[101]:" << image[idx + 1][0] << " " << image[idx + 1][1] << " " << image[idx + 11][2] << " " << image[idx + 1][3] << std::endl;
        std::cout << "image[102]:" << image[idx + width][0] << " " << image[idx + width][1] << " " << image[idx + width][2] << " " << image[idx + width][3] << std::endl;
        std::cout << "image[103]:" << image[idx + width + 1][0] << " " << image[idx + width + 1][1] << " " << image[idx + width + 1][2] << " " << image[idx + width + 1][3] << std::endl;
        if (histogram != 0) {
            std::cout << "histogram:" << histogram[0][0] << std::endl;
        }
        std::cout << "==========================================================" << std::endl;
        return "";
    }
};

#endif // RAW_IMAGE_H
