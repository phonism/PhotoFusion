// PhotoFusionAPI.cpp
#include "PhotoFusionAPI.h"
#include <iostream>
#include "ImageIO.h"
#include "RawProcessor.h"

class PhotoFusionImpl {
public:
    bool process_raw(const unsigned short* libraw_image, int width, int height, int black_level, const char* output_file) {
        ImageIO imageIO;

        try {
            RawImage raw_image = imageIO.ReadRawImage("");
            raw_image.width = width;
            raw_image.height = height;
            raw_image.black_level = black_level;
            for (int i = 0; i < height * width; ++i) {
                raw_image.raw_data.push_back(libraw_image[i]);
            }

            // Process the raw image and convert it to an RGB image
            RawProcessor raw_processor;;
            raw_processor.raw2image(raw_image);
            raw_processor.scale_colors(raw_image);
            raw_processor.pre_interpolate(raw_image);
            raw_processor.ahd_interpolate(raw_image);
            raw_processor.convert_to_rgb(raw_image);
            raw_processor.ppm_tiff_writer(raw_image, output_file);
        } catch (const std::exception& e) {
            std::cerr << "An error occurred: " << e.what() << std::endl;
            return 1;
        }
    }
};

// 实现C风格的函数
bool process_raw(const unsigned short* libraw_image, int width, int height, int black_level, const char* output_file) {
    PhotoFusionImpl impl;
    return impl.process_raw(libraw_image, width, height, black_level, output_file);
}
