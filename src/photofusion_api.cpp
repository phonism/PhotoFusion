// PhotoFusionAPI.cpp
#include <iostream>
#include "photofusion_api.h"
#include "image_io.h"
#include "raw_processor.h"
#include "stack.h"
#include "raw_image.h"
#include "denoise.h"


class PhotoFusionImpl {
public:
    void raw2tiff(const char* input_file, const char* output_file) {
        auto io = new ImageIO();
        io->read_raw_image(input_file);

        auto raw_image = process_raw(io->image, 
            io->width, io->height, 
            io->black_level, io->cam_mul, output_file);
        delete raw_image;
        delete io;
    }

    RawImage* process_raw(
            const unsigned short* libraw_image, 
            int width, int height, int black_level, 
            const float* cam_mul,
            const char* output_file) {
        RawImage* raw_image = new RawImage();
        try {
            raw_image->width = width;
            raw_image->height = height;
            raw_image->black_level = black_level;
            raw_image->raw_data = static_cast<unsigned short*>(std::malloc(width * height * sizeof(unsigned short)));
            long long sum = 0;
            for (int i = 0; i < height * width; ++i) {
                raw_image->raw_data[i] = libraw_image[i];
            }
            if (cam_mul != nullptr) {
                for (int i = 0; i < 4; ++i) {
                    raw_image->cam_mul[i] = cam_mul[i];
                }
            }
            // Process the raw image and convert it to an RGB image
            RawProcessor raw_processor;
            raw_processor.raw2image(*raw_image);
            raw_processor.scale_colors(*raw_image);
            raw_processor.pre_interpolate(*raw_image);
            raw_processor.ahd_interpolate(*raw_image);
            raw_processor.convert_to_rgb(*raw_image);
            raw_processor.gamma_adjustment(*raw_image);

            DenoiseProcessor denoise_processor;
            denoise_processor.process(raw_image, 3);

#ifdef USE_JPEG
            raw_processor.write_jpeg(*raw_image, output_file);
#else
            raw_processor.ppm_tiff_writer(*raw_image, output_file);
#endif
            return raw_image;
        } catch (const std::exception& e) {
            std::cerr << "An error occurred: " << e.what() << std::endl;
            return raw_image;
        }
        return raw_image;
    }
public:
    bool process_stack(std::vector<RawImage*>& images) {
        auto stack_proc = StackProcessor();
        auto raw_processor = RawProcessor();
        auto result = stack_proc.process(images);
        raw_processor.gamma_adjustment(*result);
        const char* output_file = "output_stack.tiff";
        raw_processor.ppm_tiff_writer(*result, output_file);
        return true;
    }
};

// 实现C风格的函数
RawImage* process_raw(
        const unsigned short* libraw_image, 
        int width, int height, int black_level, 
        const float* cam_mul, const char* output_file) {
    PhotoFusionImpl impl;
    return impl.process_raw(libraw_image, width, height, black_level, cam_mul, output_file);
}

bool process_stack(RawImage** images, int count) {
    PhotoFusionImpl impl;
    std::vector<RawImage*> image_vector(images, images + count);
    return impl.process_stack(image_vector);
}

void raw2tiff(const char* input_file, const char* output_file) {
    PhotoFusionImpl impl;
    return impl.raw2tiff(input_file, output_file);
}
