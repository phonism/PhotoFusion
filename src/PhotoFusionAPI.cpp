#include "PhotoFusionAPI.h"

bool PhotoFusionAPI::process_raw(const char* input_file, const char* output_file) {
    LibRaw LibRawProcessor;

    int ret = LibRawProcessor.open_file(input_file);
    if (ret != LIBRAW_SUCCESS) {
        std::cerr << "Cannot open file: " << input_file << std::endl;
        return 1;
    }

    // 解压RAW数据
    ret = LibRawProcessor.unpack();
    if (ret != LIBRAW_SUCCESS) {
        std::cerr << "Cannot unpack " << input_file << std::endl;
        return 1;
    }

    ImageIO imageIO;

    try {
        RawImage raw_image = imageIO.ReadRawImage(input_file);
        ushort *libraw_image = LibRawProcessor.imgdata.rawdata.raw_image;
        int width = LibRawProcessor.imgdata.sizes.width;
        int height = LibRawProcessor.imgdata.sizes.height;
        raw_image.width = width;
        raw_image.height = height;
        raw_image.black_level = LibRawProcessor.imgdata.color.black;
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

    LibRawProcessor.recycle();
    return 1;
}
