#include <iostream>
#include <fstream>
#include "PhotoFusionAPI.h"
#include "libraw/libraw.h"

int main(int argc, char** argv) {
    libraw_data_t *librawData = libraw_init(0);
    const char *rawFilePath = "./photo.nef";
    int rc = libraw_open_file(librawData, rawFilePath);
    // 调用处理函数
    if (librawData == NULL) {
        std::cout << "ERROR" << std::endl;
    }
    rc = libraw_unpack(librawData);
    if (!process_raw(librawData->rawdata.raw_image, 
                librawData->sizes.width, librawData->sizes.height, 
                librawData->color.black, "output.tiff")) {
        std::cerr << "Failed to process raw image" << std::endl;
        return 1;
    }

    libraw_close(librawData);
    return 0;
}
