#include <iostream>
#include <fstream>
#include "photofusion_api.h"

int main(int argc, char** argv) {
    /*
    const char *rawFilePath = "./photo.nef";
    auto io = new ImageIO();
    io->read_raw_image(rawFilePath);

    auto raw_image = process_raw(io->image, 
            io->width, io->height, 
            io->black_level, io->cam_mul, "output.tiff");
            */
    raw2tiff("./photo.nef", "output.tiff");

    return 0;
}
