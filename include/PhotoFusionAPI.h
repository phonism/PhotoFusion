#ifndef PHOTOFUSION_API_H
#define PHOTOFUSION_API_H

#include <iostream>
#include <fstream>
#include <libraw/libraw.h>
#include "ImageIO.h"
#include "RawProcessor.h"

class PhotoFusionAPI {
public:
    static bool process_raw(const char* input_file, const char* output_file);
};
#endif // PHOTOFUSION_API_H
