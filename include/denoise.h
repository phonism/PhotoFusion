#ifndef DENOISE_H
#define DENOISE_H

#include <vector>
#include <iostream>
#include "common.h"
#include "raw_image.h"

class DenoiseProcessor {
public:
    static void process(RawImage* raw_image, int algorithm);
};

#endif // DENOISE_H
