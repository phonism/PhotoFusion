#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include <string>
#include "RawImage.h"

class ImageIO {
public:
    static RawImage ReadRawImage(const std::string& filename);
    
};

#endif // IMAGE_IO_H
