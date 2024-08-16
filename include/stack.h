#ifndef STACK_H
#define STACK_H

#include <vector>
#include <iostream>
#include "common.h"
#include "raw_image.h"

class StackProcessor {
public:
    static RawImage* process(std::vector<RawImage*>& images);
};

#endif // STACK_H
