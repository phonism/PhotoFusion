#ifndef STACK_H
#define STACK_H

#include <vector>
#include <iostream>
#include "Common.h"
#include "RawImage.h"

class StackProcessor {
public:
    static RawImage* process(std::vector<RawImage*>& images);
};

#endif // STACK_H
