#ifndef PHOTOFUSION_API_H
#define PHOTOFUSION_API_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

// 定义一个结构体来表示PhotoFusion的上下文（如果有需要）
typedef struct PhotoFusionContext PhotoFusionContext;

// C风格的函数接口
bool process_raw(const unsigned short* libraw_image, int width, int height, int black_level, const char* output_file);

#ifdef __cplusplus
}
#endif

#endif // PHOTOFUSION_API_H
