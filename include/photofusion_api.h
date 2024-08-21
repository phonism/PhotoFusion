#ifndef PHOTOFUSION_API_H
#define PHOTOFUSION_API_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

struct RawImage;

// 定义一个结构体来表示PhotoFusion的上下文（如果有需要）
typedef struct PhotoFusionContext PhotoFusionContext;

// C风格的函数接口
void raw2tiff(const char* input_file, const char* output_file);

struct RawImage* process_raw(
        const unsigned short* image, 
        int width, int height, int black_level, 
        const float* cam_mul,
        const char* output_file);

bool process_stack(struct RawImage** images, int count);

#ifdef __cplusplus
}
#endif

#endif // PHOTOFUSION_API_H
