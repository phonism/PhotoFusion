#include <iostream>
#include <filesystem>
#include <string>
#include <algorithm>
#include <vector>
#include "ImageIO.h"
#include "PhotoFusionAPI.h"

// namespace fs = std::__fs::filesystem;
namespace fs = std::filesystem;

bool has_suffix(const std::string& str, const std::string& suffix) {
    if (suffix.size() > str.size()) {
        return false;
    }
    return str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main() {
    std::string path = "stack_images";
    std::vector<RawImage*> images;
    std::vector<fs::directory_entry> entries;

    // 收集所有目录条目
    for (const auto& entry : fs::directory_iterator(path)) {
        entries.push_back(entry);
    }

    // 按照字典序排序
    std::sort(entries.begin(), entries.end(), [](const fs::directory_entry& a, const fs::directory_entry& b) {
        return a.path().filename().string() < b.path().filename().string();
    });

    // 遍历排序后的条目
#ifdef USE_OPENMP
    // TODO: openmp is not working for now
    // #pragma omp parallel for
#endif
    for (int i = 0; i < entries.size(); ++i) {
        const auto& entry = entries[i];
        if (fs::is_regular_file(entry)) {
            auto photo_name = entry.path().filename().string();
            if (!has_suffix(photo_name, ".NEF")) {
                continue;
            }
            std::cout << "Processing photo: " << photo_name << std::endl;
            std::string fullPath = path + "/" + photo_name;
            const char* rawFilePath = fullPath.c_str();
            auto io = new ImageIO();
            io->read_raw_image(rawFilePath);
            auto raw_image = process_raw(io->image, 
                    io->width, io->height, 
                    io->black_level, io->cam_mul, "output.tiff");
            std::cout << "Processing DONE!" << std::endl;
#ifdef USE_OPENMP
            // #pragma omp critical
#endif
            {
                images.push_back(raw_image);
            }
        }
    }

    RawImage** imageArray = new RawImage*[images.size()];
    for (size_t i = 0; i < images.size(); ++i) {
        imageArray[i] = images[i];
    }
    process_stack(imageArray, images.size());
    return 0;
}
