#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include <sys/stat.h>
#include "RawImage.h"

struct tiff_data {
    int t_width, t_height, bps, comp, phint, offset, t_flip, samples, bytes, extrasamples;
    int t_tile_width, t_tile_length, sample_format, predictor;
    int rows_per_strip;
    int *strip_offsets, strip_offsets_count;
    int *strip_byte_counts, strip_byte_counts_count;
    unsigned t_filters;
    int t_vwidth, t_vheight, t_lm,t_tm;
    int t_fuji_width; 
    float t_shutter;
    int newsubfiletype;
};

static const int tagtype_bytes [20] = { 
    1, 1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8, 4, 2, 8, 8, 8, 8, 8
};

class DataStream {
public:
    DataStream(const char *fname);
    virtual ~DataStream();
    virtual int valid();
    virtual int read(void *ptr, size_t size, size_t nmemb);  
    virtual int eof();
    virtual int seek(int64_t o, int whence);
    virtual int64_t tell();
    virtual int64_t size() { return _fsize; }
    virtual char *gets(char *str, int sz);
    virtual int scanf_one(const char *fmt, void *val);
    virtual const char *fname();
    virtual int get_char() {
        return fgetc(f);
    }

protected:
    FILE *f;
    char* filename;
    int64_t _fsize;
};

class ImageIO {
public:
    RawImage* read_raw_image(const char* filename);
    int open_datastream(DataStream* stream);
    void tiff_get(DataStream* stream, unsigned base, unsigned *tag, unsigned *type, unsigned *len, unsigned *save);
    int parse_tiff_ifd(DataStream* stream, int base);
    int parse_tiff(DataStream* stream);
    int parse_exif(DataStream* stream, int base);
    void parse_makernote(DataStream* stream, int base, int uptag);
    void apply_tiff(DataStream* stream);
    void open_raw(DataStream* stream);

public:
    tiff_data tiff_ifd[8];
    int ifd = -1;
    int width;
    int height;
    int bps;
    unsigned short* image;
    int offset;
    int meta_offset;
    int black_level;
    float cam_mul[4];

    ~ImageIO() {
        if (image) {
            free(image);
            image = nullptr;
        }
    }
    
};

#endif // IMAGE_IO_H
