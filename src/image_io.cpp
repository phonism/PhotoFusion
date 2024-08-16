#include <iostream>
#include <time.h>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "image_io.h"

RawImage* ImageIO::read_raw_image(const char* filename) {
    RawImage* rawImage;
    // Add the implementation for reading a raw image from a file
    // For now, we'll just create an empty raw image.
    auto stream = new DataStream(filename);
    open_datastream(stream);
    open_raw(stream);
    delete stream;
    return rawImage;
}

unsigned short get2(DataStream* stream) {
    unsigned char s[2] = {0xff, 0xff};
    stream->read(s, 1, 2);
    return s[0] | s[1] << 8;
}

unsigned get4(DataStream* stream) {
    unsigned char s[4] = {0xff, 0xff, 0xff, 0xff};
    stream->read(s, 1, 4);
    return s[0] | s[1] << 8 | s[2] << 16 | s[3] << 24;;
}

double getreal(DataStream* stream, int type) {
    switch (type) {
        case 4:
            return (unsigned int)get4(stream);
        case 5:
            double ud = (unsigned int)get4(stream);
            double vd = (unsigned int)get4(stream);
            return ud / (vd ? vd : 1);
    }
    return 0.0;
}

int getint(DataStream* stream, int type) {
    switch (type) {
        case 3:
            return get2(stream);
        default:
            return get4(stream);
    }
    return 0.0;
}

void ImageIO::tiff_get(DataStream* stream, unsigned base, unsigned *tag, unsigned *type, unsigned *len, unsigned *save) {
    *tag = get2(stream);
    *type = get2(stream);
    *len = get4(stream);
    *save = stream->tell() + 4;
    if (*len * tagtype_bytes[*type] > 4) {
        stream->seek(get4(stream) + base, SEEK_SET);
    }
}

time_t get_timestamp(DataStream* stream) {
    struct tm t;
    char str[20];
    int i;

    str[19] = 0;
    stream->read(str, 19, 1);
    memset(&t, 0, sizeof t); 
    sscanf(str, "%d:%d:%d %d:%d:%d", &t.tm_year, &t.tm_mon, &t.tm_mday, &t.tm_hour, &t.tm_min, &t.tm_sec);
    t.tm_year -= 1900;
    t.tm_mon -= 1;
    t.tm_isdst = -1; 
    return mktime(&t);
}

void ImageIO::parse_makernote(DataStream* stream, int base, int uptag) {
    int64_t fsize = stream->size();
    char buf[10];
    stream->read(buf, 1, 10);
    base = stream->tell();
    int order = get2(stream);
    get2(stream);
    int offset = get4(stream);
    int entries = get2(stream);

    int cblack[4];
    int color_balance_version = 0;
    int ColorBalanceData_ready = 0;
    char ColorBalanceData_buf[324];
    while (entries--) {
        unsigned tag, type, len, save;
        tiff_get(stream, base, &tag, &type, &len, &save);
        tag |= uptag;
        switch (tag) {
            case 0x008c:
            case 0x0096:
                meta_offset = stream->tell();
                break;
            case 0x003d:
                for (int i = 0; i < 4; ++i) {
                    cblack[i] = get2(stream);
                }
                black_level = cblack[3];
                break;
            case 0x0097:
                for (int i = 0; i < 4; ++i) {
                    color_balance_version = color_balance_version * 10 + stream->get_char() - '0';
                }
                stream->seek(0x118L, SEEK_CUR);
                ColorBalanceData_ready = stream->read(ColorBalanceData_buf, 324, 1) == 1;
                // std::cout << "cam_mul:" << ColorBalanceData_ready << std::endl;
                break;
            case 0x000c:
                cam_mul[0] = getreal(stream, type);
                cam_mul[2] = getreal(stream, type);
                cam_mul[1] = getreal(stream, type);
                cam_mul[3] = getreal(stream, type);
                break;
        }
        stream->seek(save, SEEK_SET);
    }
}

int ImageIO::parse_exif(DataStream* stream, int base) {
    int entries = get2(stream);
    int fsize = stream->size();
    while (entries--) {
        unsigned tag, type, len, save;
        tiff_get(stream, base, &tag, &type, &len, &save);
        int64_t savepos = stream->tell();
        switch(tag) {
            case 0x829a: // 33434 TODO
                // double shutter = getreal(stream, type);
                getreal(stream, type);
                break;
            case 0x829d: // 33437, FNumber aperture
                getreal(stream, type);
                break;
            case 0x8827: // 34855 iso_speed
                get2(stream);
                break;
            case 0x8832: // 34866 iso_speed
                getreal(stream, type);
                break;
            case 0x9003: // 36867
            case 0x9004: // 36868 timestamp
                get_timestamp(stream);
                break;
            case 0x9209: // 37385 flash_used
                getreal(stream, type);
                break;
            case 0x920a: // 37386 focal_len
                getreal(stream, type);
                break;
            case 0x927c: // 37500
                parse_makernote(stream, base, 0);
                break;
        }
        stream->seek(save, SEEK_SET);
    }
    return 1;
}

int ImageIO::parse_tiff_ifd(DataStream* stream, int base) {
    unsigned short entries = get2(stream);
    if (ifd > (int)(sizeof tiff_ifd / sizeof tiff_ifd[0])) {
        return 1;
    }
    if (entries > 512) {
        return 1;
    }
    ifd++;
    int64_t fsize = stream->size();
    char make[64], model[64], software[64], artist[64];
    time_t timestamp;
    while (entries--) {
        unsigned tag, type, len, save;
        tiff_get(stream, base, &tag, &type, &len, &save);
        int64_t savepos = stream->tell();
        switch(tag) {
            case 0xf007:
                break;
            case 0x00fe: // 254
                tiff_ifd[ifd].newsubfiletype = getreal(stream, type);
                break;
            case 0x0103: // 259
                tiff_ifd[ifd].comp = getint(stream, type);
                break;
            case 0x010f: // 271 品牌
                stream->gets(make, 64);
                break;
            case 0x0110: // 272 相机型号
                stream->gets(model, 64);
                break;
            case 0x0112: // 274
                tiff_ifd[ifd].t_flip = "50132467"[get2(stream) & 7] - '0';
                break;
            case 0x0131: // 305
                stream->gets(software, 64);
                break;
            case 0x0132: // 306
                timestamp = get_timestamp(stream);
                break;
            case 0x013b: // 315 artist
                stream->read(artist, 64, 1);
                break;
            case 0x014a: // 330 subifds
                while (len--) {
                    int idx = stream->tell();
                    stream->seek(get4(stream) + base, SEEK_SET);
                    if (parse_tiff_ifd(stream, base)) 
                        break;
                    stream->seek(idx + 4, SEEK_SET);
                }
                break;
            case 0x0106: // 262
                tiff_ifd[ifd].phint = get2(stream);
                break;
            case 0x0115: /* 277, SamplesPerPixel */
                tiff_ifd[ifd].samples = getint(stream, type) & 7;
                break;
            case 0x0116: // 278
                tiff_ifd[ifd].rows_per_strip = getint(stream, type);
                break;
            case 0x828d: /* 33421, CFARepeatPatternDim */
                break;
            case 0x8769: /* 34665, EXIF tag */
                stream->seek(get4(stream) + base, SEEK_SET);
                // TODO parse exif
                parse_exif(stream, base);
                break;
            case 0x8825: /* 34853, GPSInfo tag */
                break;
            case 0x0100: // 256
                tiff_ifd[ifd].t_width = getint(stream, type);
                break;
            case 0x0101: // 257
                tiff_ifd[ifd].t_height = getint(stream, type);
                break;
            case 0x0102: /* 258, BitsPerSample */
                tiff_ifd[ifd].bps = getint(stream, type);
                break;
            case 0x0202: // 514
                tiff_ifd[ifd].bytes = get4(stream);
                break;
            case 0x0117: // 279
                tiff_ifd[ifd].bytes = get4(stream);
                break;
            case 0x0111: // 273
            case 0x0201: // 513 TODO: 可能有jpeg解析
                tiff_ifd[ifd].offset = get4(stream) + base;
                break;
        }
        stream->seek(save, SEEK_SET);
    }
    return 0;
}

int ImageIO::parse_tiff(DataStream* stream) {
    int doff;
    stream->seek(0, SEEK_SET);
    unsigned short order = get2(stream);
    get2(stream);
    int cnt = 1;
    while (doff = get4(stream)) {
        int64_t doff64 = doff;
        if (doff64 > stream->size()) {
            break;
        }
        stream->seek(doff64, SEEK_SET);
        if (parse_tiff_ifd(stream, 0)) {
            break;
        }
    }
    return 1;
}

void ImageIO::apply_tiff(DataStream* stream) {
    int shutter = 0;
    for (int i = ifd; i--;) {
        if (tiff_ifd[i].t_shutter) {
            shutter = tiff_ifd[i].t_shutter;
        }
        tiff_ifd[i].t_shutter = shutter;
    }
    int raw_idx = 0, max_height = 0;
    for (int i = 0; i <= ifd; ++i) {
        if (tiff_ifd[i].t_height > max_height) {
            max_height = tiff_ifd[i].t_height;
            raw_idx = i;
        }
    }
    width = tiff_ifd[raw_idx].t_width;
    height = tiff_ifd[raw_idx].t_height;
    offset = tiff_ifd[raw_idx].offset;
    bps = tiff_ifd[raw_idx].bps;
}

int ImageIO::open_datastream(DataStream* stream) {
    if (!stream) {
        return -1;
    }
    if (!stream->valid()) {
        return -1;
    }

    char head[64] = {0};
   
    unsigned short order = get2(stream);
    unsigned hlen = get4(stream);
    // 移到文件头
    stream->seek(0, SEEK_SET);
    if (stream->read(head, 1, 64) < 64) {
        return -1;
    }
    stream->seek(0, SEEK_END);

    int file_len = stream->tell();
    parse_tiff(stream);
    apply_tiff(stream);
    return 0;
}

unsigned short* make_decoder_ref(const unsigned char **source) {
    int max, len, h, i, j;
    const unsigned char *count;
    unsigned short *huff;

    count = (*source += 16) - 17;
    for (max = 16; max && !count[max]; max--) {
    }
    huff = (unsigned short *)calloc(1 + (1 << max), sizeof *huff);
    huff[0] = max;
    for (h = len = 1; len <= max; len++) {
        for (i = 0; i < count[len]; i++, ++*source) {
            for (j = 0; j < 1 << (max - len); j++) {
                if (h <= 1 << max) {
                    huff[h++] = len << 8 | **source;
                }
            }
        }
    }
    return huff;
}

ushort* make_decoder(const unsigned char *source) {
    return make_decoder_ref(&source);
}

unsigned short getbithuff(DataStream* stream, int nbits, unsigned short *huff) {
    static unsigned bitbuf = 0;
    static int vbits = 0, reset = 0;
    unsigned c;
    if (nbits > 25) {
        return 0;
    }
    if (nbits < 0) {
        return bitbuf = vbits = reset = 0; 
    }
    if (nbits == 0 || vbits < 0) {
        return 0;
    }
    while (!reset && vbits < nbits && (c = stream->get_char()) != (unsigned)EOF &&
         !(reset = 0 && c == 0xff && stream->get_char())) {
        bitbuf = (bitbuf << 8) + (unsigned char)c;
        vbits += 8;
    }
    c = vbits == 0 ? 0 : bitbuf << (32 - vbits) >> (32 - nbits);
    if (huff) {
        vbits -= huff[c] >> 8;
        c = (unsigned char)huff[c];
    } else { 
        vbits -= nbits;
    }
    return c;
}

void read_shorts(DataStream* stream, unsigned short* pixel, unsigned count) {
    stream->read(pixel, 2, count);
}

void ImageIO::open_raw(DataStream* stream) {
    static const unsigned char nikon_tree[][32] = {
        {0, 1, 5, 1, 1, 1, 1, 1, 1, 2, 0,  0,  0, 0, 0, 0, 
        5, 4, 3, 6, 2, 7, 1, 0, 8, 9, 11, 10, 12}, /* 12-bit lossy */
        {0,    1,    5,    1,    1,    1, 1, 1, 1, 2, 0, 0,  0,  0,
        0,    0, 0x39, 0x5a, 0x38, 0x27, 0x16, 5, 4, 3, 2, 1, 0, 11, 12, 12}, /* 12-bit lossy after split */
        {0, 1, 4, 2, 3, 1, 2, 0, 0, 0, 0,  0,  0, 0, 0, 0, /* 12-bit lossless */
        5, 4, 6, 3, 7, 2, 8, 1, 9, 0, 10, 11, 12}, 
        {0, 1, 4, 3, 1, 1, 1, 1, 1, 2, 0,  0,  0,  0,  0, 0, /* 14-bit lossy */
        5, 6, 4, 7, 8, 3, 9, 2, 1, 0, 10, 11, 12, 13, 14}, 
        {0, 1,    5,    1,    1,    1, 1, 1, 1, 1, 2, 0, 0, 0,  0,
        0, /* 14-bit lossy after split */
        8, 0x5c, 0x4b, 0x3a, 0x29, 7, 6, 5, 4, 3, 2, 1, 0, 13, 14}, 
        {0, 1, 4, 2, 2, 3, 1,  2, 0,  0,  0, 0, 0, 0,  0, 0, /* 14-bit lossless */
        7, 6, 8, 5, 9, 4, 10, 3, 11, 12, 2, 0, 1, 13, 14}};
    int tree = 0;
    int min, row, col, idx, len, shl, diff;
    unsigned short hpred[2], vpred[2][2];

    image = (unsigned short*)malloc(width * height * sizeof(unsigned short));
    stream->seek(meta_offset, SEEK_SET);
    char ver0 = stream->get_char();
    char ver1 = stream->get_char();
    if (ver0 == 0x46) {
        tree = 2;
    }
    if (bps == 14) {
        tree += 3;
    }
    read_shorts(stream, vpred[0], 4);
    int max = 1 << bps & 0x7fff;
    ushort* huff = make_decoder(nikon_tree[tree]);
    stream->seek(offset, SEEK_SET);
    
    getbithuff(stream, -1, 0);

    for (min = row = 0; row < height; ++row) {
        for (col = 0; col < width; col++) {
            idx = getbithuff(stream, *huff, huff + 1);
            len = idx & 15;
            shl = idx >> 4;
            diff = ((getbithuff(stream, len - shl, 0) << 1) + 1) << shl >> 1;
            if (len > 0 && (diff & (1 << (len - 1))) == 0) {
                diff -= (1 << len) - !shl;
            }
            if (col < 2) {
                hpred[col] = vpred[row & 1][col] += diff;
            } else {
                hpred[col & 1] += diff;
            }
            image[row * width + col] = hpred[col & 1];
        }
    }
    free(huff);
}

DataStream::DataStream(const char *fname) : filename(nullptr), f(nullptr), _fsize(0) {
    if (fname && strlen(fname) > 0) {
        // 动态分配内存给 filename
        filename = new char[strlen(fname) + 1];
        strcpy(filename, fname);
        
        struct stat st;
        if (stat(filename, &st) == 0) {
            _fsize = st.st_size;
        }
        f = fopen(filename, "rb");
    }
}

DataStream::~DataStream() {
    if (f) {
        fclose(f);
    }
    delete[] filename;
}

int DataStream::valid() {
    return f ? 1 : 0; 
}

int DataStream::read(void *ptr, size_t size, size_t nmemb) {
    return int(fread(ptr, size, nmemb, f));
}

int DataStream::eof() {
    return feof(f);
}

int DataStream::seek(int64_t o, int whence) {
    return fseeko(f, o, whence);
}

int64_t DataStream::tell() {
    return ftello(f);
}

char* DataStream::gets(char *str, int sz) {
    if (sz < 1) {
        return NULL;
    }
    return fgets(str, sz, f); 
}

int DataStream::scanf_one(const char *fmt, void *val) {
    return fscanf(f, fmt, val);
}

const char *DataStream::fname() {
    return (filename && strlen(filename) > 0) ? filename : nullptr;
}
