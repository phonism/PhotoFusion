//
//  PhotoFusionBridge.m
//  RawDecoderDemo
//
//  Created by luqi03 on 2024/8/7.
//

#import "photofusion_api.h"
#import "raw_image.h"
#import "PhotoFusion-Bridging-Header.h"


@implementation Raw2TiffWrapper

- (void)convertRawToTiff:(NSString *)rawFilePath tiffFilePath:(NSString *)tiffFilePath {
    const char *rawPath = [rawFilePath cStringUsingEncoding:NSUTF8StringEncoding];
    const char *tiffPath = [tiffFilePath cStringUsingEncoding:NSUTF8StringEncoding];
    raw2tiff(rawPath, tiffPath);
}

@end
