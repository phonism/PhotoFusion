//
//  Use this file to import your target's public headers that you would like to expose to Swift.
//

#import <Foundation/Foundation.h>

@interface Raw2TiffWrapper : NSObject
- (void)convertRawToTiff:(NSString *)rawFilePath tiffFilePath:(NSString *)tiffFilePath;
@end
