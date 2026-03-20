#import <Foundation/Foundation.h>
#import <Cocoa/Cocoa.h>
#include <pthread.h>

struct MyPlugin;
void MacInputEvent(struct MyPlugin *plugin, int32_t cursorX, int32_t cursorY, int8_t button);

@interface MainView : NSView
@property (nonatomic) struct MyPlugin *plugin;
@property (nonatomic) uint32_t *bits;
@property (nonatomic) uint32_t width;
@property (nonatomic) uint32_t height;
@property (nonatomic) BOOL hasSuperView;
@property (nonatomic, copy) NSString *debugText;
@property (nonatomic) pthread_mutex_t *bitsLock; // shared with Plugin::paint()
@end

@implementation MainView
- (void)drawRect:(NSRect)dirtyRect {
	// Lock the pixel buffer for the duration of this read so that
	// Plugin::paint() (running on the animation thread) cannot write
	// to bits[] while we are drawing from it.
	if (_bitsLock) pthread_mutex_lock(_bitsLock);

	const unsigned char *data = (const unsigned char *) _bits;
	NSDrawBitmap(self.bounds, _width, _height, 8 /* bits per channel */, 4 /* channels per pixel */, 32 /* bits per pixel */,
			4 * _width /* bytes per row */, NO /* planar */, NO /* has alpha */, 
			NSDeviceRGBColorSpace /* color space */, &data /* data */);

	// Draw debug text overlay (green, top-left).
	// debugText is written on the main thread via dispatch_async, so it cannot
	// race with this drawRect: — both run serially on the main thread.
	if (_debugText) {
		NSDictionary *attrs = @{
			NSFontAttributeName: [NSFont monospacedSystemFontOfSize:11 weight:NSFontWeightRegular],
			NSForegroundColorAttributeName: [NSColor colorWithSRGBRed:0.0 green:1.0 blue:0.0 alpha:0.85],
		};
		[_debugText drawAtPoint:NSMakePoint(4, _height - 16) withAttributes:attrs];
	}

	if (_bitsLock) pthread_mutex_unlock(_bitsLock);
}

- (BOOL)acceptsFirstMouse:(NSEvent *)event {
	return YES;
}

- (void)mouseDown:(NSEvent *)event {
	NSPoint cursor = [self convertPoint:[event locationInWindow] fromView:nil];
	MacInputEvent(_plugin, cursor.x, cursor.y, 1);
}

- (void)mouseUp:(NSEvent *)event {
	NSPoint cursor = [self convertPoint:[event locationInWindow] fromView:nil];
	MacInputEvent(_plugin, cursor.x, cursor.y, -1);
}

- (void)mouseDragged:(NSEvent *)event {
	NSPoint cursor = [self convertPoint:[event locationInWindow] fromView:nil];
	MacInputEvent(_plugin, cursor.x, cursor.y, 0);
}
@end

void *MacInitialise(struct MyPlugin *plugin, uint32_t *bits, pthread_mutex_t *bitsLock, uint32_t width, uint32_t height) {
	NSRect frame;
	frame.origin.x = 0;
	frame.origin.y = 0;
	frame.size.width = width;
	frame.size.height = height;
	MainView *mainView = [[MainView alloc] initWithFrame:frame];
	mainView.plugin   = plugin;
	mainView.bits     = bits;
	mainView.bitsLock = bitsLock;
	mainView.width    = width;
	mainView.height   = height;
	return mainView;
}

void MacDestroy(void *mainView) {
	[((MainView *) mainView) release];
}

void MacSetParent(void *_mainView, void *_parentView) {
	MainView *mainView = (MainView *) _mainView;
	NSView *parentView = (NSView *) _parentView;
	if (mainView.hasSuperView) [mainView removeFromSuperview];
	[parentView addSubview:mainView];
	mainView.hasSuperView = true;
}

void MacSetVisible(void *_mainView, bool show) {
	MainView *mainView = (MainView *) _mainView;
	[mainView setHidden:(show ? NO : YES)];
}

void MacSetDebugText(void *_mainView, const char *text) {
	MainView *mainView = (MainView *) _mainView;
	// Convert to NSString on the calling thread — this copies the C string bytes
	// immediately, before the caller's stack frame (and its char buffer) is gone.
	NSString *str = [NSString stringWithUTF8String:text];
	dispatch_async(dispatch_get_main_queue(), ^{
		mainView.debugText = str;  // block retains both mainView and str
	});
}

void MacPaint(void *_mainView) {
	MainView *mainView = (MainView *) _mainView;
	// AppKit must only be touched from the main thread; dispatch_async
	// marshals setNeedsDisplayInRect: there safely from any background thread.
	dispatch_async(dispatch_get_main_queue(), ^{
		[mainView setNeedsDisplayInRect:mainView.bounds];  // block retains mainView
	});
}
