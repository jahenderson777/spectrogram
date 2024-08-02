#define GUI_API CLAP_WINDOW_API_COCOA

struct GUI {
	void *mainView;
	uint32_t *bits;
};

extern "C" void *MacInitialise(Plugin *plugin, uint32_t *bits, uint32_t width, uint32_t height);
extern "C" void MacDestroy(void *mainView);
extern "C" void MacSetParent(void *_mainView, void *_parentView);
extern "C" void MacSetVisible(void *_mainView, bool show);
extern "C" void MacPaint(void *_mainView);

static void GUIPaint(Plugin *plugin, bool internal) {
	if (internal) plugin->paint(plugin->gui->bits);
	MacPaint(plugin->gui->mainView);
}

static void GUICreate(Plugin *plugin) {
	assert(!plugin->gui);
	plugin->gui = (GUI *) calloc(1, sizeof(GUI));
	plugin->gui->bits = (uint32_t *) calloc(1, GUI_WIDTH * GUI_HEIGHT * 4);
	plugin->paint(plugin->gui->bits);
	plugin->gui->mainView = MacInitialise(plugin, plugin->gui->bits, GUI_WIDTH, GUI_HEIGHT);
}

static void GUIDestroy(Plugin *plugin) {
	assert(plugin->gui);
	MacDestroy(plugin->gui->mainView);
	free(plugin->gui->bits);
	free(plugin->gui);
	plugin->gui = nullptr;
}

static void GUISetParent(Plugin *plugin, const clap_window_t *parent) { MacSetParent(plugin->gui->mainView, parent->cocoa); }
static void GUISetVisible(Plugin *plugin, bool visible) { MacSetVisible(plugin->gui->mainView, visible); }
static void GUIOnPOSIXFD(Plugin *) {}

extern "C" void MacInputEvent(Plugin *plugin, int32_t cursorX, int32_t cursorY, int8_t button) {
	if (button == -1) plugin->processMouseRelease();
	if (button ==  0) plugin->processMouseDrag   (cursorX, GUI_HEIGHT - 1 - cursorY);
	if (button ==  1) plugin->processMousePress  (cursorX, GUI_HEIGHT - 1 - cursorY);
	GUIPaint(plugin, true);
}
