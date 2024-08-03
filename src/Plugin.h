#pragma once
#include <clap/helpers/plugin.hh>
#include <atomic>
#include <thread>
#include <chrono>
#include <future>
#include <vector>
#include <mutex>
#include <Accelerate/Accelerate.h>
#include <array>
#include <cmath>


// GUI size.
#define GUI_WIDTH (512)
#define GUI_HEIGHT (512)


class Plugin : public clap::helpers::Plugin<clap::helpers::MisbehaviourHandler::Terminate,
                                            clap::helpers::CheckingLevel::Maximal>
{
public:
    static constexpr const char* features[5] =
    {
        CLAP_PLUGIN_FEATURE_AUDIO_EFFECT,
        CLAP_PLUGIN_FEATURE_UTILITY,
        CLAP_PLUGIN_FEATURE_STEREO,
        CLAP_PLUGIN_FEATURE_INSTRUMENT,
        nullptr
    };

    static constexpr clap_plugin_descriptor descriptor =
    {
        .clap_version = CLAP_VERSION,
        .id = "org.earthics.Spectrogram",
        .name = "Spectrogram",
        .vendor = "Alexander Henderson",
        .url = "https://earthics.org",
        .manual_url = "",
        .support_url = "",
        .version = "1.0.0",
        .description = "Spectrogram",
        .features = features
    };


    Plugin(const clap_host* host);
    ~Plugin();

    bool init() noexcept override;

    bool implementsAudioPorts() const noexcept override
    {
        return true;
    }

    uint32_t audioPortsCount(bool /*isInput*/) const noexcept override
    {
        return 1;
    }

    bool audioPortsInfo(uint32_t index, bool isInput, clap_audio_port_info* info) const noexcept override;


    bool implementsParams() const noexcept override
    {
        return true;
    }

    bool isValidParamId(clap_id paramId) const noexcept override
    {
        return paramId == gainPrmId_;
    }

    uint32_t paramsCount() const noexcept override
    {
        return 1;
    }

    bool paramsInfo(uint32_t paramIndex, clap_param_info* info) const noexcept override;

    bool paramsValue(clap_id paramId, double* value) noexcept override;

    bool paramsValueToText(clap_id paramId, double value, char* display, uint32_t size) noexcept override;
    bool paramsTextToValue(clap_id paramId, const char* display, double* value) noexcept override;

    bool activate(double sampleRate, uint32_t, uint32_t) noexcept override;

    clap_process_status process(const clap_process* process) noexcept override;

    void paintRectangle(uint32_t *bits, uint32_t l, uint32_t r, uint32_t t, uint32_t b, uint32_t border, uint32_t fill);

    void paint(uint32_t *bits);

    void processMouseDrag(int32_t x, int32_t y);

    void processMousePress(int32_t x, int32_t y);

    void processMouseRelease();

    struct GUI *gui;


    bool implementsGui() const noexcept override;
    bool guiIsApiSupported(const char *api, bool isFloating) noexcept override;
    bool guiGetPreferredApi(const char **api, bool *is_floating) noexcept override;
    bool guiCreate(const char *api, bool isFloating) noexcept override;
    void guiDestroy() noexcept override;
    bool guiSetScale(double scale) noexcept override;
    bool guiShow() noexcept override;
    bool guiHide() noexcept override;
    bool guiGetSize(uint32_t *width, uint32_t *height) noexcept override;
    bool guiCanResize() const noexcept override;
    bool guiGetResizeHints(clap_gui_resize_hints_t *hints) noexcept override;
    bool guiAdjustSize(uint32_t *width, uint32_t *height) noexcept override;
    bool guiSetSize(uint32_t width, uint32_t height) noexcept override;
    void guiSuggestTitle(const char *title) noexcept override;
    bool guiSetParent(const clap_window *window) noexcept override;
    bool guiSetTransient(const clap_window *window) noexcept override;


    //---------------------------//
    // clap_plugin_timer_support //
    //---------------------------//
    bool implementsTimerSupport() const noexcept override;
    void onTimer(clap_id timerId) noexcept override;

    void startAnimationLoop(int fps);
    void stopAnimationLoop();

    bool implementsNotePorts() const noexcept override
    {
        return true;
    }

    uint32_t notePortsCount (bool isInput) const noexcept override
    {
        return isInput ? 1 : 0;
    }

    bool notePortsInfo (uint32_t index, bool isInput, clap_note_port_info *info) const noexcept override;

private:
    static constexpr clap_id gainPrmId_ = 2345;

    double gain_ = 0.0;
    double targetGain_ = 0.0;
    double stepSizeToTargetGain_ = 0.0;
    int fadeLengthInSamples_ = 0;
    int currentFadeIndex_ = fadeLengthInSamples_;

    int x = 0;

    clap_id timerId;

    std::thread animationThread;
    std::atomic<bool> running;

    // Buffer to store audio samples
    std::vector<float> audioBuffer;
    std::mutex bufferMutex;
    size_t bufferIndex = 0;
    static const size_t bufferSize = 2048 * 16; // Adjust the size as needed

    
    const int fftSize = 2048;// 4096;
    const int fftOverlap = fftSize / 4;
    int fftX;
    std::vector<float> fftBuffer; //(fftSize, 0.0f);
    FFTSetup fftSetup;
    DSPSplitComplex fftResult;
    float *fftWindow;
    float *fftBufferWindowed;
    float *fftPrevMagnitudes;
    float *fftMagnitudes;

    void paintVerticalLine(uint32_t* bits, size_t x, const float* magnitudes, size_t numBins);
    void paintInterpolatedVerticalLines(uint32_t* bits, size_t x, const float* prevMagnitudes, const float* currMagnitudes, size_t numBins, size_t numLines);

    static constexpr std::array<std::array<float, 3>, 256> jetColormap = [] {
        std::array<std::array<float, 3>, 256> colormap = {};
        for (int i = 0; i < 256; ++i) {
            float r = std::min(1.0f, std::max(0.0f, 1.5f - std::abs(4.0f * (i / 255.0f) - 3.0f)));
            float g = std::min(1.0f, std::max(0.0f, 1.5f - std::abs(4.0f * (i / 255.0f) - 2.0f)));
            float b = std::min(1.0f, std::max(0.0f, 1.5f - std::abs(4.0f * (i / 255.0f) - 1.0f)));
            colormap[i] = {r, g, b};
        }
        return colormap;
    }();

    uint32_t mapValueToColor(float value);
    uint32_t getMatlabRgb(float ordinal);

    int filterLength = 101;
    float* filter;
    void createLowPassFIRFilter();

};
