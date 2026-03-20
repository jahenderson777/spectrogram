#include "Plugin.h"
#include <clap/helpers/host-proxy.hxx>
#include <clap/helpers/plugin.hxx>
#include "Utils.h"
#include <cstring>
#include "tinycolormap.hpp"

#if defined(_WIN32)
#include "gui_w32.cpp"
#elif defined(__linux__)
#include "gui_x11.cpp"
#elif defined(__APPLE__)
#include "gui_mac.cpp"
#endif


Plugin::Plugin(const clap_host* host) 
    : clap::helpers::Plugin<clap::helpers::MisbehaviourHandler::Terminate, clap::helpers::CheckingLevel::Maximal>(
        &descriptor, host)
        , running(false)
{
}

Plugin::~Plugin()
{
    // Defensive cleanup in case deactivate() was not called.
    if (fftSetup) vDSP_destroy_fftsetup(fftSetup);
    free(fftResult.realp);
    free(fftResult.imagp);
    free(fftWindow);
    free(fftBufferWindowed);
    free(fftPrevMagnitudes);
    free(fftMagnitudes);
    free(filter);
}

bool Plugin::init() noexcept {
    return true;
}

bool Plugin::audioPortsInfo(uint32_t index, bool /*isInput*/, clap_audio_port_info* info) const noexcept
{
    if (index != 0)
        return false;

    info->id = 0;
    info->in_place_pair = CLAP_INVALID_ID;
    strncpy(info->name, "main", sizeof(info->name));
    info->flags = CLAP_AUDIO_PORT_IS_MAIN;
    info->channel_count = 2;
    info->port_type = CLAP_PORT_STEREO;

    return true;
}

bool Plugin::paramsInfo(uint32_t paramIndex, clap_param_info* info) const noexcept
{
    if (paramIndex >= 1)
        return false;

    info->id = gainPrmId_;
    info->flags = CLAP_PARAM_IS_AUTOMATABLE;
    strncpy(info->name, "Gain", CLAP_NAME_SIZE);
    strncpy(info->module, "", CLAP_NAME_SIZE);
    info->min_value = 0.0;
    info->max_value = 1.0;
    info->default_value = utils::decibelsToGain(0.0);

    return true;
}

bool Plugin::paramsValue(clap_id paramId, double* value) noexcept
{
    if (paramId != gainPrmId_)
        return false;

    *value = utils::toLinearCurve(gain_);
    return true;
}

bool Plugin::paramsValueToText(clap_id paramId, double value, char* display, uint32_t size) noexcept
{
    if (paramId != gainPrmId_)
        return false;

    const auto valueIndB = utils::gainToDecibels(utils::toExponentialCurve(value));

    if (valueIndB <= utils::minusInfinitydB)
    {
        snprintf(display, size, "-inf dB");
    }
    else
    {
        snprintf(display, size, "%.2f dB", valueIndB);
    }

    return true;
}

bool Plugin::paramsTextToValue(clap_id paramId, const char* display, double* value) noexcept
{
    if (paramId != gainPrmId_)
        return false;

    const double value_ = strtod(display, nullptr);
    *value = utils::toLinearCurve(utils::decibelsToGain(value_));

    return true;
}

bool Plugin::activate(double sampleRate, uint32_t /*minFrameCount*/, uint32_t /*maxFrameCount*/) noexcept
{
    const auto samplesInOneMs = sampleRate / 1000;
    const auto fadeLengthInMs = 5;
    fadeLengthInSamples_ = samplesInOneMs * fadeLengthInMs;

    audioBuffer.resize(bufferSize, 0.0f);
    fftBuffer.resize(fftSize, 0.0f);
    bufferIndex = 0;
    fftX = 0;

    fftSetup = vDSP_create_fftsetup(log2f(fftSize), kFFTRadix2);
    fftResult.realp = (float*)malloc(fftSize / 2 * sizeof(float));
    fftResult.imagp = (float*)malloc(fftSize / 2 * sizeof(float));
    fftWindow = (float*)malloc(fftSize * sizeof(float));
    fftBufferWindowed = (float*)malloc(fftSize * sizeof(float));
    fftPrevMagnitudes = (float*)malloc(fftSize / 2 * sizeof(float));
    fftMagnitudes = (float*)malloc(fftSize / 2 * sizeof(float));

    memset(fftPrevMagnitudes, 0, fftSize / 2 * sizeof(float));
    memset(fftMagnitudes, 0, fftSize / 2 * sizeof(float));
    paintPrevMagnitudes_.assign(fftSize / 2, 0.0f);

    vDSP_hann_window(fftWindow, fftSize, vDSP_HANN_NORM);

    filter = (float*)malloc(filterLength * sizeof(float));
    createLowPassFIRFilter();

    return true;
}

void Plugin::deactivate() noexcept
{
    if (fftSetup) {
        vDSP_destroy_fftsetup(fftSetup);
        fftSetup = nullptr;
    }
    free(fftResult.realp);  fftResult.realp = nullptr;
    free(fftResult.imagp);  fftResult.imagp = nullptr;
    free(fftWindow);        fftWindow = nullptr;
    free(fftBufferWindowed); fftBufferWindowed = nullptr;
    free(fftPrevMagnitudes); fftPrevMagnitudes = nullptr;
    free(fftMagnitudes);    fftMagnitudes = nullptr;
    free(filter);           filter = nullptr;

    fftBuffer.clear();
    audioBuffer.clear();
}

clap_process_status Plugin::process(const clap_process* process) noexcept
{
    if (process->audio_outputs_count <= 0)
        return CLAP_PROCESS_CONTINUE;

    float** input = process->audio_inputs[0].data32;
    float** output = process->audio_outputs[0].data32;
    const auto outputChannelsCount = process->audio_outputs->channel_count;

    auto event = process->in_events;
    auto eventsSize = event->size(event);

    const clap_event_header_t* nextEvent{nullptr};
    uint32_t nextEventIndex{0};
    if (eventsSize != 0)
    {
        nextEvent = event->get(event, nextEventIndex);
    }

    for (uint32_t index = 0; index < process->frames_count; index++)
    {
        for (uint32_t channel = 0; channel < outputChannelsCount; ++channel)
        {
            output[channel][index] = input[channel][index];
        }
    }

    // Handle MIDI events (audio thread — no GUI work here)
    for (uint32_t i = 0; i < process->in_events->size(process->in_events); ++i) {
        const clap_event_header_t *event = process->in_events->get(process->in_events, i);
        if (event->type == CLAP_EVENT_MIDI) {
            const clap_event_midi_t *midiEvent = reinterpret_cast<const clap_event_midi_t*>(event);
            uint8_t status   = midiEvent->data[0] & 0xF0;
            uint8_t velocity = midiEvent->data[2];

            if (status == 0x90 && velocity > 0) {
                bufferIndex = 0; // Reset FFT accumulator on note-on
                // Tell the paint thread to clear the display and reset the draw cursor
                std::lock_guard<std::mutex> lock(fftQueueMutex_);
                fftQueue_.push({std::vector<float>(), true});
                break;
            }
        }
    }

    // Accumulate samples, run FFT, and push magnitude frames for the paint thread
    if (fftSetup && !fftBuffer.empty())
    {
        float* fftInput = process->audio_inputs[0].data32[0];
        size_t remaining = process->frames_count;
        while (remaining > 0) {
            size_t toCopy = std::min(remaining, (size_t)(fftSize - bufferIndex));
            std::copy(fftInput, fftInput + toCopy, fftBuffer.begin() + bufferIndex);
            bufferIndex += toCopy;
            fftInput    += toCopy;
            remaining   -= toCopy;

            if (bufferIndex == fftSize) {
                // Apply Hann window
                for (size_t i = 0; i < fftSize; ++i)
                    fftBufferWindowed[i] = fftBuffer.data()[i] * fftWindow[i];

                // Perform FFT
                vDSP_ctoz((DSPComplex*)fftBufferWindowed, 2, &fftResult, 1, fftSize / 2);
                vDSP_fft_zrip(fftSetup, &fftResult, 1, log2f(fftSize), FFT_FORWARD);

                // Convert to magnitude (dB)
                vDSP_zvmags(&fftResult, 1, fftMagnitudes, 1, fftSize / 2);
                float zero = 1.0f;
                vDSP_vdbcon(fftMagnitudes, 1, &zero, fftMagnitudes, 1, fftSize / 2, 0);

                // Hand the magnitude data off to the paint thread via the queue
                {
                    std::lock_guard<std::mutex> lock(fftQueueMutex_);
                    fftQueue_.push({std::vector<float>(fftMagnitudes, fftMagnitudes + fftSize / 2), false});
                }

                // Slide buffer for 75% overlap
                std::copy(fftBuffer.begin() + fftOverlap, fftBuffer.end(), fftBuffer.begin());
                bufferIndex = fftSize - fftOverlap;
            }
        }
    }

    return CLAP_PROCESS_CONTINUE;
}


/*
bandStart = 5120
bandEnd = 20480
numLines = 8

while bandEnd >= 20
    apply a 512 window to the samples
    apply a 512 fft to the samples
    plot the fft result, bandStart, bandEnd, fftResult, fftX, numLines
    low pass filter the samples
    reduce sample rate by 4
    bandStart *= 0.25
    bandEnd *= 0.25

// Slide buffer for overlap
std::copy(fftBuffer.begin() + fftOverlap, fftBuffer.end(), fftBuffer.begin());
bufferIndex = fftSize - fftOverlap;

*/

/*void Plugin::paintVerticalLine(uint32_t* bits, size_t x, const float* magnitudes, size_t numBins, uint32_t width, uint32_t height) {
    for (size_t y = 0; y < GUI_HEIGHT; ++y) {
        float mag = magnitudes[y];
        float normalizedMag = 1.0f * (mag / 100.0f) ;

        const tinycolormap::Color c = tinycolormap::GetColor(normalizedMag, tinycolormap::ColormapType::Magma);

        uint32_t color = ((uint8_t)(c.b()*256.0f) << 16) | ((uint8_t)(c.g()*256.0f) << 8) | (uint8_t)(c.r()*256.0f);
        
        // Plot the pixel
        size_t pixelIndex = (GUI_HEIGHT - y) * GUI_WIDTH + (x % GUI_WIDTH); // Invert y for correct orientation
        bits[pixelIndex] = color;
    }
}*/

/*void Plugin::paintVerticalLine(uint32_t* bits, size_t x, const float* magnitudes, size_t numBins, uint32_t width, uint32_t height) {
    // Define the frequency range for the spectrogram (e.g., 20 Hz to 20 kHz)
    const float minFrequency = 20.0f;
    const float maxFrequency = 20000.0f;

    // Compute the log scale factors
    const float logMin = std::log10(minFrequency);
    const float logMax = std::log10(maxFrequency);
    const float logRange = logMax - logMin;

    for (size_t y = 0; y < height; ++y) {
        // Map y-coordinate to a frequency in the logarithmic scale
        float logFrequency = logMin + (logRange * y / height);
        float frequency = std::pow(10.0f, logFrequency);

        // Map frequency to the corresponding bin index
        size_t binIndex = static_cast<size_t>((frequency - minFrequency) / (maxFrequency - minFrequency) * numBins);
        binIndex = std::min(std::max(binIndex, static_cast<size_t>(0)), numBins - 1); // Clamp the bin index

        // Get the magnitude for the current frequency bin
        float mag = magnitudes[binIndex];
        float normalizedMag = 1.0f * (mag / 100.0f);

        // Get the color based on the normalized magnitude
        const tinycolormap::Color c = tinycolormap::GetColor(normalizedMag, tinycolormap::ColormapType::Magma);
        uint32_t color = ((uint8_t)(c.b() * 256.0f) << 16) | ((uint8_t)(c.g() * 256.0f) << 8) | (uint8_t)(c.r() * 256.0f);

        // Plot the pixel
        size_t pixelIndex = (GUI_HEIGHT - y - 1) * GUI_WIDTH + x; // Invert y for correct orientation
        bits[pixelIndex] = color;
    }
}*/

static float hertzToMel(float hertz) {
    //return 2595.0f * std::log10(1.0f + hertz / 700.0f);
    //return 2100.0f * std::log10(1.0f + hertz / 300.0f);
    //return 1800.0f * std::log10(1.0f + hertz / 150.0f);
    return 901.0f * std::log10(1.0f + hertz);
}

static float melToHertz(float mel) {
    return 1.0f * (std::pow(10.0f, mel / 901.0f) - 1.0f);
}

void Plugin::paintVerticalLine(uint32_t* bits, size_t x, const float* magnitudes, size_t numBins) {
    const float minFrequency = 20.0f;
    const float maxFrequency = 20000.0f;

    const float minMel = hertzToMel(minFrequency);
    const float maxMel = hertzToMel(maxFrequency);
    const float melRange = maxMel - minMel;

    for (size_t y = 0; y < GUI_HEIGHT; ++y) {
        float mel = minMel + (melRange * y / GUI_HEIGHT);
        float frequency = melToHertz(mel);

        // Map frequency to the corresponding linear bin indices
        float linearBinIndex = (frequency - minFrequency) / (maxFrequency - minFrequency) * numBins;
        size_t binIndexLow = static_cast<size_t>(std::floor(linearBinIndex));
        size_t binIndexHigh = static_cast<size_t>(std::ceil(linearBinIndex));
        binIndexLow = std::min(std::max(binIndexLow, static_cast<size_t>(0)), numBins - 1);
        binIndexHigh = std::min(std::max(binIndexHigh, static_cast<size_t>(0)), numBins - 1);

        // Interpolate the magnitude between the two adjacent bins
        float lowMag = magnitudes[binIndexLow];
        float highMag = magnitudes[binIndexHigh];
        float interpolationFactor = linearBinIndex - binIndexLow;
        float mag = (1.0f - interpolationFactor) * lowMag + interpolationFactor * highMag;
        float normalizedMag = mag / 100.0f;

        const tinycolormap::Color c = tinycolormap::GetColor(normalizedMag, tinycolormap::ColormapType::Magma);
        uint32_t color = ((uint8_t)(c.b() * 256.0f) << 16) | ((uint8_t)(c.g() * 256.0f) << 8) | (uint8_t)(c.r() * 256.0f);

        size_t pixelIndex = (GUI_HEIGHT - y - 1) * GUI_WIDTH + x;
        bits[pixelIndex] = color;
    }
}

void Plugin::paintInterpolatedVerticalLines(uint32_t* bits, size_t x, const float* prevMagnitudes, const float* currMagnitudes, size_t numBins, size_t numLines) {
    const float minFrequency = 25.0f;
    const float maxFrequency = 20000.0f;

    const float minMel = hertzToMel(minFrequency);
    const float maxMel = hertzToMel(maxFrequency);
    const float melRange = maxMel - minMel;

    for (size_t i = 0; i < numLines; ++i) {
        float xInterpolationFactor = static_cast<float>(i) / static_cast<float>(numLines);

        for (size_t y = 0; y < GUI_HEIGHT; ++y) {
            float mel = minMel + (melRange * y / GUI_HEIGHT);
            float frequency = melToHertz(mel);

            // Map frequency to the corresponding linear bin indices
            float linearBinIndex = (frequency - minFrequency) / (maxFrequency - minFrequency) * numBins;
            size_t binIndexLow = static_cast<size_t>(std::floor(linearBinIndex));
            size_t binIndexHigh = static_cast<size_t>(std::ceil(linearBinIndex));
            binIndexLow = std::min(std::max(binIndexLow, static_cast<size_t>(0)), numBins - 1);
            binIndexHigh = std::min(std::max(binIndexHigh, static_cast<size_t>(0)), numBins - 1);

            // Interpolate the magnitude between the two adjacent bins
            float prevLowMag = prevMagnitudes[binIndexLow];
            float prevHighMag = prevMagnitudes[binIndexHigh];
            float currLowMag = currMagnitudes[binIndexLow];
            float currHighMag = currMagnitudes[binIndexHigh];
            float yInterpolationFactor = linearBinIndex - binIndexLow;
            if (yInterpolationFactor < 0.0f) 
                yInterpolationFactor = 0.0f;
            else if (yInterpolationFactor > 1.0f) yInterpolationFactor = 1.0f;

            float prevMag = prevLowMag;// (1.0f - yInterpolationFactor) * prevLowMag + yInterpolationFactor * prevHighMag;
            float currMag = currLowMag;// (1.0f - yInterpolationFactor) * currLowMag + yInterpolationFactor * currHighMag;
            float mag = (1.0f - xInterpolationFactor) * prevMag + xInterpolationFactor * currMag;
            float normalizedMag = (1.2f + (float)y * 0.4f / (float)GUI_HEIGHT) * mag / 140.0f;
            //normalizedMag = std::pow(normalizedMag, 0.25);
            normalizedMag += 0.4f;
            if (normalizedMag < 0.0f) 
                normalizedMag = 0.0f;
            else if (normalizedMag > 1.0f) normalizedMag = 1.0f;


            const tinycolormap::Color c = tinycolormap::GetColor(normalizedMag, tinycolormap::ColormapType::Turbo);
            uint32_t color = ((uint8_t)(c.b() * 255.0f) << 16) | ((uint8_t)(c.g() * 255.0f) << 8) | (uint8_t)(c.r() * 255.0f);

            size_t pixelIndex = (GUI_HEIGHT - y - 1) * GUI_WIDTH + x + i;
            bits[pixelIndex] = color;
        }
    }
}

void Plugin::paintRectangle(uint32_t *bits, uint32_t l, uint32_t r, uint32_t t, uint32_t b, uint32_t border, uint32_t fill) {
	for (uint32_t i = t; i < b; i++) {
		for (uint32_t j = l; j < r; j++) {
			bits[i * GUI_WIDTH + j] = (i == t || i == b - 1 || j == l || j == r - 1) ? border : fill;
		}
	}
}

void Plugin::paint(uint32_t *bits) {
    if (!bits) return;

    // Drain the FFT frame queue that was filled by the audio thread.
    // All pixel-buffer writes happen here, on the GUI/animation thread.
    std::vector<FFTFrame> frames;
    {
        std::lock_guard<std::mutex> lock(fftQueueMutex_);
        while (!fftQueue_.empty()) {
            frames.push_back(std::move(fftQueue_.front()));
            fftQueue_.pop();
        }
    }

    if (frames.empty()) return; // nothing to paint; don't bother locking

    // Lock the pixel buffer while we write so that drawRect: (main thread)
    // cannot read a partially-updated frame.
    if (gui) pthread_mutex_lock(&gui->bitsMutex);

    for (auto& frame : frames) {
        if (frame.reset) {
            // MIDI note-on: clear the display and reset the draw cursor
            paintRectangle(bits, 0, GUI_WIDTH, 0, GUI_HEIGHT, 0x000000, 0x000000);
            fftX = 0;
            std::fill(paintPrevMagnitudes_.begin(), paintPrevMagnitudes_.end(), 0.0f);
        } else if (!frame.magnitudes.empty()) {
            // Lazy-initialise the previous-frame buffer if needed
            if (paintPrevMagnitudes_.size() != frame.magnitudes.size())
                paintPrevMagnitudes_.assign(frame.magnitudes.size(), 0.0f);

            if (fftX >= GUI_WIDTH) fftX = 0;
            paintInterpolatedVerticalLines(bits, fftX,
                                           paintPrevMagnitudes_.data(),
                                           frame.magnitudes.data(),
                                           frame.magnitudes.size(), 3);
            fftX += 3;
            paintPrevMagnitudes_ = std::move(frame.magnitudes);
        }
    }

    if (gui) pthread_mutex_unlock(&gui->bitsMutex);
}

void Plugin::processMouseDrag(int32_t x, int32_t y) {
	// TODO.
}

void Plugin::processMousePress(int32_t x, int32_t y) {
	// TODO.
}

void Plugin::processMouseRelease() {
	// TODO.
}

bool Plugin::implementsGui() const noexcept { return true; }

bool Plugin::guiIsApiSupported(const char *api, bool isFloating) noexcept {
    return 0 == strcmp(api, GUI_API) && !isFloating;
}

bool Plugin::guiGetPreferredApi(const char **api, bool *is_floating) noexcept {
    *api = GUI_API;
    *is_floating = false;
    return true;
}

bool Plugin::guiCreate(const char *api, bool isFloating) noexcept { 
    if (!guiIsApiSupported(api, isFloating)) return false;
    // We'll define GUICreate in our platform specific code file.
    /*if (_host.canUseTimerSupport()) {
        _host.timerSupportRegister(200, &timerId);
    }*/
    GUICreate((Plugin *) this);
    running = true;
    startAnimationLoop(30);
    return true;
}

void Plugin::guiDestroy() noexcept {
    // We'll define GUIDestroy in our platform specific code file.
    //_host.timerSupportUnregister(timerId);

    stopAnimationLoop();

    GUIDestroy((Plugin *) this);

    
}

bool Plugin::guiSetScale(double scale) noexcept { 
    return false; 
}

bool Plugin::guiGetSize(uint32_t *width, uint32_t *height) noexcept {
    *width = GUI_WIDTH;
    *height = GUI_HEIGHT;
    return true;
}

bool Plugin::guiCanResize() const noexcept {
    return false;
}

bool Plugin::guiGetResizeHints(clap_gui_resize_hints_t *hints) noexcept { return false; }

bool Plugin::guiAdjustSize(uint32_t *width, uint32_t *height) noexcept {
    return guiGetSize(width, height);
}

bool Plugin::guiSetSize(uint32_t width, uint32_t height) noexcept { return true; }

bool Plugin::guiSetParent(const clap_window *window) noexcept { 
    assert(0 == strcmp(window->api, GUI_API));
    // We'll define GUISetParent in our platform specific code file.
    GUISetParent((Plugin *) this, window);
    return true;
}

bool Plugin::guiSetTransient(const clap_window *window) noexcept { return false; }

void Plugin::guiSuggestTitle(const char *title) noexcept {}

bool Plugin::guiShow() noexcept { 
    // We'll define GUISetVisible in our platform specific code file.
    GUISetVisible((Plugin *) this, true);
    return true;
}

bool Plugin::guiHide() noexcept { 
    GUISetVisible((Plugin *) this, false);
    return true;
}

bool Plugin::implementsTimerSupport() const noexcept { return true; }

void Plugin::onTimer(clap_id timerId) noexcept {
   /* if (gui //&& PluginSyncAudioToMain(plugin)
        ) {
        // Repaint the GUI.
        GUIPaint(this, true);
    }*/
}

void Plugin::startAnimationLoop(int fps) {
    x = 1;
    auto interval = std::chrono::milliseconds(1000 / fps);
    animationThread = std::thread([this, interval]() {
        while (running) {
            std::this_thread::sleep_for(interval);
            GUIPaint(this, true);
        }
    });
}

void Plugin::stopAnimationLoop() {
    running = false;
    if (animationThread.joinable()) {
        animationThread.join();
    }
}

bool Plugin::stateSave(const clap_ostream *stream) noexcept {
    // Write the gain value as a double.
    const auto written = stream->write(stream, &gain_, sizeof(gain_));
    return written == sizeof(gain_);
}

bool Plugin::stateLoad(const clap_istream *stream) noexcept {
    double value = 0.0;
    const auto read = stream->read(stream, &value, sizeof(value));
    if (read != sizeof(value))
        return false;
    gain_ = value;
    return true;
}

bool Plugin::notePortsInfo (uint32_t index, bool isInput, clap_note_port_info *info) const noexcept
{
    if (isInput)
    {
        info->id = 1;
        info->supported_dialects = CLAP_NOTE_DIALECT_MIDI;// | CLAP_NOTE_DIALECT_CLAP;
        info->preferred_dialect = CLAP_NOTE_DIALECT_MIDI;
        strncpy(info->name, "NoteInput", CLAP_NAME_SIZE);
        return true;
    }
    return false;
}

uint32_t Plugin::getMatlabRgb(float ordinal)
{
    uint8_t r, g, b;
    r = (ordinal < 0.0)  ? (uint8_t)0 : (ordinal >= 0.5)  ? (uint8_t)255 : (uint8_t)(ordinal / 0.5 * 255);
    g = (ordinal < -0.5) ? (uint8_t)((ordinal + 1) / 0.5 * 255) : (ordinal > 0.5) ? (uint8_t)(255 - ((ordinal - 0.5) / 0.5 * 255)) : (uint8_t)255;
    b = (ordinal > 0.0)  ? (uint8_t)0 : (ordinal <= -0.5) ? (uint8_t)255 : (uint8_t)(ordinal * -1.0 / 0.5 * 255);
    return (r << 16) | (g << 8) | b;
}

// Function to create a low-pass FIR filter using a windowed sinc function
void Plugin::createLowPassFIRFilter() {
    float normalizedCutoff = 0.25;
    float M = filterLength - 1.0;
    for (int i = 0; i < filterLength; ++i) {
        if (i == M / 2.0) {
            filter[i] = normalizedCutoff;
        } else {
            filter[i] = normalizedCutoff * (sin(M_PI * normalizedCutoff * (i - M / 2.0)) / (M_PI * normalizedCutoff * (i - M / 2.0)));
        }
        // Apply a Hann window
        filter[i] *= 0.5 * (1.0 - cos(2.0 * M_PI * i / M));
    }
}
