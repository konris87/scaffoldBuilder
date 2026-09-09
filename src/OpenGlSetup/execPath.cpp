// Resolves the running executable's own directory and makes it the current
// working directory. Kept in its own translation unit on purpose: it needs
// <windows.h> on Windows, whose ERROR / min / max / near / far macros would
// otherwise break Logger.h, Eigen and the OpenGlSetup camera headers if pulled
// into main.cpp. Nothing project-specific is included here.

#include <filesystem>
#include <string>
#include <system_error>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <cstdint>
#else
#include <unistd.h>
#endif

void set_cwd_to_executable_dir() {
    namespace fs = std::filesystem;
    fs::path exePath;

#ifdef _WIN32
    wchar_t buf[MAX_PATH];
    DWORD n = GetModuleFileNameW(NULL, buf, MAX_PATH);
    if (n > 0 && n < MAX_PATH) exePath = fs::path(std::wstring(buf, buf + n));
#elif defined(__APPLE__)
    char buf[4096];
    uint32_t size = sizeof(buf);
    if (_NSGetExecutablePath(buf, &size) == 0) exePath = fs::path(buf);
#else
    char buf[4096];
    ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf));
    if (n > 0) exePath = fs::path(std::string(buf, buf + n));
#endif

    if (!exePath.empty()) {
        std::error_code ec;
        fs::current_path(exePath.parent_path(), ec);
    }
}
