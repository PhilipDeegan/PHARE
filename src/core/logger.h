
#ifndef PHARE_CORE_LOGGER_H
#define PHARE_CORE_LOGGER_H


#include <tuple>
#include <chrono>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iostream>
#include <unordered_map>


#include "thread_queue.h"
#include "initializer/data_provider.h"

namespace PHARE::core
{
inline auto get_thread_id()
{
    std::ostringstream os;
    os << std::hex << pthread_self();
    return os.str();
}

inline auto now_in_nanos()
{
    return static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::system_clock::now().time_since_epoch())
                                        .count());
}

struct Logbject
{
    Logbject() = delete;

    unsigned line;    //    = __LINE__;
    char const* file; // = __FILE__;
    std::string stuff;
};

class LogWriter
{
public:
    virtual void operator<<(std::string const& str) = 0;
};

class StdOutLogWriter : public LogWriter
{
public:
    void operator<<(std::string const& str) override { std::cout << str << std::endl; }
};

class FileStringLogWriter : public LogWriter
{
public:
    FileStringLogWriter(std::string filePath)
        : fileWriter{filePath + "/Logger." + get_thread_id() + ".txt"}
    {
    }

    void operator<<(std::string const& str) override { fileWriter << str << std::endl; }

private:
    std::ofstream fileWriter;
};

class NanoLogWriter : public LogWriter
{
public:
    NanoLogWriter() {}

    void operator<<(std::string const& str) override
    { /*fileWriter << str << std::endl;*/
    }

private:
    // std::ofstream fileWriter;
};


class Logger
{
public:
    // Logger() = default;
    Logger(std::unique_ptr<LogWriter>&& logWriter)
        : logWriter_{std::move(logWriter)}
    {
    }

    struct StartStopLogger
    {
        StartStopLogger(std::string const& key)
            : key_{key}
        {
        }
        std::string key_;
        std::size_t start_{now_in_nanos()};
    };

    void start(std::string const& key) { nestings.emplace_back(key); }

    void stop()
    {
        auto& startStop = nestings.back();
        log(indent(nestings.size()) + startStop.key_ + " " + Logger::mtimer(startStop.start_));
        nestings.pop_back();
    }

    struct ScopeLogger
    {
        ScopeLogger(Logger& outer, std::string const& key)
            : outer_{outer}
        {
            outer_.start(key);
        }
        ~ScopeLogger() { outer_.stop(); }

        Logger& outer_;
    };

    auto scope(std::string const& key) { return ScopeLogger{*this, key}; }


private:
    Logger(Logger const&) = delete;
    Logger(Logger&&)      = delete;

    static std::string mtimer(std::size_t nanos) { return std::to_string(now_in_nanos() - nanos); }
    static std::string indent(std::size_t size) { return std::string(size, ' '); }

    void log(std::string const& s0)
    {
        queue.enqueue([this](std::string s1) { (*this->logWriter_) << s1; }, s0);
    }
    void log(std::string&& s) { log(s); }

    std::unique_ptr<LogWriter> logWriter_;
    ThreadQueue queue;
    std::vector<StartStopLogger> nestings;
};


struct LogWriterFactory
{
    static auto make(PHARE::initializer::PHAREDict const& dict)
    {
        return std::make_unique<FileStringLogWriter>(".");
    }
};



class LogMan
{
public: // static
    static LogMan& start(PHARE::initializer::PHAREDict const& dict)
    {
        assert(!self);

        return *(self = std::make_unique<LogMan>(LogWriterFactory::make(dict)));
    }
    static void kill() { self.release(); }

    static LogMan& get() { return *self; }

public:
    LogMan(std::unique_ptr<LogWriter>&& logWriter)
        : logger{std::move(logWriter)}
    {
    }

    void start(Logbject&& lob) { logger.start(lob.stuff); }
    void stop() { logger.stop(); }
    auto scope(Logbject&& lob) { return logger.scope(lob.stuff); }

private:
    Logger logger;

    static std::unique_ptr<LogMan> self;
};

} // namespace PHARE::core

#define PHARE_LOBJECT(str)                                                                         \
    PHARE::core::Logbject { __LINE__, __FILE__, str }

#define PHARE_LOG_START(str) PHARE::core::LogMan::get().start(PHARE_LOBJECT(str))
#define PHARE_LOG_STOP() PHARE::core::LogMan::get().stop()
#define PHARE_LOG_SCOPE(str) auto _scopeLog = PHARE::core::LogMan::get().scope(PHARE_LOBJECT(str))
#define PHARE_LOG_NAMED_SCOPE(name, str)                                                           \
    auto name = PHARE::core::LogMan::get().scope(PHARE_LOBJECT(str))


#endif /* PHARE_CORE_LOGGER_H */
