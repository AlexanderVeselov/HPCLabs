#ifndef TIME_PROFILER_HPP_
#define TIME_PROFILER_HPP_

#include <chrono>
#include <iostream>
#include <string>

class TimeProfiler
{
public:
    TimeProfiler(std::string const& section_name)
        : section_name_(section_name)
    {
        start_ = std::chrono::high_resolution_clock::now();
    }

    ~TimeProfiler()
    {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start_).count();

        std::cout << "(TIME PROFILER) Execution of " << section_name_ << " took " << elapsed << " s." << std::endl;
    }


private:
    std::chrono::high_resolution_clock::time_point start_;
    std::string section_name_;

};

#ifdef NDEBUG
#define PROFILE_TIME(section_name) (void)0
#else
#define PROFILE_TIME(section_name) TimeProfiler __profiler(section_name)
#endif

#endif // TIME_PROFILER_HPP_
