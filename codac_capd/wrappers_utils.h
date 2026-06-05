#include <fstream>
#include <string>
#include <iostream>
#include <tuple>

using namespace std;

std::tuple<std::string, std::string, std::string> readIMapWrappers(const std::string &filename)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    std::string f_wrap;
    std::string f_wrap_disturb;
    std::string vectorField_wrap;

    while (std::getline(infile, line))
    {
        // Skip empty lines or comment lines
        if (line.empty() || line.rfind("###", 0) == 0)
            continue;

        // Detect which wrapper it is
        if (line.find("IMapWrapper f_wrap") != std::string::npos)
        {
            size_t start = line.find("\"");
            size_t end = line.rfind("\"");
            if (start != std::string::npos && end != std::string::npos && end > start)
            {
                f_wrap = line.substr(start + 1, end - start - 1);
            }
        }
        if (line.find("IMapWrapper f_disturb_wrap") != std::string::npos)
        {
            size_t start = line.find("\"");
            size_t end = line.rfind("\"");
            if (start != std::string::npos && end != std::string::npos && end > start)
            {
                f_wrap_disturb = line.substr(start + 1, end - start - 1);
            }
        }
        else if (line.find("IMapWrapper vectorField_wrap") != std::string::npos)
        {
            size_t start = line.find("\"");
            size_t end = line.rfind("\"");
            if (start != std::string::npos && end != std::string::npos && end > start)
            {
                vectorField_wrap = line.substr(start + 1, end - start - 1);
            }
        }
    }

    infile.close();
    return {f_wrap, f_wrap_disturb, vectorField_wrap};
}