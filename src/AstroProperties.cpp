#include "../include/AstroProperties.hpp"
#include <stdexcept>
#include <unordered_map>

double getDensityKgPerM3(const std::string& astro_name) {
    static const std::unordered_map<std::string, double> density_map = {
        {"Chury",   533.0},   // kg/m^3
        {"Phobos",  1860.0},
        {"Bennu",   1190.0},
        {"Ryugu",   1190.0},
        {"Lune",    3340.0},
        {"Mars",    3930.0},
        {"Mercure", 5427.0}
    };

    auto it = density_map.find(astro_name);
    if (it == density_map.end()) {
        throw std::runtime_error("Unknown astro name: " + astro_name);
    }

    return it->second;
}