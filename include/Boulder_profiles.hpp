#pragma once
#include <limits>
#include <string>
#include "Boulder_detector.hpp"

namespace boulder_profiles {

    inline boulder::PipelineParams MakeBaseParams(const std::string& mesh_path,
                                                  const std::string& output_dir) {
        boulder::PipelineParams p;
    
        p.mesh_path = mesh_path;
        p.output_dir = output_dir;
    
        p.visualize = true;
        p.verbose = true;
        p.use_absolute_height = false;
    
        p.height_threshold = std::numeric_limits<double>::quiet_NaN();
        p.roughness_max = std::numeric_limits<double>::quiet_NaN();
    
        return p;
    }
    
    inline boulder::PipelineParams Chury(const std::string& mesh_path,
                                         const std::string& output_dir = "out_chury") {
        auto p = MakeBaseParams(mesh_path, output_dir);
    
        p.smooth_iterations = 1;
        p.radius_multiplier = 12;
        p.min_neighbors = 26;
        p.height_sigma_factor = 1.80;
        p.outer_ring_depth = 2;
        p.max_elongation = 3.2;
    
        return p;
    }
    
    inline boulder::PipelineParams Bennu(const std::string& mesh_path,
                                         const std::string& output_dir = "out_bennu") {
        auto p = MakeBaseParams(mesh_path, output_dir);
    
        p.smooth_iterations = 0;
        p.radius_multiplier = 5.4;
        p.min_neighbors = 20;
        p.height_sigma_factor = 1.30;
        p.outer_ring_depth = 1;
        p.max_elongation = 3.3;
    
        return p;
    }
    
    inline boulder::PipelineParams Ryugu(const std::string& mesh_path,
                                         const std::string& output_dir = "out_ryugu") {
        auto p = MakeBaseParams(mesh_path, output_dir);
    
        p.smooth_iterations = 0;
        p.radius_multiplier = 5.8;
        p.min_neighbors = 20;
        p.height_sigma_factor = 1.35;
        p.outer_ring_depth = 1;
        p.max_elongation = 3.2;
    
        return p;
    }
    
    inline boulder::PipelineParams Phobos(const std::string& mesh_path,
                                          const std::string& output_dir = "out_phobos") {
        auto p = MakeBaseParams(mesh_path, output_dir);
    
        p.smooth_iterations = 1;
        p.radius_multiplier = 6.8;
        p.min_neighbors = 24;
        p.height_sigma_factor = 1.60;
        p.outer_ring_depth = 2;
        p.max_elongation = 2.8;
    
        return p;
    }
    
    inline boulder::PipelineParams GetParamsForBody(const std::string& body_name,
                                                    const std::string& mesh_path,
                                                    const std::string& output_dir = "") {
        if (body_name == "chury") {
            return Chury(mesh_path, output_dir.empty() ? "out_chury" : output_dir);
        }
        if (body_name == "bennu") {
            return Bennu(mesh_path, output_dir.empty() ? "out_bennu" : output_dir);
        }
        if (body_name == "ryugu") {
            return Ryugu(mesh_path, output_dir.empty() ? "out_ryugu" : output_dir);
        }
        if (body_name == "phobos") {
            return Phobos(mesh_path, output_dir.empty() ? "out_phobos" : output_dir);
        }
    
        throw std::invalid_argument("Unknown body_name: " + body_name);
    }
}// namespace boulder_profiles