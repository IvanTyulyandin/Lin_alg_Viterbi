#pragma once

#include "HMM.h"

#include <string>

// Return HMM with probabilities stored as negative logarithm
HMM read_HMM(const std::string& HMM_file_name);
