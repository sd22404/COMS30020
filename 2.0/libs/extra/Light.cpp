#include "Light.h"

void Light::sample() {
    samples.clear();
    for (int i = 0; i < uSteps; ++i) {
        for (int j = 0; j < vSteps; ++j) {
            float u = static_cast<float>(i + 0.5f) / static_cast<float>(uSteps);
            float v = static_cast<float>(j + 0.5f) / static_cast<float>(vSteps);
            // jitter inside cell for softer shadows
            u += (rand() / static_cast<float>(RAND_MAX) - 0.5f) / static_cast<float>(uSteps);
            v += (rand() / static_cast<float>(RAND_MAX) - 0.5f) / static_cast<float>(vSteps);
            u = glm::clamp(u, 0.0f, 1.0f);
            v = glm::clamp(v, 0.0f, 1.0f);
            samples.push_back(position + u * uVec + v * vVec);
        }
    }
}