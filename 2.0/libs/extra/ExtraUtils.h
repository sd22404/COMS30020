#pragma once

enum RenderMode {
    WIREFRAME,
    RASTERISED,
    RAYTRACED
};

enum ShadingMode {
    FLAT,
    GOURAUD,
    PHONG,
    GLASS
};

enum Direction {
    UP,
    DOWN,
    LEFT,
    RIGHT,
    FORWARD,
    BACKWARD
};

static uint32_t packColour(const glm::vec3 colour) {
    return (255 << 24) + (
           (lround(colour.r * 255) << 16) +
           (lround(colour.g * 255) << 8) +
           (lround(colour.b * 255)));
}

static glm::vec3 unpackColour(const uint32_t colour) {
    return {
        static_cast<float>((colour >> 16) & 0xFF) / 255.0f,
        static_cast<float>((colour >> 8) & 0xFF) / 255.0f,
        static_cast<float>(colour & 0xFF) / 255.0f,
    };
}