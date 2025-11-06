#pragma once

enum RenderMode {
    WIREFRAME,
    RASTERISED,
    RAYTRACED
};

enum LightingMode {
    FLAT,
    GOURAUD,
    PHONG
};

enum Direction {
    UP,
    DOWN,
    LEFT,
    RIGHT,
    FORWARD,
    BACKWARD
};

static uint32_t packColour(glm::vec3 colour) {
    return (255 << 24) + (
           (int(colour.r * 255) << 16) +
           (int(colour.g * 255) << 8) +
           (int(colour.b * 255)));
}

static glm::vec3 unpackColour(uint32_t colour) {
    return {
        ((colour >> 16) & 0xFF) / 255.0f,
        ((colour >> 8) & 0xFF) / 255.0f,
        (colour & 0xFF) / 255.0f,
    };
}