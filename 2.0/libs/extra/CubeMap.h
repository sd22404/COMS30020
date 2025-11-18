#pragma once
#include "TextureMap.h"

struct CubeMap {
    TextureMap front;
    TextureMap back;
    TextureMap right;
    TextureMap left;
    TextureMap top;
    TextureMap bottom;
};