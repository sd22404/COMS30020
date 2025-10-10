#include <DrawingWindow.h>
#include <Utils.h>

class Renderer {
    private:
        DrawingWindow &window;
    public:
        Renderer(DrawingWindow &window) : window(window) {}
        void draw();
};
