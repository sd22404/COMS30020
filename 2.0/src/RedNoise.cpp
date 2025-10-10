#include <DrawingWindow.h>
#include "Renderer.h"
#include "Camera.h"

#define WIDTH 640
#define HEIGHT 480

#define SCENE std::string("./assets/scenes/scene-3.obj")

void handleEvent(SDL_Event event, DrawingWindow &window, Camera &cam) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_w) cam.move(FORWARD);
		if (event.key.keysym.sym == SDLK_s) cam.move(BACKWARD);
		if (event.key.keysym.sym == SDLK_a) cam.move(LEFT);
		if (event.key.keysym.sym == SDLK_d) cam.move(RIGHT);
		if (event.key.keysym.sym == SDLK_q) cam.move(DOWN);
		if (event.key.keysym.sym == SDLK_e) cam.move(UP);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	Camera cam = Camera(WIDTH, HEIGHT);
	Scene scene = Scene(SCENE);
	Renderer r = Renderer(window, scene, cam);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, cam);
		r.draw(scene);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
