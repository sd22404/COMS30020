#include <DrawingWindow.h>
#include "Renderer.h"
#include "Camera.h"
#include <vector>

#define WIDTH 640
#define HEIGHT 480

#define SCENE std::string("./assets/scenes/scene-2.obj")

void handleEvent(SDL_Event event, DrawingWindow &window, Camera &cam, Renderer &r) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_w) cam.move(Camera::FORWARD);
		if (event.key.keysym.sym == SDLK_s) cam.move(Camera::BACKWARD);
		if (event.key.keysym.sym == SDLK_a) cam.move(Camera::LEFT);
		if (event.key.keysym.sym == SDLK_d) cam.move(Camera::RIGHT);
		if (event.key.keysym.sym == SDLK_q) cam.move(Camera::DOWN);
		if (event.key.keysym.sym == SDLK_e) cam.move(Camera::UP);
		if (event.key.keysym.sym == SDLK_1) r.setMode(Renderer::WIREFRAME);
		if (event.key.keysym.sym == SDLK_2) r.setMode(Renderer::RASTERISED);
		if (event.key.keysym.sym == SDLK_3) r.setMode(Renderer::RAYTRACED);
		if (event.key.keysym.sym == SDLK_SPACE) cam.reset();
		if (event.key.keysym.sym == SDLK_LCTRL) cam.toggleOrbit();
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

[[noreturn]] int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	Camera cam = Camera(WIDTH, HEIGHT);
	Scene scene = Scene(SCENE, cam);
	Renderer r = Renderer(window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, cam, r);
		r.draw(scene);
		cam.orbit();
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
