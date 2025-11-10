#include <DrawingWindow.h>
#include "Renderer.h"
#include <vector>

#define WIDTH 1280
#define HEIGHT 960

#define SCENE std::string("./assets/scenes/scene-2.obj")

void handleEvent(const SDL_Event &event, const DrawingWindow &window, Camera &cam, Renderer &r, const Scene &scene) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_w) cam.move(FORWARD);
		if (event.key.keysym.sym == SDLK_s) cam.move(BACKWARD);
		if (event.key.keysym.sym == SDLK_a) cam.move(LEFT);
		if (event.key.keysym.sym == SDLK_d) cam.move(RIGHT);
		if (event.key.keysym.sym == SDLK_q) cam.move(DOWN);
		if (event.key.keysym.sym == SDLK_e) cam.move(UP);
		if (event.key.keysym.sym == SDLK_1) r.setRenderMode(WIREFRAME);
		if (event.key.keysym.sym == SDLK_2) r.setRenderMode(RASTERISED);
		if (event.key.keysym.sym == SDLK_3) r.setRenderMode(RAYTRACED);
		if (event.key.keysym.sym == SDLK_SPACE) cam.reset();
		if (event.key.keysym.sym == SDLK_LCTRL) cam.toggleOrbit();
		if (event.key.keysym.sym == SDLK_LALT) r.toggleLight();
		if (event.key.keysym.sym == SDLK_i) scene.moveLight(FORWARD);
		if (event.key.keysym.sym == SDLK_k) scene.moveLight(BACKWARD);
		if (event.key.keysym.sym == SDLK_j) scene.moveLight(LEFT);
		if (event.key.keysym.sym == SDLK_l) scene.moveLight(RIGHT);
		if (event.key.keysym.sym == SDLK_o) scene.moveLight(UP);
		if (event.key.keysym.sym == SDLK_u) scene.moveLight(DOWN);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

[[noreturn]] int main(int argc, char *argv[]) {
	auto window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	auto cam = Camera(WIDTH, HEIGHT, 3.0f, glm::vec3(0, 0, 4));
	
	auto ceiling = Light(glm::vec3(0, 0.94f, 0), 15.0f);
	auto blue = Light(glm::vec3(0.7f, -0.3f, 0), glm::vec3(0, 0, 1), 5.0f);
	std::vector<Light> lights = {ceiling};

	const std::vector<Obj> objs = {
		// Obj{"./assets/sphere/glass-sphere.obj", PHONG, {-0.6, 0.3, -0.6}, 0.25f},
		// Obj{"./assets/cornell-box/cornell-box.obj", FLAT, {0, 0, 0}, 0.35f},
		Obj{SCENE, FLAT},
		Obj{"./assets/sphere/gold-sphere.obj", PHONG},
		Obj{"./assets/bunny/glass-bunny.obj", PHONG},
		Obj{"./assets/hackspace-logo/logo.obj", FLAT, {0.1f, 0.1f, -0.93f}, 0.0015f},
	};
	const auto scene = Scene(objs, lights);
	
	auto r = Renderer(window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, cam, r, scene);
		r.draw(scene, cam);
		cam.orbit();
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
