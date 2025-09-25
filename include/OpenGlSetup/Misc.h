#ifndef MISC_H
#define MISC_H

struct FrameBuffer {
	unsigned int id = 0;
	unsigned int textureId = 0;
	unsigned int rbo = 0;
	int width{ 800 };
	int height{ 600 };
	int posx{ 0 };
	int posy{ 0 };
};

void create_frame_buffer(FrameBuffer& framebuffer);

#endif