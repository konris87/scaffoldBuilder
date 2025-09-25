#include "OpenGLsetup/Misc.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

void create_frame_buffer(FrameBuffer& framebuffer){
	
	// delete old buffers if they exist
	if (framebuffer.id) {
		glDeleteFramebuffers(1, &framebuffer.id);
		glDeleteTextures(1, &framebuffer.textureId);
		glDeleteRenderbuffers(1, &framebuffer.rbo);
	}

	glGenFramebuffers(1, &framebuffer.id);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.id); 

	// generate texture
	glGenTextures(1, &framebuffer.textureId);
	glBindTexture(GL_TEXTURE_2D, framebuffer.textureId);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, framebuffer.width, framebuffer.height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	// attach it to currently bound framebuffer object
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffer.textureId, 0);

	glGenRenderbuffers(1, &framebuffer.rbo);
	glBindRenderbuffer(GL_RENDERBUFFER, framebuffer.rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, framebuffer.width, framebuffer.height);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
	
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, framebuffer.rbo);

	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
	
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
};