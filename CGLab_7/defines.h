#pragma once

#include "RasterSurface.h"
#include <random>
#include "StoneHenge.h"
#include "StoneHenge_Texture.h"
#define farPlane 100.0f
#define nearPlane 0.1f

#define Raster_W 900
#define Raster_H 700
#define R_PIXELS ((Raster_W)*(Raster_H))

unsigned int raster[R_PIXELS] = { 0 };
float ZBuffer[R_PIXELS] = { 100 };
unsigned colorR = 0xFF0000;
unsigned colorG = 0x00FF00;
unsigned colorB = 0x0000FF;
unsigned colorY = 0xFFFF00;

float fTheta = -60;

struct VERTEXES
{
	float x, y, z;
};
struct VERTEX 
{
	float x, y, z, w, u, v;
	float nx, ny, nz;
	unsigned int color;
};
struct mat4x4
{
	float m[4][4] = { 0 };
};
struct Triangle
{
	VERTEX a, b, c;
};
VERTEX stone[1457];

void GettingtheVertices() {
	for (size_t i = 0; i < 1457; i++)
	{
		stone[i].x = StoneHenge_data[i].pos[0] * 0.1f;
		stone[i].y = StoneHenge_data[i].pos[1] * 0.1f;
		stone[i].z = StoneHenge_data[i].pos[2] * 0.1f;
		stone[i].w = 1.0f;
		stone[i].u = StoneHenge_data[i].uvw[0];
		stone[i].v = StoneHenge_data[i].uvw[1];
		stone[i].nx = StoneHenge_data[i].nrm[0];
		stone[i].ny = StoneHenge_data[i].nrm[1];
		stone[i].nz = StoneHenge_data[i].nrm[2];
	}
}

VERTEX sunlightDirection = { -0.577, -0.577, 0.577 };
unsigned int sunColor = 0xFFC0C0F0;

VERTEX pointLight = { -1, 0.5, 1 };
unsigned int pointLightColor = 0xFFFFFF00;
float LightRadius = 2.0f;
struct MATRIX {
	union {

		struct{
			float x1, y1, z1;
			float x2, y2, z2;
			float x3, y3, z3;
		};
		float arr[3][3];
	};
};
struct MATRIX_4X4 {
	union {

		struct {
			float x1, y1, z1, w1;
			float x2, y2, z2, w2;
			float x3, y3, z3, w3;
			float x4, y4, z4, w4;
		};
		float arr[4][4];
	};
};
struct CUBE
{
	union {

		struct {
			Triangle bt1;
			Triangle bt2;
			Triangle t1, t2;
			Triangle l1, l2;
			Triangle r1, r2;
			Triangle fw1, fw2;
			Triangle ba1, ba2;
		};
		Triangle arr[2][6];
	};
};
struct LINE 
{
	VERTEX a, b;
};
struct GRID 
{
	union {

		struct {
			LINE a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11; 
			LINE h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11; 
		};

		LINE arr[11][2];
	};
};
