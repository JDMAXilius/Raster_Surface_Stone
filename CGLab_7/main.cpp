#include "rasterizer.h"
#include "XTime.h"
#include <vector>
#include <assert.h>
#include <iostream>
#include <Windows.h>
#include <algorithm>

int main() {

	RS_Initialize(Raster_W, Raster_H);

	VERTEX stars[3000];
	for (size_t i = 0; i < 3000; i++)
	{
		stars[i] = { Random(-1,1)*50, Random(-1,1)*50, Random(-1,1) * 50, 1 };
	}
	MATRIX_4X4 cam = matrixmultiply(matrixtranslation(0, 0, -6.0f), rotateX(-18.0f));
	MATRIX_4X4 temp;

	unsigned int red = 0xFFFF0000;

	XTime timer;

	do
	{
		clear(Raster_W * Raster_H);
		//std::fill_n(ZBuffer, RPIXELS, 1);
		timer.Signal();
		if (GetAsyncKeyState('W') & 0x1) {
			temp = matrixtranslation(0, 0, 0.5f * timer.Delta());
			cam = matrixmultiply(temp, cam);
		}
		
		if (GetAsyncKeyState('S') & 0x1) {
			temp = matrixtranslation(0, 0, -0.5f * timer.Delta());
			cam = matrixmultiply(temp, cam);
		}
		
		if (GetAsyncKeyState('D') & 0x1) {
			temp = matrixtranslation(0.5f * timer.Delta(), 0, 0);
			cam = matrixmultiply(temp, cam);
		}
	
		if (GetAsyncKeyState('A') & 0x1) {
			temp = matrixtranslation(-0.5f * timer.Delta(), 0, 0);
			cam = matrixmultiply(temp, cam);
		}
		if (GetAsyncKeyState(VK_SPACE) & 0x1) {
			temp = matrixtranslation(0, 0.5f * timer.Delta(), 0);
			cam = matrixmultiply(temp, cam);
		}
		if (GetAsyncKeyState(VK_LSHIFT) & 0x1) {
			temp = matrixtranslation(0, -0.5f * timer.Delta(), 0);
			cam = matrixmultiply(temp, cam);
		}
		LightRadius = std::abs(cos(timer.TotalTime())) * 2 + 0.25f;
		VertexShader = VS_World;
		PixelShader = PS_White;

		SV_WorldMatrix = identity();

		//viewmatrix = matrixmultiply(matrixtranslation(0, 0, -1.0f), rotateX(-18.0f));
		viewmatrix = inverse(cam);

		projectmatrix = project(Raster_W/Raster_H, nearPlane, farPlane, 90.0f);
		for (size_t i = 0; i < 3000; i++)
		{
			DrawStar(stars[i]);
		}
		GettingtheVertices();
		DrawStoneHenge();
		

	} while (RS_Update(raster, R_PIXELS));

	RS_Shutdown();

	return 0;
}
