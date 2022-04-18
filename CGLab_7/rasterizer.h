#pragma once

#include "shader.h"
#include <algorithm>

void PerspectiveDivide(VERTEX& multiplyMe) {
	multiplyMe.x /= multiplyMe.w;
	multiplyMe.y /= multiplyMe.w;
	multiplyMe.z /= multiplyMe.w;
}
void clear(unsigned int color) {
	for (size_t i = 0; i < R_PIXELS; i++)
	{
		raster[i] = color;
		ZBuffer[i] = 1;
	}
}
void PlotPixel(int x, int y, unsigned color)
{
	(x >= 0 && x <= Raster_W);
	(y >= 0 && y <= Raster_H);
	raster[Raster_W * y + x] = color;
}
void PlotPixelZ(int x, int y, float z, unsigned color)
{
	/*if (x >= Raster_W || y >= Raster_H || x < 0 || y < 0)
		return;*/


	(x >= 0 && x <= Raster_W);
	(y >= 0 && y <= Raster_H);
	//assert(z >= 0 && y <= Raster_H);
	float Zcomparer = ZBuffer[Raster_W * y + x];

	if (z < Zcomparer)
	{
		raster[Raster_W * y + x] = color;
		ZBuffer[Raster_W * y + x] = z;
	}
	else
	{
		//raster[Raster_W * y + x] = color;
	}
}
void dpixel(int x, int y, float z, int screenWidth, unsigned int color) {
	if (x >= Raster_W || x < 0 || y >= Raster_H || y < 0)
		return;
	float comparer = ZBuffer[find(x, y, screenWidth)];
	if (z < comparer)
	{
		raster[find(x, y, screenWidth)] = color;
		ZBuffer[find(x, y, screenWidth)] = z;
	}

}\
float LinearizeDepth(float depth, float nearz, float farz)
{
	float z = depth * 2.0 - 1.0; // back to NDC 
	return (2.0 * nearz * farz) / (farz + nearz - z * (farz - nearz));
}
unsigned int BGRAtoARGB(unsigned int color) {

	//BGRA -> ARGB //separate chanells, moving around and combining all together.
	color = ((0xFF000000 & color) >> 24) |
		((0x00FF0000 & color) >> 8) |
		((0x0000FF00 & color) << 8) |
		((0x000000FF & color) << 24);

	return color;
}
void DrawStar(VERTEX vec) {

	VERTEX copy = vec;
	if (VertexShader)
		VertexShader(copy);

	if (vec.z < nearPlane)
		return;

	PerspectiveDivide(copy);

	VERTEX fVec = btopixel(copy);

	unsigned int color = 0xFFFFFFFF;
	if (PixelShader)
		PixelShader(color);

	dpixel(fVec.x, fVec.y, fVec.z, Raster_W, color);

}

//unsigned int BGRAtoARGB(unsigned int color) {
//
//	unsigned int Tgreen;
//	unsigned int Tred;
//	unsigned int Tblue;
//	unsigned int Talpha;
//
//	Tblue = (color & 0xFF000000) >> 24;
//	Tgreen = (color & 0x00FF0000) >> 8;
//	Tred = (color & 0x0000FF00) << 8;
//	Talpha = (color & 0x000000FF) << 24;
//
//	return Talpha | Tred | Tgreen | Tblue;
//}

unsigned int berpC(unsigned int A, unsigned int B, unsigned int C, VERTEX r)
{
	int Aa = (A & 0xFF000000) >> 24;
	int Ar = (A & 0x00FF0000) >> 16;
	int Ag = (A & 0x0000FF00) >> 8;
	int Ab = (A & 0x000000FF);
	int Ba = (B & 0xFF000000) >> 24;
	int Br = (B & 0x00FF0000) >> 16;
	int Bg = (B & 0x0000FF00) >> 8;
	int Bb = (B & 0x000000FF);
	int Ca = (C & 0xFF000000) >> 24;
	int Cr = (C & 0x00FF0000) >> 16;
	int Cg = (C & 0x0000FF00) >> 8;
	int Cb = (C & 0x000000FF);
	unsigned int fa = ((Aa)*r.x + (Ba)*r.y + (Ca)*r.z);
	unsigned int fr = ((Ar)*r.x + (Br)*r.y + (Cr)*r.z);
	unsigned int fg = ((Ag)*r.x + (Bg)*r.y + (Cg)*r.z);
	unsigned int fb = ((Ab)*r.x + (Bb)*r.y + (Cb)*r.z);
	fa = (fa << 24);
	fr = (fr << 16);
	fg = (fg << 8);
	return (fa | fr | fg | fb);
}



void FillTriangle(const Triangle& tri) {
	 int minX = std::min(tri.a.x, std::min(tri.b.x, tri.c.x));
	 int maxX = std::max(tri.a.x, std::max(tri.b.x, tri.c.x));

	 int minY = std::min(tri.a.y, std::min(tri.b.y, tri.c.y));
	 int maxY = std::max(tri.a.y, std::max(tri.b.y, tri.c.y));

	for (int y = minY; y <= maxY; ++y)
	{
		for (int x = minX; x <= maxX; ++x)
		{
			VERTEX p = { static_cast<float>(x), static_cast<float>(y) };
			VERTEX temp = FindBary(tri.a, tri.b, tri.c, p);
			if (temp.x >= 0 && temp.x <= 1 && temp.y >= 0 && temp.y <= 1 && temp.z >= 0 && temp.z <= 1) {


				float z = tri.a.z * temp.x + tri.b.z * temp.y + tri.c.z * temp.z;
				float u = Bary((tri.a.u / tri.a.w), (tri.b.u / tri.b.w), (tri.c.u / tri.c.w), temp);
				float v = Bary((tri.a.v / tri.a.w), (tri.b.v / tri.b.w), (tri.c.v / tri.c.w), temp);
				float bar = Bary((1 / tri.a.w), (1 / tri.b.w), (1 / tri.c.w), temp); 
				u = u / bar;
				v = v / bar;
				int mipLevel = ((z - nearPlane) / (farPlane - nearPlane)) * StoneHenge_numlevels;
				int newWidth = StoneHenge_width >> mipLevel;
				int newHeight = StoneHenge_height >> mipLevel;

				unsigned int offset = StoneHenge_leveloffsets[mipLevel];

				float pixelX = (u * newWidth);
				float pixelY = (v * newHeight);

				int index = static_cast<int>(pixelX) + static_cast<int>(pixelY) * newWidth;

				unsigned int topLeft = StoneHenge_pixels[index + offset];
				unsigned int topRight = StoneHenge_pixels[index + 1 + offset];
				unsigned int botLeft = StoneHenge_pixels[index + newWidth + offset];
				unsigned int botRight = StoneHenge_pixels[index + newWidth + 1 + offset];

				float ratioX = pixelX - floorf(pixelX);
				float ratioY = pixelY - floorf(pixelY);


				unsigned int avTop = ColorLerp(topLeft, topRight, ratioX);
				unsigned int avBot = ColorLerp(botLeft, botRight, ratioX);
				unsigned int finalPixel = ColorLerp(avTop, avBot, ratioY);

				unsigned int baryColor = berpC(tri.a.color, tri.b.color, tri.c.color, temp);
				unsigned int actuallyFinal = multiplycolors(baryColor, BGRAtoARGB(finalPixel));
				dpixel(x, y, z, Raster_W, actuallyFinal);
			}
		}
	}
}
int ClipLineNearPlane(VERTEX& a, VERTEX& b) {

	if (a.z > nearPlane && b.z > nearPlane)
		return 1;  
	if (a.z < nearPlane && b.z < nearPlane)
		return -1; 


	if (a.z < nearPlane) {
		float ratio = (nearPlane - a.z) / (b.z - a.z);
		a.x = LerpLine(a.x, b.x, ratio);
		a.y = LerpLine(a.y, b.y, ratio);
		a.z = LerpLine(a.z, b.z, ratio);
		a.w = LerpLine(a.w, b.w, ratio);
		a.u = LerpLine(a.u, b.u, ratio);
		a.v = LerpLine(a.v, b.v, ratio);
	}
	if (b.z < nearPlane) {
		float ratio = (nearPlane - b.z) / (a.z - b.z);
		b.x = LerpLine(b.x, a.x, ratio);
		b.y = LerpLine(b.y, a.y, ratio);
		b.z = LerpLine(b.z, a.z, ratio);
		b.w = LerpLine(b.w, a.w, ratio);
		b.u = LerpLine(b.u, a.u, ratio);
		b.v = LerpLine(b.v, a.v, ratio);
	}

	return 0;
}
int ClipTriangleNearPlane(VERTEX& a, VERTEX& b, VERTEX& c, VERTEX& d) {

	if (a.z > nearPlane && b.z > nearPlane && c.z > nearPlane) {
		return -1;
	}
	if (a.z < nearPlane && b.z < nearPlane && c.z < nearPlane) {
		return 0;
	}

	VERTEX lineAB[2] = { a, b };
	VERTEX lineBC[2] = { b, c };
	VERTEX lineCA[2] = { c, a };

	int results[3] = {
		ClipLineNearPlane(lineAB[0], lineAB[1]),
		ClipLineNearPlane(lineBC[0], lineBC[1]),
		ClipLineNearPlane(lineCA[0], lineCA[1])
	};

	if (results[0] == -1)
	{
		a = lineCA[1];
		b = lineBC[0];
		c = c;
		return 1;
	}
	if (results[0] == 1)
	{
		a = a;
		b = b;
		c = lineBC[1];
		d = lineCA[0];
		return 2;
	}
	if (results[1] == -1)
	{
		a = a;
		b = lineAB[1];
		c = lineCA[0];
		return 1;
	}

	if (results[1] == 1)
	{
		a = lineAB[0];
		b = b;
		c = c;
		d = lineCA[1];
		return 2;
	}

	if (results[2] == -1)
	{
		a = lineAB[0];
		b = b;
		c = lineBC[1];
		return 1;
	}

	if (results[2] == 1)
	{
		a = a;
		b = lineAB[1];
		d = c;
		c = lineBC[0];
		return 2;
	}
}
void DrawTriangle(const Triangle& triangle) {

	Triangle copyTri = triangle;

	if (VertexShader) {
		VertexShader(copyTri.a);
		VertexShader(copyTri.b);
		VertexShader(copyTri.c);
	}
	VERTEX d;
	int answer = ClipTriangleNearPlane(copyTri.a, copyTri.b, copyTri.c, d);

	if (answer == 0)
		return;


	PerspectiveDivide(copyTri.a);
	PerspectiveDivide(copyTri.b);
	PerspectiveDivide(copyTri.c);
	PerspectiveDivide(d);
	Triangle tri = triangletopixel(copyTri);
	Triangle secondTri = { copyTri.a, copyTri.c, d };
	Triangle sTri = triangletopixel(secondTri);

	FillTriangle(tri);

	if (answer >= 2) {
		FillTriangle(sTri);
	}
	
}

//void ParametricLine(VERTEXES4D inia, VERTEXES4D inib, unsigned int color)
//{
//	VERTEXES4D copya = inia;
//	VERTEXES4D copyb = inib;
//	if (VertexShader) {
//		VertexShader(copya);
//		VertexShader(copyb);
//	}
//
//	if (ClipLineNearPlane(copya, copyb) == -1)
//	{
//		return;
//	}
//
//	PerspectiveDivide(copya);
//	PerspectiveDivide(copyb);
//
//	VERTEXES4D a = btopixel(copya);
//	VERTEXES4D b = btopixel(copyb);
//	int deltaX = std::abs(b.x - a.x);
//	int deltaY = std::abs(b.y - a.y);
//	int totalPixels = (deltaX > deltaY) ? deltaX : deltaY;
//	for (int i = 0; i < totalPixels; i++)
//	{
//		float R = static_cast<float>(i) / totalPixels;
//		int pX = static_cast<int>(a.x + R * (b.x - a.x) + 0.5f);
//		int pY = static_cast<int>(a.y + R * (b.y - a.y) + 0.5f);
//		float pZ = LerpLine(a.z, b.z, R);
//
//		dpixel(pX, pY, pZ, RWIDTH, color);
//	}
//}


void DrawStoneHenge() {

	for (int i = 0; i < 2532; i += 3)
	{
		Triangle temp = { stone[StoneHenge_indicies[i]], stone[StoneHenge_indicies[i + 1]] , stone[StoneHenge_indicies[i + 2]] };
		DrawTriangle(temp);
	}
}