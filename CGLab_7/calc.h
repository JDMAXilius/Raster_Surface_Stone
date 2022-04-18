#pragma once

#include "defines.h"
//lerp(a, b, c, x) = ((a * x) + (b * x) + (c * x)); or P=(uA)+(vB)+(wC)
float barycentriclerpEq(int a, int b, int c, VERTEX r)
{
	return((a * r.x) + (b * r.y) + (c * r.z));

}
float ImplicitLineEq2(VERTEX a, VERTEX b, VERTEX p)
{
	return (a.y - b.y) * p.x + (b.x - a.x) * p.y + a.x * b.y - a.y * b.x;
}


float LineLerp(int y1, int y2, float ratio)
{
	return static_cast<float>(y1) + (ratio * static_cast<float>(y2 - y1));
}

float Bary(float a, float b, float c, VERTEX ratios) {

	return a * ratios.x + b * ratios.y + c * ratios.z;
}
float Random(float min, float max) {
	float ratio = rand() / (float)RAND_MAX;
	return (max - min) * ratio + min;
}
float ile(VERTEXES a, VERTEXES b, VERTEXES p)
{
	return (a.y - b.y) * p.x + (b.x - a.x) * p.y + (a.x * b.y) - (a.y * b.x);
}
int find(int x, int y, int screenWidth) {
	return (y * screenWidth + x);
}
float lerpline(int y1, int y2, float ratio)
{
	return static_cast<float>(y1) + (ratio * static_cast<float>(y2 - y1));
}
float ImplicitLineEquation(VERTEX a, VERTEX b, VERTEX p) {
	float answer = (a.y - b.y) * p.x + (b.x - a.x) * p.y + (a.x * b.y) - (a.y * b.x);
	return answer;
}
float LerpLine(float y1, float y2, float ratio) {
	return (y1)+(ratio * (y2 - y1));
}
float dot(VERTEX a, VERTEX b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
VERTEX cross(VERTEX a, VERTEX b)
{
	VERTEX vec = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
	return vec;
}
float vlength(VERTEX vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}
VERTEX NormalizeVector(VERTEX vec) {
	float magnitude = vlength(vec);
	if (magnitude != 0)
	{
		vec.x /= magnitude;
		vec.y /= magnitude;
		vec.z /= magnitude;
	}
	else
		vec.x = vec.y = vec.z = 0.0f;

	return vec;
}
float saturate(float v) {
	if (v > 1)
		return 1;
	else if (v < 0)
		return 0;
	else
		return v;
}
unsigned int addcolors(unsigned int A, unsigned int B )
{
	int Aa = ((A & 0xFF000000) >> 24) / 255.0f;
	int Ar = ((A & 0x00FF0000) >> 16) / 255.0f;
	int Ag = ((A & 0x0000FF00) >> 8 )/ 255.0f;
	int Ab = (A & 0x000000FF) / 255.0f;
	int Ba = ((B & 0xFF000000) >> 24) / 255.0f;
	int Br = ((B & 0x00FF0000) >> 16) / 255.0f;
	int Bg = ((B & 0x0000FF00) >> 8) / 255.0f;
	int Bb = (B & 0x000000FF) / 255.0f;
	unsigned int fa = static_cast<unsigned>(saturate(Aa + Ba) * 255) << 24;
	unsigned int fr = static_cast<unsigned>(saturate(Ar + Br) * 255) << 16;
	unsigned int fg = static_cast<unsigned>(saturate(Ag + Bg) * 255) << 8;
	unsigned int fb = static_cast<unsigned>(saturate(Ab + Bb) * 255);
	return (fa | fr | fg | fb);
}
unsigned int multiplycolors(unsigned int A, unsigned int B)
{
	float Aa = ((A & 0xFF000000) >> 24) / 255.0f;
	float Ar = ((A & 0x00FF0000) >> 16) / 255.0f;
	float Ag = ((A & 0x0000FF00) >> 8) / 255.0f;
	float Ab = ((A & 0x000000FF)) / 255.0f;
	float Ba = ((B & 0xFF000000) >> 24) / 255.0f;
	float Br = ((B & 0x00FF0000) >> 16) / 255.0f;
	float Bg = ((B & 0x0000FF00) >> 8) / 255.0f;
	float Bb = (B & 0x000000FF) / 255.0f;

	float finala = Aa + Ba;
	float finalr = Ar + Br;
	float finalg = Ag + Bg;
	float finalb = Ab + Bb;

	if (finala > 1)
		finala = 1;
	if (finala < 0)
		finala = 0;
	if (finalb > 1)
		finalb = 1;
	if (finalb < 0)
		finalb = 0;
	if (finalr > 1)
		finalr = 1;
	if (finalr < 0)
		finalr = 0;
	if (finalg > 1)
		finalg = 1;
	if (finalg < 0)
		finalg = 0;

	unsigned int fa = static_cast<unsigned>(finala * 255) << 24;
	unsigned int fr = static_cast<unsigned>(finalr * 255) << 16;
	unsigned int fg = static_cast<unsigned>(finalg * 255) << 8;
	unsigned int fb = static_cast<unsigned>(finalb * 255);
	return (fa | fr | fg | fb);
}
unsigned int ColorLerp(unsigned int A, unsigned int B, float ratio) {

	unsigned int Aalpha = (A & 0xFF000000) >> 24;
	unsigned int Ared = (A & 0x00FF0000) >> 16;
	unsigned int Agreen = (A & 0x0000FF00) >> 8;
	unsigned int Ablue = (A & 0x000000FF);


	unsigned int Balpha = (B & 0xFF000000) >> 24;
	unsigned int Bred = (B & 0x00FF0000) >> 16;
	unsigned int Bgreen = (B & 0x0000FF00) >> 8;
	unsigned int Bblue = (B & 0x000000FF);

	unsigned int Falpha = (static_cast<int>(Balpha) - static_cast<int>(Aalpha)) * ratio + Aalpha;
	unsigned int Fred = (static_cast<int>(Bred) - static_cast<int>(Ared)) * ratio + Ared;
	unsigned int Fgreen = (static_cast<int>(Bgreen) - static_cast<int>(Agreen)) * ratio + Agreen;
	unsigned int Fblue = (static_cast<int>(Bblue) - static_cast<int>(Ablue)) * ratio + Ablue;

	return ((Falpha << 24) | (Fred << 16) | (Fgreen << 8) | (Fblue));
}
//Barycentric Interpolation
VERTEX BarycentricCoordinates(VERTEX a, VERTEX b, VERTEX c, VERTEX p)
{
	/*β = ImplicitLineEquation(B, line AC)
	γ = ImplicitLineEquation(C, line BA)
	α = ImplicitLineEquation(A, line CB)
	b = ImplicitLineEquation(P, line AC)
	y = ImplicitLineEquation(P, line BA)
	a = ImplicitLineEquation(P, line CB)
	Pβγα = (b / β, y / γ, a / α)*/

	VERTEX Pβγα;

	float β = ImplicitLineEq2(b, a, c);
	float γ = ImplicitLineEq2(c, b, a);
	float α = ImplicitLineEq2(a, c, b);
	float B = ImplicitLineEq2(p, a, c);
	float Y = ImplicitLineEq2(p, b, a);
	float A = ImplicitLineEq2(p, c, b);

	Pβγα = { (B / β), (Y / γ), (A / α) };

	return Pβγα;
}

VERTEX FindBary(VERTEX a, VERTEX b, VERTEX c, VERTEX point) {

	unsigned int returnColor = 0xFF000000;

	float beta = ImplicitLineEquation(a, c, b);
	float gamma = ImplicitLineEquation(b, a, c);
	float alpha = ImplicitLineEquation(c, b, a);
	float bLine = ImplicitLineEquation(a, c, point);
	float yLine = ImplicitLineEquation(b, a, point);
	float aLine = ImplicitLineEquation(c, b, point);

	VERTEX P = { (aLine / alpha), (bLine / beta), (yLine / gamma) };

	return P;
}
VERTEXES vertexmatrixmultiply(VERTEXES v, MATRIX m)
{
	VERTEXES t;
	t.x = (v.x * m.x1) + (v.y * m.x2) + (v.z * m.x3);
	t.y = (v.x * m.y1) + (v.y * m.y2) + (v.z * m.y3);
	t.z = (v.x * m.z1) + (v.y * m.z2) + (v.z * m.z3);
	return t;
}
VERTEX vertexmatrixmultiply(VERTEX v, MATRIX_4X4 m) 
{
	VERTEX t;
	t.x = (v.x * m.x1) + (v.y * m.x2) + (v.z * m.x3) + (v.w * m.x4);
	t.y = (v.x * m.y1) + (v.y * m.y2) + (v.z * m.y3) + (v.w * m.y4);
	t.z = (v.x * m.z1) + (v.y * m.z2) + (v.z * m.z3) + (v.w * m.z4);
	t.w = (v.x * m.w1) + (v.y * m.w2) + (v.z * m.w3) + (v.w * m.w4);
	t.u = v.u;
	t.v = v.v;
	t.nx = v.nx;
	t.ny = v.ny;
	t.nz = v.nz;
	t.color = v.color;
	return t;
}
MATRIX_4X4 identity() 
{
	MATRIX_4X4 t = { 0 };
	t.x1 = 1;
	t.y2 = 1;
	t.z3 = 1;
	t.w4 = 1;
	return t;
}
MATRIX matrixmultiply(MATRIX m1, MATRIX m2)
{
	MATRIX t = { 0 };
	t.x1 = (m1.x1 * m2.x1) + (m1.y1 * m2.x2) + (m1.z1 * m2.x3);
	t.y1 = (m1.x1 * m2.y1) + (m1.y1 * m2.y2) + (m1.z1 * m2.y3);
	t.z1 = (m1.x1 * m2.z1) + (m1.y1 * m2.z2) + (m1.z1 * m2.z3);
	t.x2 = (m1.x2 * m2.x1) + (m1.y2 * m2.x2) + (m1.z2 * m2.x3);
	t.y2 = (m1.x2 * m2.y1) + (m1.y2 * m2.y2) + (m1.z2 * m2.y3);
	t.z2 = (m1.x2 * m2.z1) + (m1.y2 * m2.z2) + (m1.z2 * m2.z3);
	t.x3 = (m1.x3 * m2.x1) + (m1.y3 * m2.x2) + (m1.z3 * m2.x3);
	t.y3 = (m1.x3 * m2.y1) + (m1.y3 * m2.y2) + (m1.z3 * m2.y3);
	t.z3 = (m1.x3 * m2.z1) + (m1.y3 * m2.z2) + (m1.z3 * m2.z3);
	return t;
}
MATRIX_4X4 matrixmultiply(MATRIX_4X4 m1, MATRIX_4X4 m2) 
{
	MATRIX_4X4 t = { 0 };
	t.x1 = (m1.x1 * m2.x1) + (m1.y1 * m2.x2) + (m1.z1 * m2.x3) + (m1.w1 * m2.x4);
	t.y1 = (m1.x1 * m2.y1) + (m1.y1 * m2.y2) + (m1.z1 * m2.y3) + (m1.w1 * m2.y4);
	t.z1 = (m1.x1 * m2.z1) + (m1.y1 * m2.z2) + (m1.z1 * m2.z3) + (m1.w1 * m2.z4);
	t.w1 = (m1.x1 * m2.w1) + (m1.y1 * m2.w2) + (m1.z1 * m2.w3) + (m1.w1 * m2.w4);
	t.x2 = (m1.x2 * m2.x1) + (m1.y2 * m2.x2) + (m1.z2 * m2.x3) + (m1.w2 * m2.x4);
	t.y2 = (m1.x2 * m2.y1) + (m1.y2 * m2.y2) + (m1.z2 * m2.y3) + (m1.w2 * m2.y4);
	t.z2 = (m1.x2 * m2.z1) + (m1.y2 * m2.z2) + (m1.z2 * m2.z3) + (m1.w2 * m2.z4);
	t.w2 = (m1.x2 * m2.w1) + (m1.y2 * m2.w2) + (m1.z2 * m2.w3) + (m1.w2 * m2.w4);
	t.x3 = (m1.x3 * m2.x1) + (m1.y3 * m2.x2) + (m1.z3 * m2.x3) + (m1.w3 * m2.x4);
	t.y3 = (m1.x3 * m2.y1) + (m1.y3 * m2.y2) + (m1.z3 * m2.y3) + (m1.w3 * m2.y4);
	t.z3 = (m1.x3 * m2.z1) + (m1.y3 * m2.z2) + (m1.z3 * m2.z3) + (m1.w3 * m2.z4);
	t.w3 = (m1.x3 * m2.w1) + (m1.y3 * m2.w2) + (m1.z3 * m2.w3) + (m1.w3 * m2.w4);
	t.x4 = (m1.x4 * m2.x1) + (m1.y4 * m2.x2) + (m1.z4 * m2.x3) + (m1.w4 * m2.x4);
	t.y4 = (m1.x4 * m2.y1) + (m1.y4 * m2.y2) + (m1.z4 * m2.y3) + (m1.w4 * m2.y4);
	t.z4 = (m1.x4 * m2.z1) + (m1.y4 * m2.z2) + (m1.z4 * m2.z3) + (m1.w4 * m2.z4);
	t.w4 = (m1.x4 * m2.w1) + (m1.y4 * m2.w2) + (m1.z4 * m2.w3) + (m1.w4 * m2.w4);
	return t;																	 
}
MATRIX_4X4 matrixtranslation(float x, float y, float z) 
{
	MATRIX_4X4 t = identity();
	t.x4 = x;
	t.y4 = y;
	t.z4 = z;
	return t;
}
MATRIX_4X4 rotateX(float d) 
{
	d = d * 3.14f / 180.0f;
	MATRIX_4X4 t = { 0 };
	t.y2 = cos(d);
	t.z2 = -sin(d);
	t.y3 = sin(d);
	t.z3 = cos(d);
	t.x1 = 1;
	t.w4 = 1;
	return t;
}
MATRIX_4X4 rotateY(float d) 
{
	d = d * 3.14f / 180.0f;
	MATRIX_4X4 t = { 0 };
	t.x1 = cos(d);
	t.z1 = sin(d);
	t.x3 = -sin(d);
	t.z3 = cos(d);
	t.y2 = 1;
	t.w4 = 1;
	return t;
}
MATRIX_4X4 rotateZ(float d)
{
	d = d * 3.14f / 180.0f;
	MATRIX_4X4 t = { 0 };
	t.x1 = cos(d);
	t.y1 = -sin(d);
	t.x2 = sin(d);
	t.y2 = cos(d);
	t.z3 = 1;
	return t;
}
MATRIX_4X4 inverse(MATRIX_4X4 m) {
	float t = m.x2;
	m.x2 = m.y1;
	m.y1 = t;
	t = m.x3;
	m.x3 = m.z1;
	m.z1 = t;
	t = m.y3;
	m.y3 = m.z2;
	m.z2 = t;
	VERTEXES tv = { m.x4, m.y4, m.z4 };
	MATRIX tm = { m.x1, m.y1, m.z1, m.x2, m.y2, m.z2, m.x3, m.y3, m.z3 };
	tv = vertexmatrixmultiply(tv, tm);
	m.x4 = tv.x * -1;
	m.y4 = tv.y * -1;
	m.z4 = tv.z * -1;
	return m;
}
MATRIX_4X4 project(float ar, float np, float fp, float fov) {
	float Ys = 1 / tan((fov / 2) * 3.141592f / 180);
	float Xs = Ys * ar;
	MATRIX_4X4 t =
	{
	   Xs, 0, 0, 0	, 0, Ys, 0, 0,	 0, 0, fp / (fp - np), 1,   0, 0, -(fp * np) / (fp - np), 0
	};
	return t;
}
VERTEX NDC2(VERTEX ndc)
{
	//NDC -> Screen
	VERTEX screen;

	screen.x = ((ndc.x + 1) * (Raster_W / 2));
	screen.y = ((1 - ndc.y) * (Raster_H / 2));
	//result
	return screen;
}
VERTEXES btopixel(VERTEXES v)
{
	v.x = floor((v.x + 1) * (Raster_W / 2));
	v.y = floor((1 - v.y) * (Raster_H / 2));
	return v;
}
VERTEX btopixel(VERTEX v) 
{
	v.x = floor((v.x + 1) * (Raster_W / 2));
	v.y = floor((1 - v.y) * (Raster_H / 2));
	return v;
}
Triangle triangletopixel(Triangle& tri) {

	tri.a = btopixel(tri.a);

	tri.b = btopixel(tri.b);
	tri.c = btopixel(tri.c);
	return tri;
}
//Triangle triangletopixel(Triangle& tri) {
//
//	tri.a = NDC2(tri.a);
//	tri.b = NDC2(tri.b);
//	tri.c = NDC2(tri.c);
//	return tri;
//}

#pragma region Matrix Calculations


VERTEX vertexMultMatrix4x4(VERTEX v, mat4x4 m)
{

	VERTEX vex;

	vex.x = (v.x * m.m[0][0]) + (v.y * m.m[1][0]) + (v.z * m.m[2][0]) + (v.w * m.m[3][0]);
	vex.y = (v.x * m.m[0][1]) + (v.y * m.m[1][1]) + (v.z * m.m[2][1]) + (v.w * m.m[3][1]);
	vex.z = (v.x * m.m[0][2]) + (v.y * m.m[1][2]) + (v.z * m.m[2][2]) + (v.w * m.m[3][2]);
	vex.w = (v.x * m.m[0][3]) + (v.y * m.m[1][3]) + (v.z * m.m[2][3]) + (v.w * m.m[3][3]);

	float w = v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + m.m[3][3];

	/*if (w != 0.0f)
	{
		vex.x /= w; vex.y /= w; vex.z /= w;
	}*/

	/*for (int i = 0; i < 3; i++)
	{
		t.x = (v.x * m.m[0][0]) + (v.y * m.m[1][0]) + (v.z * m.m[2][0]) +
			(v.w * m.m[3][0]); t.x = (v.x * m.m[0][0]) + (v.y * m.m[1][0]) +
			(v.z * m.m[2][0]) + (v.w * m.m[3][0]);
	}*/

	return vex;
}

mat4x4 matrixMultMatrix4x4b(mat4x4 m, mat4x4 l)
{
	mat4x4 results;
	float product = 0;
	int i;
	int j;

	for (i = 0; i <= 3; i++)
	{
		for (j = 0; j <= 3; j++)
		{
			results.m[i][j] += m.m[i][j] * l.m[j][i];
		}
	}
	return results;
}
mat4x4 matrixMultMatrix4x4(mat4x4 m1, mat4x4 m2)
{
	mat4x4 mat = { 0 };

	int i = 0;
	int j;
	while (i <= 3)
	{
		mat.m[0][i] = (m1.m[0][0] * m2.m[0][i]) + (m1.m[0][1] * m2.m[1][i]) + (m1.m[0][2] * m2.m[2][i]) + (m1.m[0][3] * m2.m[3][i]);

		mat.m[1][i] = (m1.m[1][0] * m2.m[0][i]) + (m1.m[1][1] * m2.m[1][i]) + (m1.m[1][2] * m2.m[2][i]) + (m1.m[1][3] * m2.m[3][i]);

		mat.m[2][i] = (m1.m[2][0] * m2.m[0][i]) + (m1.m[2][1] * m2.m[1][i]) + (m1.m[2][2] * m2.m[2][i]) + (m1.m[2][3] * m2.m[3][i]);

		mat.m[3][i] = (m1.m[3][0] * m2.m[0][i]) + (m1.m[3][1] * m2.m[1][i]) + (m1.m[3][2] * m2.m[2][i]) + (m1.m[3][3] * m2.m[3][i]);
		i++;
	};

	return mat;
}

mat4x4 matrixIdentity()
{
	mat4x4 m =
	{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	};

	/*Identity
		[
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1,
		]*/

		/*mat4x4 m = { 0 };

		m.m[0][0] = 1;
		m.m[1][1] = 1;
		m.m[2][2] = 1;
		m.m[3][3] = 1;*/

	return m;
}
mat4x4 matrixRotationZ4x4(float fTheta)
{
	mat4x4 matRotZ = { 0 };
	//time of rotation
	//fTheta = fTheta / 180.0f * M_PI;
	fTheta = fTheta / 180.0f * 3.14;

	// Rotation Z
	matRotZ.m[0][0] = cosf(fTheta * 0.5f);
	matRotZ.m[0][1] = sinf(fTheta * 0.5f);
	matRotZ.m[1][0] = -sinf(fTheta * 0.5f);
	matRotZ.m[1][1] = cosf(fTheta * 0.5f);
	matRotZ.m[2][2] = 1;
	matRotZ.m[3][3] = 1;

	return matRotZ;
}

mat4x4 matrixRotationX4x4b(float fTheta)
{
	//MATRIX_4X4 matRotX = { 0 };
	//time of rotation
	//fTheta = fTheta / 180.0f * M_PI;
	fTheta = fTheta = fTheta * (3.14 / 180.0f);
	mat4x4 matRotX =
	{
		1, 0, 0, 0,
		0, cos(fTheta), -sin(fTheta), 0,
		0, sin(fTheta), cos(fTheta), 0,
		0, 0, 0, 1
	};

	//// Rotation X
	//matRotX.m[0][0] = 1;
	//matRotX.m[1][1] = cosf(fTheta * 0.5f);
	//matRotX.m[1][2] = sinf(fTheta * 0.5f);
	//matRotX.m[2][1] = -sinf(fTheta * 0.5f);
	//matRotX.m[2][2] = cosf(fTheta * 0.5f);
	//matRotX.m[3][3] = 1;

	return matRotX;
}

mat4x4 matrixRotationY4x4(float fTheta)
{
	fTheta = fTheta * (3.14 / 180.0f);
	//Y Rotation
	mat4x4 matRotY =
	{
		cosf(fTheta), 0, sinf(fTheta), 0,
		0, 1, 0, 0,
		-sinf(fTheta), 0, cosf(fTheta), 0,
		0, 0, 0, 1
	};

	//mat4x4 matRotY = { 0 };
	////time of rotation
	////fTheta = fTheta / 180.0f * M_PI;
	//fTheta = fTheta / 180.0f * 3.14;
	//// Rotation X
	//matRotY.m[0][0] = cosf(fTheta * 0.5f);
	//matRotY.m[0][2] = sinf(fTheta * 0.5f);
	//matRotY.m[2][1] = -sinf(fTheta * 0.5f);
	//matRotY.m[2][1] = cosf(fTheta * 0.5f);
	//matRotY.m[2][2] = 1;
	//matRotY.m[3][3] = 1;

	return matRotY;
}

mat4x4 matrixTranslationb(VERTEX t)
{
	mat4x4 m =
	{
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			t.x, t.y, t.z, 1,
	};

	/*mat4x4 m = matrixIdentity();

	m.m[3][0] = x;
	m.m[3][1] = y;
	m.m[3][2] = z;*/

	return m;
}

mat4x4 inverseb(mat4x4 m) {

	mat4x4 mat;
	float t;
	int i;
	VERTEX verx;

	/*t = m.m[1][0];
	t = m.m[2][0];
	t = m.m[2][1];*/
	for (i = 1; i <= 1; i++)
	{
		t = m.m[i][0] = m.m[0][i];
		//m.m[i][0] = m.m[0][i];
	}
	/*t = m.m[1][0];
	m.m[1][0] = m.m[0][1];
	m.m[0][1] = t;
	t = m.m[2][0];
	m.m[2][0] = m.m[0][2];
	m.m[0][2] = t;*/
	t = m.m[2][1];
	m.m[2][1] = m.m[1][2];
	m.m[1][2] = t;

	verx = { m.m[3][0],  m.m[3][1],  m.m[3][2],  m.m[3][3] };
	mat =
	{ m.m[0][0],  m.m[0][1], m.m[0][2], m.m[0][3],
	  m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
	  m.m[2][0],  m.m[2][1],  m.m[2][2],  m.m[2][3]
		/*, m.m[3][0],  m.m[3][1],  m.m[3][2],  m.m[3][3]*/
	};

	verx = vertexMultMatrix4x4(verx, mat);

	m.m[3][0] = verx.x * -1;
	m.m[3][1] = verx.y * -1;
	m.m[3][2] = verx.z * -1;
	//m.m[3][3] = tv.w * -1;
	return m;
}

mat4x4 ProjectionMatrixb()
{

	// Projection Matrix
	float fNear = 0.1f;
	float fFar = 10.0f;
	float fFov = 90.0f;
	float fAspectRatio = (float)Raster_H / (float)Raster_W;
	float fFovRad = 1.0f / tanf((fFov * 0.5f) * 3.14159f / 180.0f);
	float fFovRad2 = 1.0f / tanf((fFov * 0.5f));
	float YScale = fFovRad;
	float XScale = fAspectRatio * YScale;

	mat4x4 matProj =
	{
			XScale, 0, 0, 0,
			0, YScale, 0, 0,
			0, 0, fFar / (fFar - fNear), 1,
			0, 0, (-fFar * fNear) / (fFar - fNear), 0
	};

	/*matProj.m[0][0] = fAspectRatio * fFovRad;
	matProj.m[1][1] = fFovRad;
	matProj.m[2][2] = fFar / (fFar - fNear);
	matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matProj.m[2][3] = 1.0f;
	matProj.m[3][3] = 0.0f;*/

	return matProj;
}

mat4x4 ProjectionMatrixc()
{

	// Projection Matrix
	float fNear = nearPlane;
	float fFar = farPlane;
	float fFov = 90.0f;
	float fAspectRatio = (float)Raster_H / (float)Raster_W;
	float fFovRad = 1.0f / tanf((fFov * 0.5f) * 3.14159f / 180.0f);
	float fFovRad2 = 1.0f / tanf((fFov * 0.5f));
	float YScale = fFovRad;
	float XScale = fAspectRatio * YScale;

	mat4x4 matProj =
	{
			XScale, 0, 0, 0,
			0, YScale, 0, 0,
			0, 0, fFar / (fFar - fNear), 1,
			0, 0, (-fFar * fNear) / (fFar - fNear), 0
	};

	/*matProj.m[0][0] = fAspectRatio * fFovRad;
	matProj.m[1][1] = fFovRad;
	matProj.m[2][2] = fFar / (fFar - fNear);
	matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matProj.m[2][3] = 1.0f;
	matProj.m[3][3] = 0.0f;*/

	return matProj;
}

#pragma endregion Matrix Calculations
