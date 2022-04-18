#pragma once

#include "calc.h"


void (*VertexShader)(VERTEX&) = 0;
void (*PixelShader)(unsigned int&) = 0;
MATRIX_4X4 SV_WorldMatrix;
MATRIX_4X4 viewmatrix;
MATRIX_4X4 projectmatrix;

void VS_World(VERTEX& multiplyMe)
{
   
    multiplyMe = vertexmatrixmultiply(multiplyMe, SV_WorldMatrix);
    VERTEX copy = multiplyMe;
    multiplyMe = vertexmatrixmultiply(multiplyMe, viewmatrix);
    multiplyMe = vertexmatrixmultiply(multiplyMe, projectmatrix);
    VERTEX normal = { multiplyMe.nx, multiplyMe.ny, multiplyMe.nz };
        normal = NormalizeVector(normal);

        VERTEX dir = { -sunlightDirection.x, -sunlightDirection.y, -sunlightDirection.z };
        dir = NormalizeVector(dir);

        float dotdot = dot(dir, normal);
        float lightRatio = saturate(dotdot);

        unsigned int result = ColorLerp(0xFF000000, sunColor, lightRatio);


        VERTEX pointdir = { (pointLight.x - copy.x), (pointLight.y - copy.y), (pointLight.z - copy.z) };

        float attenuation = 1 - saturate(vlength(pointdir) / LightRadius);

        pointdir = NormalizeVector(pointdir);

        dotdot = dot(pointdir, normal);
        lightRatio = saturate(dotdot) * attenuation;


        unsigned int pointresult = ColorLerp(0xFF000000, pointLightColor, lightRatio);


        multiplyMe.color = addcolors(result, pointresult);
    /*multiplyMe.x = multiplyMe.x / multiplyMe.w;
    multiplyMe.y = multiplyMe.y / multiplyMe.w;
    multiplyMe.z = multiplyMe.z / multiplyMe.w;*/
}
void PS_White(unsigned int& color) 
{
    color = 0xFFFFFFFF;
}
void PS_Green(unsigned int& color) 
{
    color = 0xFF00FF00;
}