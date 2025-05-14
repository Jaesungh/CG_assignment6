#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <limits>
#include <cstring>
#include <cstdio>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const int WIDTH = 512;
const int HEIGHT = 512;

enum ShadingMode { FLAT, GOURAUD, PHONG };

// 3D vector
struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator*(const Vec3& o) const { return Vec3(x * o.x, y * o.y, z * o.z); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    friend Vec3 operator*(float s, const Vec3& v) { return v * s; }
};

inline float dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}
inline float length(const Vec3& v) {
    return std::sqrt(dot(v, v));
}
inline Vec3 normalize(const Vec3& v) {
    float len = length(v);
    return (len > 1e-6f ? v / len : Vec3(0, 0, 0));
}

// 4D vector for pipeline
struct Vec4 { float x, y, z, w; };
// raw vertex
struct Vertex { float x, y, z; };
// RGB pixel
struct Color { unsigned char r, g, b; };


Vertex* gVertexBuffer = nullptr;
Vec3* gNormalBuffer = nullptr;
int* gIndexBuffer = nullptr;
int     gNumVertices = 0;
int     gNumTriangles = 0;
float   zBuffer[HEIGHT][WIDTH];
Color   image[HEIGHT][WIDTH];

ShadingMode gMode = PHONG;


Vec4 modelTransform(const Vertex& v) {
    return { 2 * v.x, 2 * v.y, 2 * v.z - 7.0f, 1.0f };
}
Vec4 cameraTransform(const Vec4& v) {
    return v;
}
Vec4 perspectiveProject(const Vec4& v) {
    float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, n = -0.1f, f = -1000.0f;
    if (std::fabs(v.z) < 1e-5f) return { 0,0,0,1 };
    float x = (2 * n * v.x) / ((r - l) * v.z);
    float y = (2 * n * v.y) / ((t - b) * v.z);
    float z = (f + n) / (n - f) + (2 * f * n) / ((n - f) * v.z);
    return { x,y,z,1 };
}

inline int viewportX(float x_ndc) { return int((x_ndc + 1) * 0.5f * WIDTH); }
inline int viewportY(float y_ndc) { return int((1 - (y_ndc + 1) * 0.5f) * HEIGHT); }


Color shade(const Vec3& N, const Vec3& V, const Vec3& P) {
    Vec3 ka{ 0,1,0 }, kd{ 0,0.5f,0 }, ks{ 0.5f,0.5f,0.5f };
    float p = 32;
    Vec3 Ia{ 0.2f,0.2f,0.2f }, Il{ 1,1,1 };
    Vec3 Lpos{ -4,4,-3 };
    Vec3 L = normalize(Lpos - P);
    Vec3 H = normalize(L + V);
    float NdotL = std::max(dot(N, L), 0.0f);
    float NdotH = std::max(dot(N, H), 0.0f);
    Vec3 I = Ia * ka + Il * (kd * NdotL + ks * std::pow(NdotH, p));
    float gamma = 2.2f;
    I.x = std::pow(I.x, 1.0f / gamma);
    I.y = std::pow(I.y, 1.0f / gamma);
    I.z = std::pow(I.z, 1.0f / gamma);
    Color c;
    c.r = (unsigned char)(255 * std::min(I.x, 1.0f));
    c.g = (unsigned char)(255 * std::min(I.y, 1.0f));
    c.b = (unsigned char)(255 * std::min(I.z, 1.0f));
    return c;
}

void setPixel(int x, int y, float z, const Color& c) {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return;
    if (z < zBuffer[y][x]) {
        zBuffer[y][x] = z;
        image[y][x] = c;
    }
}

void drawTriangle(const Vec4& v0, const Vec4& v1, const Vec4& v2,
    const Vec3& p0, const Vec3& p1, const Vec3& p2,
    const Vec3& n0, const Vec3& n1, const Vec3& n2) {
    int x0 = viewportX(v0.x), y0 = viewportY(v0.y);
    int x1 = viewportX(v1.x), y1 = viewportY(v1.y);
    int x2 = viewportX(v2.x), y2 = viewportY(v2.y);
    int minX = std::max(0, std::min({ x0,x1,x2 }));
    int maxX = std::min(WIDTH - 1, std::max({ x0,x1,x2 }));
    int minY = std::max(0, std::min({ y0,y1,y2 }));
    int maxY = std::min(HEIGHT - 1, std::max({ y0,y1,y2 }));
    float denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
    if (std::fabs(denom) < 1e-6f) return;

    Vec3 faceN = normalize(cross(p1 - p0, p2 - p0));
    Vec3 cent = (p0 + p1 + p2) / 3;
    Vec3 Vc = normalize(-cent);
    Color flatC = shade(faceN, Vc, cent);
   
    Vec3 V0 = normalize(-p0), V1 = normalize(-p1), V2 = normalize(-p2);
    Color c0 = shade(n0, V0, p0), c1 = shade(n1, V1, p1), c2 = shade(n2, V2, p2);
    for (int y = minY;y <= maxY;++y) for (int x = minX;x <= maxX;++x) {
        float l1 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / denom;
        float l2 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / denom;
        float l3 = 1 - l1 - l2;
        if (l1 < 0 || l2 < 0 || l3 < 0) continue;
        float z = l1 * v0.z + l2 * v1.z + l3 * v2.z;
        Color out;
        switch (gMode) {
        case FLAT:     out = flatC; break;
        case GOURAUD:
            out.r = (unsigned char)(l1 * c0.r + l2 * c1.r + l3 * c2.r);
            out.g = (unsigned char)(l1 * c0.g + l2 * c1.g + l3 * c2.g);
            out.b = (unsigned char)(l1 * c0.b + l2 * c1.b + l3 * c2.b);
            break;
        case PHONG: {
            Vec3 Np = normalize(n0 * l1 + n1 * l2 + n2 * l3);
            Vec3 Pp = p0 * l1 + p1 * l2 + p2 * l3;
            out = shade(Np, normalize(-Pp), Pp);
            break;
        }
        }
        setPixel(x, y, z, out);
    }
}

void create_scene() {
    int w = 32, h = 16;
    gNumVertices = (h - 2) * w + 2;
    gNumTriangles = (h - 3) * (w - 1) * 2 + 2 * (w - 1);
    gVertexBuffer = new Vertex[gNumVertices];
    gNormalBuffer = new Vec3[gNumVertices];
    gIndexBuffer = new int[3 * gNumTriangles];
    int t = 0;
    for (int j = 1;j < h - 1;++j) for (int i = 0;i < w;++i) {
        float theta = j / (float)(h - 1) * M_PI;
        float phi = i / (float)(w - 1) * 2 * M_PI;
        float x = std::sinf(theta) * std::cosf(phi);
        float y = std::cosf(theta);
        float z = -std::sinf(theta) * std::sinf(phi);
        gVertexBuffer[t] = { x,y,z };
        gNormalBuffer[t++] = normalize(Vec3{ x,y,z });
    }
    gVertexBuffer[t] = { 0,1,0 };   gNormalBuffer[t++] = Vec3{ 0,1,0 };
    gVertexBuffer[t] = { 0,-1,0 };  gNormalBuffer[t++] = Vec3{ 0,-1,0 };
    t = 0;
    for (int j = 0;j < h - 3;++j) for (int i = 0;i < w - 1;++i) {
        gIndexBuffer[t++] = j * w + i;
        gIndexBuffer[t++] = (j + 1) * w + (i + 1);
        gIndexBuffer[t++] = j * w + (i + 1);
        gIndexBuffer[t++] = j * w + i;
        gIndexBuffer[t++] = (j + 1) * w + i;
        gIndexBuffer[t++] = (j + 1) * w + (i + 1);
    }
    int top = (h - 2) * w, bot = top + 1;
    for (int i = 0;i < w - 1;++i) {
        gIndexBuffer[t++] = top;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = bot;
        gIndexBuffer[t++] = (h - 3) * w + (i + 1);
        gIndexBuffer[t++] = (h - 3) * w + i;
    }
}


int main() {
    create_scene();

    
    ShadingMode modes[3] = { FLAT, GOURAUD, PHONG };
    const char* names[3] = { "flat", "gouraud", "phong" };

    for (int m = 0; m < 3; ++m) {
        gMode = modes[m];
        
        std::memset(image, 0, sizeof(image));
        for (int y = 0; y < HEIGHT; ++y)
            for (int x = 0; x < WIDTH; ++x)
                zBuffer[y][x] = std::numeric_limits<float>::infinity();

        
        for (int i = 0; i < gNumTriangles; ++i) {
            int k0 = gIndexBuffer[3 * i + 0];
            int k1 = gIndexBuffer[3 * i + 1];
            int k2 = gIndexBuffer[3 * i + 2];
            Vec4 v0 = perspectiveProject(cameraTransform(modelTransform(gVertexBuffer[k0])));
            Vec4 v1 = perspectiveProject(cameraTransform(modelTransform(gVertexBuffer[k1])));
            Vec4 v2 = perspectiveProject(cameraTransform(modelTransform(gVertexBuffer[k2])));
            Vec4 e0 = cameraTransform(modelTransform(gVertexBuffer[k0]));
            Vec4 e1 = cameraTransform(modelTransform(gVertexBuffer[k1]));
            Vec4 e2 = cameraTransform(modelTransform(gVertexBuffer[k2]));
            Vec3 p0{ e0.x, e0.y, e0.z };
            Vec3 p1{ e1.x, e1.y, e1.z };
            Vec3 p2{ e2.x, e2.y, e2.z };
            drawTriangle(v0, v1, v2,
                p0, p1, p2,
                gNormalBuffer[k0], gNormalBuffer[k1], gNormalBuffer[k2]);
        }

        
        char filename[64];
        std::snprintf(filename, sizeof(filename), "output_%s.ppm", names[m]);
        FILE* f = std::fopen(filename, "wb");
        std::fprintf(f, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
        std::fwrite(image, sizeof(Color), WIDTH * HEIGHT, f);
        std::fclose(f);
        std::cout << "Rendered " << filename << "\n";
    }

    delete[] gVertexBuffer;
    delete[] gNormalBuffer;
    delete[] gIndexBuffer;
    return 0;
}
