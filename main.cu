#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <zlib.h>
#include <cuda_runtime.h>
#include <map>
#include <filesystem>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " \
                      << cudaGetErrorString(err) << std::endl; \
            exit(1); \
        } \
    } while(0)

struct Vec3 {
    double x, y, z;
};

struct Vec2 {
    double u, v;
};

struct Vec3i {
    int x, y, z;
};

struct Color {
    unsigned char r, g, b;
};

struct Face {
    int v1, v2, v3;
    int t1, t2, t3;
    int materialId;
};

struct Material {
    std::string name;
    std::string textureFile;
    int textureId;
    int texWidth;
    int texHeight;
};

struct BlockTexture {
    std::string name;
    Color avgColor;
    int paletteId;
};

struct OBJ {
    std::vector<Vec3> verts;
    std::vector<Vec2> texCoords;
    std::vector<Face> faces;
    std::vector<Material> materials;
    double minX = 0, minY = 0, minZ = 0;
    double maxX = 0, maxY = 0, maxZ = 0;
    std::string mtlFile;

    void calcBounds() {
        if (verts.empty()) return;
        minX = maxX = verts[0].x;
        minY = maxY = verts[0].y;
        minZ = maxZ = verts[0].z;

        for (const auto& v : verts) {
            minX = std::min(minX, v.x); maxX = std::max(maxX, v.x);
            minY = std::min(minY, v.y); maxY = std::max(maxY, v.y);
            minZ = std::min(minZ, v.z); maxZ = std::max(maxZ, v.z);
        }
    }
};

Color getAverageColor(const std::string& filepath) {
    int w, h, channels;
    unsigned char* img = stbi_load(filepath.c_str(), &w, &h, &channels, 3);
    
    if (!img) {
        std::cerr << "Failed to load texture: " << filepath << std::endl;
        return {128, 128, 128};
    }

    long long r = 0, g = 0, b = 0;
    int count = 0;

    for (int i = 0; i < w * h * 3; i += 3) {
        r += img[i];
        g += img[i + 1];
        b += img[i + 2];
        count++;
    }

    stbi_image_free(img);

    return {
        (unsigned char)(r / count),
        (unsigned char)(g / count),
        (unsigned char)(b / count)
    };
}

std::vector<BlockTexture> loadBlockTextures(const std::string& texturesDir) {
    std::vector<BlockTexture> blocks;
    int paletteId = 1;

    std::cout << "Loading block textures from: " << texturesDir << std::endl;

    for (const auto& entry : std::filesystem::directory_iterator(texturesDir)) {
        if (entry.path().extension() == ".png") {
            std::string filename = entry.path().stem().string();
            Color avgColor = getAverageColor(entry.path().string());
            
            blocks.push_back({
                "minecraft:" + filename,
                avgColor,
                paletteId++
            });
        }
    }

    std::cout << "Loaded " << blocks.size() << " block textures" << std::endl;
    return blocks;
}

OBJ loadOBJ(const std::string& file) {
    OBJ obj;
    std::ifstream in(file);
    std::string line;
    int currentMaterial = -1;

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        
        if (type == "mtllib") {
            iss >> obj.mtlFile;
            break;
        }
    }
    
    in.clear();
    in.seekg(0);

    if (!obj.mtlFile.empty()) {
        std::string mtlPath;
        std::string objDir = std::filesystem::path(file).parent_path().string();
        
        if (objDir.empty()) {
            mtlPath = obj.mtlFile;
        } else {
            mtlPath = objDir + "/" + obj.mtlFile;
        }
        
        std::cout << "Looking for MTL file: " << mtlPath << std::endl;
        
        std::ifstream mtlIn(mtlPath);
        if (!mtlIn.is_open()) {
            std::cerr << "Failed to open MTL file: " << mtlPath << std::endl;
            mtlIn.open(obj.mtlFile);
            if (mtlIn.is_open()) {
                std::cout << "Found MTL in current directory: " << obj.mtlFile << std::endl;
                mtlPath = obj.mtlFile;
            }
        }
        
        std::string mtlLine;
        Material currentMtl;
        bool hasMaterial = false;
        
        while (std::getline(mtlIn, mtlLine)) {
            std::istringstream mtlIss(mtlLine);
            std::string mtlType;
            mtlIss >> mtlType;
            
            if (mtlType == "newmtl") {
                if (hasMaterial) {
                    obj.materials.push_back(currentMtl);
                    std::cout << "  Added material: " << currentMtl.name << std::endl;
                }
                currentMtl = Material();
                mtlIss >> currentMtl.name;
                currentMtl.textureFile = "";
                hasMaterial = true;
                std::cout << "Found material: " << currentMtl.name << std::endl;
            } else if (mtlType == "map_Kd") {
                std::string texPath;
                std::getline(mtlIss, texPath);
                size_t start = texPath.find_first_not_of(" \t");
                if (start != std::string::npos) {
                    texPath = texPath.substr(start);
                }
                
                std::cout << "  Texture path from MTL: " << texPath << std::endl;
                
                std::string filename = texPath;
                size_t lastSlash = filename.find_last_of("/\\");
                if (lastSlash != std::string::npos) {
                    filename = filename.substr(lastSlash + 1);
                }
                
                std::cout << "  Extracted filename: " << filename << std::endl;
                
                std::string objDir = std::filesystem::path(file).parent_path().string();
                if (objDir.empty()) objDir = ".";
                
                std::string mtlDir = std::filesystem::path(mtlPath).parent_path().string();
                if (mtlDir.empty()) mtlDir = ".";
                
                std::vector<std::string> searchPaths = {
                    texPath,
                    objDir + "/" + filename,
                    mtlDir + "/" + filename,
                    filename
                };
                
                bool found = false;
                for (const auto& searchPath : searchPaths) {
                    if (std::filesystem::exists(searchPath)) {
                        currentMtl.textureFile = searchPath;
                        std::cout << "  Found texture: " << searchPath << std::endl;
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    std::cerr << "  Could not find texture: " << filename << std::endl;
                    std::cerr << "    Searched in:" << std::endl;
                    for (const auto& searchPath : searchPaths) {
                        std::cerr << "      " << searchPath << std::endl;
                    }
                }
            }
        }
        
        if (hasMaterial) {
            obj.materials.push_back(currentMtl);
            std::cout << "  Added material: " << currentMtl.name << std::endl;
        }
        
        if (mtlIn.is_open()) {
            std::cout << "Loaded " << obj.materials.size() << " materials from MTL" << std::endl;
        } else {
            std::cout << "Could not open MTL file" << std::endl;
        }
    }

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v") {
            double x, y, z;
            iss >> x >> y >> z;
            obj.verts.push_back({x, y, z});
        } else if (type == "vt") {
            double u, v;
            iss >> u >> v;
            obj.texCoords.push_back({u, v});
        } else if (type == "usemtl") {
            std::string mtlName;
            iss >> mtlName;
            currentMaterial = -1;
            for (size_t i = 0; i < obj.materials.size(); i++) {
                if (obj.materials[i].name == mtlName) {
                    currentMaterial = i;
                    break;
                }
            }
        } else if (type == "f") {
            std::string v1s, v2s, v3s;
            iss >> v1s >> v2s >> v3s;

            auto parseVertex = [](const std::string& s, int& v, int& t) {
                size_t slash1 = s.find('/');
                if (slash1 == std::string::npos) {
                    v = std::stoi(s) - 1;
                    t = -1;
                } else {
                    v = std::stoi(s.substr(0, slash1)) - 1;
                    size_t slash2 = s.find('/', slash1 + 1);
                    std::string tstr = s.substr(slash1 + 1, slash2 - slash1 - 1);
                    t = tstr.empty() ? -1 : std::stoi(tstr) - 1;
                }
            };

            Face f;
            parseVertex(v1s, f.v1, f.t1);
            parseVertex(v2s, f.v2, f.t2);
            parseVertex(v3s, f.v3, f.t3);
            f.materialId = currentMaterial;
            obj.faces.push_back(f);
        }
    }

    obj.calcBounds();
    std::cout << "Loaded " << obj.verts.size() << " vertices, " 
              << obj.texCoords.size() << " tex coords, " 
              << obj.faces.size() << " faces" << std::endl;
    return obj;
}

__device__ inline int max3(int a, int b, int c) {
    int m = a;
    if (b > m) m = b;
    if (c > m) m = c;
    return m;
}

__device__ inline int colorDistance(unsigned char r1, unsigned char g1, unsigned char b1,
                                   unsigned char r2, unsigned char g2, unsigned char b2) {
    int dr = (int)r1 - (int)r2;
    int dg = (int)g1 - (int)g2;
    int db = (int)b1 - (int)b2;
    return dr*dr + dg*dg + db*db;
}

__device__ int findBestBlock(unsigned char r, unsigned char g, unsigned char b,
                            const Color* blockColors, int numBlocks) {
    int bestIdx = 0;
    int bestDist = 999999;

    for (int i = 0; i < numBlocks; i++) {
        int dist = colorDistance(r, g, b, blockColors[i].r, blockColors[i].g, blockColors[i].b);
        if (dist < bestDist) {
            bestDist = dist;
            bestIdx = i;
        }
    }

    return bestIdx;
}

__device__ void getTexColor(const unsigned char* textures, int texOffset, int texW, int texH,
                           double u, double v, unsigned char& r, unsigned char& g, unsigned char& b) {
    u = u - floor(u);
    v = v - floor(v);
    
    int x = (int)(u * texW) % texW;
    int y = (int)((1.0 - v) * texH) % texH;
    
    if (x < 0) x += texW;
    if (y < 0) y += texH;
    
    int idx = texOffset + (y * texW + x) * 3;
    r = textures[idx];
    g = textures[idx + 1];
    b = textures[idx + 2];
}

__device__ void drawLine3D(unsigned int* voxels, int* blockIds, int w, int h, int l,
                          Vec3i p1, Vec3i p2, int blockId) {
    int dx = abs(p2.x - p1.x);
    int dy = abs(p2.y - p1.y);
    int dz = abs(p2.z - p1.z);

    int sx = (p1.x < p2.x) ? 1 : -1;
    int sy = (p1.y < p2.y) ? 1 : -1;
    int sz = (p1.z < p2.z) ? 1 : -1;

    int dm = max3(dx, dy, dz);
    int i = dm;

    int x = p1.x, y = p1.y, z = p1.z;

    int x1 = dx - dm/2;
    int y1 = dy - dm/2;
    int z1 = dz - dm/2;

    while (i-- >= 0) {
        if (x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l) {
            int idx = x + z * w + y * w * l;
            atomicOr(&voxels[idx], 1);
            atomicMax(&blockIds[idx], blockId);
        }

        x1 += dx;
        if (x1 >= dm) { x += sx; x1 -= dm; }

        y1 += dy;
        if (y1 >= dm) { y += sy; y1 -= dm; }

        z1 += dz;
        if (z1 >= dm) { z += sz; z1 -= dm; }
    }
}

__device__ void rasterizeTriangle(unsigned int* voxels, int* blockIds, int w, int h, int l,
                                 Vec3i v1, Vec3i v2, Vec3i v3,
                                 double u1, double v1t, double u2, double v2t, double u3, double v3t,
                                 const unsigned char* textures, int texOffset, int texW, int texH,
                                 const Color* blockColors, const int* blockPaletteIds, int numBlocks) {
    
    unsigned char r1, g1, b1, r2, g2, b2, r3, g3, b3;
    getTexColor(textures, texOffset, texW, texH, u1, v1t, r1, g1, b1);
    getTexColor(textures, texOffset, texW, texH, u2, v2t, r2, g2, b2);
    getTexColor(textures, texOffset, texW, texH, u3, v3t, r3, g3, b3);
    
    unsigned char avgR = (r1 + r2 + r3) / 3;
    unsigned char avgG = (g1 + g2 + g3) / 3;
    unsigned char avgB = (b1 + b2 + b3) / 3;
    
    int blockIdx = findBestBlock(avgR, avgG, avgB, blockColors, numBlocks);
    int blockId = blockPaletteIds[blockIdx];

    drawLine3D(voxels, blockIds, w, h, l, v1, v2, blockId);
    drawLine3D(voxels, blockIds, w, h, l, v2, v3, blockId);
    drawLine3D(voxels, blockIds, w, h, l, v3, v1, blockId);

    if (v1.y > v2.y) {
        Vec3i tv = v1; v1 = v2; v2 = tv;
        double tu = u1; u1 = u2; u2 = tu;
        double tvt = v1t; v1t = v2t; v2t = tvt;
    }
    if (v2.y > v3.y) {
        Vec3i tv = v2; v2 = v3; v3 = tv;
        double tu = u2; u2 = u3; u3 = tu;
        double tvt = v2t; v2t = v3t; v3t = tvt;
    }
    if (v1.y > v2.y) {
        Vec3i tv = v1; v1 = v2; v2 = tv;
        double tu = u1; u1 = u2; u2 = tu;
        double tvt = v1t; v1t = v2t; v2t = tvt;
    }

    for (int y = v1.y; y <= v3.y; y++) {
        Vec3i intersections[3];
        int count = 0;

        auto addIntersection = [&](Vec3i a, Vec3i b, int cy) {
            if (a.y == b.y) return;
            if ((a.y <= cy && cy < b.y) || (b.y <= cy && cy < a.y)) {
                double t = (double)(cy - a.y) / (b.y - a.y);
                int x = a.x + t * (b.x - a.x);
                int z = a.z + t * (b.z - a.z);
                if (count < 3) {
                    intersections[count++] = {x, cy, z};
                }
            }
        };

        addIntersection(v1, v2, y);
        addIntersection(v2, v3, y);
        addIntersection(v3, v1, y);

        if (count > 1) {
             if (intersections[0].x > intersections[1].x) {
                 Vec3i temp = intersections[0];
                 intersections[0] = intersections[1];
                 intersections[1] = temp;
             }
             drawLine3D(voxels, blockIds, w, h, l, intersections[0], intersections[1], blockId);
        }
    }
}

__global__ void rasterizeFaces(const Vec3i* vertices, const Vec2* texCoords, const Face* faces,
                              int numFaces, unsigned int* voxels, int* blockIds,
                              int w, int h, int l,
                              const unsigned char* textures, const int* texOffsets,
                              const int* texWidths, const int* texHeights,
                              const Color* blockColors, const int* blockPaletteIds, int numBlocks) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numFaces) return;

    Face f = faces[idx];
    Vec3i v1 = vertices[f.v1];
    Vec3i v2 = vertices[f.v2];
    Vec3i v3 = vertices[f.v3];

    double u1 = 0.5, vt1 = 0.5;
    double u2 = 0.5, vt2 = 0.5;
    double u3 = 0.5, vt3 = 0.5;

    if (f.t1 >= 0) {
        Vec2 tc = texCoords[f.t1];
        u1 = tc.u; vt1 = tc.v;
    }
    if (f.t2 >= 0) {
        Vec2 tc = texCoords[f.t2];
        u2 = tc.u; vt2 = tc.v;
    }
    if (f.t3 >= 0) {
        Vec2 tc = texCoords[f.t3];
        u3 = tc.u; vt3 = tc.v;
    }

    int matId = f.materialId >= 0 ? f.materialId : 0;
    int texOffset = texOffsets[matId];
    int texW = texWidths[matId];
    int texH = texHeights[matId];

    rasterizeTriangle(voxels, blockIds, w, h, l, v1, v2, v3,
                     u1, vt1, u2, vt2, u3, vt3,
                     textures, texOffset, texW, texH,
                     blockColors, blockPaletteIds, numBlocks);
}

void writeShort(std::vector<uint8_t>& buf, uint16_t val) {
    buf.push_back((val >> 8) & 0xff);
    buf.push_back(val & 0xff);
}

void writeInt(std::vector<uint8_t>& buf, int32_t val) {
    buf.push_back((val >> 24) & 0xff);
    buf.push_back((val >> 16) & 0xff);
    buf.push_back((val >> 8) & 0xff);
    buf.push_back(val & 0xff);
}

void writeString(std::vector<uint8_t>& buf, const std::string& str) {
    writeShort(buf, str.size());
    buf.insert(buf.end(), str.begin(), str.end());
}

void writeVarint(std::vector<uint8_t>& buf, int val) {
    uint32_t uval = val;
    while ((uval & 0xFFFFFF80) != 0) {
        buf.push_back((uval & 0x7F) | 0x80);
        uval >>= 7;
    }
    buf.push_back(uval & 0x7F);
}
//Sponge V3 see https://github.com/SpongePowered/Schematic-Specification/blob/master/versions/schematic-3.md
void writeSchematic(const std::string& file, const std::vector<std::vector<std::vector<int>>>& blocks_data,
                   int w, int h, int l, const std::vector<BlockTexture>& blockTextures) {
    std::vector<uint8_t> nbt;

    nbt.push_back(10);
    writeString(nbt, "");

    nbt.push_back(10);
    writeString(nbt, "Schematic");

    nbt.push_back(3);
    writeString(nbt, "Version");
    writeInt(nbt, 3);

    nbt.push_back(3);
    writeString(nbt, "DataVersion");
    writeInt(nbt, 2975);

    nbt.push_back(2);
    writeString(nbt, "Width");
    writeShort(nbt, w);

    nbt.push_back(2);
    writeString(nbt, "Height");
    writeShort(nbt, h);

    nbt.push_back(2);
    writeString(nbt, "Length");
    writeShort(nbt, l);

    nbt.push_back(10);
    writeString(nbt, "Metadata");
    nbt.push_back(0);

    nbt.push_back(10);
    writeString(nbt, "Blocks");

    nbt.push_back(10);
    writeString(nbt, "Palette");

    nbt.push_back(3);
    writeString(nbt, "minecraft:air");
    writeInt(nbt, 0);

    for (const auto& block : blockTextures) {
        nbt.push_back(3);
        writeString(nbt, block.name);
        writeInt(nbt, block.paletteId);
    }

    nbt.push_back(0);

    std::vector<uint8_t> blockData;
    for (int y = 0; y < h; y++) {
        for (int z = 0; z < l; z++) {
            for (int x = 0; x < w; x++) {
                writeVarint(blockData, blocks_data[x][y][z]);
            }
        }
    }

    nbt.push_back(7);
    writeString(nbt, "Data");
    writeInt(nbt, blockData.size());
    nbt.insert(nbt.end(), blockData.begin(), blockData.end());

    nbt.push_back(0);

    nbt.push_back(0);
    nbt.push_back(0);

    gzFile gz = gzopen(file.c_str(), "wb9");
    if (!gz) {
        std::cerr << "Failed to open output file" << std::endl;
        return;
    }
    gzwrite(gz, nbt.data(), nbt.size());
    gzclose(gz);
}

void convert(const OBJ& obj, const std::string& output, int target, 
            const std::vector<BlockTexture>& blockTextures) {
    double sizeX = obj.maxX - obj.minX;
    double sizeY = obj.maxY - obj.minY;
    double sizeZ = obj.maxZ - obj.minZ;
    double maxSize = std::max({sizeX, sizeY, sizeZ});
    double scale = target / maxSize;

    int w = (int)(sizeX * scale) + 1;
    int h = (int)(sizeY * scale) + 1;
    int l = (int)(sizeZ * scale) + 1;

    std::cout << "Dimensions: " << w << "x" << h << "x" << l << std::endl;

    std::vector<Vec3i> vertices(obj.verts.size());
    for (size_t i = 0; i < obj.verts.size(); i++) {
        vertices[i].x = (int)((obj.verts[i].x - obj.minX) * scale);
        vertices[i].y = (int)((obj.verts[i].y - obj.minY) * scale);
        vertices[i].z = (int)((obj.verts[i].z - obj.minZ) * scale);
    }

    std::vector<unsigned char> allTextures;
    std::vector<int> texOffsets(obj.materials.size());
    std::vector<int> texWidths(obj.materials.size());
    std::vector<int> texHeights(obj.materials.size());
    
    for (size_t i = 0; i < obj.materials.size(); i++) {
        texOffsets[i] = allTextures.size();
        
        if (!obj.materials[i].textureFile.empty()) {
            int w, h, channels;
            unsigned char* img = stbi_load(obj.materials[i].textureFile.c_str(), &w, &h, &channels, 3);
            
            if (img) {
                std::cout << "Loaded material " << obj.materials[i].name << " texture: " 
                         << obj.materials[i].textureFile << " (" << w << "x" << h << ")" << std::endl;
                texWidths[i] = w;
                texHeights[i] = h;
                allTextures.insert(allTextures.end(), img, img + w * h * 3);
                stbi_image_free(img);
            } else {
                std::cerr << "Failed to load texture: " << obj.materials[i].textureFile << std::endl;
                texWidths[i] = 1;
                texHeights[i] = 1;
                allTextures.push_back(128);
                allTextures.push_back(128);
                allTextures.push_back(128);
            }
        } else {
            texWidths[i] = 1;
            texHeights[i] = 1;
            allTextures.push_back(128);
            allTextures.push_back(128);
            allTextures.push_back(128);
        }
    }

    if (obj.materials.empty()) {
        std::cout << "No materials found, using default color" << std::endl;
        texOffsets.push_back(0);
        texWidths.push_back(1);
        texHeights.push_back(1);
        allTextures = {128, 128, 128};
    }

    std::vector<Color> blockColors(blockTextures.size());
    std::vector<int> blockPaletteIds(blockTextures.size());
    for (size_t i = 0; i < blockTextures.size(); i++) {
        blockColors[i] = blockTextures[i].avgColor;
        blockPaletteIds[i] = blockTextures[i].paletteId;
    }

    Vec3i* d_vertices;
    Vec2* d_texCoords;
    Face* d_faces;
    unsigned int* d_voxels;
    int* d_blockIds;
    unsigned char* d_textures;
    int* d_texOffsets;
    int* d_texWidths;
    int* d_texHeights;
    Color* d_blockColors;
    int* d_blockPaletteIds;

    size_t voxelSize = w * h * l;
    size_t voxelMemSize = voxelSize * sizeof(unsigned int);

    CUDA_CHECK(cudaMalloc(&d_vertices, vertices.size() * sizeof(Vec3i)));
    CUDA_CHECK(cudaMalloc(&d_texCoords, obj.texCoords.size() * sizeof(Vec2)));
    CUDA_CHECK(cudaMalloc(&d_faces, obj.faces.size() * sizeof(Face)));
    CUDA_CHECK(cudaMalloc(&d_voxels, voxelMemSize));
    CUDA_CHECK(cudaMalloc(&d_blockIds, voxelSize * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_textures, allTextures.size()));
    CUDA_CHECK(cudaMalloc(&d_texOffsets, texOffsets.size() * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_texWidths, texWidths.size() * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_texHeights, texHeights.size() * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_blockColors, blockColors.size() * sizeof(Color)));
    CUDA_CHECK(cudaMalloc(&d_blockPaletteIds, blockPaletteIds.size() * sizeof(int)));

    CUDA_CHECK(cudaMemcpy(d_vertices, vertices.data(), vertices.size() * sizeof(Vec3i), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_texCoords, obj.texCoords.data(), obj.texCoords.size() * sizeof(Vec2), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_faces, obj.faces.data(), obj.faces.size() * sizeof(Face), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_textures, allTextures.data(), allTextures.size(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_texOffsets, texOffsets.data(), texOffsets.size() * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_texWidths, texWidths.data(), texWidths.size() * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_texHeights, texHeights.data(), texHeights.size() * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_blockColors, blockColors.data(), blockColors.size() * sizeof(Color), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_blockPaletteIds, blockPaletteIds.data(), blockPaletteIds.size() * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_voxels, 0, voxelMemSize));
    CUDA_CHECK(cudaMemset(d_blockIds, 0, voxelSize * sizeof(int)));

    std::cout << "Rasterizing " << obj.faces.size() << " triangles with " << obj.materials.size() << " materials..." << std::endl;
    int threadsPerBlock = 256;
    int numBlocks = (obj.faces.size() + threadsPerBlock - 1) / threadsPerBlock;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    rasterizeFaces<<<numBlocks, threadsPerBlock>>>(
        d_vertices, d_texCoords, d_faces, obj.faces.size(),
        d_voxels, d_blockIds, w, h, l,
        d_textures, d_texOffsets, d_texWidths, d_texHeights,
        d_blockColors, d_blockPaletteIds, blockTextures.size()
    );

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "GPU took: " << milliseconds << " ms" << std::endl;

    CUDA_CHECK(cudaGetLastError());

    std::vector<unsigned int> h_voxels(voxelSize);
    std::vector<int> h_blockIds(voxelSize);
    CUDA_CHECK(cudaMemcpy(h_voxels.data(), d_voxels, voxelMemSize, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_blockIds.data(), d_blockIds, voxelSize * sizeof(int), cudaMemcpyDeviceToHost));

    int count = 0;
    for (auto v : h_voxels) if (v) count++;
    std::cout << "Generated " << count << " voxels" << std::endl;

    std::vector<std::vector<std::vector<int>>> blocks_data(
        w, std::vector<std::vector<int>>(h, std::vector<int>(l, 0)));

    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            for (int z = 0; z < l; z++) {
                int idx = x + z * w + y * w * l;
                blocks_data[x][y][z] = h_voxels[idx] ? h_blockIds[idx] : 0;
            }
        }
    }

    cudaFree(d_vertices);
    cudaFree(d_texCoords);
    cudaFree(d_faces);
    cudaFree(d_voxels);
    cudaFree(d_blockIds);
    cudaFree(d_textures);
    cudaFree(d_texOffsets);
    cudaFree(d_texWidths);
    cudaFree(d_texHeights);
    cudaFree(d_blockColors);
    cudaFree(d_blockPaletteIds);

    writeSchematic(output, blocks_data, w, h, l, blockTextures);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input.obj> <output.schem> [scale] [textures_dir]" << std::endl;
        return 1;
    }

    std::string input = argv[1];
    std::string output = argv[2];
    int target = argc > 3 ? std::stoi(argv[3]) : 100;
    std::string texturesDir = argc > 4 ? argv[4] : "textures";

    std::cout << "Loading block textures..." << std::endl;
    std::vector<BlockTexture> blockTextures = loadBlockTextures(texturesDir);

    if (blockTextures.empty()) {
        std::cerr << "No block textures found in " << texturesDir << std::endl;
        return 1;
    }

    std::cout << "Loading: " << input << std::endl;
    OBJ obj = loadOBJ(input);

    std::cout << "Converting with texture mapping..." << std::endl;
    convert(obj, output, target, blockTextures);

    std::cout << "Saved: " << output << std::endl;
    return 0;
}
