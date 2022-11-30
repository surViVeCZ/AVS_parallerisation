/**
 * @file    tree_mesh_builder.h
 *
 * @author  Petr Pouč <xpoucp01@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    30.11.2022
 **/

#ifndef TREE_MESH_BUILDER_H
#define TREE_MESH_BUILDER_H
#include "base_mesh_builder.h"

class TreeMeshBuilder : public BaseMeshBuilder
{
public:
    TreeMeshBuilder(unsigned gridEdgeSize);

protected:
    unsigned MIN_GRID = 1;
    unsigned marchCubes(const ParametricScalarField &field);
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data();  }
     std::vector<Triangle_t> mTriangles; ///< Temporary array of triangles

private:
    unsigned decompose(const ParametricScalarField &field, const Vec3_t<float> &offset, unsigned gridSize);
};

#endif // TREE_MESH_BUILDER_H
