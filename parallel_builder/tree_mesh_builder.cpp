/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  Petr Pouč <xpoucp01@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    30.11.2022
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{   
    size_t totalCubesCount = mGridSize*mGridSize*mGridSize;
    unsigned totalTriangles = 0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            totalTriangles = triangular_cutter(field, Vec3_t<float>(), mGridSize);
        }
    }
    return totalTriangles;
}

unsigned TreeMeshBuilder::triangular_cutter(const ParametricScalarField &field, const Vec3_t<float> &offset, unsigned gridSize)
{
    unsigned totalTriangles = 0;
    //cut actual grid to half
    const auto grid_resized = gridSize / 2;
	
    //too small grid - isosurface is not in this grid
	if (isBlockEmpty(field, gridSize, offset))
	{
        //computing isosurface value in the middle of the grid
        const Vec3_t<float> midpoint(
            offset.x * mGridResolution + (gridSize * mGridResolution) / 2,
            offset.y * mGridResolution + (gridSize * mGridResolution) / 2,
            offset.z * mGridResolution + (gridSize * mGridResolution) / 2
        );
        const float midpoint_value = evaluateFieldAt(midpoint, field);
        const float empty_block = mIsoLevel + (sqrtf(3.F) / 2.F) * gridSize * mGridResolution;
        if (midpoint_value > empty_block){
            return 0;
        }
	}

    //grid is smaller than minimal grid size, computing triangles
	if (gridSize <= MIN_GRID)
	{
		return buildCube(offset, field);
	}

    //cutting grid into 8 parts
	for (const Vec3_t<float> vertex_pos : sc_vertexNormPos)
	{
        //computing offset for each children
        const Vec3_t<float> newoffset(offset.x + vertex_pos.x * grid_resized, offset.y + vertex_pos.y * grid_resized, offset.z + vertex_pos.z * grid_resized);
        //computing triangles for each children
        const unsigned trianglesCount = triangular_cutter(field, newoffset, grid_resized);
        totalTriangles += trianglesCount;
	}

	return totalTriangles;
}

//zde není potřeba nic upravovat
float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

   
   //#pragma omp parallel for schedule(static) reduction(min:value)
    for(unsigned i = 0; i < count; ++i)

    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        value = std::min(value, distanceSquared);
    }
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    #pragma omp critical(tree_emitTriangle)
    mTriangles.push_back(triangle);
}
