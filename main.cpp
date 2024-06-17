#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <cmath>

using namespace std;
using namespace parser;


typedef unsigned char RGB[3];

typedef struct Ray
{
    Vec3f origin;
    Vec3f direction;
    bool shadowRay;
} ray;


enum ObjectType
{
	SPHERE, TRIANGLE, MESH
};


typedef struct Hit
{
	bool hitHappened;
	Vec3f intersectionPoint;
	Vec3f surfaceNormal;
	int materialId;
	float t;
	ObjectType objectType;
	int objectId;
} hit;


Vec3f subtract(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.x-b.x;
	result.y = a.y-b.y;
	result.z = a.z-b.z;
	return result;
}

Vec3f add(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.x+b.x;
	result.y = a.y+b.y;
	result.z = a.z+b.z;
	return result;
}

float dotProduct(const Vec3f &a, const Vec3f &b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    Vec3f result;
    result.x = a.y*b.z-a.z*b.y;
    result.y = a.z*b.x-a.x*b.z;
    result.z = a.x*b.y-a.y*b.x;

    return result;
}

float findVectorLength(const Vec3f &a)
{
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

float findDistance(const Vec3f &a, const Vec3f &b)
{
	return sqrt(pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow(a.z-b.z,2));
}

Vec3f normalizeVector(const Vec3f &a)
{
	Vec3f result;
	result.x = a.x / findVectorLength(a);
	result.y = a.y / findVectorLength(a);
	result.z = a.z / findVectorLength(a);
	return result;
}

Vec3f findIntersectionPoint(const Ray &ray, float t)
{
	Vec3f result;
	result.x = ray.origin.x + t*ray.direction.x;
	result.y = ray.origin.y + t*ray.direction.y;
	result.z = ray.origin.z + t*ray.direction.z;

	return result;
}

/*
Hit findHitPoint(vector<Hit> &hitDataVector)
{
	Hit closestHitPoint;
	closestHitPoint.hitHappened = false;

	if (!hitDataVector.empty())
	{
		closestHitPoint = hitDataVector[0];

		for (int i = 1; i < hitDataVector.size(); i++)
		{
			if (hitDataVector[i].p < closestHitPoint.p)
			{
				closestHitPoint = hitDataVector[i];
			}
		}
		closestHitPoint.hitHappened = true;
	}

	return closestHitPoint;
}

*/

Hit findHit(vector<Hit> &hitDataVector)
{
	Hit result;
	result.hitHappened = false;

	if(hitDataVector.size() != 0)
	{
		result = hitDataVector[0];

		for(int i=1;i<hitDataVector.size();i++)
		{
			if(hitDataVector[i].t<result.t)
			{
				result = hitDataVector[i];
			}
		}
    	result.hitHappened = true;
	}
	return result;
}


/*
Ray generateRay(const Camera &camera, int i, int j)
{
    float left = camera.near_plane.x;
    float right = camera.near_plane.y;
    float bottom = camera.near_plane.z;
    float top = camera.near_plane.w;

    Vec3f gaze = normalizeVector(camera.gaze);

    float su = (right - left) * (j + 0.5) / camera.image_width;
    float sv = (top - bottom) * (i + 0.5) / camera.image_height;

    Vec3f m = camera.position + gaze * camera.near_distance;
    Vec3f u = normalizeVector(crossProduct(gaze, camera.up));
    Vec3f v = normalizeVector(crossProduct(u, gaze));
    Vec3f q = m + (u * left) + (v * top);
    Vec3f s = q + (u * su) - (v * sv);

    Ray ray;
    ray.origin = camera.position;
    ray.direction = normalizeVector(s - camera.position);
    ray.shadowRay = false;

    return ray;
}

*/

Ray generateRay(const Camera &camera,int screenWidth,int screenHeight)
{
    float left = camera.near_plane.x;
    float right = camera.near_plane.y;
    float bottom = camera.near_plane.z;
    float top = camera.near_plane.w;

    Vec3f gazeVector = normalizeVector(camera.gaze);

    float su = (right-left)*(screenHeight+0.5)/camera.image_width;
    float sv = (top-bottom)*(screenWidth+0.5)/camera.image_height;

    Vec3f m, q, u, v;
    m.x = camera.position.x + (gazeVector.x * camera.near_distance);
    m.y = camera.position.y + (gazeVector.y * camera.near_distance);
    m.z = camera.position.z + (gazeVector.z * camera.near_distance);
    u = crossProduct(gazeVector, camera.up);
    u = normalizeVector(u);
    v = crossProduct(u,gazeVector);
    q.x = m.x + (u.x*left) + (v.x*top);
    q.y = m.y + (u.y*left) + (v.y*top);
    q.z = m.z + (u.z*left) + (v.z*top);
    Vec3f s;
    s.x = q.x + (u.x*su) - (v.x * sv);
    s.y = q.y + (u.y*su) - (v.y * sv);
    s.z = q.z + (u.z*su) - (v.z * sv);
    Ray ray;
    ray.origin = camera.position;
    ray.direction = subtract(s,camera.position);
    ray.direction = normalizeVector(ray.direction);
    ray.shadowRay = false;

    return ray;
}

/*
Hit sphereIntersection(const Ray &ray, const Vec3f &center, float radius, int material_id, int obj_id)
{
    Hit hit;
    hit.material_id = material_id;
    hit.objType = SPHERE;
    hit.obj_id = obj_id;

    Vec3f oc = subtract(ray.origin, center);
    float a = dotProduct(ray.direction, ray.direction);
    float b = 2.0f * dotProduct(oc, ray.direction);
    float c = dotProduct(oc, oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {
        float temp = (-b - sqrtf(discriminant)) / (2.0f * a);
        if (temp > 0) {
            hit.t = temp;
            hit.hitHappened = true;
            hit.intersectionPoint = findIntersectionPoint(ray, temp);
            hit.surfaceNormal = subtract(hit.intersectionPoint, center) / radius;
            return hit;
        }
    }

    hit.hitHappened = false;
    return hit;
}
*/

Hit sphereIntersection(const Ray &ray, const Vec3f &sphereCenter, float sphereRadius, int materialId, int objectId)
{
	Hit hitResult;

	const float a = dotProduct(ray.direction,ray.direction);
	Vec3f originToCenter = subtract(ray.origin,sphereCenter);
	const float b = 2 * dotProduct(ray.direction,originToCenter);
	const float c = dotProduct(originToCenter,originToCenter)-sphereRadius*sphereRadius;

	const float discriminant = b*b - 4*a*c;

	if(discriminant < 0)			
	{
		hitResult.hitHappened = false;
	}
	else							 
	{	
		const float t1 = (-1 * b + sqrtf(discriminant))/2*a;
		const float t2 = (-1 * b - sqrtf(discriminant))/2*a;

		hitResult.materialId = materialId;
		hitResult.hitHappened = true;
		hitResult.objectType = SPHERE;
		hitResult.objectId = objectId;

		const float t = fmin(t1, t2);
		hitResult.intersectionPoint = findIntersectionPoint(ray, t);	
		hitResult.surfaceNormal = subtract(hitResult.intersectionPoint,sphereCenter);	
		hitResult.surfaceNormal.x /= sphereRadius;
		hitResult.surfaceNormal.y /= sphereRadius;
		hitResult.surfaceNormal.z /= sphereRadius;

		hitResult.t = t;
	}
	return hitResult;
}


float determinant(const Vec3f &v0,const Vec3f &v1,const Vec3f &v2)
{
	return v0.x*(v1.y*v2.z-v2.y*v1.z)+v0.y*(v2.x*v1.z-v1.x*v2.z)+v0.z*(v1.x*v2.y-v1.y*v2.x);
}

/*
Hit triangleIntersection(const Ray &ray, const Vec3f &vertexA, const Vec3f &vertexB, const Vec3f &vertexC, int material_id, int obj_id)
{
    Hit hit;
    hit.hitHappened = false;

    Vec3f rayOrigin = ray.origin;
    Vec3f rayDirection = ray.direction;

    Vec3f edgeAB = subtract(vertexA, vertexB);
    Vec3f edgeAC = subtract(vertexA, vertexC);

    Vec3f h = crossProduct(rayDirection, edgeAC);
    float detA = dotProduct(edgeAB, h);

    if (detA == 0.0)
    {
        return hit;
    }

    float invDetA = 1.0f / detA;

    Vec3f s = subtract(rayOrigin, vertexA);
    float u = dotProduct(s, h) * invDetA;

    if (u < 0.0 || u > 1.0)
    {
        return hit;
    }

    Vec3f q = crossProduct(s, edgeAB);
    float v = dotProduct(rayDirection, q) * invDetA;

    if (v < 0.0 || (u + v) > 1.0)
    {
        return hit;
    }

    float t = dotProduct(edgeAC, q) * invDetA;

    if (t <= 0.0)
    {
        return hit;
    }

    hit.hitHappened = true;
    hit.objType = TRIANGLE;
    hit.obj_id = obj_id;
    hit.material_id = material_id;
    hit.t = t;
    hit.intersectionPoint = findIntersectionPoint(ray, t);
    hit.surfaceNormal = normalizeVector(crossProduct(edgeAB, edgeAC));

    return hit;
}
*/

Hit triangleIntersection(const Ray &ray, const Vec3f &vertexA, const Vec3f &vertexB, const Vec3f &vertexC, int materialId, int objectId)
{
	Hit hitResult;
	hitResult.hitHappened = false;

	Vec3f o = ray.origin;
	Vec3f direction = ray.direction;

	Vec3f edge1 = subtract(vertexA,vertexB);
	Vec3f edge2 = subtract(vertexA,vertexC);
	Vec3f aToOrigin = subtract(vertexA,o);

	float detA = determinant(edge1,edge2,direction);
	if(detA == 0.0)
	{
		return hitResult;
	}

	float t = (determinant(edge1, edge2, aToOrigin))/detA;
	if(t <= 0.0) {
		return hitResult;
	}

	float gamma = (determinant(edge1,aToOrigin, direction))/detA;
	if(gamma < 0 || gamma > 1) {
		return hitResult;
	}

	float beta = (determinant(aToOrigin, edge2, direction))/detA;
	if(beta < 0 || beta > (1 - gamma)) {
		return hitResult;
	}

	hitResult.hitHappened = true;
	hitResult.objectType = TRIANGLE;
	hitResult.objectId = objectId;
	hitResult.materialId = materialId;
	hitResult.t = t;
	hitResult.intersectionPoint = findIntersectionPoint(ray, t);
	hitResult.surfaceNormal = crossProduct(subtract(vertexB, vertexA), subtract(vertexC, vertexA));
	hitResult.surfaceNormal = normalizeVector(hitResult.surfaceNormal);

	return hitResult;
}

/*Hit meshIntersection(const Ray &ray, const Mesh &mesh, const Scene &scene, int material_id, int obj_id)
{
    Hit hit;
    hit.hitHappened = false;
    vector<Hit> hitInfoVector;

    for (int faceNumber = 0; faceNumber < mesh.faces.size(); faceNumber++)
    {
        Vec3f v0 = scene.vertex_data[mesh.faces[faceNumber].v0_id - 1];
        Vec3f v1 = scene.vertex_data[mesh.faces[faceNumber].v1_id - 1];
        Vec3f v2 = scene.vertex_data[mesh.faces[faceNumber].v2_id - 1];

        Hit triangleHit = triangleIntersection(ray, v0, v1, v2, material_id, obj_id);

        if (triangleHit.hitHappened && triangleHit.t >= 0)
        {
            triangleHit.material_id = material_id;
            triangleHit.objType = MESH;
            triangleHit.obj_id = obj_id;
            triangleHit.intersectionPoint = findIntersectionPoint(ray, triangleHit.t);
            triangleHit.surfaceNormal = crossProduct(subtract(v1, v0), subtract(v2, v0));
            triangleHit.surfaceNormal = normalizeVector(triangleHit.surfaceNormal);

            hitInfoVector.push_back(triangleHit);
        }
    }

    hit = findClosestHit(hitInfoVector); // Değişen isim, "findHit" yerine "findClosestHit"
    return hit;
}
*/

Hit meshIntersection(const Ray &ray, const Mesh &mesh, const Scene &scene, int materialId, int objectId)
{
	Hit hit;
	hit.hitHappened = false;
	vector<Hit> hitDataVector;

	for(int faceIndex = 0; faceIndex < mesh.faces.size(); faceIndex++)
	{
		Vec3f v0 = scene.vertex_data[mesh.faces[faceIndex].v0_id - 1];		
		Vec3f v1 = scene.vertex_data[mesh.faces[faceIndex].v1_id - 1];
		Vec3f v2 = scene.vertex_data[mesh.faces[faceIndex].v2_id - 1];

		hit = triangleIntersection(ray, v0, v1, v2, materialId, objectId);
		if(hit.hitHappened && hit.t >= 0)
		{
			hit.materialId = materialId;
			hit.objectType = MESH;
			hit.objectId = objectId;
			hit.intersectionPoint = findIntersectionPoint(ray, hit.t);
			hit.surfaceNormal = crossProduct(subtract(v1, v0), subtract(v2, v0));
            hit.surfaceNormal = normalizeVector(hit.surfaceNormal);

            hitDataVector.push_back(hit);
		}
	}

	hit = findHit(hitDataVector);
	return hit;
}

/*
Vec3f findIrradiance(const PointLight &pointLight, const Vec3f &intersectionPoint)
{
    Vec3f irradiance;
    Vec3f d = subtract(pointLight.position, intersectionPoint);
    float d_square = dotProduct(d, d);

    if (d_square != 0.0)
    {
        float inv_d_square = 1.0f / d_square;
        irradiance.x = pointLight.intensity.x * inv_d_square;
        irradiance.y = pointLight.intensity.y * inv_d_square;
        irradiance.z = pointLight.intensity.z * inv_d_square;
    }
    
    return irradiance;
}
*/

Vec3f findIrradiance(const PointLight &pointLight, const Vec3f &intersectionPoint)
{
    Vec3f irradiance;
    Vec3f lightDirection = subtract(pointLight.position,intersectionPoint);
    float squareOfDistance = dotProduct(lightDirection,lightDirection);

    if(squareOfDistance != 0.0)
    {
	    irradiance.x = pointLight.intensity.x/squareOfDistance;
	    irradiance.y = pointLight.intensity.y/squareOfDistance;
	    irradiance.z = pointLight.intensity.z/squareOfDistance;
    }
    return irradiance;
}

/*
Vec3f findDiffuse(const PointLight &currentLight, const Scene &scene, int material_id, const Vec3f &normal, const Vec3f &intersectionPoint)
{
    Vec3f diffuse;
    const Material &material = scene.materials[material_id - 1];

    Vec3f irradiance = findIrradiance(currentLight, intersectionPoint);
    Vec3f l = normalizeVector(subtract(currentLight.position, intersectionPoint));

    float dotPro = dotProduct(l, normal);
    dotPro = std::max(dotPro, 0.0f);  // Ensure dotPro is non-negative

    diffuse.x = material.diffuse.x * dotPro * irradiance.x;
    diffuse.y = material.diffuse.y * dotPro * irradiance.y;
    diffuse.z = material.diffuse.z * dotPro * irradiance.z;

    return diffuse;
}
*/

const Vec3f findDiffuse(const PointLight &light, const Scene &scene, int materialId, const Vec3f &surfaceNormal, const Vec3f &intersectionPoint)
{
	Vec3f diffuse;

	Vec3f irradiance = findIrradiance(light,intersectionPoint);

	Vec3f lightDirection = subtract(light.position,intersectionPoint);
	lightDirection = normalizeVector(lightDirection);


	float dotproduct = dotProduct(lightDirection, surfaceNormal);
	if(dotproduct < 0)
	{
		dotproduct = 0;
	}

	diffuse.x = scene.materials[materialId-1].diffuse.x*dotproduct*irradiance.x;
	diffuse.y = scene.materials[materialId-1].diffuse.y*dotproduct*irradiance.y;
	diffuse.z = scene.materials[materialId-1].diffuse.z*dotproduct*irradiance.z;

	return diffuse;
}

/*
Vec3f calculateSpecular(const PointLight &light, const Scene &scene, const Ray &ray, int material_id, const Vec3f &normal, const Vec3f &intersectionPoint)
{
    Vec3f specular;

    Material material = scene.materials[material_id - 1];

    Vec3f irradiance = calculateIrradiance(light, intersectionPoint);

    Vec3f lightDirection = subtract(light.position, intersectionPoint);
    lightDirection = normalizeVector(lightDirection);

    Vec3f halfVector = subtract(lightDirection, ray.direction);
    halfVector = normalizeVector(halfVector);

    float dotProduct = dotProduct(normal, halfVector);
    if (dotProduct < 0)
    {
        dotProduct = 0;
    }

    specular.x = material.specular.x * pow(dotProduct, material.phong_exponent) * irradiance.x;
    specular.y = material.specular.y * pow(dotProduct, material.phong_exponent) * irradiance.y;
    specular.z = material.specular.z * pow(dotProduct, material.phong_exponent) * irradiance.z;

    return specular;
}
*/////////////

Vec3f findSpecular(const PointLight &currentLight,const Scene &scene, const Ray &ray,int material_id,const Vec3f &normal,const Vec3f &intersectionPoint)
{
	Vec3f specular;

	Material material = scene.materials[material_id - 1];

	Vec3f irradiance = findIrradiance(currentLight,intersectionPoint);

	Vec3f lightDirection = subtract(currentLight.position, intersectionPoint);
	lightDirection = normalizeVector(lightDirection);

	Vec3f h = subtract(lightDirection, ray.direction);
	h = normalizeVector(h);

	float dotproduct = dotProduct(normal, h);
	if(dotproduct < 0)
	{
		dotproduct = 0;
	}

	specular.x = material.specular.x * pow(dotproduct, material.phong_exponent) * irradiance.x;
	specular.y = material.specular.y * pow(dotproduct, material.phong_exponent) * irradiance.y;
	specular.z = material.specular.z * pow(dotproduct, material.phong_exponent) * irradiance.z;

	return specular;
}

bool isMirror(const Scene &scene, int materialId)
{
    const Material &material = scene.materials[materialId - 1];
    return (material.mirror.x>0||material.mirror.y > 0 || material.mirror.z>0);
}



Hit findHitResult(const Scene &scene, const Ray &ray)
{
    vector<Hit> hitDataVector;

    for (const Sphere &sphere : scene.spheres)
    {
        Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
        Hit hit = sphereIntersection(ray, center, sphere.radius, sphere.material_id, 0);
        if (hit.hitHappened && hit.t >= 0)
        {
            hitDataVector.push_back(hit);
        }
    }

    for (const Triangle &triangle : scene.triangles)
    {
        Vec3f vertex0 = scene.vertex_data[triangle.indices.v0_id - 1];
        Vec3f vertex1 = scene.vertex_data[triangle.indices.v1_id - 1];
        Vec3f vertex2 = scene.vertex_data[triangle.indices.v2_id - 1];

        Hit hit = triangleIntersection(ray, vertex0, vertex1, vertex2, triangle.material_id, 0);

        if (hit.hitHappened && hit.t >= 0)
        {
            hitDataVector.push_back(hit);
        }
    }

    for (int meshNumber = 0; meshNumber < scene.meshes.size(); meshNumber++)
    {
        const Mesh &currentMesh = scene.meshes[meshNumber];
        Hit hit = meshIntersection(ray, currentMesh, scene, currentMesh.material_id, meshNumber);

        if (hit.hitHappened && hit.t >= 0)
        {
            hitDataVector.push_back(hit);
        }
    }

    Hit hitResult = findHit(hitDataVector);

    return hitResult;
}



Vec3f findPixelColor(const Scene &scene, const Hit &hitResult, const Camera &currentCamera, const Ray &ray, int maxDepth)
{
    int numberOfLights = scene.point_lights.size();
    int nSpheres = scene.spheres.size();
    int nTriangles = scene.triangles.size();
    int nMeshes = scene.meshes.size();

    float pixel1 = 0;
    float pixel2 = 0;
    float pixel3 = 0;

    Vec3f pixelColor;

    if(hitResult.hitHappened)
    {
    	int materialId = hitResult.materialId;

		pixel1 = scene.materials[materialId - 1].ambient.x * scene.ambient_light.x;
		pixel2 = scene.materials[materialId - 1].ambient.y * scene.ambient_light.y;
		pixel3 = scene.materials[materialId - 1].ambient.z * scene.ambient_light.z;
      	
		for(int lightNumber = 0; lightNumber < numberOfLights; lightNumber++)
		{
		    bool shadowFlag = false;

		    PointLight currentLight = scene.point_lights[lightNumber];
		    float lightToCam = findDistance(currentLight.position, currentCamera.position);

			Vec3f lightDirection = subtract(currentLight.position,hitResult.intersectionPoint);
			lightDirection = normalizeVector(lightDirection);


			Vec3f lightDirectionEpsilon;
			lightDirectionEpsilon.x = lightDirection.x * scene.shadow_ray_epsilon;
			lightDirectionEpsilon.y = lightDirection.y * scene.shadow_ray_epsilon;
			lightDirectionEpsilon.z = lightDirection.z * scene.shadow_ray_epsilon;

			Ray shadowRay = { add(hitResult.intersectionPoint, lightDirectionEpsilon), lightDirection, true };

			Hit shadowHit;
			vector<Hit> shadowHitVector;

			float tLight = subtract(currentLight.position, shadowRay.origin).x / shadowRay.direction.x;

        	for(int sphereNumber = 0; sphereNumber < nSpheres; sphereNumber++)
        	{
          		Sphere currentSphere = scene.spheres[sphereNumber];
          		Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];
          		float radius = currentSphere.radius;

          		shadowHit = sphereIntersection(shadowRay, center, radius, currentSphere.material_id, sphereNumber);

          		if(shadowHit.hitHappened)
          		{
            		if(tLight > shadowHit.t && shadowHit.t >= 0)
            		{
              			shadowFlag = true;
            		}
          		}
        	}


	        for(int triangleNumber = 0; triangleNumber < nTriangles; triangleNumber++)
	        {
	          	Triangle currentTriangle = scene.triangles[triangleNumber];
	          	Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
	          	Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
	          	Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];

	          	shadowHit = triangleIntersection(shadowRay, v0, v1, v2, currentTriangle.material_id, triangleNumber);

	          	if(shadowHit.hitHappened)
	          	{
	            	if(tLight > shadowHit.t && shadowHit.t >= 0)
	            	{
	              		shadowFlag = true;
	            	}
	          	}
	        }



	        if(!shadowFlag)
	        {
	          	for(int meshNumber = 0; meshNumber < nMeshes; meshNumber++)
	          	{
	            	Mesh currentMesh = scene.meshes[meshNumber];

	          		shadowHit = meshIntersection(shadowRay, currentMesh, scene, currentMesh.material_id, meshNumber);

	          		if(shadowHit.hitHappened)
	            	{
	              		if(tLight > shadowHit.t && shadowHit.t >= 0)
	              		{
	                		shadowFlag = true;
	              		}
	            	}
	        	}
        	}


	        if(!shadowFlag || (shadowFlag && lightToCam == 0))
	        {
		        int materialId = hitResult.materialId;

		        Vec3f diffuse = findDiffuse(currentLight, scene, materialId, hitResult.surfaceNormal, hitResult.intersectionPoint);

		        Vec3f specular = findSpecular(currentLight, scene, ray, materialId, hitResult.surfaceNormal, hitResult.intersectionPoint);
		                    
              	pixel1 += diffuse.x + specular.x;
  		        pixel2 += diffuse.y + specular.y;
  		        pixel3 += diffuse.z + specular.z;
	        }    
    	}

        bool mirrorness = isMirror(scene,materialId);
        Vec3f reflection;
        reflection.x = 0;
        reflection.y = 0;
        reflection.z = 0;

        if(maxDepth > 0 && mirrorness)
        {
        	float lightDirection = -2 * dotProduct(ray.direction, hitResult.surfaceNormal);
        	Vec3f normal_lightDirection;
        	normal_lightDirection.x = hitResult.surfaceNormal.x * lightDirection + ray.direction.x;
        	normal_lightDirection.y = hitResult.surfaceNormal.y * lightDirection + ray.direction.y;
        	normal_lightDirection.z = hitResult.surfaceNormal.z * lightDirection + ray.direction.z;

      		normal_lightDirection = normalizeVector(normal_lightDirection);

      		Vec3f lightDirectionEpsilon;
        	lightDirectionEpsilon.x = normal_lightDirection.x * scene.shadow_ray_epsilon;
        	lightDirectionEpsilon.y = normal_lightDirection.y * scene.shadow_ray_epsilon;
        	lightDirectionEpsilon.z = normal_lightDirection.z * scene.shadow_ray_epsilon;

        	
            Ray reflectionRay = { add(hitResult.intersectionPoint, lightDirectionEpsilon), normal_lightDirection, false};

            Hit hitresultt = findHitResult(scene, reflectionRay);


            if(!(hitresultt.objectId == hitResult.objectId && hitresultt.objectType == hitResult.objectType))
            {
                reflection = findPixelColor(scene, hitresultt, currentCamera, reflectionRay, (maxDepth-1));
            }

      		pixel1+=reflection.x * scene.materials[materialId - 1].mirror.x;
      		pixel2+=reflection.y * scene.materials[materialId - 1].mirror.y;
      		pixel3+=reflection.z * scene.materials[materialId - 1].mirror.z;

   		}
  	}
  	else        
 	{
      	pixel1 = scene.background_color.x;
      	pixel2 = scene.background_color.y;
      	pixel3 = scene.background_color.z;
  	}

  	pixelColor.x = pixel1;
  	pixelColor.y = pixel2;
  	pixelColor.z = pixel3;

  	return pixelColor;
}

int main(int argc, char* argv[])
{
    Scene scene;

    scene.loadFromXml(argv[1]);

    int nCameras = scene.cameras.size();

    for(int cameraN = 0;cameraN<nCameras;cameraN++)
    {
        Camera currentCamera = scene.cameras[cameraN];
        int width = currentCamera.image_width;
        int height = currentCamera.image_height;
        int numberOfLights = scene.point_lights.size();

        unsigned char* image = new unsigned char [width*height*3];
        int pixelNumber = 0;

        for(int i=0;i<height;i++)
        {
            for(int j=0;j<width;j++)
            {
                Ray ray = generateRay(currentCamera,i,j);

                Hit hitResult = findHitResult(scene,ray);


             	Vec3f pixelColor = findPixelColor(scene,hitResult,currentCamera,ray,(scene.max_recursion_depth));

            	if(pixelColor.x > 255)
            		image[pixelNumber] = 255;
              	else
                	image[pixelNumber] = round(pixelColor.x);

            	if(pixelColor.y > 255)
            		image[pixelNumber + 1] = 255;
             	else
                	image[pixelNumber + 1] = round(pixelColor.y);

            	if(pixelColor.z > 255)
            		image[pixelNumber + 2] = 255;
              	else
                	image[pixelNumber + 2] = round(pixelColor.z);

              	pixelNumber += 3;

            }
        }

        write_ppm(currentCamera.image_name.c_str(),image,width,height);
    }
}
