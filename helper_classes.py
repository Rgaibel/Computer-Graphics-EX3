import numpy as np


# This function gets a vector and returns its normalized form.
def normalize(vector):
    return vector / np.linalg.norm(vector)


# This function gets a vector and the normal of the surface it hit
# This function returns the vector that reflects from the surface
# lecture 5.10
def reflected(vector, normal):
    v = vector - 2 * np.dot(vector, normal) * normal
    return normalize(v)


def refracted(vector, normal, theta2):
    theta1 = np.arccos(np.dot(normal, normalize(vector)))
    n1 = np.sin(theta1)
    n2 = np.sin(theta2)
    medium_division = n1/n2
    return np.dot(medium_division*np.cos(theta1)-np.cos(theta2), normal) + medium_division*vector


# Lights


class LightSource:

    def __init__(self, intensity):
        self.intensity = intensity


class DirectionalLight(LightSource):

    def __init__(self, intensity, direction):
        super().__init__(intensity)
        self.direction = np.array(direction)

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection_point):
        return Ray(intersection_point, self.direction)

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self, intersection):
        return np.inf

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        return self.intensity


class PointLight(LightSource):

    def __init__(self, intensity, position, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.kc = kc
        self.kl = kl
        self.kq = kq

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection):
        return Ray(intersection, normalize(self.position - intersection))

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self, intersection):
        return np.linalg.norm(intersection - self.position)

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        d = self.get_distance_from_light(intersection)
        return self.intensity / (self.kc + self.kl*d + self.kq * (d**2))


class SpotLight(LightSource):

    def __init__(self, intensity, position, direction, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.direction = np.array(direction)
        self.kc = kc
        self.kl = kl
        self.kq = kq

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection):
        return Ray(intersection, normalize(self.direction))

    def get_distance_from_light(self, intersection):
        return np.linalg.norm(intersection - self.position)

        # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        cosAngle = np.dot(normalize(intersection - self.position), normalize(self.direction))
        d = self.get_distance_from_light(intersection)
        return (self.intensity * cosAngle) / (self.kc + self.kl * d + self.kq * (d ** 2))


class Ray:
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction

    # The function is getting the collection of objects in the scene and looks for the one with minimum distance.
    # The function should return the nearest object and its distance (in two different arguments)
    def nearest_intersected_object(self, objects):
        object_distances = [object.intersect(self) for object in objects]
        nearest_object = None
        min_distance = np.inf
        for i, distance in enumerate(object_distances):
            # print(distance)
            if (distance) and (distance[0] < min_distance):
                min_distance = distance[0]
                nearest_object = objects[i]

        return nearest_object, min_distance


class Object3D:

    def set_material(self, ambient, diffuse, specular, shininess, reflection):
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.reflection = reflection


class Plane(Object3D):
    def __init__(self, normal, point):
        self.normal = np.array(normal)
        self.point = np.array(point)

    def intersect(self, ray: Ray):
        v = self.point - ray.origin
        t = (np.dot(v, self.normal) / np.dot(self.normal, ray.direction))
        if t > 0:
            return t, self
        else:
            return None

    def compute_normal(self, intersection):
        return self.normal


class Triangle(Object3D):
    # Triangle gets 3 points as arguments
    def __init__(self, a, b, c):
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.ab = self.b - self.a
        self.ac = self.c - self.a
        self.bc = self.c - self.b
        self.ca = self.a - self.c
        self.normal = self.compute_normal(None)
        self.plane = Plane(self.normal, a)

    def compute_normal(self, intersection):
        n = np.cross(self.ab, self.ac)
        return n

    # Hint: First find the intersection on the plane
    # Later, find if the point is in the triangle using barycentric coordinates
    def intersect(self, ray: Ray):
        plane_intersect = self.plane.intersect(ray)
        if not plane_intersect:
            return None

        distance = plane_intersect[0]
        ray_plane_intersect = ray.origin + ray.direction * distance # (- .00005)?
        ia = ray_plane_intersect - self.a
        ib = ray_plane_intersect - self.b
        ic = ray_plane_intersect - self.c

        x = np.cross(self.ab, ia)
        y = np.cross(self.bc, ib)
        z = np.cross(self.ca, ic)

        if np.dot(x, self.normal) > 0 and np.dot(y, self.normal) > 0 and np.dot(z, self.normal) > 0:
            return distance, self
        else:
            return None


class Sphere(Object3D):
    def __init__(self, center, radius: float):
        self.center = center
        self.radius = radius

    def intersect(self, ray: Ray):
        b = 2 * np.dot(ray.direction, ray.origin - self.center)
        c = np.linalg.norm(ray.origin - self.center) ** 2 - self.radius ** 2
        delta = b ** 2 - 4*c
        if delta > 0:
            t1 = (-b + np.sqrt(delta)) / 2
            t2 = (-b - np.sqrt(delta)) / 2
            if t1 > 0 and t2 > 0:
                return min(t1, t2), self
        return None

    def compute_normal(self, intersection):
        return normalize(intersection - self.center)





class Mesh(Object3D):
    # Mesh are defined by a list of vertices, and a list of faces.
    # The faces are triplets of vertices by their index number.
    def __init__(self, v_list, f_list):
        self.v_list = v_list
        self.f_list = f_list
        self.triangle_list = self.create_triangle_list()

    def create_triangle_list(self):
        l = []
        # TODO
        return l

    def apply_materials_to_triangles(self):
        for t in self.triangle_list:
            t.set_material(self.ambient, self.diffuse, self.specular, self.shininess, self.reflection)

    # Hint: Intersect returns both distance and nearest object.
    # Keep track of both.
    def intersect(self, ray: Ray):
        #TODO
        pass

    def compute_normal(self, intersection):
        #TODO
        pass


class Mesh(Object3D):
    # Mesh are defined by a list of vertices, and a list of faces.
    # The faces are triplets of vertices by their index number.
    def __init__(self, v_list, f_list):
        self.v_list = v_list
        self.f_list = f_list
        self.triangle_list = self.create_triangle_list()

    def create_triangle_list(self):
        triangle_list = []
        for face in self.f_list:
            triangle = Triangle(self.v_list[face[0]], self.v_list[face[1]], self.v_list[face[2]])
            triangle_list.append(triangle)
        return triangle_list

    def apply_materials_to_triangles(self):
        for t in self.triangle_list:
            t.set_material(self.ambient, self.diffuse, self.specular, self.shininess, self.reflection)

    # Hint: Intersect returns both distance and nearest object.
    # Keep track of both.
    def intersect(self, ray: Ray):
        nearest_object = None
        min_distance = np.inf
        for i, triangle in enumerate(self.triangle_list):
            t = triangle.intersect(ray)
            if t and (t[0] < min_distance):
                min_distance = t[0]
                nearest_object = t[1]

        return nearest_object, min_distance

    def compute_normal(self, intersection):
        intersect_triangle = self.intersect(Ray(1e-5*intersection, normalize(intersection-1e-5*intersection )))
        return intersect_triangle[0].compute_normal(intersection)
