from helper_classes import *
import matplotlib.pyplot as plt


def render_scene(camera, ambient, lights, objects, screen_size, max_depth):
    width, height = screen_size
    ratio = float(width) / height
    screen = (-1, 1 / ratio, 1, -1 / ratio)  # left, top, right, bottom
    image = np.zeros((height, width, 3))

    for i, y in enumerate(np.linspace(screen[1], screen[3], height)):
        for j, x in enumerate(np.linspace(screen[0], screen[2], width)):
            pixel = np.array([x, y, 0])
            color = np.zeros(3)

            for light in lights:

                direction_vector = normalize(pixel - camera)
                ray = Ray(camera, direction_vector)
                reflection = 1

                for k in range(max_depth):

                    # check for intersections
                    nearest_object, min_distance = ray.nearest_intersected_object(objects)
                    if nearest_object is None:
                        break

                    intersection_point = ray.origin + min_distance * ray.direction
                    normal_to_object = nearest_object.compute_normal(intersection_point)
                    shifted_point = intersection_point + 1e-5 * normal_to_object
                    light_intensity = light.get_intensity(intersection_point)
                    intersection_to_light_direction_vector = light.get_light_ray(intersection_point).direction

                    _, min_distance = Ray(shifted_point, intersection_to_light_direction_vector).nearest_intersected_object(objects)
                    intersection_to_light_distance = light.get_distance_from_light(intersection_point)
                    is_shadowed = min_distance < intersection_to_light_distance

                    if is_shadowed:
                        break

                    phong_reflectance_model_equation = np.zeros(3)

                    # ambiant reflectance
                    # note I need to add something here, in the case that there are multiple lights, so this is only added once
                    if not color.any():
                        phong_reflectance_model_equation += np.array(nearest_object.ambient) * ambient
                    # diffuse reflectance
                    phong_reflectance_model_equation += light_intensity * nearest_object.diffuse * np.dot(normal_to_object, intersection_to_light_direction_vector)

                    # specular reflectance
                    l_hat = reflected(np.negative(intersection_to_light_direction_vector), normal_to_object)
                    v = np.negative(direction_vector)
                    phong_reflectance_model_equation += np.array(nearest_object.specular) * light_intensity * np.dot(v, l_hat) ** (nearest_object.shininess)

                    # reflection
                    color += reflection * phong_reflectance_model_equation
                    reflection *= nearest_object.reflection

                    ray = Ray(shifted_point, reflected(ray.direction, normal_to_object))
            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color, 0, 1)

    return image


# Write your own objects and lights
# TODO
def your_own_scene():
    camera = np.array([0,0,1])
    lights = []
    objects = []
    return camera, lights, objects

