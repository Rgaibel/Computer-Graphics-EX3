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
            # This is the main loop where each pixel color is computed.
            # TODO
            # ray = constructRayThroughPixel(camera, pixel)

            origin = camera
            direction = normalize(pixel - origin)

            color = np.zeros(3)
            reflection = 1

            for k in range(max_depth):
                # check for intersections
                nearest_object, min_distance = Ray(origin, direction).nearest_intersected_object(objects)
                if nearest_object is None:
                    break

                intersection = origin + min_distance * direction

                normal_to_object = nearest_object.compute_normal(intersection)
                shifted_point = intersection + 1e-5 * normal_to_object
                for light in lights:
                    intersection_to_light = normalize(light.position - shifted_point)

                    _, min_distance = Ray(shifted_point, intersection_to_light).nearest_intersected_object(objects)
                    intersection_to_light_distance = np.linalg.norm(light.position - intersection)
                    is_shadowed = min_distance < intersection_to_light_distance

                    if is_shadowed:
                        break

                    illumination = np.zeros(3)

                    # ambiant
                    illumination += np.array(nearest_object.ambient) * light.kc

                    # diffuse
                    illumination += np.array(nearest_object.diffuse) * light.kl * np.dot(intersection_to_light, normal_to_object)

                    # specular
                    intersection_to_camera = normalize(camera - intersection)
                    H = normalize(intersection_to_light + intersection_to_camera)
                    illumination += np.array(nearest_object.specular) * light.kq * np.dot(normal_to_object, H) ** (nearest_object.shininess / 4)

                    # reflection
                    color += reflection * illumination
                    reflection *= nearest_object.reflection

                    origin = shifted_point
                    direction = reflected(direction, normal_to_object)
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

