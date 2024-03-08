import numpy as np
import matplotlib.pyplot as plt

class Circle():
  def __init__(self, center, radius):
    self.center = np.clongdouble(center)
    self.radius = np.clongdouble(radius)

  def contains_point(self, z):
    # on circle itself
    dist_from_center = np.absolute(z - self.center)

    if isclose(dist_from_center, self.radius):
      return True
    else:
      return False

  def contains_interior_point(self, z):

    dist_from_center = np.absolute(z - self.center)
    diff = dist_from_center - self.radius

    if isclose(diff, 0): # On circle itself
      return False
    elif diff > 0: # On exterior of circle
      return False
    else:
      return True


  def get_intersection(self, circle):
    # Based on http://paulbourke.net/geometry/circlesphere/

    z1 = self.center
    r1 = self.radius

    z2 = circle.center
    r2 = circle.radius

    d = np.absolute(self.center - circle.center) # Distance between centers

    if d > r1 + r2:
      return []

    if d < np.abs(r1-r2):
      return []

    if isclose(d,0):
      return []

    a = (r1**2 - r2**2 + d**2) / (2*d)
    h = np.sqrt(r1**2 - a**2)
    u = (z2 - z1)/d

    z3 = z1 + a*u

    int1 = z3 + h*1.0j*u
    int2 = z3 - h*1.0j*u

    return np.array([int1, int2], dtype=np.clongdouble)


  def intersection_angle_in_unit_circle(self, circle):

    intersection = self.get_intersection(circle)

    if len(intersection) == 0:
      return

    in_unit_circle = False

    for p in intersection:
      if np.sqrt(p[0]**2 + p[1]**2) < 1:
        in_unit_circle = True
        x_int = p[0]
        y_int = p[1]

    if not in_unit_circle:
      return

    t1 = np.arccos((x_int - self.center_x)/self.radius)
    t2 = np.arccos((x_int - circle.center_x)/circle.radius)

    arg_diff = get_arg(self.center_x, self.center_y) - get_arg(circle.center_x, circle.center_y)

    if (arg_diff > np.pi and arg_diff < 2*np.pi) or (arg_diff < 0 and arg_diff > -np.pi):
      # Positively oriented circles
      if y_int - self.center_y < 0:
        t1 = 2*np.pi - t1

      if y_int - circle.center_y > 0:
        t2 = 2*np.pi - t2

      angle = np.arccos(-np.cos(t1 + t2))

    else:

      if y_int - self.center_y > 0:
        t1 = 2*np.pi - t1

      if y_int - circle.center_y < 0:
        t2 = 2*np.pi - t2

      angle = 2*np.pi - np.arccos(-np.cos(t1 + t2))

    return angle


  def plot(self):
    ts = np.linspace(0, 2*pi, 100)
    xs = self.radius*np.cos(ts) + np.real(self.center)
    ys = self.radius*np.sin(ts) + np.imag(self.center)

    plt.plot(xs, ys)


class Segment():
  def __init__(self, z1, z2):

    self.z1 = np.clongdouble(z1)
    self.z2 = np.clongdouble(z2)

    self.x1 = np.real(z1)
    self.y1 = np.imag(z1)
    self.x2 = np.real(z2)
    self.y2 = np.imag(z2)

    self.min_x = self.x1 if self.x1 < self.x2 else self.x2
    self.max_x = self.x2 if self.x1 < self.x2 else self.x1

    self.min_y = self.y1 if self.y1 < self.y2 else self.y2
    self.max_y = self.y2 if self.y1 < self.y2 else self.y1

    self.slope = (self.y2 - self.y1) / (self.x2 - self.x1)


  def contains_point(self, z):
    z = np.clongdouble(z)

    x = np.real(z)
    y = np.imag(z)

    if x < self.min_x or x > self.max_x:
      return False

    if y < self.min_y or y > self.max_y:
      return False

    if isclose(self.slope*(x-self.x1) + self.y1, y):
      return True
    else:
      return False

  def get_intersection(self, circle):
    # Intersection with a circle

    x1 = self.x1
    x2 = self.x2
    x3 = np.real(circle.center)

    y1 = self.y1
    y2 = self.y2
    y3 = np.imag(circle.center)

    r = circle.radius

    a = (x2 - x1)**2 + (y2-y1)**2
    b = 2*((x2-x1)*(x1-x3) + (y2-y1)*(y1-y3))
    c = x3**2 + y3**2 + x1**2 + y1**2 - 2*(x3*x1 + y3*y1) - r**2

    disc = b**2 - 4*a*c
    u1 = (-b + np.sqrt(disc))/(2*a)
    u2 = (-b - np.sqrt(disc))/(2*a)

    if isclose(disc, 0):

      x_int = x1 + u1*(x2-x1)
      y_int = y1 + u1*(y2-y1)

      if self.contains_point(x_int + 1j*y_int):
        return (x_int, y_int)
      else:
        return []

    elif disc < 0:
      return []


    else:

      x_int1 = x1 + u1*(x2-x1)
      y_int1 = y1 + u1*(y2-y1)

      x_int2 = x1 + u2*(x2-x1)
      y_int2 = y1 + u2*(y2-y1)

      points = []

      if self.contains_point(x_int1 + 1j*y_int1):
        points.append(x_int1 + 1j*y_int1)

      if self.contains_point(x_int2 + 1j*y_int2):
        points.append(x_int2 + 1j*y_int2)

      if len(points) == 0:
        return []
      else:
        return np.array(points, dtype = np.clongdouble)

  def intersection_angle_in_unit_circle(self, circle):

    center_x = np.real(circle.center)
    center_y = np.imag(circle.center)

    intersection = self.get_intersection(circle)

    if len(intersection) == 0:
      return

    in_unit_circle = False


    for p in intersection:

      if np.absolute(p) < 1:
        in_unit_circle = True
        x_int = np.real(p)
        y_int = np.imag(p)

    if not in_unit_circle:
      return

    t = np.arccos((x_int - center_x)/circle.radius)

    if (self.x1 - center_x)**2 + (self.y1 - center_y)**2 > circle.radius**2:
      # Orient circle clockwise

      if y_int - center_y > 0:
        t = 2*np.pi - t

      angle = np.arccos(np.sin(t))

    else:

      if y_int - center_y < 0:
        t = 2*np.pi - t

      angle = 2*np.pi - np.arccos(np.sin(t))

    return angle

  def plot(self):
    plt.plot([self.x1, self.x2], [self.y1, self.y2])