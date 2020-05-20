"""
    Geometry routines for phylo_tree plugin

    Isaac Stead 2020
"""
class Line(object):
    """Quick class to represent a line for geometry operations"""
    
    def __init__(self, point_a, point_b):
        self.a = point_a
        self.b = point_b

    @property
    def slope(self):
        """The line's slope"""
        xa, ya = self.a
        xb, yb = self.b
        return (yb - ya) / (xb - xa)

    def slope_xy(self):
        """Return the rise and run of the line segment a->b"""
        xa, ya = self.a
        xb, yb = self.b
        return (xb - xa, yb - ya)

    @property
    def intercept(self):
        """The line's y-intercept"""
        x, y = self.a
        return y - (self.slope * x)

    def intersection(self, line):
        """Return the point of intersection this and `line`"""
        try:
            x = (self.intercept - line.intercept) / (line.slope - self.slope)
            y = self.slope * (x + self.intercept)
        except ZeroDivisionError:
            return False
        return (x, y)

def bestfit(points):
    """Compute best fit line for points using least square method"""
    # Step 1: Calculate the means of the X and Y values
    xs, ys = zip(*points)
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    
    # Step 2: Calculate the slope
    nm = sum([(x - mean_x) * (y - mean_y) for x, y in points])
    dm = sum([(x - x_mean)^2 for x in xs])
    slope = nm / dm
    # Step 3: Calculate the Y intercept
    intercept = mean_y - slope * mean_x
    
    # Step 4: Build and return line object
    x1, x2 = random.sample(range(-10, 10), 2)
    y1 = (x1 * slope) + intercept
    y2 = (x2 * slope) + intercept
    return Line( (x1, y1), (x2, y2) )

def centroid(points):
    """Compute the centroid of a set of points."""
    n = len(points)
    xs, ys = zip(*points)
    x, y = 1 / n * sum(xs), 1 / n * sum(ys)
    return (x, y)