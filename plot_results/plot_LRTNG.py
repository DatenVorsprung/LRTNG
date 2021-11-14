# plotting the outputs from LRT-NG (ellipse, circle, intersections etc.)

import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import svd
from numpy.linalg import inv
from numpy.linalg import norm

import argparse


def draw_ellipse_circle_intersection(ellipse, radius_circle, color, alph, axis_3d):
    c = ellipse[1:3]
    r = ellipse[3]
    M = np.reshape(ellipse[4:8], (2, 2)) / r ** 2

    # "singular value decomposition" to extract the orientation and the
    # axes of the ellipsoid
    _, D, V = svd(M)

    plot_grid = 50

    # get the major and minor axes
    a = 1 / np.sqrt(D[0])
    b = 1 / np.sqrt(D[1])

    theta = np.arange(0, 2 * np.pi + 1 / plot_grid, 1 / plot_grid)

    # parametric equation of the ellipse
    state = np.zeros((2, np.size(theta)))
    state[0, :] = a * np.cos(theta)
    state[1, :] = b * np.sin(theta)

    # coordinate transform
    X = V @ state

    # intersection with circle
    for idx, (x, y) in enumerate(zip(X[0], X[1])):
        dist_ellipse = norm([x, y])  # radius ellipse in this direction in CC
        if dist_ellipse > radius_circle:
            X[:, idx] = [x, y] / dist_ellipse * radius_circle

    X[0] += c[0]
    X[1] += c[1]

    axis_3d.plot(xs=X[0], ys=X[1], zs=ellipse[0], color=color, alpha=alph)


def plot(time_horizon, dim, axis1, axis2, file_ellipse, file_circle, color, alph, axis_3d, skip_reachsets=1, start=0):
    data_ellipse = np.loadtxt(file_ellipse)
    data_circle = np.loadtxt(file_circle)

    # permutation matrix to project on axis1 and axis2
    P = np.eye(dim)
    P[:, [0, axis1]] = P[:, [axis1, 0]]
    P[:, [1, axis2]] = P[:, [axis2, 1]]

    count = 1

    for ellipse, circle in zip(data_ellipse[start + 1:], data_circle[start:]):
        # if (ellipse[0]/time_step)%skip_reachsets != 0:
        if count != skip_reachsets:
            count += 1
            continue

        count = 1

        ellipse2 = ellipse
        radiusCircle = circle[dim + 1]  # needed to create intersection

        # create ellipse plotting values for 2d projection
        # https://math.stackexchange.com/questions/2438495/showing-positive-definiteness-in-the-projection-of-ellipsoid
        # construct ellipse2 to have a 2-dimensional ellipse as an input to
        # ellipse_plot

        if dim > 2:
            center = ellipse[1:dim + 1]
            ellipse2[1] = center[axis1]
            ellipse2[2] = center[axis2]
            radius_ellipse = ellipse[dim + 1]
            m1 = np.reshape(ellipse[dim + 2:], (dim, dim))
            m1 = m1 / radius_ellipse ** 2
            m1 = P.transpose() @ m1 @ P  # permutation to project on chosen axes
            ellipse2[3] = 1  # because radius is already in m1

            # plot ellipse onto axis1-axis2 plane
            J = m1[0:2, 0:2]
            K = m1[2:, 2:]
            L = m1[2:, 0:2]
            m2 = J - L.transpose() @ inv(K) @ L
            ellipse2[4:8] = m2.reshape(1, -1)

        draw_ellipse_circle_intersection(ellipse2[0:8], radiusCircle, color, alph, axis_3d)

        if ellipse[0] >= time_horizon:
            break

    # plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--time_step", type=float, required=True)
    parser.add_argument("--time_horizon", type=float, required=True)
    parser.add_argument("--dim", type=int, required=True)
    parser.add_argument("--output_directory", default="../saved_outputs/")
    parser.add_argument("--circle_file", required=True)
    parser.add_argument("--ellipse_file", required=True)
    parser.add_argument("--axis1", default=0, type=int)
    parser.add_argument("--axis2", default=1, type=int)
    # initial radius
    parser.add_argument("--color", default="orange")
    parser.add_argument("--alpha", default=0.6, type=float)
    parser.add_argument("--skip_reachsets", default=1, type=int)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    plot(
        args.time_horizon,
        args.dim,
        args.axis1,
        args.axis2,
        args.output_directory + args.ellipse_file,
        args.output_directory + args.circle_file,
        color=args.color,
        alph=args.alpha,
        axis_3d=ax,
        skip_reachsets=args.skip_reachsets)

    ax.view_init(elev=-10.)
    ax.tick_params(axis='both', labelsize=8)
    ax.set_xlabel('x'+str(args.axis1), fontsize=8)
    ax.set_ylabel('x'+str(args.axis2), fontsize=8)
    ax.set_zlabel('Time (s)', fontsize=8)

    if args.save:
        plt.savefig('LRTNG_plot.pdf')
    else:
        plt.show()
