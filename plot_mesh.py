import matplotlib.pyplot as plt
import pandas as pd
import sys
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import math


def func0(x, R):
    return R**3 * math.cos(3 * x) + 4 * R**5 * math.cos(5 * x)
def func1(x, R):
    return R**3 * math.sin(3 * x) - 4 * R**5 * math.sin(5 * x)
def func2(x, R):
    return -math.sqrt(15) * R**4 * math.sin(4 * x)

#def func0(x, R):
#    return (1 + 0.25*math.cos(3*x))*math.cos(2*x)
#def func1(x, R):
#    return (1 + 0.25*math.cos(3*x))*math.sin(2*x)
#def func2(x, R):
#    return 0.25*math.sin(3*x)

#def func0(x,R):
#    return R*math.cos(x) + R**3 * math.cos(3*x)/3
#def func1(x,R):
#    return R*math.sin(x) - R**3 * math.sin(3*x)/3
#def func2(x,R):
#    return R*R*math.cos(2*x)



def plot_mesh(ax, csv_filename, show_colors, title):
    # Read the data from the CSV file
    df = pd.read_csv(csv_filename)
    # Separate points and edges
    vertex_df = df[df["type"] == "vertex"]
    edge_df = df[df["type"] == "edge"]

    if show_colors:
        # Check if there are enough points to plot the surface
        if len(vertex_df) >= 3:
            # Plot the surface with colors
            trisurf = ax.plot_trisurf(vertex_df["x"], vertex_df["y"], vertex_df["z"], cmap="jet", edgecolor='none')
            # Add a color bar
            fig.colorbar(trisurf, ax=ax, label='Solution')
        else:
            print("Not enough points to plot the surface.")
    else:
        # Plot vertices (solution values) with colors
        sc = ax.scatter(vertex_df["x"], vertex_df["y"], vertex_df["z"], c=vertex_df["z"], cmap="jet", s=20, label="Solution")
        # Add a color bar
        fig.colorbar(sc, ax=ax, label='Solution')

    # Plot edges
    edges = []
    for i in range(0, len(edge_df), 2):
        edges.append([(edge_df.iloc[i]["x"], edge_df.iloc[i]["y"], edge_df.iloc[i]["z"]),
                      (edge_df.iloc[i+1]["x"], edge_df.iloc[i+1]["y"], edge_df.iloc[i+1]["z"])])

    ax.add_collection3d(Line3DCollection(edges, colors="black", linewidths=1))

    # Set labels and title
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Solution")
    ax.set_title(title)
    ax.view_init(elev=30, azim=60)  # Adjustable rotation

def plot_mesh_par(ax, csv_filename,  title):
    # Read the data from the CSV file
    df = pd.read_csv(csv_filename)
    # Separate points and edges
    vertex_df = df[df["type"] == "vertex"]
    edge_df = df[df["type"] == "edge"]

    # Plot vertices (solution values) with colors
    sc = ax.scatter(vertex_df["u1"], vertex_df["u2"], vertex_df["u3"], c=vertex_df["u3"], cmap="jet", s=20, label="Solution")
    # Add a color bar
    fig.colorbar(sc, ax=ax, label='Solution')

    # Plot edges
    edges = []
    for i in range(0, len(edge_df), 2):
        edges.append([(edge_df.iloc[i]["u1"], edge_df.iloc[i]["u2"], edge_df.iloc[i]["u3"]),
                      (edge_df.iloc[i+1]["u1"], edge_df.iloc[i+1]["u2"], edge_df.iloc[i+1]["u3"])])

    ax.add_collection3d(Line3DCollection(edges, colors="black", linewidths=1))

    # Add the parametric curve

    R = 0.5
    n = [i * 2 * math.pi / 1000 for i in range(1001)]
    x = [func0(n[i], R) for i in range(1001)]
    y = [func1(n[i], R) for i in range(1001)]
    z = [func2(n[i], R) for i in range(1001)]

    ax.plot(x, y, z, color='red', lw = 4)

    # Set labels and title
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Solution")
    ax.set_title(title)
    ax.view_init(elev=30, azim=60)  

# Run the plot
if __name__ == "__main__":
    parametric = len(sys.argv) > 1 and sys.argv[1] == "parametric"

    show_colors = len(sys.argv) > 1 + parametric and sys.argv[1 + parametric] == "color"

    if not parametric:

        fig = plt.figure(figsize=(15, 6))

        # Create subplots
        ax1 = fig.add_subplot(131, projection='3d')
        plot_mesh(ax1, "postProcessing/start.csv", show_colors, "Initial Mesh")

        ax2 = fig.add_subplot(132, projection='3d')
        plot_mesh(ax2, "postProcessing/solutionNewton.csv", show_colors, "Final mesh with Newton's method")

        ax3 = fig.add_subplot(133, projection='3d')
        plot_mesh(ax3, "postProcessing/solutionPicardi.csv", show_colors, "Final mesh with Picardi's method")

        plt.show()
    
    else:
        fig = plt.figure(figsize=(15, 6))

        # Create subplots
        ax1 = fig.add_subplot(121, projection='3d')
        plot_mesh_par(ax1, "postProcessing/startParametric.csv", "Initial Mesh")

        ax2 = fig.add_subplot(122, projection='3d')
        plot_mesh_par(ax2, "postProcessing/solutionParametric.csv", "Final mesh with Newton's method")

        plt.show()