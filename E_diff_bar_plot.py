import sys
import matplotlib.pyplot as plt
import numpy as np

def read_data(file_path):
    labels = []
    values = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            labels.append(parts[0])
            values.append(float(parts[1]))
    return labels, values

def plot_data(labels, values, wt=15.6):
    plt.figure(figsize=(12, 6))

    # Generate the positions of the bars
    x_positions = np.arange(len(labels)) * 1.5  # Increase the spacing between bars

    # Determine the color of each bar
    bar_colors = ['salmon' if value < wt else 'deepskyblue' for value in values]

    # Create the bar plot
    plt.bar(x_positions, values, width=0.4, color=bar_colors)  # Apply bar colors

    # Add horizontal line at y=10
    plt.axhline(y=wt, color='grey', linestyle='--', linewidth=1, alpha=0.5)

    # Set the x-ticks to the middle of each bar
    plt.xticks(x_positions, labels, rotation=90)

    plt.ylabel('Δ$E$')
    plt.tight_layout()  # Adjust layout to make room for label rotation
    plt.savefig("bar_plot.png", format='png')
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    labels, values = read_data(file_path)
    plot_data(labels, values)

