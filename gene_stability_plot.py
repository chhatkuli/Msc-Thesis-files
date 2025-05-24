import matplotlib.pyplot as plt
from typing import List

# gene numbers
x = [5, 4, 3, 2]

# stability values
y = [3.048, 1.827, 0.273, 0.038]

gene_names: List[str] = [
    "18srRNA, OsUBI, UBQ5, OsActin, eEF1a",
    "UBQ5, 18srRNA, OsUBI, OsActin",
    "18srRNA, OsUBI, UBQ5",
    "OsUBI, 18srRNA",
]


def plot_stability_values_for_genes(x: List[int], y: List[float], labels: List[str]):
    """
    plot gene stability values against number of genes used

    params:
    x: num of genes used (int)
    y: avg stability values (M) (int or float)
    labels: gene names (str)

    returns:
    a plot with x, y and labels

    raises:
    -ValueError if inputs are not proper types or lengths dont match
    """

    # validate input types
    if not (isinstance(x, list) and all(isinstance(i, int) for i in x)):
        raise ValueError("x should be a list of integers")

    if not (isinstance(y, list) and all(isinstance(i, (float, int)) for i in y)):
        raise ValueError("y should be a list of floats or integers")

    if not (isinstance(labels, list) and all(isinstance(i, str) for i in labels)):
        raise ValueError("labels should be a list of strings")

    if not (len(x) == len(y) == len(labels)):
        raise ValueError("inputs should have same length")

    plt.scatter(x, y)

    # labels
    for px, py, pl in zip(x, y, labels):
        plt.text(px, py + 0.1, f"{pl}", fontsize=9)

    # padding so labels dont get clipped out
    plt.xlim(min(x) - 1, max(x) + 1)
    plt.ylim(min(y) - 0.5, max(y) + 0.5)

    plt.gca().invert_xaxis()
    plt.xticks(ticks=sorted(x))
    plt.xlabel("Genes (least stable to most stable)")
    plt.ylabel("M value")
    # plt.title("gene stability ranking")
    plt.tight_layout()

    plt.savefig("gene_stability_plot.png", dpi=300, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    plot_stability_values_for_genes(x, y, gene_names)
