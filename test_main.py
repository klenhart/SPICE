import numpy as np


def main():
    expression_vector = np.array([[0.9], [0.1]])
    distance_matrix = np.array([[0.0, 0.25], [0.5, 0.0]])
    print(np.dot(distance_matrix, expression_vector))


if __name__ == "__main__":
    main()
