import numpy as np

def knn(data: np.ndarray, k: int):
  length = np.array([[np.sum((x - y) ** 2) for y in data] for x in data])
  print(length[0])