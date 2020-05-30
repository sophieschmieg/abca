package fields;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class MatrixAlgebra<T extends Element> extends AbstractAlgebra<T, Matrix<T>> {
	private Ring<T> ring;
	private int dimension;
	private Matrix<T> zero;
	private Matrix<T> id;
	private FreeModule<T> free;

	public MatrixAlgebra(Ring<T> ring, int dimension) {
		this.ring = ring;
		this.dimension = dimension;
		this.free = new FreeModule<T>(ring, dimension);
		List<List<T>> zeroes = new ArrayList<>();
		List<T> row = new ArrayList<>();
		for (int j = 0; j < dimension; j++) {
			row.add(ring.zero());
		}
		for (int i = 0; i < dimension; i++) {
			zeroes.add(row);
		}
		this.zero = new Matrix<T>(zeroes);
		List<List<T>> id = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> idRow = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					idRow.add(ring.one());
				} else {
					idRow.add(ring.zero());
				}
			}
			id.add(idRow);
		}
		this.id = new Matrix<T>(id);
	}

	@Override
	public Ring<T> getRing() {
		return ring;
	}

	@Override
	public Matrix<T> zero() {
		return zero;
	}

	@Override
	public Matrix<T> add(Matrix<T> s1, Matrix<T> s2) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				row.add(ring.add(s1.entry(i + 1, j + 1), s2.entry(i + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public Matrix<T> negative(Matrix<T> s) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				row.add(ring.negative(s.entry(i + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public Matrix<T> getEmbedding(T t) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					row.add(t);
				} else {
					row.add(ring.zero());
				}
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public Matrix<T> scalarMultiply(T t, Matrix<T> s) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				row.add(ring.multiply(t, s.entry(i + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public Matrix<T> getRandomElement() {
		List<List<T>> rand = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				row.add(ring.getRandomElement());
			}
			rand.add(row);
		}
		return new Matrix<T>(rand);
	}

	@Override
	public boolean isFinite() {
		return ring.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return ring.getNumberOfElements().pow(dimension * dimension);
	}

	@Override
	public Iterator<Matrix<T>> iterator() {
		return new Iterator<Matrix<T>>() {
			private List<List<Iterator<T>>> it = null;
			private List<List<T>> list = new ArrayList<List<T>>();

			private void init() {
				if (this.it != null)
					return;
				this.it = new ArrayList<>();
				for (int i = 0; i < dimension; i++) {
					List<Iterator<T>> rowIt = new ArrayList<>();
					List<T> rowList = new ArrayList<>();
					this.it.add(rowIt);
					this.list.add(rowList);
					for (int j = 0; j < dimension; j++) {
						rowIt.add(ring.iterator());
						if (i == 0 && j == 0)
							rowList.add(null);
						else
							rowList.add(this.it.get(i).get(j).next());
					}
				}
			}

			@Override
			public boolean hasNext() {
				init();
				for (int i = 0; i < dimension; i++) {
					for (int j = 0; j < dimension; j++) {
						if (this.it.get(i).get(j).hasNext()) {
							return true;
						}
					}
				}
				return false;
			}

			@Override
			public Matrix<T> next() {
				init();
				boolean broken = false;
				outerLoop: for (int i = 0; i < dimension; i++) {
					for (int j = 0; j < dimension; j++) {
						if (this.it.get(i).get(j).hasNext()) {
							this.list.get(i).set(j, this.it.get(i).get(j).next());
							for (int l = 0; l < i; l++) {
								for (int k = 0; k < dimension; k++) {
									this.it.get(l).set(k, ring.iterator());
									this.list.get(l).set(k, this.it.get(l).get(k).next());
								}
							}
							for (int k = 0; k < j; k++) {
								this.it.get(i).set(k, ring.iterator());
								this.list.get(i).set(k, this.it.get(i).get(k).next());
							}
							broken = true;
							break outerLoop;
						}
					}
				}
				if (!broken)
					throw new RuntimeException();
				return new Matrix<T>(list);
			}
		};
	}

	@Override
	public Matrix<T> one() {
		return id;
	}

	@Override
	public BigInteger characteristic() {
		return ring.characteristic();
	}

	@Override
	public Matrix<T> multiply(Matrix<T> t1, Matrix<T> t2) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			result.add(new ArrayList<>());
			for (int j = 0; j < dimension; j++) {
				T r = ring.zero();
				for (int k = 0; k < dimension; k++) {
					r = ring.add(r, ring.multiply(t1.entry(i + 1, k + 1), t2.entry(k + 1, j + 1)));
				}
				result.get(i).add(r);
			}
		}
		return new Matrix<T>(result);
	}

	public Matrix<T> swapOperation(int iIndex, int jIndex) {
		if (iIndex == jIndex) {
			return one();
		}
		iIndex--;
		jIndex--;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(one().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if ((i == iIndex && j == iIndex) || (i == jIndex && j == jIndex)) {
					matrix[i][j] = ring.zero();
				} else if ((i == iIndex && j == jIndex) || (i == jIndex && j == iIndex)) {
					matrix[i][j] = ring.one();
				} else if (i == j) {
					matrix[i][j] = ring.one();
				} else {
					matrix[i][j] = ring.zero();
				}
			}
		}
		return new Matrix<T>(matrix);
	}

	public Matrix<T> rowShearOperation(int iIndex, int jIndex, T value) {
		if (iIndex == jIndex || value.equals(ring.zero())) {
			return one();
		}
		iIndex--;
		jIndex--;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(one().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					matrix[i][j] = ring.one();
				} else if (i == jIndex && j == iIndex) {
					matrix[i][j] = value;
				} else {
					matrix[i][j] = ring.zero();
				}
			}
		}
		return new Matrix<T>(matrix);
	}

	public Matrix<T> colShearOperation(int iIndex, int jIndex, T value) {
		if (iIndex == jIndex || value.equals(ring.zero())) {
			return one();
		}
		iIndex--;
		jIndex--;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(one().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					matrix[i][j] = ring.one();
				} else if (i == iIndex && j == jIndex) {
					matrix[i][j] = value;
				} else {
					matrix[i][j] = ring.zero();
				}
			}
		}
		return new Matrix<T>(matrix);
	}

	public class GaussianEliminationResult {
		private Matrix<T> rowOperations = one();
		private Matrix<T> colOperations = one();
		private Matrix<T> rowOperationsInverse = one();
		private Matrix<T> colOperationsInverse = one();
		private int swaps = 0;
		private Matrix<T> diagonalMatrix;

		public Matrix<T> getRowOperations() {
			return rowOperations;
		}

		public Matrix<T> getColOperations() {
			return colOperations;
		}

		public Matrix<T> getRowOperationsInverse() {
			return rowOperationsInverse;
		}

		public Matrix<T> getColOperationsInverse() {
			return colOperationsInverse;
		}

		public Matrix<T> getDiagonalMatrix() {
			return diagonalMatrix;
		}

		public int getSwaps() {
			return swaps;
		}
	}

	public GaussianEliminationResult gaussianElimination(Matrix<T> t) {
		if (!ring.isEuclidean()) {
			throw new RuntimeException();
		}
		if (t.gaussianEliminationResult != null) {
			return t.gaussianEliminationResult;
		}
		GaussianEliminationResult r = new GaussianEliminationResult();
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				matrix[i][j] = t.entry(i + 1, j + 1);
			}
		}
		int pivotIndex = 0;
		while (pivotIndex < dimension) {
			BigInteger minMeasure = matrix[pivotIndex][pivotIndex].equals(ring.zero()) ? null
					: ring.euclidMeasure(matrix[pivotIndex][pivotIndex]);
			int minRow = pivotIndex;
			int minCol = pivotIndex;
			for (int i = pivotIndex; i < dimension; i++) {
				for (int j = pivotIndex; j < dimension; j++) {
					if (i == pivotIndex && j == pivotIndex) {
						continue;
					}
					if (matrix[i][j].equals(ring.zero())) {
						continue;
					}
					BigInteger measure = ring.euclidMeasure(matrix[i][j]);
					if (minMeasure == null || measure.compareTo(minMeasure) < 0) {
						minMeasure = measure;
						minRow = i;
						minCol = j;
					}
				}
			}
			if (minMeasure == null) {
				break;
			}
			if (pivotIndex == minRow && pivotIndex == minCol) {
				boolean found = false;
				for (int i = pivotIndex + 1; i < dimension; i++) {
					if (!matrix[i][pivotIndex].equals(ring.zero())) {
						found = true;
						break;
					}
				}
				for (int j = pivotIndex + 1; j < dimension; j++) {
					if (found || !matrix[pivotIndex][j].equals(ring.zero())) {
						found = true;
						break;
					}
				}
				if (!found) {
					pivotIndex++;
					continue;
				}
			}
			if (minRow != pivotIndex) {
				r.rowOperations = multiply(r.rowOperations, swapOperation(pivotIndex + 1, minRow + 1));
				r.rowOperationsInverse = multiply(swapOperation(pivotIndex + 1, minRow + 1), r.rowOperationsInverse);
				r.swaps++;
				for (int j = 0; j < dimension; j++) {
					T tmp = matrix[pivotIndex][j];
					matrix[pivotIndex][j] = matrix[minRow][j];
					matrix[minRow][j] = tmp;
				}
			}
			if (minCol != pivotIndex) {
				r.colOperations = multiply(swapOperation(minCol + 1, pivotIndex + 1), r.colOperations);
				r.colOperationsInverse = multiply(r.colOperationsInverse, swapOperation(minCol + 1, pivotIndex + 1));
				r.swaps++;
				for (int i = 0; i < dimension; i++) {
					T tmp = matrix[i][pivotIndex];
					matrix[i][pivotIndex] = matrix[i][minCol];
					matrix[i][minCol] = tmp;
				}
			}

			T pivot = matrix[pivotIndex][pivotIndex];
			for (int i = pivotIndex + 1; i < dimension; i++) {
				T shear = ring.negative(ring.quotientAndRemainder(matrix[i][pivotIndex], pivot).get(0));
				r.rowOperations = multiply(r.rowOperations,
						rowShearOperation(pivotIndex + 1, i + 1, ring.negative(shear)));
				r.rowOperationsInverse = multiply(rowShearOperation(pivotIndex + 1, i + 1, shear),
						r.rowOperationsInverse);
				for (int j = 0; j < dimension; j++) {
					matrix[i][j] = ring.add(matrix[i][j], ring.multiply(shear, matrix[pivotIndex][j]));
				}
			}
			for (int j = pivotIndex + 1; j < dimension; j++) {
				T shear = ring.negative(ring.quotientAndRemainder(matrix[pivotIndex][j], pivot).get(0));
				r.colOperations = multiply(colShearOperation(pivotIndex + 1, j + 1, ring.negative(shear)),
						r.colOperations);
				r.colOperationsInverse = multiply(r.colOperationsInverse,
						colShearOperation(pivotIndex + 1, j + 1, shear));
				for (int i = 0; i < dimension; i++) {
					matrix[i][j] = ring.add(matrix[i][j], ring.multiply(shear, matrix[i][pivotIndex]));
				}
			}
		}

		r.diagonalMatrix = new Matrix<T>(matrix);
		t.gaussianEliminationResult = r;
		return r;
	}

	public Polynomial<T> characteristicPolynomial(Matrix<T> t) {
		if (!(ring instanceof Field)) {
			throw new UnsupportedOperationException();
		}
		PolynomialRing<T> r = new PolynomialRing<T>((Field<T>) ring, 1, Polynomial.GREVLEX);
		MatrixAlgebra<Polynomial<T>> a = new MatrixAlgebra<>(r, dimension);
		List<List<Polynomial<T>>> asPolynomial = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			asPolynomial.add(new ArrayList<>());
			for (int j = 0; j < dimension; j++) {
				asPolynomial.get(i).add(r.getEmbedding(t.entry(i + 1, j + 1)));
			}
		}
		return a.determinant(a.subtract(a.getEmbedding(r.getVar(1)), new Matrix<>(asPolynomial)));
	}

	public T determinant(Matrix<T> t) {
		GaussianEliminationResult gauss = gaussianElimination(t);
		T result = ring.one();
		if (gauss.swaps % 2 == 1) {
			result = ring.negative(result);
		}
		for (int i = 0; i < dimension; i++) {
			result = ring.multiply(result, gauss.diagonalMatrix.entry(i + 1, i + 1));
		}
		return result;
	}

	public T trace(Matrix<T> t) {
		T tr = ring.zero();
		for (int i = 0; i < dimension; i++) {
			tr = ring.add(tr, t.entry(i + 1, i + 1));
		}
		return tr;
	}

	public Vector<T> multiply(Matrix<T> t, Vector<T> x) {
		List<T> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			result.add(free.innerProduct(t.row(i + 1), x));
		}
		return new Vector<T>(result);
	}

	public List<Vector<T>> kernelBasis(Matrix<T> t) {
		List<Vector<T>> result = new ArrayList<>();
		GaussianEliminationResult gauss = gaussianElimination(t);
		for (int i = 0; i < dimension; i++) {
			if (!gauss.diagonalMatrix.entry(i + 1, i + 1).equals(ring.zero())) {
				continue;
			}
			result.add(gauss.colOperationsInverse.column(i + 1));
		}
		return result;
	}

	public boolean hasSolution(Matrix<T> t, Vector<T> b) {
		GaussianEliminationResult gauss = gaussianElimination(t);
		Vector<T> rhs = multiply(gauss.rowOperationsInverse, b);
		for (int i = 0; i < dimension; i++) {
			if (!ring.isDivisible(rhs.get(i + 1), gauss.diagonalMatrix.entry(i + 1, i + 1))) {
				return false;
			}
		}
		return true;
	}

	public Vector<T> solve(Matrix<T> t, Vector<T> b) {
		List<T> result = new ArrayList<>();
		GaussianEliminationResult gauss = gaussianElimination(t);
		Vector<T> rhs = multiply(gauss.rowOperationsInverse, b);
		for (int i = 0; i < dimension; i++) {
			T r = ring.quotientAndRemainder(rhs.get(i + 1), gauss.diagonalMatrix.entry(i + 1, i + 1)).get(0);
			result.add(r);
		}
		return multiply(gauss.colOperationsInverse, new Vector<>(result));
	}

	@Override
	public boolean isUnit(Matrix<T> t) {
		return ring.isUnit(determinant(t));
	}

	@Override
	public Matrix<T> inverse(Matrix<T> t) {
		GaussianEliminationResult gauss = gaussianElimination(t);
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					matrix[i][j] = ring.inverse(gauss.diagonalMatrix.entry(i + 1, j + 1));
				} else {
					matrix[i][j] = ring.zero();
				}
			}
		}
		return multiply(gauss.colOperationsInverse, new Matrix<T>(matrix), gauss.rowOperationsInverse);
	}

	@Override
	public boolean isIntegral() {
		return false;
	}

	@Override
	public boolean isZeroDivisor(Matrix<T> t) {
		return kernelBasis(t).size() == 0;
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isDivisible(Matrix<T> dividend, Matrix<T> divisor) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<Matrix<T>> quotientAndRemainder(Matrix<T> dividend, Matrix<T> divisor) {
		throw new RuntimeException();
	}

	@Override
	public BigInteger euclidMeasure(Matrix<T> t) {
		return null;
	}

	@Override
	public Iterable<Matrix<T>> getUnits() {
		return new Iterable<Matrix<T>>() {

			@Override
			public Iterator<Matrix<T>> iterator() {
				return new Iterator<Matrix<T>>() {
					private Matrix<T> next = null;
					private Iterator<Matrix<T>> it = MatrixAlgebra.this.iterator();

					private void setNext() {
						if (next == null) {
							while (it.hasNext()) {
								Matrix<T> c = it.next();
								if (MatrixAlgebra.this.isUnit(c)) {
									next = c;
									break;
								}
							}
						}
					}

					@Override
					public boolean hasNext() {
						setNext();
						return next != null;
					}

					@Override
					public Matrix<T> next() {
						setNext();
						Matrix<T> result = next;
						next = null;
						setNext();
						return result;
					}
				};
			}
		};
	}

	@Override
	public List<Matrix<T>> getGenerators() {
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(one().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				matrix[i][j] = ring.zero();
			}
		}
		List<Matrix<T>> generators = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				matrix[i][j] = ring.one();
				generators.add(new Matrix<T>(matrix));
				matrix[i][j] = ring.zero();
			}
		}
		return generators;
	}
}
