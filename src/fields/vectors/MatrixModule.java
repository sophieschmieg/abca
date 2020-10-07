package fields.vectors;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Ring;

public class MatrixModule<T extends Element<T>> extends AbstractModule<T, Matrix<T>> {
	private Ring<T> ring;
	private int rows;
	private int columns;
	private Matrix<T> zero;
	private FreeModule<T> domain;
	private MatrixAlgebra<T> domainAlgebra;
	private FreeModule<T> codomain;
	private MatrixAlgebra<T> codomainAlgebra;

	MatrixModule(Ring<T> ring, int dimension, FreeModule<T> free, MatrixAlgebra<T> algebra) {
		this.ring = ring;
		this.rows = dimension;
		this.columns = dimension;
		this.domain = free;
		this.domainAlgebra = algebra;
		this.codomain = free;
		this.codomainAlgebra = algebra;
		List<List<T>> zeroes = new ArrayList<>();
		List<T> row = new ArrayList<>();
		for (int j = 0; j < columns; j++) {
			row.add(ring.zero());
		}
		for (int i = 0; i < rows; i++) {
			zeroes.add(row);
		}
		this.zero = new Matrix<T>(zeroes);
	}
	
	public MatrixModule(Ring<T> ring, int rows, int columns) {
		this.ring = ring;
		this.rows = rows;
		this.columns = columns;
		this.domain = new FreeModule<T>(ring, columns);
		this.domainAlgebra = this.domain.matrixAlgebra();
		this.codomain = new FreeModule<T>(ring, rows);
		this.codomainAlgebra = this.codomain.matrixAlgebra();
		List<List<T>> zeroes = new ArrayList<>();
		List<T> row = new ArrayList<>();
		for (int j = 0; j < columns; j++) {
			row.add(ring.zero());
		}
		for (int i = 0; i < rows; i++) {
			zeroes.add(row);
		}
		this.zero = new Matrix<T>(zeroes);
	}

	public FreeModule<T> domain() {
		return domain;
	}

	public FreeModule<T> codomain() {
		return codomain;
	}

	public MatrixAlgebra<T> domainAlgebra() {
		return domainAlgebra;
	}

	public MatrixAlgebra<T> codomainAlgebra() {
		return codomainAlgebra;
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
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < columns; j++) {
				row.add(ring.add(s1.entry(i + 1, j + 1), s2.entry(i + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public Matrix<T> negative(Matrix<T> s) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < columns; j++) {
				row.add(ring.negative(s.entry(i + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public Matrix<T> scalarMultiply(T t, Matrix<T> s) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < columns; j++) {
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
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < columns; j++) {
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
		return ring.getNumberOfElements().pow(rows * columns);
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
				for (int i = 0; i < rows; i++) {
					List<Iterator<T>> rowIt = new ArrayList<>();
					List<T> rowList = new ArrayList<>();
					this.it.add(rowIt);
					this.list.add(rowList);
					for (int j = 0; j < columns; j++) {
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
				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < columns; j++) {
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
				outerLoop: for (int i = 0; i < rows; i++) {
					for (int j = 0; j < columns; j++) {
						if (this.it.get(i).get(j).hasNext()) {
							this.list.get(i).set(j, this.it.get(i).get(j).next());
							for (int l = 0; l < i; l++) {
								for (int k = 0; k < columns; k++) {
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

	public Matrix<T> multiply(Matrix<T> t1, Matrix<T> t2) {
		if (t1.rows() != rows || t2.columns() != columns || t1.columns() != t2.rows()) {
			System.err.println(t1);
			System.err.println(t2);
			throw new ArithmeticException("Matricies have wrong dimension");
		}
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			result.add(new ArrayList<>());
			for (int j = 0; j < columns; j++) {
				T r = ring.zero();
				for (int k = 0; k < t1.columns(); k++) {
					r = ring.add(r, ring.multiply(t1.entry(i + 1, k + 1), t2.entry(k + 1, j + 1)));
				}
				result.get(i).add(r);
			}
		}
		return new Matrix<T>(result);
	}

	public Matrix<T> rowSwapOperation(int iIndex, int jIndex) {
		return codomainAlgebra.swapOperation(iIndex, jIndex);
	}

	public Matrix<T> colSwapOperation(int iIndex, int jIndex) {
		return domainAlgebra.swapOperation(iIndex, jIndex);
	}

	public Matrix<T> rowShearOperation(int iIndex, int jIndex, T value) {
		if (iIndex == jIndex || value.equals(ring.zero())) {
			return codomainAlgebra.one();
		}
		iIndex--;
		jIndex--;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(codomainAlgebra.one().entry(1, 1).getClass(), rows, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
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
			return domainAlgebra.one();
		}
		iIndex--;
		jIndex--;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(domainAlgebra.one().entry(1, 1).getClass(), columns, columns);
		for (int i = 0; i < columns; i++) {
			for (int j = 0; j < columns; j++) {
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
		private Matrix<T> rowOperations = codomainAlgebra.one();
		private Matrix<T> colOperations = domainAlgebra.one();
		private Matrix<T> rowOperationsInverse = codomainAlgebra.one();
		private Matrix<T> colOperationsInverse = domainAlgebra.one();
		private int swaps = 0;
		private Matrix<T> diagonalMatrix;
		private int rank = 0;

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
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = t.entry(i + 1, j + 1);
			}
		}
		int pivotIndex = 0;
		while (pivotIndex < Math.min(rows, columns)) {
			BigInteger minMeasure = matrix[pivotIndex][pivotIndex].equals(ring.zero()) ? null
					: ring.euclidMeasure(matrix[pivotIndex][pivotIndex]);
			int minRow = pivotIndex;
			int minCol = pivotIndex;
			for (int i = pivotIndex; i < rows; i++) {
				for (int j = pivotIndex; j < columns; j++) {
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
				for (int i = pivotIndex + 1; i < rows; i++) {
					if (!matrix[i][pivotIndex].equals(ring.zero())) {
						found = true;
						break;
					}
				}
				for (int j = pivotIndex + 1; j < columns; j++) {
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
				r.rowOperations = codomainAlgebra.multiply(r.rowOperations,
						rowSwapOperation(pivotIndex + 1, minRow + 1));
				r.rowOperationsInverse = codomainAlgebra.multiply(rowSwapOperation(pivotIndex + 1, minRow + 1),
						r.rowOperationsInverse);
				r.swaps++;
				for (int j = 0; j < columns; j++) {
					T tmp = matrix[pivotIndex][j];
					matrix[pivotIndex][j] = matrix[minRow][j];
					matrix[minRow][j] = tmp;
				}
			}
			if (minCol != pivotIndex) {
				r.colOperations = domainAlgebra.multiply(colSwapOperation(minCol + 1, pivotIndex + 1), r.colOperations);
				r.colOperationsInverse = domainAlgebra.multiply(r.colOperationsInverse,
						colSwapOperation(minCol + 1, pivotIndex + 1));
				r.swaps++;
				for (int i = 0; i < rows; i++) {
					T tmp = matrix[i][pivotIndex];
					matrix[i][pivotIndex] = matrix[i][minCol];
					matrix[i][minCol] = tmp;
				}
			}

			T pivot = matrix[pivotIndex][pivotIndex];
			for (int i = pivotIndex + 1; i < rows; i++) {
				T shear = ring.negative(ring.quotientAndRemainder(matrix[i][pivotIndex], pivot).get(0));
				r.rowOperations = codomainAlgebra.multiply(r.rowOperations,
						rowShearOperation(pivotIndex + 1, i + 1, ring.negative(shear)));
				r.rowOperationsInverse = codomainAlgebra.multiply(rowShearOperation(pivotIndex + 1, i + 1, shear),
						r.rowOperationsInverse);
				for (int j = 0; j < columns; j++) {
					matrix[i][j] = ring.add(matrix[i][j], ring.multiply(shear, matrix[pivotIndex][j]));
				}
			}
			for (int j = pivotIndex + 1; j < columns; j++) {
				T shear = ring.negative(ring.quotientAndRemainder(matrix[pivotIndex][j], pivot).get(0));
				r.colOperations = domainAlgebra.multiply(colShearOperation(pivotIndex + 1, j + 1, ring.negative(shear)),
						r.colOperations);
				r.colOperationsInverse = domainAlgebra.multiply(r.colOperationsInverse,
						domainAlgebra.colShearOperation(pivotIndex + 1, j + 1, shear));
				for (int i = 0; i < rows; i++) {
					matrix[i][j] = ring.add(matrix[i][j], ring.multiply(shear, matrix[i][pivotIndex]));
				}
			}
		}

		r.rank = pivotIndex;
		r.diagonalMatrix = new Matrix<T>(matrix);
		t.gaussianEliminationResult = r;
		return r;
	}

	public Vector<T> multiply(Matrix<T> t, Vector<T> x) {
		if (t.rows() != rows || x.dimension() != columns || t.columns() != columns) {
			System.err.println(t);
			System.err.println(x);
			throw new ArithmeticException("Matrix and Vector have wrong dimension");
		}
		List<T> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			T rowValue = ring.zero();
			for (int j = 0; j < columns; j++) {
				rowValue = ring.add(rowValue, ring.multiply(t.entry(i+1, j+1), x.get(j+1)));
			}
			result.add(rowValue);
		}
		return new Vector<>(result);
	}

	public List<Vector<T>> kernelBasis(Matrix<T> t) {
		List<Vector<T>> result = new ArrayList<>();
		GaussianEliminationResult gauss = gaussianElimination(t);
		for (int j = 0; j < columns; j++) {
			if (j < rows && !gauss.diagonalMatrix.entry(j + 1, j + 1).equals(ring.zero())) {
				continue;
			}
			result.add(gauss.colOperationsInverse.column(j + 1));
		}
		return result;
	}

	public boolean hasSolution(Matrix<T> t, Vector<T> b) {
		GaussianEliminationResult gauss = gaussianElimination(t);
		Vector<T> rhs = codomainAlgebra.multiply(gauss.rowOperationsInverse, b);
		for (int i = 0; i < rows; i++) {
			if (i < columns && !ring.isDivisible(rhs.get(i + 1), gauss.diagonalMatrix.entry(i + 1, i + 1))) {
				return false;
			}
			if (i >= columns && !rhs.get(i + 1).equals(ring.zero())) {
				return false;
			}
		}
		return true;
	}

	public Vector<T> solve(Matrix<T> t, Vector<T> b) {
		List<T> result = new ArrayList<>();
		GaussianEliminationResult gauss = gaussianElimination(t);
		Vector<T> rhs = codomainAlgebra.multiply(gauss.rowOperationsInverse, b);
		for (int i = 0; i < columns; i++) {
			if (i < rows && (!rhs.get(i + 1).equals(ring.zero()) || !gauss.diagonalMatrix.entry(i + 1, i + 1).equals(ring.zero()))) {
				T r = ring.divide(rhs.get(i + 1), gauss.diagonalMatrix.entry(i + 1, i + 1));
				result.add(r);
			} else {
				result.add(ring.zero());
			}
		}
		return domainAlgebra.multiply(gauss.colOperationsInverse, new Vector<>(result));
	}
	
	public int rank(Matrix<T> m) {
		return gaussianElimination(m).rank;
	}

	@Override
	public List<Matrix<T>> getModuleGenerators() {
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(zero().entry(1, 1).getClass(), rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = ring.zero();
			}
		}
		List<Matrix<T>> generators = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = ring.one();
				generators.add(new Matrix<T>(matrix));
				matrix[i][j] = ring.zero();
			}
		}
		return generators;
	}

	public Vector<T> asVector(Matrix<T> m) {
		List<T> asList = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				asList.add(m.entry(i + 1, j + 1));
			}
		}
		return new Vector<>(asList);
	}
	
	public FreeModule<T> asFreeModule() {
		return new FreeModule<>(ring, rows * columns);
	}

	@Override
	public boolean isLinearIndependent(List<Matrix<T>> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Matrix<T> m : s) {
			asVectors.add(asVector(m));
		}
		return asFreeModule().isLinearIndependent(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<Matrix<T>> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Matrix<T> m : s) {
			asVectors.add(asVector(m));
		}
		return asFreeModule().isGeneratingModule(asVectors);
	}

}
