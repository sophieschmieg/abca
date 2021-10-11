package fields.vectors;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Field;
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

	@Override
	public Exactness exactness() {
		return ring.exactness();
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

	public Matrix<T> multiply(int[] permutation, Matrix<T> t) {
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(codomainAlgebra.one().entry(1, 1).getClass(), rows, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				matrix[i][j] = t.entry(permutation[i] + 1, j + 1);
			}
		}
		return new Matrix<>(matrix);
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

	// LD^-1U=PA
	public class LDUPResult {
		private int[] permutation;
		private Matrix<T> lowerTriangle;
		private Matrix<T> inverseDiagonal;
		private Matrix<T> upperTriangle;
		private int signum;
		private T determinant;
		private Set<Integer> slips;

		public int[] getPermutation() {
			return permutation;
		}

		public Matrix<T> getLowerTriangle() {
			return lowerTriangle;
		}

		public Matrix<T> getUpperTriangle() {
			return upperTriangle;
		}

		public int getSignum() {
			return signum;
		}

		public T getDeterminant() {
			return determinant;
		}

		public Matrix<T> getInverseDiagonal() {
			return inverseDiagonal;
		}

		public Set<Integer> getSlips() {
			return slips;
		}
	}

	public LDUPResult ldup(Matrix<T> t) {
		if (t.rows() != rows || t.columns() != columns) {
			throw new RuntimeException("Wrong Algebra");
		}
		if (!ring.isIntegral()) {
			throw new RuntimeException();
		}
		if (t.ldupResult != null) {
			return t.ldupResult;
		}
		LDUPResult r = new LDUPResult();
		r.slips = new TreeSet<>();
		r.determinant = ring.one();
		r.signum = 1;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, columns);
		int[] permutation = new int[rows];
		@SuppressWarnings("unchecked")
		T[][] inverseDiagonal = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, rows);
		@SuppressWarnings("unchecked")
		T[][] lower = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				lower[i][j] = i == j ? ring.one() : ring.zero();
				inverseDiagonal[i][j] = i == j ? ring.one() : ring.zero();
			}
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = t.entry(i + 1, j + 1);
			}
			permutation[i] = i;
		}
		int pivotColumn = -1;
		T prevPivot = ring.one();
		for (int pivotRow = 0; pivotRow < rows; pivotRow++) {
			pivotColumn++;
			if (pivotColumn >= columns) {
				break;
			}
			List<T> column = new ArrayList<>();
			for (int i = pivotRow; i < rows; i++) {
				column.add(matrix[i][pivotColumn]);
			}
			int minRow = ring.preferredPivotStrategy().pivot(column);
			if (minRow < 0) {
				pivotRow--;
				r.determinant = ring.zero();
				r.slips.add(pivotColumn + 1);
				continue;
			}
			minRow += pivotRow;
			if (minRow != pivotRow) {
				r.signum *= -1;
				int tmpIndex = permutation[pivotRow];
				permutation[pivotRow] = permutation[minRow];
				permutation[minRow] = tmpIndex;
				for (int j = pivotColumn; j < columns; j++) {
					T tmp = matrix[pivotRow][j];
					matrix[pivotRow][j] = matrix[minRow][j];
					matrix[minRow][j] = tmp;
				}
				for (int j = 0; j < pivotRow; j++) {
					T tmp = lower[pivotRow][j];
					lower[pivotRow][j] = lower[minRow][j];
					lower[minRow][j] = tmp;
				}
			}
			T pivot = matrix[pivotRow][pivotColumn];
			lower[pivotRow][pivotRow] = pivot;
			inverseDiagonal[pivotRow][pivotRow] = ring.multiply(prevPivot, pivot);
			r.determinant = ring.divideChecked(ring.multiply(r.determinant, pivot), prevPivot);
			for (int i = pivotRow + 1; i < rows; i++) {
				T negativeShear = matrix[i][pivotColumn];
				T shear = ring.negative(negativeShear);
				lower[i][pivotRow] = negativeShear;
				for (int j = pivotColumn; j < columns; j++) {
					matrix[i][j] = ring.divideChecked(
							ring.add(ring.multiply(pivot, matrix[i][j]), ring.multiply(shear, matrix[pivotRow][j])),
							prevPivot);
				}
			}
			prevPivot = pivot;
		}
		for (int j = pivotColumn + 1; j < columns; j++) {
			r.slips.add(j);
		}
		r.upperTriangle = new Matrix<T>(matrix);
		r.inverseDiagonal = new Matrix<T>(inverseDiagonal);
		r.lowerTriangle = new Matrix<T>(lower);
		r.permutation = invertPermutation(permutation);
		t.ldupResult = r;
		return r;
	}

	public Matrix<T> lupUpperTriangle(Matrix<T> t) {
		LDUPResult ldup = ldup(t);
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = ring.divideChecked(
						ring.multiply(ldup.getLowerTriangle().entry(i + 1, i + 1),
								ldup.getUpperTriangle().entry(i + 1, j + 1)),
						ldup.getInverseDiagonal().entry(i + 1, i + 1));
			}
		}
		return new Matrix<>(matrix);
	}

	public Matrix<T> permuteRows(int[] permutation, Matrix<T> t) {
		if (t.rows() != rows || t.columns() != columns || permutation.length != rows) {
			throw new RuntimeException("Wrong Algebra");
		}
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < columns; j++) {
				row.add(t.entry(permutation[i] + 1, j + 1));
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	public Matrix<T> permuteColumns(int[] permutation, Matrix<T> t) {
		if (t.rows() != rows || t.columns() != columns || permutation.length != rows) {
			throw new RuntimeException("Wrong Algebra");
		}
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < columns; j++) {
				row.add(t.entry(i + 1, permutation[j] + 1));
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	public Vector<T> permuteVector(int[] permution, Vector<T> t) {
		List<T> result = new ArrayList<>();
		for (int i = 0; i < permution.length; i++) {
			result.add(t.get(permution[i] + 1));
		}
		return new Vector<>(result);
	}

	// R^-1.A.C^-1 = D
	// R.D.C = A
	public class SmithNormalFormResult {
		private Matrix<T> rowOperations;
		private Matrix<T> colOperations;
		private Matrix<T> rowOperationsInverse;
		private Matrix<T> colOperationsInverse;
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

		public int getRank() {
			return rank;
		}
	}

	public SmithNormalFormResult smithNormalForm(Matrix<T> t) {
		if (!ring.isEuclidean()) {
			throw new RuntimeException();
		}
		if (t.smithNormalFormResult != null) {
			return t.smithNormalFormResult;
		}
		if (t.rows() != rows || t.columns() != columns) {
			throw new RuntimeException("Wrong Algebra");
		}
		SmithNormalFormResult r = new SmithNormalFormResult();
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, columns);
		@SuppressWarnings("unchecked")
		T[][] rowOps = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, rows);
		@SuppressWarnings("unchecked")
		T[][] colOps = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), columns, columns);
		@SuppressWarnings("unchecked")
		T[][] rowOpsInv = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), rows, rows);
		@SuppressWarnings("unchecked")
		T[][] colOpsInv = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), columns, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = t.entry(i + 1, j + 1);
			}
		}
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				rowOps[i][j] = i == j ? ring.one() : ring.zero();
				rowOpsInv[i][j] = i == j ? ring.one() : ring.zero();
			}
		}
		for (int i = 0; i < columns; i++) {
			for (int j = 0; j < columns; j++) {
				colOps[i][j] = i == j ? ring.one() : ring.zero();
				colOpsInv[i][j] = i == j ? ring.one() : ring.zero();
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
				r.swaps++;
				for (int j = 0; j < rows; j++) {
					T tmp = rowOps[j][pivotIndex];
					rowOps[j][pivotIndex] = rowOps[j][minRow];
					rowOps[j][minRow] = tmp;
					tmp = rowOpsInv[pivotIndex][j];
					rowOpsInv[pivotIndex][j] = rowOpsInv[minRow][j];
					rowOpsInv[minRow][j] = tmp;
				}
				for (int j = 0; j < columns; j++) {
					T tmp = matrix[pivotIndex][j];
					matrix[pivotIndex][j] = matrix[minRow][j];
					matrix[minRow][j] = tmp;
				}
			}
			if (minCol != pivotIndex) {
				r.swaps++;
				for (int i = 0; i < columns; i++) {
					T tmp = colOps[pivotIndex][i];
					colOps[pivotIndex][i] = colOps[minCol][i];
					colOps[minCol][i] = tmp;
					tmp = colOpsInv[i][pivotIndex];
					colOpsInv[i][pivotIndex] = colOpsInv[i][minCol];
					colOpsInv[i][minCol] = tmp;
				}
				for (int i = 0; i < rows; i++) {
					T tmp = matrix[i][pivotIndex];
					matrix[i][pivotIndex] = matrix[i][minCol];
					matrix[i][minCol] = tmp;
				}
			}

			T pivot = matrix[pivotIndex][pivotIndex];
			for (int i = pivotIndex + 1; i < rows; i++) {
				if (matrix[i][pivotIndex].equals(ring.zero())) {
					continue;
				}
				T shear = ring.negative(ring.divide(matrix[i][pivotIndex], pivot));

				for (int j = 0; j < rows; j++) {
					rowOps[j][pivotIndex] = ring.add(rowOps[j][pivotIndex],
							ring.multiply(ring.negative(shear), rowOps[j][i]));
					rowOpsInv[i][j] = ring.add(rowOpsInv[i][j], ring.multiply(shear, rowOpsInv[pivotIndex][j]));
				}
				for (int j = 0; j < columns; j++) {
					matrix[i][j] = ring.add(matrix[i][j], ring.multiply(shear, matrix[pivotIndex][j]));
				}
			}
			for (int j = pivotIndex + 1; j < columns; j++) {
				if (matrix[pivotIndex][j].equals(ring.zero())) {
					continue;
				}
				T shear = ring.negative(ring.divide(matrix[pivotIndex][j], pivot));

				for (int i = 0; i < columns; i++) {
					colOps[pivotIndex][i] = ring.add(colOps[pivotIndex][i],
							ring.multiply(ring.negative(shear), colOps[j][i]));
					colOpsInv[i][j] = ring.add(colOpsInv[i][j], ring.multiply(shear, colOpsInv[i][pivotIndex]));
				}
				for (int i = 0; i < rows; i++) {
					matrix[i][j] = ring.add(matrix[i][j], ring.multiply(shear, matrix[i][pivotIndex]));
				}
			}
		}

		r.rank = pivotIndex;
		r.diagonalMatrix = new Matrix<T>(matrix);
		r.rowOperations = new Matrix<T>(rowOps);
		r.rowOperationsInverse = new Matrix<T>(rowOpsInv);
		r.colOperations = new Matrix<T>(colOps);
		r.colOperationsInverse = new Matrix<T>(colOpsInv);
		t.smithNormalFormResult = r;
		return r;
	}

	public Matrix<T> invertUpperTriangleMatrix(Matrix<T> t) {
		if (t.rows() != rows || t.columns() != columns) {
			throw new ArithmeticException("Matrix and Vector have wrong dimension");
		}
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < i; j++) {
				row.add(ring.zero());
			}
			row.add(ring.inverse(t.entry(i + 1, i + 1)));
			for (int j = i + 1; j < columns; j++) {
				T value = ring.zero();
				for (int k = i; k < j; k++) {
					value = ring.add(ring.multiply(row.get(k), t.entry(k + 1, j + 1)), value);
				}
				row.add(ring.divideChecked(ring.negative(value), t.entry(j + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	public Vector<T> multiply(Matrix<T> t, Vector<T> x) {
		if (t.rows() != rows || x.dimension() != columns || t.columns() != columns) {
			throw new ArithmeticException("Matrix and Vector have wrong dimension");
		}
		List<T> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			T rowValue = ring.zero();
			for (int j = 0; j < columns; j++) {
				rowValue = ring.add(rowValue, ring.multiply(t.entry(i + 1, j + 1), x.get(j + 1)));
			}
			result.add(rowValue);
		}
		return new Vector<>(result);
	}

	public Vector<T> multiply(int[] rowPermutation, Vector<T> x) {
		if (rowPermutation.length != rows || x.dimension() != rows) {
			throw new ArithmeticException("Matrix and Vector have wrong dimension");
		}
		List<T> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			result.add(x.get(rowPermutation[i] + 1));
		}
		return new Vector<>(result);
	}

	public int[] invertPermutation(int[] permutation) {
		int[] result = new int[permutation.length];
		for (int i = 0; i < permutation.length; i++) {
			result[permutation[i]] = i;
		}
		return result;
	}

	public List<Vector<T>> kernelBasis(Matrix<T> t) {
		if (t.equals(zero())) {
			return codomain.getBasis();
		}
		List<Vector<T>> result = new ArrayList<>();
		if (ring.isEuclidean()) {
			SmithNormalFormResult smith = smithNormalForm(t);
			for (int i = 0; i < columns; i++) {
				if (i >= rows || smith.getDiagonalMatrix().entry(i + 1, i + 1).equals(ring.zero())) {
					result.add(smith.getColOperationsInverse().column(i + 1));
				}
			}
			return result;
		}
		LDUPResult gauss = ldup(t);
		// LD^-1U=PA
		// LD^-1Ux=PAx=0
		// LD^-1y = 0 <=> y = 0
		// x in ker A <=> Ux = 0
		int numPrevSlips = 0;
		for (int i : gauss.getSlips()) {
			List<T> vector = new ArrayList<>();
			for (int k = 0; k < columns; k++) {
				vector.add(ring.zero());
			}
			int row = i - numPrevSlips;
			for (int j = i - 1; j >= 0; j--) {
				if (gauss.getSlips().contains(j)) {
					continue;
				}
				vector.set(i, ring.one());
				row--;
			}
			result.add(new Vector<>(vector));
			numPrevSlips++;
		}
		return result;
	}

	public List<Vector<T>> imageBasis(Matrix<T> t) {
		if (ring instanceof Field<?>) {
			List<Vector<T>> asVectors = new ArrayList<>();
			LDUPResult gauss = ldup(t);
			for (int j = 0; j < columns; j++) {
				if (gauss.getSlips().contains(j)) {
					continue;
				}
				asVectors.add(t.column(j + 1));
			}
			return asVectors;
		}
		SmithNormalFormResult gauss = smithNormalForm(t);
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < gauss.rank; j++) {
				row.add(gauss.getDiagonalMatrix().entry(i + 1, j + 1));
			}
			result.add(row);
		}
		Matrix<T> asMatrix = new MatrixModule<>(ring, rows, gauss.rank).multiply(gauss.rowOperations,
				new Matrix<>(result));
		List<Vector<T>> asVectors = new ArrayList<>();
		for (int i = 0; i < gauss.rank; i++) {
			asVectors.add(asMatrix.column(i + 1));
		}
		return asVectors;
	}

	private boolean hasSolutionField(Matrix<T> t, Vector<T> b) {
		LDUPResult gauss = ldup(t);
		return hasSolutionEuclid(t, b);
	}

	private boolean hasSolutionEuclid(Matrix<T> t, Vector<T> b) {
		SmithNormalFormResult gauss = smithNormalForm(t);
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

	public boolean hasSolution(Matrix<T> t, Vector<T> b) {
		if (ring instanceof Field<?>) {
			return hasSolutionField(t, b);
		} else if (ring.isEuclidean()) {
			return hasSolutionEuclid(t, b);
		}
		throw new ArithmeticException("Cannot solve!");
	}

	private Vector<T> solveField(Matrix<T> t, Vector<T> b) {
//		List<T> result = new ArrayList<>();
//		LUPResult lup = lup(t);
//		b = multiply(lup.getPermutation(), b);
//		b = codomain.scalarMultiply(lup.multiplier, b);
		return solveEuclid(t, b);
	}

	private Vector<T> solveEuclid(Matrix<T> t, Vector<T> b) {
		List<T> result = new ArrayList<>();
		SmithNormalFormResult gauss = smithNormalForm(t);
		Vector<T> rhs = codomainAlgebra.multiply(gauss.rowOperationsInverse, b);
		for (int i = 0; i < columns; i++) {
			if (i < rows && (!rhs.get(i + 1).equals(ring.zero())
					|| !gauss.diagonalMatrix.entry(i + 1, i + 1).equals(ring.zero()))) {
				T r = ring.divideChecked(rhs.get(i + 1), gauss.diagonalMatrix.entry(i + 1, i + 1));
				result.add(r);
			} else {
				result.add(ring.zero());
			}
		}
		return domainAlgebra.multiply(gauss.colOperationsInverse, new Vector<>(result));
	}

	public Vector<T> solve(Matrix<T> t, Vector<T> b) {
		if (ring instanceof Field<?>) {
			return solveField(t, b);
		} else if (ring.isEuclidean()) {
			return solveEuclid(t, b);
		}
		throw new ArithmeticException("Cannot solve!");
	}

	public int rank(Matrix<T> m) {
		return columns - ldup(m).getSlips().size();
	}

	public Matrix<T> transpose(Matrix<T> t) {
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), t.rows(), t.columns());
		for (int i = 0; i < t.rows(); i++) {
			for (int j = 0; j < t.columns(); j++) {
				matrix[j][i] = t.entry(i + 1, j + 1);
			}
		}
		return new Matrix<>(matrix);
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
	public List<List<T>> nonTrivialCombinations(List<Matrix<T>> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Matrix<T> m : s) {
			asVectors.add(asVector(m));
		}
		return asFreeModule().nonTrivialCombinations(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<Matrix<T>> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (Matrix<T> m : s) {
			asVectors.add(asVector(m));
		}
		return asFreeModule().isGeneratingModule(asVectors);
	}

	@Override
	public String toString() {
		return "M(" + ring.toString() + ", " + rows + "x" + columns + ")";
	}
}
