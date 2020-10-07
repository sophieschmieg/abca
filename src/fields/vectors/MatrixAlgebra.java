package fields.vectors;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.polynomials.UnivariatePolynomial;
import fields.polynomials.UnivariatePolynomialRing;

public class MatrixAlgebra<T extends Element<T>> extends AbstractAlgebra<T, Matrix<T>> {
	private Ring<T> ring;
	private int dimension;
	private Matrix<T> id;
	private MatrixModule<T> module;

	MatrixAlgebra(Ring<T> ring, int dimension, FreeModule<T> free) {
		this.ring = ring;
		this.dimension = dimension;
		this.module = new MatrixModule<T>(ring, dimension, free, this);
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
		return module.zero();
	}

	@Override
	public Matrix<T> add(Matrix<T> s1, Matrix<T> s2) {
		return module.add(s1, s2);
	}

	@Override
	public Matrix<T> negative(Matrix<T> s) {
		return module.negative(s);
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
		return module.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public Matrix<T> getRandomElement() {
		return module.getRandomElement();
	}

	@Override
	public boolean isFinite() {
		return module.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return module.getNumberOfElements();
	}

	@Override
	public Iterator<Matrix<T>> iterator() {
		return module.iterator();
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
		return module.multiply(t1, t2);
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
		return module.rowShearOperation(iIndex, jIndex, value);
	}

	public Matrix<T> colShearOperation(int iIndex, int jIndex, T value) {
		return module.colShearOperation(iIndex, jIndex, value);
	}

	public MatrixModule<T>.GaussianEliminationResult gaussianElimination(Matrix<T> t) {
		return module.gaussianElimination(t);
	}

	public UnivariatePolynomial<T> characteristicPolynomial(Matrix<T> t) {
		UnivariatePolynomialRing<T> r = ring.getUnivariatePolynomialRing();
		MatrixAlgebra<Polynomial<T>> a = new FreeModule<>(r, dimension).matrixAlgebra();
		List<List<Polynomial<T>>> asPolynomial = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			asPolynomial.add(new ArrayList<>());
			for (int j = 0; j < dimension; j++) {
				asPolynomial.get(i).add(r.getEmbedding(t.entry(i + 1, j + 1)));
			}
		}
		return r.toUnivariate(a.determinant(a.subtract(a.getEmbedding(r.getVar(1)), new Matrix<>(asPolynomial))));
	}

	public T determinant(Matrix<T> t) {
		MatrixModule<T>.GaussianEliminationResult gauss = gaussianElimination(t);
		T result = ring.one();
		if (gauss.getSwaps() % 2 == 1) {
			result = ring.negative(result);
		}
		for (int i = 0; i < dimension; i++) {
			result = ring.multiply(result, gauss.getDiagonalMatrix().entry(i + 1, i + 1));
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
		return module.multiply(t, x);
	}

	public List<Vector<T>> kernelBasis(Matrix<T> t) {
		return module.kernelBasis(t);
	}

	public boolean hasSolution(Matrix<T> t, Vector<T> b) {
		return module.hasSolution(t, b);
	}

	public Vector<T> solve(Matrix<T> t, Vector<T> b) {
		return module.solve(t, b);
	}

	@Override
	public boolean isUnit(Matrix<T> t) {
		return ring.isUnit(determinant(t));
	}

	@Override
	public Matrix<T> inverse(Matrix<T> t) {
		MatrixModule<T>.GaussianEliminationResult gauss = gaussianElimination(t);
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					matrix[i][j] = ring.inverse(gauss.getDiagonalMatrix().entry(i + 1, j + 1));
				} else {
					matrix[i][j] = ring.zero();
				}
			}
		}
		return multiply(gauss.getColOperationsInverse(), new Matrix<T>(matrix), gauss.getRowOperationsInverse());
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
	public Matrix<T> projectToUnit(Matrix<T> t) {
		return one();
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

	public int rank(Matrix<T> m) {
		return module.rank(m);
	}

	public MatrixModule<T> asModule() {
		return module;
	}

	@Override
	public List<Matrix<T>> getModuleGenerators() {
		return module.getModuleGenerators();
	}

	@Override
	public Vector<T> asVector(Matrix<T> m) {
		return module.asVector(m);
	}

	public FreeModule<T> asFreeModule() {
		return module.asFreeModule();
	}

	@Override
	public boolean isLinearIndependent(List<Matrix<T>> s) {
		return module.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Matrix<T>> s) {
		return module.isGeneratingModule(s);
	}

	@Override
	public Ideal<Matrix<T>> getIdeal(List<Matrix<T>> generators) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public Ideal<Matrix<T>> intersect(Ideal<Matrix<T>> t1, Ideal<Matrix<T>> t2) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public int krullDimension() {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public boolean isGeneratingAlgebra(List<Matrix<T>> s) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public List<Matrix<T>> getAlgebraGenerators() {
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(zero().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				matrix[i][j] = ring.zero();
			}
		}
		List<Matrix<T>> generators = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			matrix[i][i] = ring.one();
			generators.add(new Matrix<T>(matrix));
			matrix[i][i] = ring.zero();
			if (i > 0) {
				matrix[i - 1][i] = ring.one();
				generators.add(new Matrix<T>(matrix));
				matrix[i - 1][i] = ring.zero();
				matrix[i][i - 1] = ring.one();
				generators.add(new Matrix<T>(matrix));
				matrix[i][i - 1] = ring.zero();
			}
		}
		return generators;
	}
}
